from astropy.io import ascii, fits
from astropy import constants, units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import os
import numpy as np
import re
import sys
sys.path.append('/disk/bifrost/yuanze/KBSS/CubEx_run/scripts')
#import run_cubetools_MUSE as ctools
import run_cubetools_v1 as ctools
# Variables


# Path to the KBSS data

#tab=qsos[(qsos['contam']=="False")&(qsos['Field']!="Q1623")]


def get_wcs_pixel_position(fits_file, ra, dec):
    with fits.open(fits_file) as hdul:
        wcs = WCS(hdul[0].header)
    
    sky_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    pixel_coord = sky_coord.to_pixel(wcs)
    return pixel_coord
def determine_mask(wc,line_update_list, redshift,channel="blue"):
    wave = {"Lya": 1215.67, "CIV": 1549.06, "MgII": 2799.12, "HeII": 1640.4, "OVI": 1031, "SiIV": 1402.77, "CII": 1334.53,"CIII]":1909, "NV": 1240, "OI": 1304, "continuum_blue": 1280,"continuum_red":1800}  # Angstroms
    sigma_v = {"Lya": 1.e3, "CIV": 1.5e3,"CIII]": 1e3, "SiIV": 3e3, "CII": 3e3, "NV": 3e3, "OI": 3e3, "OVI": 3e3, "HeII": 1.5e3, "continuum_blue":1.5e3,"continuum_red":1.5e3}  # km/s turbulent velocity
    inda = np.arange(len(wc))
    lmin=[]#np.full(len(line_update_list),None)
    lmax=[]#np.full(len(line_update_list),None)
    for linen,line in enumerate(line_update_list):
        if line == "continuum":
            if channel == "blue":
                line = "continuum_blue"
            elif channel == "red":
                line = "continuum_red"
        vc = (wc / (1 + redshift) - wave[line]) / wave[line] * constants.c.to(u.km/u.s).value
        if vc[-1]<sigma_v[line] or vc[0]>-sigma_v[line]:
            continue # Skip if the line is not in the wavelength range
            if channel == "blue":
                line = "Lya"# Default to Lyman-alpha if the line is not found in wavelength array
            elif channel == "red":
                line = "CIV"# Default to CIV if the line is not found in wavelength array
            vc = (wc / (1 + redshift) - wave[line]) / wave[line] * constants.c.to(u.km/u.s).value
        index = (vc > (-sigma_v[line])) & (vc < sigma_v[line])
        lmin.append(inda[index][0])  # Lyman-alpha wavelength range
        lmax.append(inda[index][-1])  # Lyman-alpha wavelength range, optical
        print(f"Line: {line}, Min index: {lmin[-1]}, Max index: {lmax[-1]}")
#    except:
#        if line=="continuum":
#            vc = (wc / (1 + redshift) - wave["Lya"]) / wave["Lya"] * constants.c.to(u.km/u.s).value
#            index = (vc < sigma_v[line])#vc > sigma_v[line]
#            lmin = inda[index][0]  # Lyman-alpha wavelength range
#            lmax = inda[index][-1]  # Lyman-alpha wavelength range, optical
#            print(inda[index])
#        else:
#            IndexError
    return np.array(lmin),np.array(lmax)
def update_mask(subdir,tab,redshiftref="z_Lyaneb",linelist = ["Lya"],\
    subdir_data=None,dtype="MUSE",scriptname="runcubetool_QSO.sh",**kwargs):
    if dtype=="MUSE":
        adp_file = ctools.find_adp_fits_file(subdir)
    elif dtype=="KBSS":
        adp_file = ctools.find_adp_fits_file(subdir_data, dtype="KBSS")
    if adp_file:
        if dtype=="KBSS":
            redshift = tab[redshiftref]
            hdu = fits.open(adp_file+".fits")[0]
        else:
            hdu = fits.open(adp_file+".fits")[1]
            redshift = tab[redshiftref][0]
        hdr = hdu.header
        wc = (np.arange(hdr['NAXIS3']) - hdr['CRPIX3'] + 1) * hdr['CD3_3'] + hdr['CRVAL3']
        lmin,lmax = determine_mask(wc,linelist, redshift,**kwargs)
        print("with redshift", redshift, "min wavelength index", lmin, "max wavelength index", lmax)
        print("Done updating wavelength mask for data in", subdir)
        print("----------------------------------------------------------------------------------------------------------------")
        # Read and modify the mask file
        mask_file = os.path.join(subdir, "zmask.txt")
        with open(mask_file, "w") as file:
            #file.write("93 143\n")
            for i in range(len(lmin)):
                if lmin[i] is not None and lmax[i] is not None:
                    file.write(f"{lmin[i]} {lmax[i]}\n")
            # Add the fixed lines for Lyman-alpha and C IV
            #file.write("2307 2357\n")
        # Modify the runcubetool.sh script

        runcubetool_script = os.path.join(subdir, scriptname)
        
        if os.path.exists(runcubetool_script):
            parameter={"maskpix":f"'{np.nanmin(lmin)} {np.nanmin(lmax)}'"}
            ctools.replace_parameters_in_file(runcubetool_script,parameter)
        else:
            print(f"Warning: {runcubetool_script} does not exist.")

if __name__ == "__main__":
    root_directory = '/disk/bifrost/yuanze/KBSS/MUSEQSO'  # Set this to the root directory you want to start the search from
    ascii_file_path = root_directory+'/meta/MUSEQSO_machine_readable_updated2.list'  # Path to the ASCII file
    qsos = ascii.read(ascii_file_path,format="ipac")
    filters = ["table['file_count'] < 2","table['M_i'] < -29.6","table['z_sys']>3.5"]
    all_directories,tab=ctools.find_directories_from_ascii(qsos,root_directory,filters=filters)
    for subdir in all_directories:
        update_mask(subdir, qsos)
        