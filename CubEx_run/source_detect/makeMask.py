from astropy.io import ascii, fits
from astropy import constants, units as u
import os
import numpy as np
# Variables

# Path to the KBSS data
KBSSpath="/disk/bifrost/yuanze/KBSS"
qsos = ascii.read(KBSSpath+"/KCWI/qsos_bright.kcwi",format="ipac")
tab=qsos[(qsos['contam']=="False")&(qsos['Field']!="Q1623")]

for i, element in enumerate(tab):
    #field = element['FILENAME'].split('-')[0]
    field=element['Field']
    if field != "Q1623":
        objname="qso"
        subdapath=KBSSpath+"/"+field+"/"+objname.upper()
        fn=subdapath+"/{}-{}_icubes_wcs.cyli.fits".format(field.lower(),objname)
    else:
        continue

    hdu = fits.open(fn)[0]
    hdr = hdu.header

    #wc = (np.arange(hdr['NAXIS1']) - hdr['CRPIX1'] + 1) * hdr['CD1_1'] + hdr['CRVAL1']
    wc = (np.arange(hdr['NAXIS3']) - hdr['CRPIX3'] + 1) * hdr['CD3_3'] + hdr['CRVAL3']
    vc = (wc / (1 + element['z_sys']) - 1215.67) / 1215.67 *3e5
    wdir=f"{KBSSpath}/CubEx_run/{field}"
    sentry=qsos[qsos["Field"]==field]
    redshift=sentry["z_sys"]


    line="Lya"
    wave={"Lya":1215.67,"CIV":1549.06,"MgII":2799.12,"HeII":1640.4,"OVI":1031,"SiIV":1402.77,"CII":1334.53,"NV":1240,"OI":1304} #Angstroms
    sigma_v={"Lya":2.e3,"CIV":3e3,"SiIV":3e3,"CII":3e3,"NV":3e3,"OI":3e3,"OVI":3e3,"HeII":3e3} #km/s turbulent velocity
    index=(vc>(-sigma_v[line])) & (vc<sigma_v[line])
    inda=np.arange(len(wc))
    lmin=inda[index][0] #Lyman-alpha wavelength range
    lmax=inda[index][-1] #Lyman-alpha wavelength range, optical
    print("with redshift",redshift,"min wavelength index",lmin,"max wavelength index",lmax)
    print("doing",field)
    wdir = f"{KBSSpath}/CubEx_run/{field}"

    # Read and modify the mask file
    mask_file = os.path.join(wdir, "zmask.txt")
    if os.path.exists(mask_file):
        with open(mask_file, "r") as file:
            lines = file.readlines()

        # Ensure the first two lines are the original ones
        if lines[0].strip() != "144 193" or lines[2].strip() != "2343 2392":
            print(f"Warning: Unexpected mask file format in {mask_file}")
        else:
            # Modify the mask file
            with open(mask_file, "w") as file:
                file.write("144 193\n")
                file.write(f"{lmin} {lmax}\n")
                file.write("2343 2392\n")
                
                for line in lines[3:]:
                    file.write(line)
    else:
        # Create the mask file if it does not exist
        with open(mask_file, "w") as file:
            file.write("144 193\n")
            file.write(f"{lmin} {lmax}\n")
            file.write("2343 2392\n")

    # Modify the runcubetool.sh script
    runcubetool_script = os.path.join(wdir, "runcubtool_QSO.sh")
    if os.path.exists(runcubetool_script):
        with open(runcubetool_script, "r") as file:
            lines = file.readlines()

        with open(runcubetool_script, "w") as file:
            for line in lines:
                if "-bpsize '1 1 150' -bfrad '0 0 3'" in line:
                    # Ensure maskpix is not added twice
                    if "-maskpix" in line:
                        line = line.split("-maskpix")[0] + f'-maskpix "{lmin} {lmax}"\n'
                    else:
                        line = line.strip() + f'-maskpix "{lmin} {lmax}"\n'
                if ("-out" in line) and ("CubeBKGSub" in line):
                    # Update the output file name to include rmax value
                    parts = line.split(" ")
                    for i, part in enumerate(parts):
                        if part.startswith("-out"):
                            parts[i + 1] = '"/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs_${rmax}.PSFCONTSub.fits"'
                            break
                    line = " ".join(parts)
                if ("-out" in line) and ("CubePSFSub" in line):
                    # Update the output file name to include rmax value
                    parts = line.split(" ")
                    for i, part in enumerate(parts):
                        if part.startswith("-out"):
                            parts[i + 1] = '"/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs.PSFSub.fits"'
                            break
                    line = " ".join(parts)
                file.write(line)
    else:
        print(f"Warning: {runcubetool_script} does not exist.")