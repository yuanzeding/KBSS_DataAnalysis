import sys
sys.path.append('/disk/bifrost/yuanze/software/KcwiKit/py')
sys.path.append('/disk/bifrost/yuanze/KBSS/CubEx_run/scripts')
import run_cubetools_v1 as ctools
import os
from astropy.io import fits, ascii
import kcwi_tools
from astropy.cosmology import Planck18
# This script reprojects the KCWI cubes of faint AGN to cylindrical coordinates
# and saves the results in the specified directory.

root_directory="/disk/bifrost/yuanze/KBSS"
ascii_file_path = root_directory+'/KCWI/KBSS_faint_AGN.list'  # Path to the ASCII file
save_to_path=root_directory+"/faint_qsos/result"  # Path to save the output files
source_table = ascii.read(ascii_file_path, format='ipac')
dtype="KBSS"
filters = ["table['KCWI'] == 'yes'","table['Name'] != 'FSzP1170'","table['Name'] != 'Lab5'"]

all_directories,qsotab,data_dir = ctools.find_directories_from_ascii(source_table, root_directory,filters=filters,KBSS=(dtype=="KBSS"))
all_directories_red,_,data_dir_red = ctools.find_directories_from_ascii(source_table, root_directory,filters=filters,KBSS=(dtype=="KBSS"),channel="red")
overwrite=True
dryrun=False
for ind,sentry in enumerate(qsotab):
    #sentry=qsos_bright[qsos_bright["Field"]==field]
    #cubename=sentry["Field"].value[0]
    Field=sentry["Field"]
    Type=sentry["Type"]
    psname=sentry["Name"]
    csname=sentry["CName"]
    dapath=data_dir[ind]
    dapath_red=data_dir_red[ind]
    subdapath=all_directories[ind]
    subdapath_red=all_directories_red[ind]
    xpix = sentry["x"]
    ypix = sentry["y"]
    if Type<1.9: 
        Subfile = subdapath+"/{}-{}_icubes.PSFCONTSub.fits".format(Field ,psname)
        Subfile_red = subdapath_red+"/{}-{}-red_icubes.PSFCONTSub.fits".format(Field ,psname)
        maskfn=subdapath+"/{}-{}_icubes.PSFSub.mask.fits".format(Field ,psname)
        #maskfn_red=subdapath+"/{}-{}-red_icubes.PSFSub.mask.fits".format(Field ,psname)
        if not os.path.exists(Subfile):
            Subfile = subdapath+"/{}-{}_icubes_wcs.PSFCONTSub.fits".format(Field ,psname)
            Subfile_red = subdapath_red+"/{}-{}-red_icubes.PSFCONTSub.fits".format(Field ,psname)
    else:
        Subfile = subdapath+"/{}-{}_icubes.CONTSub.fits".format(Field ,psname)
        Subfile_red = subdapath_red+"/{}-{}-red_icubes_wcs.CONTSub.fits".format(Field ,psname)
        maskfn=subdapath+"/{}-{}_icubes.CONTSub.mask.fits".format(Field ,psname)
        #maskfn_red=subdapath+"/{}-{}-red_icubes.CONTSub.mask.fits".format(Field ,psname)
        if not os.path.exists(Subfile):
            Subfile = subdapath+"/{}-{}_icubes_wcs.CONTSub.fits".format(Field ,psname)
            Subfile_red = subdapath_red+"/{}-{}-red_icubes.CONTSub.fits".format(Field ,psname)
    writefn=subdapath+"/{}-{}_icubes.cyli.fits".format(Field ,psname)
    writefn_red=subdapath_red+"/{}-{}-red_icubes.cyli.fits".format(Field ,psname)
    if not os.path.exists(writefn) or overwrite == True:
        print("reprojecting",Field,psname,"to cylindrical system...")
        print("reading:",Subfile)
        print("writing result to ",writefn)
        if dryrun:
            continue
        hdu=fits.open(Subfile)
        hdu[0].header["CNAME3"]='KCWI Wavelength'
        hdr=hdu[0].header
        status=kcwi_tools.cart2cyli(Subfile,[xpix,ypix],hdr=hdr,masktype='direct',cos=Planck18,r_range=[0,10],dr=0.3,maskfn=maskfn,writefn=writefn,montage=True)
        hdu.close()
    else:
        print("reprojected blue cube existed and overwrite == False: no stacking performed ")
        
    if os.path.exists(Subfile_red):#and (not os.path.exists(writefn_red)):
        print("reading:",Subfile_red)
        print("writing result to ",writefn_red)
        hdu=fits.open(Subfile_red)
        hdu[0].header["CNAME3"]='KCWI Wavelength'
        hdr=hdu[0].header
    #status=kcwi_tools.cart2cyli(Subfile,[xpix,ypix],hdr=hdr,masktype='direct',cos=Planck18,r_range=[0,10],dr=binwidth/kpc_per_arcsec,maskfn=maskfn,writefn=writefn,montage=True)
        status=kcwi_tools.cart2cyli(Subfile_red,[xpix,ypix],hdr=hdr,masktype='direct',cos=Planck18,r_range=[0,10],dr=0.3,maskfn=maskfn,writefn=writefn_red,montage=True)
        hdu.close()
        #nhdu=kcwi_tools.cart2cyli(hdu[0],[xpix,ypix],r_range=[0,15],writefn=writefn_red)
    else:
        print("No red channel data for",Field,psname, " or file exists and overwrite == False")
    
