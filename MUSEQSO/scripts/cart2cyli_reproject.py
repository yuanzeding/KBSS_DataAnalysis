import sys
sys.path.append('/disk/bifrost/yuanze/software/KcwiKit/py')
sys.path.append('/disk/bifrost/yuanze/KBSS/MUSEQSO')
import run_cubetools_MUSE as ctools
import os

from astropy import cosmology


import numpy as np
import gc



from astropy.io import fits,ascii


keep_vars = ['keep_vars','cosmology','KBSSpath','dovariance','clean','root_directory','filters','gc', 'os', 'np', 'fits', 'all_directories', 'source_table', 'ctools', 'overwrite']
    


KBSSpath="/disk/bifrost/yuanze/KBSS"



root_directory = KBSSpath+"/MUSEQSO"
source_table = ascii.read(root_directory+"/MUSEQSO_machine_readable_updated2.list",format="ipac")#QSOtab=qsos[(qsos['contam']=="False")&(qsos['Field']!="Q1623")]
filters = ["table['file_count'] < 2", "table['M_i'] > -29.6"]
#QSOtab=qsos[(qsos['Field']!="Q0142")&(qsos['Field']!="Q1623")]
all_directories,tab = ctools.find_directories_from_ascii(source_table,root_directory,filters=filters)
print("Number of directories found:",len(all_directories))
#gc.set_debug(gc.DEBUG_LEAK)

dovariance=False
overwrite=True
clean=True
for subdir in all_directories[60:]:
    import kcwi_tools
    adp_prefix = ctools.find_adp_fits_file(subdir)
    quasar_name = os.path.basename(subdir)
    #overwrite=True
    sentry = source_table[source_table['Quasar'] == quasar_name]
    Subfile = adp_prefix+".PSFCONTSub.fits"
    NSubfile = adp_prefix+".PSFSub.fits"
    psfcen = np.loadtxt(subdir+"/psfcen.txt")
    xpix = psfcen[0,1]
    ypix = psfcen[0,2]
    dr = 0.2 # radial increment of the cylindrical projection, arcsec
    writefn=adp_prefix+".PSFCONTSub.cyli.fits"
    maskfn=f"{adp_prefix}.mask.fits"
    if os.path.exists(maskfn):
        pass
    else:
        maskfn=''
        print("Warning: mask file not found for",quasar_name)
    if overwrite == True or (not os.path.exists(writefn)):
        print("reprojecting",subdir.split("/")[-1],"to cylindrical system...")
        print("writing:",writefn)
        hdu=fits.open(Subfile)
        #hdu2=fits.open(cubefile)
        hdu[0].header["CNAME3"]="Wavelength"
        hdr=hdu[0].header
        status=kcwi_tools.cart2cyli(Subfile,[xpix,ypix],hdr=hdr,masktype='direct',cos=cosmology.Planck18,ncomp=0,clean=clean,r_range=[0,10],dr=dr,maskfn=maskfn,writefn=writefn,montage=True)
        hdu.close()
     #   del hdu
    else:
        print(f"reprojected cube existed for {quasar_name:s}: no stacking performed ")
    if dovariance:
        writefn=adp_prefix+".PSFCONTSubvar.cyli.fits"
        if overwrite == True or (not os.path.exists(writefn)):
            print("reprojecting variance cube of",quasar_name,"to cylindrical system...")
            print("writing:",writefn)
            hdu=fits.open(Subfile)
            hdu[0].header["CNAME3"]="Wavelength"
            hdr=hdu[0].header
            status=kcwi_tools.cart2cyli(NSubfile,[xpix,ypix],hdr=hdr,masktype='direct',cos=cosmology.Planck18,ncomp=2,clean=clean,r_range=[0,10],dr=dr,maskfn=maskfn,writefn=writefn,montage=True)
            hdu.close()
        else:
            print(f"reprojected variance cube existed for {quasar_name:s}: no stacking performed ")

    for var in list(locals()):
        if var not in keep_vars:
            del locals()[var]
    gc.collect()
