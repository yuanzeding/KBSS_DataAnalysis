import sys
sys.path.append('/disk/bifrost/yuanze/software/KcwiKit/py')
sys.path.append('/disk/bifrost/yuanze/KBSS/MUSEQSO/scripts')
import run_cubetools_MUSE as ctools
import os

from astropy import cosmology


import numpy as np
import gc



from astropy.io import fits,ascii


keep_vars = ['tab','keep_vars','cosmology','KBSSpath','dovariance','clean','root_directory','doQSO','filters','gc', 'os', 'np', 'fits', 'all_directories', 'source_table', 'ctools', 'overwrite']
    


KBSSpath="/disk/bifrost/yuanze/KBSS"



root_directory = KBSSpath+"/MUSEQSO"
source_table = ascii.read(root_directory+"/meta/MUSEQSO_machine_readable_updated2_withMi2.list",format="ipac")#QSOtab=qsos[(qsos['contam']=="False")&(qsos['Field']!="Q1623")]
filters = ["table['file_count'] < 2","table['M_i_z2']<-29.2","table['z_sys']>3.5"]
#QSOtab=qsos[(qsos['Field']!="Q0142")&(qsos['Field']!="Q1623")]
all_directories,tab = ctools.find_directories_from_ascii(source_table,root_directory,filters=filters)
print("Number of directories found:",len(all_directories))
#gc.set_debug(gc.DEBUG_LEAK)

dovariance=True
overwrite=True
clean=True
doQSO=False
#print(tab)
for n_dir,subdir in enumerate(all_directories): 
    import kcwi_tools
    adp_prefix = ctools.find_adp_fits_file(subdir)
    source_row = tab[n_dir]#source_table[source_table['Quasar'] == quasar_name]
    print("Processing directory:",subdir)
    print(adp_prefix)
    print(source_row)
    quasar_name = source_row["Quasar"]
    if doQSO:
        Subfile=adp_prefix+".fits"
        writefn=adp_prefix+".cyli.fits"
        dovariance=False
        subncomp=1
    else:
        Subfile = adp_prefix+".PSFCONTSub.fits"
        subncomp=0
        writefn=adp_prefix+".PSFCONTSub.cyli.fits"
    NSubfile = adp_prefix+".PSFSub.fits"
    psfcen = np.loadtxt(subdir+"/psfcen.txt")
    xpix = psfcen[0,1]
    ypix = psfcen[0,2]
    dr = 0.2 # radial increment of the cylindrical projection, arcsec
    maskfn=f"{adp_prefix}.mask.fits"
    if os.path.exists(maskfn):
        pass
    else:
        maskfn=''
        print("Warning: mask file not found for",quasar_name)
    if overwrite == True or (not os.path.exists(writefn)):
        print("reprojecting",quasar_name,"to cylindrical system...")
        print("writing:",writefn)
        hdu=fits.open(Subfile)
        #hdu2=fits.open(cubefile)
        hdu[subncomp].header["CNAME3"]="Wavelength"
        hdr=hdu[subncomp].header
        status=kcwi_tools.cart2cyli(Subfile,[xpix,ypix],hdr=hdr,masktype='direct',cos=cosmology.Planck18,ncomp=subncomp,clean=clean,r_range=[0,10],dr=dr,maskfn=maskfn,writefn=writefn,montage=True)
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
