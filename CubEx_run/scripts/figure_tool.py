import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import datetime
import astropy.units as u
from astropy import constants
from astropy.io import fits,ascii
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import astropy.visualization as vis
import os,sys
from cwitools import reduction, extraction, synthesis, utils, coordinates
from astropy.cosmology import Planck18
sys.path.append('/disk/bifrost/yuanze/KBSS/CubEx_run/scripts')
import run_cubetools_v1 as ctools
from regions.shapes.circle import CircleSkyRegion,CirclePixelRegion
from regions import PixCoord
import scienceplots
plt.style.use(["science","no-latex"])
lines={'NV':[1238.821,1242.804],'CIII]':[1909], 'CIV':[1548.20,1550.77],"OI":[1304,1306] ,"SiIII":[1206.5],"SiIV":[1393.76,1402.77],"SiII":[1260.42],"CII":[1334.53],"OVI":[1031,1037],"HeII":[1640],"Lya":[1215.67]}

def make_1D_spec(wvl,spec,var,sentry,MakeQSOspec=False,redshift=2.5, line_model=None
                , linename=None,postfix=""
                , save_to_path="./"
                , Org_spec=None, popt=None, pcov=None
                , useVelocity=False, colorlst=["#CC6677","#DDCC77","#6699CC"]):
    
    regu = r"$f_{\lambda}\rm (10^{-16}~erg~s^{-1}~cm^{-2}~\AA^{-1})$"
    plt.rcParams.update({
    "font.family": "dejavu sans",   # specify font family here
    "font.serif": ["Times"],  # specify font here
    "font.size":8,
    "legend.fancybox":True})
    fig, ax = plt.subplots(figsize=(5,3),dpi=300)
    sname= sentry["Name"]
    field = sentry["Field"]

    #xlow=Sub_cube.spectral_extrema[0].value/(1+redshift)
    #xhigh=Sub_cube.spectral_extrema[1].value/(1+redshift)
    cw=lines[linename][0]  # rest wavelength of the targeted line
    if useVelocity:
        x = (wvl.value / (1 + redshift) - cw) / cw *3e5
        xlow=-3000
        xhigh=3000
        ax.axvline(x=0,ls="--",lw=1,c="k")
    else:
        x = wvl.value/(1+redshift)
        #if postfix=="_red":
         #   xlow = 5680/(1+redshift)
        #    xhigh = 8640/(1+redshift)
        if MakeQSOspec:
            xlow = 850
            xhigh = 1600
        else:
            xlow=(cw-20)
            xhigh=(cw+20)
    keys=lines.keys()
    colors=matplotlib.cm.hsv(np.arange(len(keys))/len(keys))


    # redshift of the CGM being probed. QSO sight-line
    zcgm = redshift
    rlinec= (zcgm+1)/(redshift+1)
    mask=(x>xlow) & (x<xhigh)

    if MakeQSOspec:
        yup=np.max(Org_spec[mask])#np.max(Org_spec[mask])
    else:
        yup=np.max(spec[mask])
    #yup=np.max(Org_spec[~mask][:1000])
    ax.set_ylim(-yup*1e-1, yup*1.3)
    #Marking the important lines
    #ax.axvline(rlinec*1215.67,ls=":")
    #ax.text(rlinec*1215.67, yup/1.3,r"Ly$\alpha$",rotation=90,verticalalignment='bottom',horizontalalignment='center',fontsize=10)
    dlam=7.5/2800


    #xlow=1500
    #xhigh=1580
    if useVelocity:
        ax.axvline(x=0,ls="--",lw=1,c="k")
    else:
        for ind,key in enumerate(keys):
            for line in lines[key]:
                if line < xhigh and line > xlow:
                    ax.axvline(rlinec*line,ls=":",c=colors[ind])
        #            masktmp = (x.value < ((1+redshift)*line*(1+dlam))) & (x.value > ((1+redshift)*line*(1-dlam)))
        #            mask= mask | masktmp
                # print(x[mask])
            if line < xhigh and line > xlow:
                ax.text(rlinec*line, yup*1.1,key,rotation=90,verticalalignment='bottom',horizontalalignment='center')

    if MakeQSOspec:
        ax.plot(x[mask], Org_spec[mask],"-", label = "Flux") 
    else:
        ax.plot(x[mask], spec[mask],"-", label = "Flux")
    
    if line_model is not None and popt is not None:
        plt.plot(x , line_model(wvl.value, *popt), linestyle='dashed')
        peak=np.max(line_model(wvl.value, *popt))
        snr=peak/np.sqrt(np.mean(var[mask]))
        factor=2.355*3e5/lines[linename][0]  # convert to km/s
        ax.text(0.01,0.90\
                ,"FWHM={:.1f}$\pm${:.1f} km/s".format(np.abs(popt[2])*factor,np.sqrt(pcov[2,2])*factor),fontsize=12,color=colorlst[0],verticalalignment='top',horizontalalignment='left',transform=ax.transAxes)
        ax.text(0.01,0.83\
                ,"SNR={:.1f}".format(snr),fontsize=12,color=colorlst[0],verticalalignment='top',horizontalalignment='left',transform=ax.transAxes)
    #ax.plot(x[~mask]/(1+redshift), spec1[~mask],"b:", label = "Flux from peak")
    #ax.fill_between(x.value/(1+redshift), spec.value - sigma, spec.value + sigma, color='cyan', alpha=0.2, label='1-sigma region')
    ax.plot(x[mask], np.sqrt(var[mask]), label = "Std.dev.")
    ax.axhline(0,color="r")
    ax.legend(loc="upper right")
    if useVelocity:
        ax.set_xlabel(r"Velocity ($\rm km~s^{-1}$)")
    else:
        ax.set_xlabel(r"Rest-frame wavelength $(\rm \AA)$")
    ax.set_ylabel(regu)
    if ("aka" not in sentry.keys()):
        ax.text(0.01,0.98\
            ,"{}-{}".format(field,sname),fontsize=12,color=colorlst[0],verticalalignment='top',horizontalalignment='left',transform=ax.transAxes)
    elif sentry["aka"][0] != "None":
        ax.text(0.01,0.98\
        ,"{}-{} ({})".format(field,sname,sentry["aka"][0]),color=colorlst[0],fontsize=12,verticalalignment='top',horizontalalignment='left',transform=ax.transAxes)
    #ax.set_title("{}".format(sname))
    
    ax.set_xlim([xlow,xhigh])

    if MakeQSOspec:
        fig.savefig(save_to_path+"/{}-QSO_{}-{}_1d{}.pdf".format(linename,field,sname,postfix))
    else:
        fig.savefig(save_to_path+"/{}neb_{}-{}_1d{}.pdf".format(linename,field,sname,postfix))
    return fig,ax

def spec2fits_1D(savetopath,sentry,spec,var,mask2d,header,postfix=""):
    sname= sentry["Name"]
    field = sentry["Field"]
    fn_spec = savetopath+'/'+field+"-"+sname+'_kcwi_oned'+postfix+'.fits'
    
    hdr1d = header.copy()
    #add hdr keywords:
    #hdr1d['segID'] = num
    #hdr1d['segname'] = directory+'/'+prefix+'_seg.fits'
    #hdr1d['sclfctr'] = scale
    #hdr1d['date'] = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
    hdr1d['totalpix'] = np.nansum(mask2d)
    hdr1d['ra'] = sentry["RA"]
    hdr1d['dec'] = sentry["Decl"]
    hdr1d['xpos'] = sentry["x"]
    hdr1d['ypos'] = sentry["y"]
    hdr1d['kbssnm'] = field+"-"+sname
    
    #save files:
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header=hdr1d
    primary_hdu.header["AUTHOR"] = "Yuanze Ding"
    primary_hdu.header["COMMENT"] = (
    "Primary HDU deliberately left with no data; "
    "see EXTNAME=SPEC, VARIANCE, MASK."
)
    hdu1d = fits.PrimaryHDU(spec,header=hdr1d)         
    # Extension 1 – the 1-D spectrum
    spec_hdu = fits.ImageHDU(
        data=spec.astype("float32"),  # keep file size modest
        name="SPEC"
    )
    spec_hdu.header["BUNIT"] = "10^(-8)erg/cm3/s"
    spec_hdu.header["COMMENT"] = "Extracted 1-D spectrum"
    spec_hdu.header["CRPIX1"] = hdr1d["CRPIX3"]
    spec_hdu.header["CRVAL1"] = hdr1d["CRVAL3"]
    spec_hdu.header["CD1_1"] = hdr1d["CDELT3"]
    # Extension 2 – variance array
    var_hdu = fits.ImageHDU(
        data=var.astype("float32"),
        name="VARIANCE"
    )
    var_hdu.header["BUNIT"] = "10^(-16)erg2/cm6/s2"
    # Extension 3 – 2-D mask that selected the spaxels used above
    mask_hdu = fits.ImageHDU(
        data=mask2d,               # uint8 → 0 or 1
        name="MASK"
    )
    mask_hdu.header["BUNIT"] = "dimensionless"
    mask_hdu.header["COMMENT"] = (
        "2-D pixel mask applied to the original 3-D IFU cube while extracting SPEC."
    )
    mask_hdu.header["MASKVAL"] = "1: included, 0: excluded"
    # ------------------------------------------------------------------
    # 3. Assemble and write to disk
    # ------------------------------------------------------------------
    hdul = fits.HDUList([primary_hdu, spec_hdu, var_hdu, mask_hdu])
    hdul.writeto(fn_spec, overwrite=True)
    print(f"Wrote {fn_spec} with {np.nansum(mask2d)} pixels included.")
    return 0