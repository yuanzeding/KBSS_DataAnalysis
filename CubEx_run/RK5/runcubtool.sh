sourcename="RK5"
csourcename="RK5"
cubename="0821"
xloc=37.18
yloc=52.53

#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 300 -withvar .false. -rmax 5 -zPSF_min 750 -zPSF_max 2340 -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/"$sourcename"/zmask_"$sourcename".txt"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFCONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3' -maskpix '629 650'

#Continuum subtraction only

$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.CONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3' -maskpix '629 650'

#SNR cube generation
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
