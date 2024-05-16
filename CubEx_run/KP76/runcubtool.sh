sourcename="KP76"
csourcename="KP76"
cubename="1623"
xloc=41.17
yloc=52.57

#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 300 -withvar .false. -rmax 10 -zPSF_min 750 -zPSF_max 2340 -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/"$sourcename"/zmask_"$sourcename".txt"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFCONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3'

#Continuum subtraction only

#$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.CONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3'

#SNR cube generation
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
