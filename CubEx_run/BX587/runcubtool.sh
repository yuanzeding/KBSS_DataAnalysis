sourcename="BX587"
cubename="2343"
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.fits" -x 74.4 -y 31.7 -nbins -1 -zPSFsize 300 -withvar .false. -rmax 14 -zPSF_min 750 -zPSF_max 2340 -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/"$sourcename"/zmask_"$sourcename".txt"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.PSFCONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3' -maskpix '1101 1141'

#SNR cube generation
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
