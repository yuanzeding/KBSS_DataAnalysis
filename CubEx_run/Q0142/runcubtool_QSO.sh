sourcename="QSO"
csourcename="QSO"
cubename="0142"
xloc=55.86
yloc=54.40
rmax=20
rmin=2
zPSFsize=150
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/q"$cubename"-qso_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/q"$cubename"-qso_icubes_wcs.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize $zPSFsize -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/Q"$cubename"/zmask.txt"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/q"$cubename"-qso_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/q"$cubename"-qso_icubes_wcs.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 3'

#Continuum subtraction only

#$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/kcwi_oned/q"$cubename"-"$sourcename"_icubes_wcs.CONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3'

#SNR cube generation
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
