sourcename="QSO"
csourcename="QSO"
cubename="2206"
xloc=47.76
yloc=50.55
rmax=20
rmin=2
zPSFsize=150
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/q"$cubename"-qso_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize $zPSFsize -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/Q"$cubename"/zmask.txt"




#Continuum subtraction

$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/q"$cubename"-qso_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs_${rmax}.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 3'   -maskpix "1073 1159"




#Continuum subtraction only
#Continuum subtraction only


#$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs_${rmax}.PSFCONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3'




#$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/kcwi_oned/q"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs_${rmax}.PSFCONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 3'








#SNR cube generation
#SNR cube generation
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
#$CubEx/CubEx -cube "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_icubes_wcs.PSFCONTSub.fits"  -var  "/disk/bifrost/yuanze/KBSS/Q"cubename"/"sourcename"/kcwi_oned/q"cubename"-"sourcename"_vcubes_wcs.fits" -f .true. -fv .true. -fsr 2 -sn 1.e9 -cct "SNR_F"
