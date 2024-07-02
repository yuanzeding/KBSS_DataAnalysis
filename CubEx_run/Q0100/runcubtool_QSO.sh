sourcename="QSO"
csourcename="QSO"
cubename="0100"
xloc=50.46
yloc=50.41
rmax=20
rmin=2
zPSFsize=150
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$csourcename"/q"$cubename"-qso_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize $zPSFsize -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/Q"$cubename"/zmask.txt"

#Continuum subtraction

$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/Q"$cubename"/"$sourcename"/q"$cubename"-qso_icubes_wcs.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/Q${cubename}/${sourcename}/q${cubename}-qso_icubes_wcs_${rmax}.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 3'  -maskpix "1250 1339"
















