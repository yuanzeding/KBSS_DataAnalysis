sourcename="M15"
csourcename="M15"
cubename="LyC22"
xloc=51.0
yloc=49.0

#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwi_oned/"$cubename"-"$csourcename"_icubes.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 150 -withvar .false. -rmax 15 -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/${cubename}_${sourcename}/zmask.txt"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 3' -maskpix '1941 1968'

