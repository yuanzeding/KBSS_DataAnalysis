sourcename=MD28
csourcename=MD28
cubename=Q1549
xloc=50.0
yloc=53.0
rmin=1
rmax=8
noisethr=3
maskpix='149 265'
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwir/"$cubename"-"$csourcename"-red_icubes_wcs_cleaned.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwir/"$cubename"-"$sourcename"-red_icubes.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 150 -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwir/zmask.txt" -noisethr $noisethr -psfimout "psfout.fits"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwir/"$cubename"-"$sourcename"-red_icubes.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwir/"$cubename"-"$sourcename"-red_icubes.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 2' -maskpix "$maskpix"
