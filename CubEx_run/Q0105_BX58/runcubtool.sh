sourcename="BX58"
csourcename="BX58"
cubename="Q0105"
xloc=48.03
yloc=51.82
rmin=2
rmax=8
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwi_oned/"$cubename"-"$csourcename"_icubes.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 150 -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/CubEx_run/${cubename}_${sourcename}/zmask.txt" -psfimout "psfout.fits"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwi_oned/"$cubename"-"$sourcename"_icubes.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 3' -maskpix '909 933'

