sourcename="BX58"
csourcename="BX58"
cubename="Q0105"
xloc=47.35
yloc=52.20
rmin=1
rmax=8
noisethr=3
maskpix='10 15'
#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwib/"$cubename"-"$csourcename"_icubes.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize 150 -withvar .false. -rmin $rmin -rmax $rmax -recenter .false. -masklist "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/zmask.txt" -noisethr $noisethr -psfimout "psfout.fits"

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/"$cubename"-"$sourcename"_icubes.PSFSub.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/"$cubename"-"$sourcename"_icubes.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 2' -maskpix "$maskpix"
