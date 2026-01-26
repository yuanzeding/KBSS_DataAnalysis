sourcename=M15
csourcename=M15
cubename=LyC22
maskpix='1116 1251'
#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwib/"$cubename"-"$csourcename"_icubes_wcs.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/"$cubename"-"$sourcename"_icubes.CONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 2' -maskpix "$maskpix"
