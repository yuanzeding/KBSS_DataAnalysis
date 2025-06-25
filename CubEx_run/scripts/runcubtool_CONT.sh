sourcename="BX58"
csourcename="BX58"
cubename="Q0105"
maskpix='10 15'

#Continuum subtraction
$CubEx/Tools/CubeBKGSub -cube  "/disk/bifrost/yuanze/KBSS/"$cubename"/"$csourcename"/kcwib/"$cubename"-"$csourcename"_icubes.fits" -out "/disk/bifrost/yuanze/KBSS/"$cubename"/"$sourcename"/kcwib/"$cubename"-"$sourcename"_icubes.CONTSub.fits" -bpsize '1 1 50' -bfrad '0 0 2' -maskpix "$maskpix"
