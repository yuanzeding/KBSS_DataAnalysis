inPrefix="ADP.2017-10-16T15:00:00.000"
xloc=50.0
yloc=50.0
rmax=30
rmin=5
zPSFsize=150
maskpix='10 15'
masklist='zmask.txt'
outpsfcen='psfcen.txt'
#Sky subtraction residual correction (skipped for this time)

#$CubEx/Tools/CubeFix 

#Quasar PSF subtraction
$CubEx/Tools/CubePSFSub -cube "${inPrefix}.fits" -out "${inPrefix}.PSFSub.fits" -x $xloc -y $yloc -nbins -1 -zPSFsize $zPSFsize -withvar .true. -rmin $rmin -rmax $rmax -recenter .true. -masklist "$masklist" -outpsfcen "$outpsfcen"

#Continuum subtraction

$CubEx/Tools/CubeBKGSub -cube  "${inPrefix}.PSFSub.fits" -out "${inPrefix}.PSFCONTSub.fits" -bpsize '1 1 150' -bfrad '0 0 2'  -maskpix "$maskpix" -outvarcube "${inPrefix}.PSFCONTSubvar.fits"