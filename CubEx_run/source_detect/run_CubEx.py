import subprocess
from astropy.io import ascii
from astropy import constants, units as u
# Variables
sourcename = "BX58"
#csourcename = "BX432"
cubename = "Q0105"

# Path to the KBSS data
KBSSpath="/disk/bifrost/yuanze/KBSS"
wdir=f"{KBSSpath}/CubEx_run/{cubename}_{sourcename}"
qsos = ascii.read(KBSSpath+"/KCWI/qsos.kcwi",format="ipac")
sentry=qsos[(qsos["Name"]==sourcename) & (qsos["Field"]==cubename)]
csourcename = sentry["Cube"][0]
if sentry["zneb"]>0:
    redshift=sentry["zneb"]
elif sentry["zabs"]>0:
    redshift=sentry["zabs"]
else:
    redshift=sentry["zlya"]

line="Lya"
wave={"Lya":1215.67,"CIV":1549.06,"MgII":2799.12,"HeII":1640.4,"OVI":1031,"SiIV":1402.77,"CII":1334.53,"NV":1240,"OI":1304} #Angstroms
sigma_v={"Lya":1.1e4,"CIV":3e3,"SiIV":3e3,"CII":3e3,"NV":3e3,"OI":3e3,"OVI":3e3,"HeII":3e3} #km/s turbulent velocity
c=constants.c.to(u.km/u.s).value #km/s speed of light
lmin=int((redshift+1.)*wave[line]*(1-sigma_v[line]/c)) #Lyman-alpha wavelength range
lmax=int((redshift+1.)*wave[line]*(1+sigma_v[line]/c)) #Lyman-alpha wavelength range, optical
print("with redshift",redshift,"min wavelength",lmin,"max wavelength",lmax)

MinNVox=20
SN_Threshold = 3.
ApplyFilter = ".true."
ApplyFilterVar = ".true."
FilterXYRad = 1
# Path to the CubEx executable
cubex_path = "/disk/bifrost/yuanze/software/CubEx/exe_files"

# Construct the new InpFile line

if sentry["Type"][0]=="QSO":
    Inp_file = f'InpFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/kcwi_oned/{cubename}-{sourcename}_icubes.PSFCONTSub.fits"'
else:
    Inp_file = f'InpFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/kcwi_oned/{cubename}-{sourcename}_icubes.CONTSub.fits"'
Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/kcwi_oned/{cubename}-{sourcename}_{line}.cat"'
Var_file = f'VarFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{csourcename}/kcwi_oned/{cubename}-{csourcename}_vcubes.fits"'
CheckCube = f'CheckCube = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/kcwi_oned/{cubename}-{sourcename}_icubes.PSFCONTSub.Objects_Id_{line}.fits"'
#Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/Q{cubename}/{sourcename}/kcwi_oned/q{cubename}-{sourcename}.cat"'
# Read in the current contents of the file
with open(f'{wdir}/par.in', 'r') as file:
    lines = file.readlines()

# Replace the line containing InpFile
with open(f'{wdir}/par.in', 'w') as file:
    for line in lines:
        if line.startswith('ApplyFilterVar ='):
            line = f'ApplyFilterVar = {ApplyFilterVar}\n'
        if line.startswith('FilterXYRad ='):
            line = f'FilterXYRad = {FilterXYRad}\n'
        if line.startswith('ApplyFilter ='):
            line = f'ApplyFilter = {ApplyFilter}\n'
        if line.startswith('SN_Threshold ='):
            line = f'SN_Threshold = {SN_Threshold}\n'
        if line.startswith('InpFile ='):
            line = Inp_file + '\n'
        if line.startswith('Catalogue ='):
            line = Catalogue + '\n'
        if line.startswith('VarFile ='):
            line = Var_file + '\n'
        if line.startswith('lmin ='):
            line = f'lmin = {lmin}\n'
        if line.startswith('lmax ='):
            line = f'lmax = {lmax}\n'
        if line.startswith('MinNVox ='):
            line = f'MinNVox = {MinNVox}\n'
        if line.startswith('CheckCube ='):
            line = CheckCube + '\n'
        file.write(line)

# Execute the CubEx command
#subprocess.run([f"{cubex_path}/CubEx", "par.in"])