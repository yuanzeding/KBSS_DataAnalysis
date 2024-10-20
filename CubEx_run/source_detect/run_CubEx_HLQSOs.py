import subprocess
from astropy.io import ascii
from astropy import constants, units as u
# Variables
field = "0449"
radius = 6
# Path to the KBSS data
KBSSpath="/disk/bifrost/yuanze/KBSS"
wdir=f"{KBSSpath}/CubEx_run/Q{field}"
qsos = ascii.read(KBSSpath+"/KCWI/qsos_bright_updated.list",format="ipac")
sentry=qsos[qsos["Field"]=="Q"+field]

redshift=sentry["z_sys"]


line="Lya"
wave={"Lya":1215.67,"CIV":1549.06,"MgII":2799.12,"HeII":1640.4,"OVI":1031,"SiIV":1402.77,"CII":1334.53,"NV":1240,"OI":1304} #Angstroms
sigma_v={"Lya":3.e3,"CIV":3e3,"SiIV":3e3,"CII":3e3,"NV":3e3,"OI":3e3,"OVI":3e3,"HeII":3e3} #km/s turbulent velocity
c=constants.c.to(u.km/u.s).value #km/s speed of light
lmin=int((redshift+1.)*wave[line]*(1-sigma_v[line]/c)) #Lyman-alpha wavelength range
lmax=int((redshift+1.)*wave[line]*(1+sigma_v[line]/c)) #Lyman-alpha wavelength range, optical
print("with redshift",redshift,"min wavelength",lmin,"max wavelength",lmax)

RescaleVar = ".true."
RescaleVarArea = '"25 75 25 75"'
RescaleVarFR = 5
RescalingVarOutFile = '"Revar.out"'

MinNVox=20
MinArea = 5
MinDz = 5
SN_Threshold = 3.0
ApplyFilter = ".true."
ApplyFilterVar = ".true."
FilterXYRad = 1
FilterZRad = 0
XYedge = 5
ReliabCheck = ".true."
# Path to the CubEx executable
cubex_path = "/disk/bifrost/yuanze/software/CubEx/exe_files"

# Construct the new InpFile line


Inp_file = f'InpFile = "/disk/bifrost/yuanze/KBSS/Q{field}/QSO/q{field}-qso_icubes_wcs_{int(radius/0.3)}.PSFCONTSub.fits"'
Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/Q{field}/QSO/q{field}-qso_{line}.cat"'
Var_file = f'VarFile = "/disk/bifrost/yuanze/KBSS/Q{field}/QSO/q{field}-qso_vcubes.fits"'
CheckCube = f'CheckCube = "/disk/bifrost/yuanze/KBSS/Q{field}/QSO/q{field}-qso_icubes_wcs.PSFCONTSub.Objects_Id_{line}.fits"'
#Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/Q{field}/QSO/kcwi_oned/q{field}-QSO.cat"'
# Read in the current contents of the file

with open(f'{wdir}/par.in', 'r') as file:
    lines = file.readlines()

# Replace the line containing InpFile
with open(f'{wdir}/par.in', 'w') as file:
    for line in lines:
        if line.startswith('RescaleVar ='):
            line = f'RescaleVar = {RescaleVar}\n'
        if line.startswith('RescaleVarArea ='):
            line = f'RescaleVarArea = {RescaleVarArea}\n'
        if line.startswith('RescaleVarFR ='):
            line = f'RescaleVarFR = {RescaleVarFR}\n'
        if line.startswith('RescalingVarOutFile ='):
            line = f'RescalingVarOutFile = {RescalingVarOutFile}\n'
        if line.startswith('FilterZRad ='):
            line = f'FilterZRad = {FilterZRad}\n'
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
        if line.startswith('MinArea ='):
            line = f'MinArea = {MinArea}\n'
        if line.startswith('MinDz ='):
            line = f'MinDz = {MinDz}\n'
        if line.startswith('CheckCube ='):
            line = CheckCube + '\n'
        if line.startswith('XYedge ='):
            line = f'XYedge = {XYedge}\n'
        if line.startswith('ReliabCheck ='):
            line = f'ReliabCheck = {ReliabCheck}\n'
        file.write(line)

# Execute the CubEx command
#subprocess.run([f"{cubex_path}/CubEx", "par.in"])