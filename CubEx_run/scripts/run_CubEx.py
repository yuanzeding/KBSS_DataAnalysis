import subprocess
from astropy.io import ascii
from astropy import constants, units as u
import sys
import numpy as np
sys.path.append('/disk/bifrost/yuanze/KBSS/CubeEx_run/scripts')
#import run_cubetools_MUSE as ctools
import run_cubetools_v1 as ctools
# Variables


# Path to the KBSS data
root_directory = '/disk/bifrost/yuanze/KBSS'  # Set this to the root directory you want to start the search from
ascii_file_path = root_directory+'/KCWI/KBSS_faint_AGN.list'  # Path to the ASCII file
standard_script_path = '/disk/bifrost/yuanze/KBSS/CubEx_run/scripts/par.in'  # Path to the standard script
filters = ["table['KCWI'] == 'yes'"]
source_table = ascii.read(ascii_file_path, format='ipac')
dtype="KBSS"
all_directories,tab,all_data_dir = ctools.find_directories_from_ascii(source_table, root_directory,filters=filters,KBSS=(dtype=="KBSS"))
KBSSpath="/disk/bifrost/yuanze/KBSS"

#qsos = ascii.read(KBSSpath+"/KCWI/KBSS_faint_AGN.list",format="ipac")
#sentry=qsos[(qsos["Name"]==sourcename) & (qsos["Field"]==cubename)]

emline="Lya"
wave={"Lya":1215.67,"CIV":1549.06,"MgII":2799.12,"HeII":1640.4,"OVI":1031,"SiIV":1402.77,"CII":1334.53,"NV":1240,"OI":1304} #Angstroms
sigma_v={"Lya":1.1e4,"CIV":3e3,"SiIV":3e3,"CII":3e3,"NV":3e3,"OI":3e3,"OVI":3e3,"HeII":3e3} #km/s turbulent velocity
c=constants.c.to(u.km/u.s).value #km/s speed of light


RescaleVar = ".true."
RescaleVarArea = '"35 65 35 65"'
RescaleVarFR = 5
RescalingVarOutFile = '"Revar.out"'

MinNVox=50
MinArea = 5
MinDz = 5
SN_Threshold = 3.0
ApplyFilter = ".true."
ApplyFilterVar = ".true."
FilterXYRad = 1
FilterZRad = 0
XYedge = 5
ReliabCheck = ".false."
# Path to the CubEx executable
cubex_path = "/disk/bifrost/yuanze/software/CubEx/exe_files"
instrument="kcwib"
# Construct the new InpFile line
ctools.copy_standard_script_to_directories(all_directories,standard_script_path,scriptname="par.in")
for ns,source_row in enumerate(tab):
    print(source_row)
    redshift=source_row["zlya"]
    sourcename = source_row["Name"]
    csourcename = source_row["CName"]
    #csourcename = "BX432"
    cubename = source_row["Field"]
    wdir=all_directories[ns]#f"{KBSSpath}/CubEx_run/{cubename}_{sourcename}"
    #if sentry["Type"][0]=="QSO":
    Inp_file = f'InpFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/{instrument}/{cubename}-{sourcename}_icubes.PSFCONTSub.fits"'
    #else:
    #    Inp_file = f'InpFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/{instrument}/{cubename}-{sourcename}_icubes.CONTSub.fits"'
    Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/{instrument}/{cubename}-{sourcename}_{emline}.cat"'
    Var_file = f'VarFile = "/disk/bifrost/yuanze/KBSS/{cubename}/{csourcename}/{instrument}/{cubename}-{csourcename}_vcubes.fits"'
    CheckCube = f'CheckCube = "/disk/bifrost/yuanze/KBSS/{cubename}/{sourcename}/{instrument}/{cubename}-{sourcename}_icubes.PSFCONTSub.Objects_Id_{emline}.fits"'
    #Catalogue = f'Catalogue = "/disk/bifrost/yuanze/KBSS/Q{cubename}/{sourcename}/kcwi_oned/q{cubename}-{sourcename}.cat"'
    # Read in the current contents of the file
    lmin=np.rint((redshift+1.)*wave[emline]*(1-sigma_v[emline]/c)) #Lyman-alpha observed frame wavelength range
    lmax=np.rint((redshift+1.)*wave[emline]*(1+sigma_v[emline]/c)) 
    print("with redshift",redshift,"min wavelength",lmin,"max wavelength",lmax)
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
    subprocess.run([f"{cubex_path}/CubEx", f"{wdir}/par.in"])