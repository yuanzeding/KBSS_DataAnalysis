import os
import re
import shlex
import subprocess
from mpi4py import MPI
import shutil
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy import constants, units as u
from regions import CircleSkyRegion, PixCoord

import sys
import numpy as np
sys.path.append('/disk/bifrost/yuanze/KBSS/MUSEQSO/scripts')
import makeMask_MUSE as masktools

def read_regions(file_path,pix_scale=0.3):
    # Read the region file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Regex to extract information
    pattern = re.compile(r'fk5;circle\((.*?),(.*?),(.*?)i\) #color=(.*?)  width=2 text={(.*?)}')

    regions = []

    for line in lines:
        match = pattern.match(line)
        if match:
            ra, dec, radius, color, text = match.groups()
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='fk5')
            radius = float(radius) * pix_scale * u.arcsec  # assuming the radius is in arcseconds
            
            # Create a CircleSkyRegion
            region = CircleSkyRegion(center=coord, radius=radius)
            region.meta["label"]=text
            region.visual["color"]=color
            region.visual["linewidth"]=0.5
            regions.append(region)

    # Now `regions` contains all the parsed CircleSkyRegions
    # You can manipulate them or plot them using regions library

    # Example: Printing region details
    return regions

def write_fits_cube(data,header, outputname):
    primary = fits.PrimaryHDU(data, header=header)
    hdul = fits.HDUList([primary])
    hdul.writeto(outputname, overwrite=True)
    return 0
def maskplot(data,mask,output=None,center=None,artist=None, meta=None,atype="circle"):
    from astropy.visualization.mpl_normalize import simple_norm
    sky_mean, sky_median, sky_std = sigma_clipped_stats(data, sigma=3, maxiters=5)
    norm = simple_norm([sky_std, 10*sky_std], 'linear', percent=99.5)
    fig, ax = plt.subplots(figsize=(5,5),dpi=300)
    ax.imshow(data-sky_median,origin='lower',cmap='gray_r',norm=norm)
    #plt.imshow(img.cut_sigma_image,origin='lower',cmap='gray_r')
    ax.imshow(mask,origin='lower',cmap='Blues',alpha=mask.astype(float)*0.7)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if center is not None:
        ax.plot(center[1],center[0],"+",c="r",ms=10,lw=0.01)
    if (artist is not None) and (atype == "circle"):
        if len(artist)>1:
            for ind,a in enumerate(artist):
                if meta is not None:
                    ax.text(a.center[0], a.center[1]+3, meta[ind]["label"], fontsize=meta[ind]["fontsize"], color=a.get_edgecolor(), ha='center')
                ax.add_artist(a)
        else:
            ax.add_artist(artist)

    if output is not None:
        fig.savefig(output)
    plt.close()
    return 0
def create_circular_mask(shape, center, radius):
    """
    Create a circular mask of given shape, center, and radius.

    Parameters:
    - shape (tuple): The shape of the mask (height, width).
    - center (tuple): The center coordinates of the circle (y, x).
    - radius (float): The radius of the circle. unit: pixels

    Returns:
    - mask (ndarray): The circular mask.

    """
    # Create grid of coordinates
    Y, X = np.ogrid[:shape[0], :shape[1]]
    
    # Calculate distance from the center
    dist_from_center = np.sqrt((X - center[1])**2 + (Y - center[0])**2)
    
    # Create mask
    mask = dist_from_center <= radius
    return mask

def replace_parameters_in_file(file_path, parameters):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            for parameter, new_value in parameters.items():
                if re.match(rf'{parameter}=.+', line.strip()):
                    line = re.sub(rf'{parameter}=.+', f'{parameter}={new_value}', line)
            file.write(line)

def execute_script_with_env(file_path, env_setup_path):
    quoted_file_path = shlex.quote(file_path)
    try:
        command = f"source {env_setup_path} && bash {quoted_file_path}"
        subprocess.run(command, shell=True, executable='/bin/bash', check=True)
        print(f'Successfully executed: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error executing {file_path}: {e}')

def process_directory(subdir, parameters, env_setup_path,dryrun=False):
    file_path = os.path.join(subdir, 'runcubtool_QSO.sh')
    if os.path.isfile(file_path):
        print(f'Processing file: {file_path}')
        replace_parameters_in_file(file_path, parameters)
        if dryrun == False:
            execute_script_with_env(file_path, env_setup_path)

def find_directories_from_ascii(table, root_dir, filters=None,keyword="Quasar"):
    if filters:
        for condition in filters:
            table = table[eval(condition)]
    
    directories = []
    for quasar in table[keyword]:
        subdir = os.path.join(root_dir, quasar)
        if os.path.isdir(subdir):
            directories.append(subdir)
    return directories, table

def copy_standard_script_to_directories(directories, standard_script_path):
    for subdir in directories:
        target_path = os.path.join(subdir, 'runcubtool_QSO.sh')
        shutil.copy(standard_script_path, target_path)
        print(f'Copied {standard_script_path} to {target_path}')

def find_adp_fits_file(subdir):
    for file in os.listdir(subdir):
        if file.startswith('ADP') and file.endswith('.fits') and not any(suffix in file for suffix in ['PSFCONTSub', 'PSFSub','cutout','white']):
            # Ensure it matches the pattern ADP.yyyy-mm-ddThh:mm:ss.sss.fits exactly
            name_without_extension, _ = os.path.splitext(os.path.join(subdir, file))
            return name_without_extension
    return None

def get_wcs_pixel_position(fits_file, ra, dec):
    with fits.open(fits_file) as hdul:
        wcs = WCS(hdul[1].header).celestial
    
    sky_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    pixel_coord = sky_coord.to_pixel(wcs)
    return pixel_coord

if __name__ == "__main__":
    root_directory = '/disk/bifrost/yuanze/KBSS/MUSEQSO'  # Set this to the root directory you want to start the search from
    ascii_file_path = root_directory+'/meta/MUSEQSO_machine_readable_updated2.list'  # Path to the ASCII file
    standard_script_path = root_directory+'/scripts/runcubtool_MUSEQSO.sh'  # Path to the standard script
    env_setup_file = os.environ.get('HEADAS') + '/headas-init.sh'  # Path to the environment setup file
    # Example filters
    source_table = ascii.read(ascii_file_path, format='ipac')
    filters = ["table['file_count'] < 2","table['z_sys']>3.5","table['M_i']<-29.6"]
    overwrite = True
    copy_script = False
    dryrun = False
    update_param = False
    if dryrun:
        print("Dry run: parameters updated but no files will be executed.")
    if (~update_param) and copy_script:
        print("Warning: copy_script is True but update_param is False. The script will be copied but not updated.")
    # update parameters in the running file
    parameters_to_update = {
        'rmax': 30,
        'rmin': 2  # Add more parameters as needed
    }
    unprocessed=[]
    # MPI initialization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Find all relevant directories
    if rank == 0:
        all_directories,tab = find_directories_from_ascii(source_table, root_directory,filters=filters)
        if copy_script == True:
            copy_standard_script_to_directories(all_directories, standard_script_path)
    else:
        all_directories = None

    # Distribute directories among processes
    all_directories = comm.bcast(all_directories, root=0)
    directories_per_process = len(all_directories) // size
    extra = len(all_directories) % size

    if rank < extra:
        start = rank * (directories_per_process + 1)
        end = start + directories_per_process + 1
    else:
        start = rank * directories_per_process + extra
        end = start + directories_per_process

    local_directories = all_directories[start:end]
    
    # Process directories assigned to this rank
    for subdir in local_directories:
        adp_prefix = find_adp_fits_file(subdir)
        print("reading ADP file at",adp_prefix)
        if os.path.exists(adp_prefix+".PSFCONTSub.fits") and overwrite == False:
            print("PSFCONTSubbed file exists and overwrite==False!")
            continue
        #else:
            #unprocessed.append(subdir)
        if adp_prefix:
            # Find the corresponding row in the source table
            quasar_name = os.path.basename(subdir)
            source_row = source_table[source_table['Quasar'] == quasar_name]
            if len(source_row) < 0:
                print(f'Could not find the source row for {quasar_name}')
                break
            if update_param:
                ra = source_row['RA'][0]
                dec = source_row['Decl'][0]
                if source_row['cutout'] == 'True':
                    xloc,yloc=50.0,50.0
                else:
                    xloc, yloc = get_wcs_pixel_position(adp_prefix+".fits", ra, dec)
                parameters_to_update['xloc'] = int(xloc)
                parameters_to_update['yloc'] = int(yloc)
                parameters_to_update['inPrefix'] = shlex.quote(adp_prefix)
                parameters_to_update['masklist'] = shlex.quote(os.path.join(subdir,"zmask.txt"))
                parameters_to_update['outpsfcen'] = shlex.quote(os.path.join(subdir,"psfcen.txt"))
            else:
                print(f'Skipping parameter update for {quasar_name}')
            masktools.update_mask(subdir,source_table)
        process_directory(subdir, parameters_to_update, env_setup_file,dryrun=dryrun)
    #print(unprocessed)
        
