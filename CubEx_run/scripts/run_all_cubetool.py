import os
import re
import subprocess
from mpi4py import MPI

def replace_parameters_in_file(file_path, parameters):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            for parameter, new_value in parameters.items():
                if re.match(rf'{parameter}=\d+', line.strip()):
                    line = re.sub(rf'{parameter}=\d+', f'{parameter}={new_value}', line)
            file.write(line)

def execute_script_with_env(file_path, env_setup_path):
    try:
        command = f"source {env_setup_path} && bash {file_path}"
        subprocess.run(command, shell=True, executable='/bin/bash', check=True)
        print(f'Successfully executed: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error executing {file_path}: {e}')

def process_directory(subdir, parameters, env_setup_path):
    file_path = os.path.join(subdir, 'runcubtool_QSO.sh')
    if os.path.isfile(file_path):
        print(f'Processing file: {file_path}')
        replace_parameters_in_file(file_path, parameters)
        execute_script_with_env(file_path, env_setup_path)

def find_directories(root_dir):
    pattern = re.compile(r'^Q\d{4}$')
    directories = []
    for subdir, _, _ in os.walk(root_dir):
        if pattern.match(os.path.basename(subdir)):
            directories.append(subdir)
    return directories

def process_directories_QSO(root_dir, parameters, env_setup_path):
    pattern = re.compile(r'^Q\d{4}$')
    #pattern = re.compile(r'^Q0100')
    for subdir, _, _ in os.walk(root_dir):
        if pattern.match(os.path.basename(subdir)):
            file_path = os.path.join(subdir, 'runcubtool_QSO.sh')
            if os.path.isfile(file_path):
                print(f'Processing file: {file_path}')
                replace_parameters_in_file(file_path, parameters)
                execute_script_with_env(file_path, env_setup_path)

if __name__ == "__main__":
    root_directory = '/disk/bifrost/yuanze/KBSS/CubEx_run'  # Set this to the root directory you want to start the search from
    env_setup_file = os.environ.get('HEADAS') + '/headas-init.sh'  # Path to the environment setup file
    parameters_to_update = {
        'rmax': 20,
        'rmin': 2  # Add more parameters as needed
    }
    # MPI initialization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Find all relevant directories
    if rank == 0:
        all_directories = find_directories(root_directory)
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
        process_directory(subdir, parameters_to_update, env_setup_file)
