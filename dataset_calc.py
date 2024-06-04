#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:46:58 2023

@author: akshat
"""

import os 
import uuid
import time
import subprocess
import argparse
import multiprocessing    

def read_config_file(filename):
    """
    Reads configuration settings from a file and returns them as a dictionary.
    
    This function opens a specified configuration file and reads its contents line by line.
    Lines starting with a '#' are considered comments and are ignored, as are empty lines.
    Each valid line is expected to be in the format 'key=value'. If the 'value' is a digit,
    it is converted to an integer; otherwise, it is kept as a string.

    Parameters:
        filename (str): The path to the configuration file to be read.

    Returns:
        dict: A dictionary containing key-value pairs of configuration settings.
              Keys are the configuration parameter names and values are the corresponding settings.

    Example:
        Assume we have a config.txt with the following content:
            # Configuration for system
            NUM_CPUs=10
            USE_LOAD_BALANCER=True

        Calling read_config_file("config.txt") will return:
            {'NUM_CPUs': 10, 'USE_LOAD_BALANCER': 'True'}
    """
    params = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line and not line.startswith("#"):  # skip comments and empty lines
                key, value = line.split("=", 1)
                if value.isdigit():
                    params[key] = int(value)
                else:
                    params[key] = value
    return params

# Read the all.ctrl file to get parameters
config_params = read_config_file("all.ctrl")

RECEPTOR_LOCATION  = config_params.get("RECEPTOR_LOCATION", "./DATA/docking_receptor.pdbqt")
EXHAUSTIVENESS     = str( config_params.get("EXHAUSTIVENESS", "1"))
CENTER_X           = str( config_params.get("CENTER_X"))
CENTER_Y           = str( config_params.get("CENTER_Y"))
CENTER_Z           = str( config_params.get("CENTER_Z"))
SIZE_X             = str( config_params.get("SIZE_X"))
SIZE_Y             = str( config_params.get("SIZE_Y"))
SIZE_Z             = str( config_params.get("SIZE_Z"))
DOCKING_SCORE_THRS = float(config_params.get("DOCKING_SCORE_THRESHOLD"))

def generate_unique_file_name(base_name, extension):
    """
    Generates a unique file name using a base name, a timestamp, and a UUID.

    Parameters:
        base_name (str): The base name for the file, which will precede the unique identifiers.
        extension (str): The file extension to append to the file name.

    Returns:
        str: A uniquely named file with the format "base_name_timestamp_uuid.extension".

    Example:
        >>> generate_unique_file_name('pose', 'txt')
        'pose_1657892345678_abcdef1234567890abcdef.txt'
    """
    timestamp = int(time.time() * 1000)
    unique_id = uuid.uuid4().hex
    file_name = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    return file_name


def check_energy(lig_): 
    """
    Computes the total energy of a ligand using the Open Babel 'obenergy' tool.

    Parameters:
        lig_ (str): The file name of the ligand in PDBQT format.

    Returns:
        float: The computed total energy of the ligand in kcal/mol. Returns 10000 if the calculation fails.

    Example:
        >>> check_energy('ligand.pdbqt')
        -7.5
    """
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', lig_]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy


def run_docking(lig_location, out_location, method='qvina'):
    """
    Perform molecular docking using QuickVina or Smina specified by the method parameter.

    Parameters:
        lig_location (str): The file path to the ligand file.
        out_location (str): The file path where the output should be written.
        method (str, optional): Specifies the docking method, defaults to 'qvina'.

    Returns:
        float: The best docking score obtained from the docking simulation.

    Raises:
        ValueError: If an unsupported method is specified.
        subprocess.CalledProcessError: If the docking process fails.

    Example:
        >>> run_docking('ligand.pdbqt', 'output.pdbqt')
        -9.2
    """
    if method == 'qvina':
        command_run = subprocess.run(["./DATA/qvina", "--receptor", RECEPTOR_LOCATION, "--ligand", lig_location, "--center_x", CENTER_X, "--center_y", CENTER_Y, "--center_z", CENTER_Z, "--size_x", SIZE_X, "--size_y", SIZE_Y, "--size_z", SIZE_Z, "--exhaustiveness", EXHAUSTIVENESS, "--out", out_location], capture_output=True)
    else:
        raise ValueError('Unsupported docking method specified. Only "qvina" or "smina" are supported.')

    # Check the quality of the output pose:
    pose_energy = check_energy(out_location)
    if pose_energy == 10000:
        return 10000  # Indicates a failed docking process with a placeholder error value.

    # Extract and return the best docking score:
    docking_score = float('inf')
    for line in command_run.stdout.decode("utf-8").split('\n'):
        parts = line.strip().split()
        if len(parts) == 4 and parts[0].isdigit():
            score = float(parts[1])
            docking_score = min(docking_score, score)
    
    return docking_score


def perform_calc_single(args): 
    """
    Processes a single SMILES string to generate a 3D structure, check its stability, and perform docking.

    Parameters:
        args (tuple): A tuple containing a SMILES string and its index.

    Example:
        >>> perform_calc_single(('CCO', 0))
        Writes results to output files and performs cleanup.

    Note:
        This function assumes that necessary global variables and paths are correctly set.
    """
    out_location    = generate_unique_file_name('pose', 'pdbqt') # For the docking pose
    output_filename = generate_unique_file_name('lig', 'pdbqt')  # For the 3D ligand (obabel converted smi)
    
    try: 
        smi, job_idx = args
        cmd = ["obabel", "-ismi","-:" + smi,"-O", output_filename, "--gen3d", "fastest"]
        # subprocess.run(cmd, timeout=120)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(cmd, stdout=devnull, stderr=devnull, timeout=120)

        # Ensure a stable molecule: 
        lig_energy = check_energy(output_filename)
    
        # Specifying docking input file & output file: 
        lig_location = output_filename
        
        # Perform docking: 
        if lig_energy < 10000: 
            docking_score = run_docking(lig_location, out_location, method='qvina')
                
        if docking_score > DOCKING_SCORE_THRS: 
            if os.path.exists(out_location):
                os.system('rm {}'.format(out_location))
            if os.path.exists(output_filename):
                os.system('rm {}'.format(output_filename))
            with open('./OUTPUT_{}.txt'.format(job_idx), 'a+') as f: 
                f.writelines(['{}, {}\n'.format(smi, docking_score)]) 
        else: 
            with open('./OUTPUT_{}.txt'.format(job_idx), 'a+') as f: 
                f.writelines(['{}, {}\n'.format(smi, docking_score)]) 
                
            # Make the directory if it does not exist 
            os.makedirs("OUTPUTS", exist_ok=True)
            
            # Move the files lig_location, out_location to the directory OUTPUTS/ there is no need to keep the file in the original location anymore
            os.system('mv {} OUTPUTS/'.format(lig_location))
            os.system('mv {} OUTPUTS/'.format(out_location))

    except: # For docking failure cases 

        if os.path.exists(out_location):
            os.system('rm {}'.format(out_location))
        if os.path.exists(output_filename):
            os.system('rm {}'.format(output_filename))
        
        with open('./OUTPUT_{}.txt'.format(job_idx), 'a+') as f: 
            f.writelines(['{}, {}\n'.format(smi, 10000)]) 


def main(filename, job_idx):
    """
    Processes a batch of SMILES strings for molecular docking in parallel. By default 
    the function parallizes the docking calculations across the number of CPU available. 
    This function invokes perform_calc_single to perform the docking calculation. 

    Parameters:
        filename (str): Path to the file containing SMILES strings.
        job_idx (int): The index of the job, used for output differentiation.

    Example:
        >>> main('smiles.txt', 1)
        Begins processing of SMILES strings and outputs results to designated files.
    """
    smiles_all     = []
    with open(filename, 'r') as f:
        smiles_all = f.readlines()
    
    print('Num smiles:', len(smiles_all))
    data = [(smiles, job_idx) for smiles in enumerate(smiles_all)]
    
    # pool object with number of element
    pool = multiprocessing.Pool()
    
    # Parallel time: 
    start_time = time.time()
    pool.map(perform_calc_single, data)
    total_time = time.time() - start_time
    print('Total Time: {}, Jobs Idx: {}'.format(total_time, job_idx))
    
        
    
    
parser = argparse.ArgumentParser()
parser.add_argument("job_id", help="Array ID: SLURM_ARRAY_TASK_ID")
args = parser.parse_args()
job_idx = int(args.job_id)    

SMILES_FILE = './DATA/partition_{}.smi'.forma(job_idx)
print('Operating on partition file: ', SMILES_FILE)

main(SMILES_FILE, job_idx)

    
    
    
