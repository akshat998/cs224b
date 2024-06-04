# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 04:43:56 2024

@author: aksha
"""
import os 
import sys
import argparse
import subprocess
import concurrent.futures

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

def update_config_file(filename, num_mols, smiles_file_path):
    """
    Updates the configuration file with new numerical and filepath settings for molecule processing.

    This function modifies specific entries in a given configuration file, specifically the number of molecules 
    (`NUM_MOLS`) and the file path for SMILES files (`SMILES_FILES`). It reads the existing configuration file line by line, 
    replaces lines containing these entries with updated values, and writes the modified lines back to the file.

    Parameters:
        filename (str): The path to the configuration file to be updated.
        num_mols (int): The new total number of molecules to be processed, updating the 'NUM_MOLS' entry.
        smiles_file_path (str): The new file path for the SMILES file, updating the 'SMILES_FILES' entry.

    Example Usage:
        Suppose 'all.ctrl' originally contains:
            NUM_MOLS=105338
            SMILES_FILES=./DATA/docking.smi

        After calling update_config_file('all.ctrl', 50000, './DATA/missing_smiles.smi'), the file 'all.ctrl' will be updated to:
            NUM_MOLS=50000
            SMILES_FILES=./DATA/missing_smiles.smi

    This function assumes the configuration file is well-formed with entries in the format 'key=value'. Lines that do not match 
    the entries being updated are written back to the file unchanged.
    """
    updated_lines = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "NUM_MOLS=" in line:
                updated_lines.append(f"NUM_MOLS={num_mols}\n")
            elif "SMILES_FILES=" in line:
                updated_lines.append(f"SMILES_FILES={smiles_file_path}\n")
            else:
                updated_lines.append(line)
    with open(filename, "w") as f:
        f.writelines(updated_lines)
    print(f"Configuration file {filename} updated.")
    
def delete_file(file_path):
    """Helper function to delete a single file. 
    
    Used for cleanup of intermediate files. 
    
    """
    os.remove(file_path)
    
def create_and_submit_job(script_name, job_part):
    """
    Creates a Slurm batch job script tailored for a specific job part and submits it for execution.

    This function dynamically generates a batch job script with custom configurations per job part,
    including dedicated output and error logs for each run. The script is written to the file system,
    submitted to Slurm with the `sbatch` command, and then deleted to clean up and prevent clutter.

    Parameters:
        script_name (str): The filename for the Slurm batch job script. This script is temporary
                           and will be deleted after submission.
        job_part (int): The specific job part or task ID to be processed. This ID is used to
                        specify which part of the dataset the `dataset_calc.py` script should
                        work on, and it also influences the naming of the job and log files.

    The generated script sets up the required environment by loading necessary modules and
    specifies resource allocations (e.g., number of tasks per node, memory per node, etc.).
    It configures job output settings such as standard output and error files which include the
    job part in their names for easy identification and troubleshooting.

    The script directly executes a Python script (`dataset_calc.py`) with the job part as an argument,
    which is intended to perform computations or data processing specific to that part of the overall task.

    Examples:
        - create_and_submit_job("resubmit_job123_part1.sh", 1)
        This would create a job script named "resubmit_job123_part1.sh" that, when executed by Slurm,
        runs `dataset_calc.py 1`, targeting the dataset segment partition_1.smi

    After submission, the function ensures that the script file is immediately deleted to maintain
    a clean workspace and avoid the accumulation of temporary files.
    """
    with open(script_name, 'w') as script_file:
        script_file.write(
            f"""#!/bin/bash
            #SBATCH --account=akshat998
            #SBATCH --ntasks-per-node=40
            #SBATCH --mem=7000M               # memory (per node)
            #SBATCH --time=12:0:00
            #SBATCH --job-name nodeTask
            #SBATCH -e stderr_{job_part}.txt
            #SBATCH -o stdout_{job_part}.txt
            #SBATCH --open-mode=append
            #SBATCH --export=NONE
            
            module --force purge
            module load nixpkgs/16.09
            module load gcc/7.3.0
            module load rdkit/2019.03.4
            module load scipy-stack/2019b
            module load openbabel
            
            python3 dataset_calc.py {job_part}
            """)
    subprocess.run(['sbatch', script_name])
    os.remove(script_name)
    print(f"Job part {job_part} resubmitted and script {script_name} deleted.")
    

# Read the all.ctrl file to get parameters
config_params = read_config_file("all.ctrl")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Monitor and resubmit script")
    parser.add_argument(
        "monitor_type", 
        choices=["check_progress", "finish_and_resubmit"],
        help="Type of monitoring to perform. 'check_progress' for checking the current progress of submitted jobs, "
             "'finish_and_resubmit' to collect outputs, clean up, and prepare and resubmit jobs for any missing molecules."
    )
    parser.add_argument(
        "job_id",
        type=int,
        help="The ID of the job or batch to monitor or manage."
    )
    args = parser.parse_args()
    
    
    if args.monitor_type == "check_progress":
        # Check job progress logic
        job_info = subprocess.run(['squeue', '-j', str(args.job_id)], capture_output=True, text=True, check=True)
        job_info = job_info.stdout.split('\n')
        if len(job_info) <= 1:  # No jobs pending/running
            print('All jobs have completed, progress is completely done. Please consider running finish_and_resubmit')
            sys.exit()
        else:
            print('Proceeding with job analysis')

        completed_job_ids = {line.split()[0].split('_')[1] for line in job_info if line.strip()}
        crashed_job_partitions = [i for i in range(1, config_params['MAX_NUM_JOBS'] + 1) if str(i) not in completed_job_ids]

        if not crashed_job_partitions:
            print('No crashed jobs detected. Please check again later. Exiting.')
            sys.exit()

        print('Detected Crashed Jobs:', [f"{args.job_id}_{idx}" for idx in crashed_job_partitions])
        for idx in crashed_job_partitions:
            script_name = f"resubmit_{args.job_id}_{idx}.sh"
            create_and_submit_job(script_name, idx)
                
    elif args.monitor_type == "finish_and_resubmit":
        # Check for currently running jobs using the job ID.
        job_info = subprocess.run(['squeue', '-j', str(args.job_id)], capture_output=True, text=True, check=True)
        job_info = job_info.stdout.split('\n')
        if len(job_info) > 1:  # If there's more than one line, jobs are still running
            print("Not all jobs have finished, skipping resubmission steps. Please consider using 'check_progress'.")
            sys.exit() 

        print("All jobs have finished, proceeding with analysis.")

        # Combine output files into a single file for analysis and remove the originals.
        output_files = [file for file in os.listdir(".") if file.startswith("OUTPUT_") and file.endswith(".txt")]
        if output_files:
            with open("DATA/combined_output.txt", "w") as combined_file:
                for file in output_files:
                    with open(file, "r") as f:
                        combined_file.write(f.read())
                    os.remove(file)
            print("Output files combined and original files deleted.")

        # Delete all intermediate .pdbqt files in parallel to clean up the workspace.
        pdbqt_files = [os.path.join(".", file) for file in os.listdir(".") if file.endswith(".pdbqt")]
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(delete_file, pdbqt_files)
        print("Cleanup completed. Intermediate PDBQT files deleted.")

        # Similarly, delete all SMI partition files in the ./DATA directory in parallel.
        smi_files = [os.path.join("./DATA", file) for file in os.listdir("./DATA") if file.startswith("partition_") and file.endswith(".smi")]
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(delete_file, smi_files)
        print("SMI partition files deleted.")

        # Identify incomplete molecules from the combined output.
        with open("DATA/combined_output.txt", 'r') as f:
            all_molecules = f.readlines()
        with open(config_params['SMILES_FILES'], 'r') as f:
            completed_mols = f.readlines()
        completed_mols = {molecule.split(',')[0].strip(): i for i, molecule in enumerate(completed_mols)}

        incomplete_calculations = [item+'\n' for item in all_molecules if item not in completed_mols]
        NUM_MISSING_MOLS = len(incomplete_calculations)
        print(f'Total number of incomplete molecules: {NUM_MISSING_MOLS}')

        # Write missing molecules to a new SMILES file for further processing.
        with open('./DATA/missing_smiles.smi', 'w') as f:
            f.writelines(incomplete_calculations)

        # Update the control file to reflect the new number of molecules and their data file.
        update_config_file("all.ctrl", NUM_MISSING_MOLS, './DATA/missing_smiles.smi')

        # Run load balancing on the new collection of molecules, then submit the jobs.
        os.system('python3 load_balancer.py')
        os.system('sbatch submit.sh')    
        
        
 