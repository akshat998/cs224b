# Distributed System for Enhanced Drug Discovery Through Molecular Docking
This application establishes a distributed system for drug discovery by utilizing load balancing and fault tolerance to identify promising molecules that interact with a specified protein of interest.

The core workflow allows users to input a SMILES text file, the receptor of interest, and the docking parameters to facilitate large-scale docking simulations. Out-of-the-box, the system employs QuickVina 2.0 for docking.

Execution within SLURM environments is highly optimized, with computations distributed in parallel across multiple CPUs and nodes. This architecture ensures efficient linear scaling relative to the number of molecules processed. Moreover, the application supports robust load balancing and enhanced fault tolerance mechanisms.
## Prerequisites
Please clone the repository using: 
```
git clone https://github.com/akshat998/cs224b.git
```
Please ensure that the following packages are installed: 
- [RDKit version 2021.09.5](https://www.rdkit.org/docs/Install.html)
- [Open Babel 3.1.0](https://openbabel.org/docs/dev/Installation/install.html)
- [Python 3.7.13](https://www.python.org/downloads/) (or higher)

## File Navigator
* `DATA`: This directory is where users can place the receptor file and the corresponding executables for running the docking process.
* `OUTPUTS`: This directory is designated for storing the results of the docking simulations (created on the fly).
* `all.ctrl`: Contains all user-specifiable parameters required for the screening process, including the docking parameterization.
* `dataset_calc.py`: A Python script for running the docking on specified ligands.
* `submit.sh`:  A Slurm submission script for submitting an array of jobs for processing.
* `load_balancer.py`:  Script that evenly distributes the workload across different nodes based on the number of molecules in a user-provided file and the number of nodes specified by the user.
* `monitor_and_resubmit.py`:  A script that can monitor ongoing jobs and resubmit tasks for nodes that have unexpectedly crashed.

## Quick Start Guide

To get started with the docking simulations, follow the steps outlined below. These steps ensure that your configuration is correctly set up for your specific docking scenario:

1. **Configure Receptor Location:**
   - Open `all.ctrl` and specify the exact location of your receptor in the designated section.

2. **Set Docking Parameters:**
   - Within `all.ctrl`, enter the appropriate `CENTER-X/Y/Z` and `SIZE-X/Y/Z` coordinates to define your docking area.

3. **Specify SMILES List Path:**
   - In `all.ctrl`, input the path to your SMILES list file. This file is crucial for defining the molecular inputs for the simulation.
   - Ensure the file adheres to the format: each line contains a SMILES followed by a newline character (e.g., `C[C@@H](N)C(=O)O \n`).

4. **Slurm Cluster Account:**
   - In `submit.sh`, replace `TODO` in `#SBATCH --account=TODO` with your actual Slurm cluster account name to ensure proper job submission.

5. **Job Submission Configuration:**
   - Adjust the number of jobs to submit for your docking calculation in `submit.sh` by modifying `#SBATCH --array=1-999` accordingly. Ensure this number matches the `MAX_NUM_JOBS` parameter set in `all.ctrl`.

6. **Executable Permissions:**
   - Make sure the docking executables have the correct executable permissions by running `chmod 777 ./DATA/qvina`.
   
7. **Assign Subtasks to Each Node with Load Balancing:**
   - Within `all.ctrl`, please ensure that the appropriate flag is set for the parameter `USE_LOAD_BALANCER`. If set to False, molecules are distributed randomly across the nodes. If set to True, molecules are divided based on the number of atoms, which leads to a more balanced compute load per node, enhancing efficiency.
   - Following this setup, run `python3 load_balancer.py`. This will generate various `partition_i.smi` files in the DATA directory, where each i-th node will process all the molecules listed in its respective file.

8. **Submit Your Job:**
   - Finally, submit your job to the Slurm cluster with the command: `sbatch submit.sh`.

By following these steps, you'll be properly set up to conduct your docking simulations. Ensure all paths and parameters are double-checked for accuracy before submitting your job.

## Analyzing Your Jobs

The `monitor_and_resubmit.py` script provides robust monitoring and management functionalities for the distributed molecular docking processes. This script helps in overseeing the progress of ongoing jobs and managing the resubmission of jobs for nodes that may have unexpectedly crashed, ensuring your simulations run efficiently and effectively.

### Features
- **Job Monitoring**: Check the progress of running jobs to ensure they are proceeding as expected.
- **Automatic Resubmission**: Automatically detects and resubmits failed or crashed jobs to maintain continuous operation without manual intervention.

### Usage
To utilize the `monitor_and_resubmit.py` script, execute it with specific arguments based on your monitoring and resubmission needs:

1. **Check Progress**:
   - Use this mode to check the current status of submitted jobs. If a crash is detected in any node, the affected node's job will be resubmitted. If all jobs have completed, you will be advised to consider using the finish-and-resubmit mode. If all jobs are in progress without any crashes, no action will be taken.
   - Command: `python3 monitor_and_resubmit.py check_progress [job_id]`

2. **Finish and Resubmit**:
   - This mode should be used once all jobs have finished. It performs output collection, cleanup, and preparation for resubmission of jobs for any missing molecules.
   - It consolidates output files, cleans up intermediate files, and prepares a new batch of submissions if there are incomplete calculations.
   - Command: `python3 monitor_and_resubmit.py finish_and_resubmit [job_id]`
   - The final output file, which contains the SMILES string and the completed docking score, will be stored in `DATA/combined_output.txt` in the format `SMILES,Docking Score`.

### Operational Details
- The script first checks job status using SLURM's `squeue` command to identify if jobs are still running or have crashed.
- For crashed jobs, the script generates a new SLURM script for each crashed job part, resubmits it, and ensures that the workspace remains clean by deleting temporary files.
- In the finish-and-resubmit mode, the script combines output files for analysis, deletes intermediate files to free up space, and identifies any incomplete molecule calculations. It then updates the configuration file to reflect the new number of molecules and their file paths, and prepares the system for a new round of submissions.


## Experiments


## Questions, problems?
Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (akshat98m[AT]stanford[DOT]edu)
