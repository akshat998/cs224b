# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 04:43:44 2024

@author: aksha
"""
import os
import random
from rdkit import Chem
from typing import List



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

def count_atoms(smiles: str) -> int:
    """
    Calculates and returns the number of atoms in a molecule based on its SMILES string.

    This function uses RDKit to convert a SMILES string into a molecule object and then
    returns the number of atoms in the molecule. If the SMILES string is invalid and cannot
    be converted into a molecule, it raises an exception.

    Parameters:
        smiles (str): The SMILES string representation of the molecule.

    Returns:
        int: The number of atoms in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid and the molecule cannot be created.

    Example:
        If smiles = "CCO", then count_atoms("CCO") will return 3.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol.GetNumAtoms()
    else:
        raise ValueError(f"Invalid SMILES string provided: '{smiles}'. Unable to create molecule.")

def load_balance_smiles(smiles_list: List[str], m: int) -> List[List[str]]:
    """
    Distributes SMILES strings into 'm' partitions to achieve a balanced computational load based on the number of atoms in each molecule.
    
    The function first calculates the number of atoms for each molecule represented by a SMILES string. It then sorts these molecules
    based on their atom counts in descending order. This sorted list ensures that molecules with the greatest computational demand 
    are considered first. Using a greedy algorithm, the function allocates each molecule to the partition with the current lowest 
    cumulative atomic count, thereby aiming for an even distribution of computational load across all partitions.
    
    Parameters:
        smiles_list (List[str]): A list of SMILES strings representing the molecules to be distributed.
        m (int): The number of partitions into which the molecules should be distributed.

    Returns:
        List[List[str]]: A list of 'm' sub-lists, where each sub-list contains the SMILES strings allocated to that partition.
    
    Example:
        smiles_data = ["CCO", "Oc1ccccc1O", "CCCC", "NN", "c1ccccc1"]
        partitions = load_balance_smiles(smiles_data, 3)
        for i, partition in enumerate(partitions):
            print(f"Partition {i+1}: {partition}")
    """
    # Calculate atoms for each molecule
    molecules = [(smiles, count_atoms(smiles)) for smiles in smiles_list]
    
    # Sort molecules by number of atoms in descending order
    molecules.sort(key=lambda x: x[1], reverse=True)
    
    # Initialize partitions and their load counts
    partitions = [[] for _ in range(m)]
    partition_loads = [0] * m
    
    # Distribute molecules using a greedy approach
    for smiles, atom_count in molecules:
        # Find the partition with the minimum load
        min_index = partition_loads.index(min(partition_loads))
        partitions[min_index].append(smiles)
        partition_loads[min_index] += atom_count
    
    return partitions


def random_load_balance_smiles(smiles_list: List[str], m: int) -> List[List[str]]:
    """
    Evenly distributes SMILES strings into 'm' partitions randomly. This function shuffles the list of SMILES strings
    and then evenly assigns them to the specified number of partitions, ensuring each partition approximately receives
    the same number of molecules.

    Parameters:
        smiles_list (List[str]): A list of SMILES strings representing the molecules to be distributed.
        m (int): The number of partitions into which the molecules should be evenly distributed.

    Returns:
        List[List[str]]: A list of 'm' sub-lists, where each sub-list contains the SMILES strings allocated to that partition.
        The distribution aims for an equal number of SMILES strings in each partition, subject to rounding when the total number
        of molecules isn't perfectly divisible by 'm'.

    Example:
        smiles_data = ["CCO", "Oc1ccccc1O", "CCCC", "NN", "c1ccccc1"]
        partitions = random_load_balance_smiles(smiles_data, 3)
        for i, partition in enumerate(partitions):
            print(f"Partition {i+1}: {partition}")
    """
    random.shuffle(smiles_list)
    partitions = [[] for _ in range(m)]    # Initialize partitions

    # Evenly distribute SMILES strings to each partition
    for index, smiles in enumerate(smiles_list):
        partition_index = index % m
        partitions[partition_index].append(smiles)
    
    return partitions


if __name__ == '__main__':
    # Read configuration parameters from the control file
    config_params = read_config_file("all.ctrl")

    # Extract configuration settings
    SMILES_FILE = str(config_params.get("SMILES_FILES"))
    MAX_NUM_JOBS = int(config_params.get("MAX_NUM_JOBS"))
    USE_LOAD_BALANCER = bool(config_params.get("USE_LOAD_BALANCER"))

    # Read the SMILES data from the file specified in the configuration
    with open(SMILES_FILE, 'r') as f:
        smiles_data = f.readlines()

    if USE_LOAD_BALANCER:
        # When load balancing is enabled, distribute SMILES strings to optimize computational load
        print('Partitioning based on the load balancer')
        partitions = load_balance_smiles(smiles_data, MAX_NUM_JOBS)
    else:
        # When load balancing is not enabled, distribute SMILES strings randomly
        print('Partitioning WITHOUT a load balancer')
        partitions = random_load_balance_smiles(smiles_data, MAX_NUM_JOBS)

    # Output and save the results for each partition
    for i, partition in enumerate(partitions):
        partition_filename = f'./DATA/partition_{i+1}.smi'
        print(f"Partition {i+1}: len = {len(partition)}")
        with open(partition_filename, 'w') as f:
            f.writelines(partition)