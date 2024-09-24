# MUTADOCK

## Introduction
MUTADOCK is a comprehensive library designed for mutation studies and multiple receptor-ligand docking. It provides tools and methods to analyze and predict the effects of mutations on receptor-ligand interactions, enabling researchers to study protein function and drug binding affinity in a detailed manner.

### Description
Our software is designed to facilitate protein mutation analysis and molecular docking. It integrates automated protein mutation using PyRosetta and a docking library capable of docking multiple proteins with multiple ligands. 

### Key Features
#### Automated Protein Mutation:
- Utilizes PyRosetta for systematic protein mutations.
- Supports various mutation strategies (e.g., single-point mutations, double-point mutations and triple-point mutations).
- Allows customization of mutation and docking parameters.
#### Docking Library:
- Capable of docking a list of proteins against a list of ligands.
- Employs AutoDock Vina to predict binding affinities and best poses.
- Provides detailed output files with docking scores and poses.
#### User Interface:
- Command-Line Interface: Simple CLI for both beginner and expert users.
- Python Bindings: The library can be imported in other codes for increased customizability by expert users


## Try it Now
The basic codes to perform mutation studies as used by us for our project can be found in a Jupyter Notebook here. The same notebook can be found on collab here.

## How To Guide

### Installation
MutaDock has been deployed on PyPi, making installation quick and simple
```pip install mutadock```

The Pyrosetta Installer will be automatically installed but Pyrosetta should be installed using
```python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'```

- Currently there is a problem with the vina on PyPi, so vina needs to be installed separately, the installation guide can be found at https://autodock-vina.readthedocs.io/en/latest/installation.html


### Mutation Studies
Mutation Studies for a protein is a very fast process with just a PDB file of the protein as the input. (We assume for the tutorial that the name of the PDB file is “protein.pdb”)

```md_mutate -i protein.pdb```

Other optional arguments can be changed as required, to check the usage run 
```md_mutate -h```

The md_mutate will output CSV files and one text file, their description is in the table below:
| No. | File Name                           | Description                                                                                       |
|-----|-------------------------------------|---------------------------------------------------------------------------------------------------|
| 1.  | protein_modified_mutations_all.csv  | Contains all possible mutations for the given protein                                             |
| 2.  | protein_modified_mutations.csv      | Contains mutations that are possible according to the PAM matrix for the given protein            |
| 3.  | protein_modified_mutations_ddG.csv  | The single mutation ddG values for the mutation in the File-2                                      |
| 4.  | protein_modified_mutations_ddG_sorted.csv | Sorted File-3 from lowest to highest ddG values                                               |
| 5.  | protein_modified_double_ddg.csv     | The ddG values of the double mutation for all the combinations of the most negative single ddG compounds |
| 6.  | protein_modified_double_ddg_sorted.csv | Sorted File-5 from lowest to highest ddG values                                               |
| 7.  | protein_modified_triple_ddg.csv     | The ddG values of the triple mutation for all the combinations of the most negative double ddG compounds |
| 8.  | protein_modified_triple_ddg_sorted.csv | Sorted File-7 from lowest to highest ddG values                                               |
| 9.  | protein_modified_mutants.txt        | Generates a list of all the mutated PDB files created. Can be directly used as input for the md_dock command in our mutadock library |

