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
- **Command-Line Interface:** Simple CLI for both beginner and expert users.
- **Python Bindings:** The library can be imported in other codes for increased customizability by expert users


## Try it Now
The basic codes to perform mutation studies as used by us for our project can be found in a Jupyter Notebook here. The same notebook can be found on collab here.

## How To Guide

### Installation
MutaDock has been deployed on PyPi, making installation quick and simple
```
pip install mutadock
```

The Pyrosetta Installer will be automatically installed but Pyrosetta should be installed using
```
python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```

- Currently there is a problem with the vina on PyPi, so vina needs to be installed separately, the installation guide can be found at https://autodock-vina.readthedocs.io/en/latest/installation.html


### Mutation Studies
Mutation Studies for a protein is a very fast process with just a PDB file of the protein as the input. (We assume for the tutorial that the name of the PDB file is “protein.pdb”)

```
md_mutate -i protein.pdb
```

Other optional arguments can be changed as required, to check the usage run 
```
md_mutate -h
```

The md_mutate will output CSV files and one text file, their description is in the table below:
| **No.** | **File Name**                           | **Description**                                                                                       |
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


### Docking Studies
Docking for multiple receptors and ligands is made simple and efficient by mutadock. The text files containing the names of the receptors and ligands need to be given as input, after that everything is automated. (If md_mutate is used, the text file for receptor is generated automatically)
Every receptor in the receptor file will be docked with every ligand in the ligand file. A standard Vina configuration file or an AutoSIte prediction output is required.
Example:
```
md_dock -r receptors.txt -l ligands.txt -c config.txt
```
Other optional arguments can be changed as required, to check the usage run 
```
md_dock -h
```
The output of md_dock with their description is in the table below:
| **No.** | **Output**              | **Description**                                                                                       |
|-----|---------------------|---------------------------------------------------------------------------------------------------|
| 1.  | PDBQT files         | The receptors and ligands will be converted to PDBQT files for AutoDock Vina.                     |
| 2.  | Output Log          | The output of AutoDock Vina with the docking scores will be stored in a log file for each combination. |
| 3.  | Output PDB          | The output of AutoDock Vina with the 5 best docking poses will be stored in a PDB file for each combination. |
| 4.  | Output PDBQT        | The output of AutoDock Vina Split with the best pose will be stored in a PDBQT file for each combination. |
| 5.  | Output SDF          | The best pose after docking will be stored in a SDF file for visualization and better usability.  |
| 6.  | Docking Results CSV | All the docking affinities are tabulated in a CSV to make analysis trivial.                       |


### All CLI Scripts
| **No.** | **Command**            | **Description**                                                                           |
|-----|--------------------|---------------------------------------------------------------------------------------|
| 1.  | md_mutate          | Predicts the best mutation of the given protein                                       |
| 2.  | md_dock            | Docked all combinations from a list of receptors and ligands                          |
| 3.  | md_vina_dock       | CLI for AutoDock Vina                                                                 |
| 4.  | md_csv_generator   | Generates all possible mutations for a protein and also the mutations possible according to PAM Matrix |
| 5.  | md_csv_sort        | Can sort any CSV file according to the column name or number chosen                   |
| 6.  | md_ddg_single      | Calculates single ddG values for a given CSV of mutations                             |
| 7.  | md_ddg_double      | Calculates double ddG values for all combinations using a given CSV of mutations      |
| 8.  | md_ddg_triple      | Calculates triple ddG values for all combinations using a given CSV of mutations      |


## Applications
- **Protein Engineering:** Designing mutated proteins with enhanced stability or new functionalities.
- **Drug Discovery:** Screening potential drug candidates by predicting binding affinities.
- **Biochemical Research:** Studying protein-ligand interactions to understand biological processes.


## Documentation
- README is included in the repository to serve as a comprehensive guide
- ReadtheDocs Page for updated documentation can be found [here](https://mutadock.readthedocs.io/en/latest/#)


## Future Developments
- **Developing a Graphical User Interface (GUI):** Enhancing user experience by providing a user-friendly interface for easier interaction with the software.
- **Creating a Web Server:** Allowing remote access and usage of the software through a web-based platform, making it accessible from anywhere.
- **Increasing Parameter Customizability:** Offering more options for users to fine-tune mutation and docking parameters to suit specific research needs and conditions.


## Acknowledgements
- Open-source tools and libraries used in the development.


## Contact
- For questions, suggestions, or collaboration, please contact [Naisarg Patel](mailto:naisarg.patel14@hotmail.com).
