################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2024 Naisarg Patel (https://github.com/naisarg14)        #
#                                                                              #
#          Project: https://github.com/naisarg14/mutadock                      #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
################################################################################


def backup(file_path):
    import os
    from datetime import datetime

    if not os.path.exists(file_path):
        return False
    master_folder, file = os.path.split(os.path.abspath(file_path))
    target_directory = os.path.join(master_folder, 'backups')
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    modified_time = os.path.getmtime(file_path)
    timestamp = datetime.fromtimestamp(modified_time).strftime("%b-%d-%Y_%H.%M")
    name, ext = os.path.splitext(file)
    target_file = os.path.join(target_directory, f'{name}_{timestamp}{ext}')
    os.rename(file_path, target_file)
    return True

def read_pdb_file(file_path):
    import re
    try:
        atoms = []
        pattern = r'^ATOM\s+(\d+)\s+([A-Z]+)\s+([A-Z]{2,3})\s+([A-Z]?)\s*(\d+)\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(\d+\.\d{1,2})\s+(\d+\.\d{1,3})(?:\s+(\d+\.\d{1,3}))?\s+[A-Z]$'
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith(('ATOM', 'HETATM')):
                    match = re.match(pattern, line)
                if match:
                    parts = {
                        "atom_serial_number": int(match.group(1)),  # Atom Serial Number
                        "atom_name": match.group(2).strip(),        # Atom Name
                        "residue_name": match.group(3).strip(),     # Residue Name
                        "chain_id": match.group(4).strip() or None, # Chain Identifier (optional)
                        "residue_sequence_number": int(match.group(5)),  # Residue Sequence Number
                        "x": float(match.group(6)),                 # X Coordinate
                        "y": float(match.group(7)),                 # Y Coordinate
                        "z": float(match.group(8)),                 # Z Coordinate
                        "occupancy": float(match.group(9)),         # Occupancy
                        "temp_factor": float(match.group(10)),      # Temperature Factor
                        "extra_factor": float(match.group(11)) if match.group(11) else None,  # Extra Factor (optional)
                    }
                    atoms.append(parts)
        return (atoms)
    except Exception as e:
        return (False, e)

def calculate_geometric_center(pdb_file):
    atoms = read_pdb_file(pdb_file)
    num_atoms = len(atoms)
    x_sum = sum(atom['x'] for atom in atoms)
    y_sum = sum(atom['y'] for atom in atoms)
    z_sum = sum(atom['z'] for atom in atoms)
    
    return (x_sum / num_atoms, y_sum / num_atoms, z_sum / num_atoms)

def calculate_radius(pdb_file):
    import math
    atoms = read_pdb_file(pdb_file)
    if isinstance(atoms, tuple):
        return atoms

    center = calculate_geometric_center(pdb_file)
    max_distance = 0

    for atom in atoms:
        distance = math.sqrt(
            (atom['x'] - center[0]) ** 2 +
            (atom['y'] - center[1]) ** 2 +
            (atom['z'] - center[2]) ** 2
        )
        if distance > max_distance:
            max_distance = distance

    return max_distance

def vina_split(input_file, output_file=None):
    from meeko import PDBQTMolecule, RDKitMolCreate

    if output_file is None:
        output_file = input_file.replace('.pdbqt', '_ligand_1.sdf')

    pdbqt_string = ""
    with open(input_file, 'r') as infile:
        for line in infile:
            if "vina result" in line.lower():
                score = [float(x) for x in line.split() if x.replace('.', '', 1).replace('-', '', 1).isdigit()][0]
            pdbqt_string += line
            if line.startswith('ENDMDL'):
                break
    molecule = PDBQTMolecule(pdbqt_string)
    sdf_string, failures = RDKitMolCreate.write_sd_string(molecule)

    if len(failures) > 0:
        msg = "\nCould not convert to RDKit. Maybe this library was not used for preparing\n"
        msg += "the input PDBQT for docking, and the SMILES string is missing?\n"
        msg += "Except for standard protein sidechains, all ligands and flexible residues\n"
        msg += "require a REMARK SMILES line in the PDBQT, which is added automatically by meeko."
        raise RuntimeError(msg)

    footer_string = f"> <Docking Score>\n{score}\n> <Credits>\nCreated using a script written by Naisarg Patel (Github:@naisarg14) for VIT's iGEM team.\n$$$$\n"
    with open(output_file, 'w') as outfile:
        outfile.write(sdf_string.replace('$$$$', footer_string))
    
    return (score, output_file)

def prepare_ligand(in_file, out_file=None):
    try:
        import sys
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing ligand files for Docking.\n"
        msg += "Easaies way to fix this is to install meeko and rdkit using the following command:\n\n"
        msg += "python -m pip install meeko rdkit\n"
        msg += "If you already have meeko and rdkit installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        print(msg)
        sys.exit(2)

    try:
        #Add support for PDB files
        if out_file is None:
            if in_file.endswith(".sdf"):
                out_file = f"{in_file.removesuffix(".sdf")}.pdbqt"
            elif in_file.endswith(".mol2"):
                out_file = f"{in_file.removesuffix(".mol2")}.pdbqt"
            else:
                return (False, "Input file is not in SDF or MOL2 format.")
        
        if in_file.endswith(".sdf"):
            mol = Chem.SDMolSupplier(in_file)[0]
        if in_file.endswith(".mol2"):
            mol = Chem.MolFromMol2File(in_file)
        
        mol = Chem.AddHs(mol)
        mp = MoleculePreparation()
        molecule_setups = mp.prepare(mol)
        pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(molecule_setups[0])

        if not success:
            raise RuntimeError(f"Could not convert to PDBQT: {error_msg}")
        with open(out_file, "w") as output_file:
            output_file.write(pdbqt_string)
    except Exception as e:
        return (False, e)

    return (True, pdbqt_string)

def prepare_receptor(receptor_filename, outputfilename="None"):
    #python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
    try:
        import sys, os
        from MolKit import Read
        from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing receptor files for Docking.\n"
        msg += "Easaies way to fix this is to install AutoDockTools_py3 using the following command:\n\n"
        msg += "python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3\n"
        msg += "If you already have AutoDockTools_py3 installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        print(msg)
        sys.exit(2)
    finally:
        original_stdout = os.dup(1)
        original_stderr = os.dup(2)
    with open(os.devnull, 'w') as fnull:
        os.dup2(fnull.fileno(), 1)
        os.dup2(fnull.fileno(), 2)

        if outputfilename is None:
            outputfilename = f"{receptor_filename}qt"
        
        # initialize required parameters
        repairs = 'hydrogens'
        charges_to_add = 'gasteiger'
        cleanup  = "waters"

        mode = 'automatic'
        delete_single_nonstd_residues = None
        dictionary = None
        unique_atom_names = False

        try:
            mols = Read(receptor_filename)
            mol = mols[0]
            if unique_atom_names:
                for at in mol.allAtoms:
                    if mol.allAtoms.get(at.name) >1:
                        at.name = at.name + str(at._uniqIndex +1)

            if len(mols)>1:
                #use the molecule with the most atoms
                ctr = 1
                for m in mols[1:]:
                    ctr += 1
                    if len(m.allAtoms)>len(mol.allAtoms):
                        mol = m

            mol.buildBondsByDistance()
            alt_loc_ats = mol.allAtoms.get(lambda x: "@" in x.name)
            len_alt_loc_ats = len(alt_loc_ats)
            if len_alt_loc_ats:
                print("WARNING!", mol.name, "has",len_alt_loc_ats, ' alternate location atoms!\nUse prepare_pdb_split_alt_confs.py to create pdb files containing a single conformation.\n')

            RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                                cleanup, outputfilename=outputfilename,
                                delete_single_nonstd_residues=delete_single_nonstd_residues,
                                dict=dictionary)
        except Exception as e:
            os.dup2(original_stdout, 1)
            os.dup2(original_stderr, 2)
            os.close(original_stdout)
            os.close(original_stderr)
            return (False, e)
        finally:
            os.dup2(original_stdout, 1)
            os.dup2(original_stderr, 2)
            os.close(original_stdout)
            os.close(original_stderr)
     
    return (True, "")

def add_score_to_csv(out_pdb, csv_file, score):
    import csv, os
    try:
        with open(csv_file, "r") as lc:
            final_line = lc.readlines()[-1]
            count = int(final_line.split(",")[0]) + 1
        name = f"{os.path.basename(out_pdb).removesuffix('_out.pdb')}"
        with open(csv_file, "a+") as out:
            writer = csv.DictWriter(out, fieldnames=["sr", "name", "affinity"])
            writer.writerow({"sr": count, "name": name, "affinity": score})
        return (True, name)
    except FileNotFoundError:
        count = 1
    except ValueError:
        count = 1
    except Exception as e:
        return (False, e)
    
def read_config(config):
    try:
        config = {}
        with open(config, 'r') as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith('#') and "=" in line:
                    key, value = line.split('=')
                    config[key.strip()] = value.strip()

        center_x = config.get('center_x', '0.0')
        center_y = config.get('center_y', '0.0')
        center_z = config.get('center_z', '0.0')
        size_x = config.get('size_x', '30.0')
        size_y = config.get('size_y', '30.0')
        size_z = config.get('size_z', '30.0')
        exhaustiveness = config.get('exhaustiveness', '32')
        n_poses = config.get('n_poses', '20')
        n_poses_write = config.get('n_poses_write', '5')
        overwrite = config.get('overwrite', 'True')

        center = (center_x, center_y, center_z)
        box_size = (size_x, size_y, size_z)

        return (True, center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite)
            
    except Exception as e:
        return (False, e)

def dock_vina(receptor, ligand, output, log_file, config=None, autosite=None, center=[0, 0, 0], box_size=[30, 30, 30], exhaustiveness=32, n_poses=20, n_poses_write=5, overwrite=True):
    if config is not None:
        values = read_config(config)
        if values[0] is False:
            return values
        _, center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite = values

    if autosite is not None:
        center = calculate_geometric_center(autosite)
    
    import subprocess
    commands = [
        "python3", "vina_dock.py", 
        "--receptor", receptor, 
        "--ligand", ligand, 
        "--output", output, 
        "--center", str(center[0]), str(center[1]), str(center[2]), 
        "--box_size", str(box_size[0]), str(box_size[1]), str(box_size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--n_poses", str(n_poses),
        "--n_poses_write", str(n_poses_write),
        ]
    if not overwrite: commands.append("--nooverwrite")
    with open(log_file, 'w') as log_file:
        result = subprocess.run(commands, stdout=log_file, stderr=log_file, text=True) 

    if result.returncode != 0:
        return (False, f"Check the error in {log_file}")

    return (True, "")


def main():
    print("This is a dependency file for mutadock library's docking module.")



if __name__ == "__main__":
    main()





#prepare_ligand("antheraxanthin.sdf")
#prepare_receptor("BCH.pdbqt")
#vina_dock("BCH.pdbqt", "antheraxanthin.pdbqt", "log.txt",exhaustiveness=4)
#vina_split("dock_try.pdbqt")
