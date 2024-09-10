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

from itertools import product
import os, time, sys
import argparse
from contextlib import contextmanager
from tqdm import tqdm
from vina_helper import prepare_receptor, prepare_ligand, vina_split, add_score_to_csv, dock_vina, read_config, calculate_geometric_center, backup


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def naisarg():
    ignore_existing = False

    start_time = time.time()
    receptors, ligands, config, autosite, quiet = prepare_inputs()

    if config is None and autosite is None:
        print("Both config and autosite not provided. Assuming center as [0,0,0] and box_size as [30,30,30].")

    if config is not None:
        values = read_config
        values = read_config(config)
        if values[0] is False:
            return values
        _, center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite = values

    if autosite is not None:
        site = calculate_geometric_center(autosite)
        if site[0] is not False:
            center = site[0]
        else:
            sys.exit(f"Error while calculating the geometric center of the autosite file {autosite} \nError: {site[1]}")


    completed_name = f"{receptor.removesuffix(".txt")}_{ligand.removesuffix(".txt")}_completed.txt"

    with open(completed_name, "r") as f:
        completed = f.readlines().split("\n")

    combinations = list(product(receptors, ligands))

    combinations = [x for x in combinations if x not in completed]

    for combination in tqdm(combinations):
        receptor = combination[0]
        ligand = combination[1]
        prepared_receptor = f"{receptor}qt"
        
        if ligand.endswith(".sdf"):
            prepared_ligand = f"{ligand.removesuffix(".sdf")}.pdbqt"
        elif ligand.endswith(".mol2"):
            prepared_ligand = f"{ligand.removesuffix(".mol2")}.pdbqt"

        if not quiet: print(f"Docking the receptor {receptor} to the ligand {ligand}")
        master_folder_receptor, receptor_file = os.path.split(os.path.abspath(receptor))
        master_folder_ligand, ligand_file = os.path.split(os.path.abspath(ligand))
        
        output_dir = os.path.join(master_folder_receptor, 'out')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        out_pdb = os.path.join(output_dir, f'{receptor_file.removesuffix(".pdb")}_{ligand_file.removesuffix(".sdf")}_out.pdb')
        log_file = os.path.join(output_dir, f'{receptor_file.removesuffix(".pdb")}_{ligand_file.removesuffix(".sdf")}_log.txt')
        csv_file = os.path.join(output_dir, "docking_results.csv")

        backup(out_pdb)
        backup(log_file)

        try:
            if not quiet: print(f"Docking for {ligand} with {receptor}")
            if not quiet: print("Press Ctrl+D (EOFE Error) to skip this receptor-ligand combination.")
            if not os.path.exists(prepared_receptor) or ignore_existing:
                if not quiet: print(f"Preparing receptor {receptor}")
                with suppress_stdout(): rec_out = prepare_receptor(receptor_filename=receptor, outputfilename=prepared_receptor)
                if not rec_out[0]:
                    print(f"Error while preparing receptor {receptor} \nError: {rec_out[1]} \nSkipping this receptor-ligand combination.")
                    continue

            if not os.path.exists(prepared_ligand) or ignore_existing:
                if not quiet: print(f"Preparing ligand {ligand}")
                with suppress_stdout(): lig_out = prepare_ligand(in_file=ligand, out_file=prepared_ligand)
                if not lig_out[0]:
                    print(f"Error while preparing ligand {ligand} \nError: {lig_out[1]} \nSkipping this receptor-ligand combination.")
                    continue

            if not quiet: print(f"Starting docking")
            with suppress_stdout(): vina_out = dock_vina(prepared_receptor, prepared_ligand, out_pdb, log_file, center=center, box_size=box_size, exhaustiveness=exhaustiveness, n_poses=n_poses, n_poses_write=n_poses_write, overwrite=overwrite)
            if not vina_out[0]:
                print(f"Error while docking {ligand} to {receptor} \nError: {vina_out[1]} \nSkipping this receptor-ligand combination.")
                continue
            if not quiet: print("Docking Completed, writing log file")

            if not quiet: print("Getting the ligand 1 after docking")
            ligand_1 = out_pdb.replace('.pdbqt', '_ligand_1.sdf')
            with suppress_stdout(): score, _ = vina_split(input_file=out_pdb, output_file=out_pdb)
            if not quiet: print("Adding affinity to CSV")
            _ = add_score_to_csv(out_pdb, csv_file, score)
            
            with open(completed_name, "a+") as file: file.write(f"{combination}\n")

            if not quiet: print(f"Docking completed, log file is {log_file}, ligand_1 is {ligand_1.replace("pdbqt", "sdf")}, docking affinity is {score}. \n")

        except EOFError:
            continue

    end_time = time.time()
    elapsed_time = (end_time - start_time)/60
    print(f"Completed in {elapsed_time:.2f} minutes!")


def prepare_inputs():
    #ADD AUTOSITE
    parser = argparse.ArgumentParser(prog="np_dock", description=None, epilog="Part of mutadock library. Written by Naisarg Patel (https://github.com/naisarg14)")
    parser.add_argument('receptor_txt', help="Text File with all receptors", metavar="RECEPTOR")
    parser.add_argument("ligand_txt", help="Text File with all ligands", metavar="LIGAND")
    parser.add_argument("--config",default=None, help="Text File with all Vina Configuration Settings", metavar="CONFIG")
    parser.add_argument("--autosite",default=None, help="PDB generated by autosite for binding site of protein.", metavar="AUTOSITE")
    parser.add_argument("--quiet", action="store_true", help="Run the Docking in quiet mode (default: False).")

    args = parser.parse_args()

    receptor_txt = args.receptor_txt
    ligand_txt = args.ligand_txt
    config = args.config

    print("Receptor file:", receptor_txt)
    print("Ligands file:", ligand_txt)
    print("Config file:", config)
    for f in [receptor_txt, ligand_txt]:
        if not f.endswith(".txt"):
            f = f"{f}.txt"
        if not os.path.exists(f):
            print(f"{f} not found or cannot be opened.")
    try:
        with open(receptor_txt, "r") as rec:
            receptors = rec.readlines()
    except IOError as err:
        sys.exit(f"Error reading the file {receptor_txt}: ".format(receptor_txt, err))
    try:
        with open(ligand_txt, "r") as lig:
            ligands = lig.readlines()
    except IOError as err:
        sys.exit(f"Error reading the file {ligand_txt}: ".format(ligand_txt, err))

    for i in range(len(receptors)):
        receptors[i] = receptors[i].replace("\n", "")
        if not os.path.isabs(receptors[i]):
            receptors[i] = os.path.join(os.getcwd(), receptors[i])
    for i in range(len(ligands)):
        ligands[i] = ligands[i].replace("\n", "")
        if not os.path.isabs(ligands[i]):
            ligands[i] = os.path.join(os.getcwd(), ligands[i])


    return receptors, ligands, config, args.autosite, args.quiet


if __name__ == "__main__":
    naisarg()
