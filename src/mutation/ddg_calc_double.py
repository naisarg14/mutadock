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


import mutation.predict_ddG as predict_ddG
import csv
from itertools import combinations
from mutation.helpers import backup, Mutation
import os, sys, argparse

try:
    from pyrosetta import *
except ImportError:
    msg = "Error with importing pyrosetta module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install pyrosetta using the following command:\n\n"
    msg += "python -m pip install pyrosetta_installer && python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'\n"
    msg += "Alternative is to install from PyRosetta's official website.\n"
    msg += "If you already have pyrosetta installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    print(msg)
    sys.exit(2)

try:
    from tqdm import tqdm
except ImportError:
    msg = "Error with importing tqdm module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install tqdm using the following command:\n\n"
    msg += "python -m pip install tqdm\n"
    msg += "If you already have tqdm installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    print(msg)
    sys.exit(2)


def main():
    full_pdb_path, full_csv_path, total, out_csv = get_inputs()
    calc_double_ddg(pdb_file=full_pdb_path, single_csv=full_csv_path, total=total, out_csv=out_csv)


def calc_double_ddg(pdb_file, single_csv=None, out_csv=None, total=0):
    if not out_csv:
        out_csv = f"{pdb_file.removesuffix('.pdb')}_double_ddg.csv"
    backup(out_csv)

    title = get_mut_csv(single_csv, int(total))
    
    sys.stdout = open(os.devnull, 'w')
    init('-mute all')
    sys.stdout = sys.__stdout__
    
    with open(out_csv, "w+") as out:
        writer = csv.DictWriter(out, fieldnames=["sr", "pdb", "combination", "mut1_wtAA", "chain1", "mut1_position", "mut1_prAA", "mut2_wtAA", "chain2", "mut2_position", "mut2_prAA", "double_ddG_value"])
        writer.writeheader()
        count = 1
        pose = pose_from_pdb(pdb_file)
        sfxn = get_fa_scorefxn()
        score_1 = sfxn.score(pose)
        for x in tqdm(combinations(title, 2)):
            mut1 = x[0]
            mut2 = x[1]

            mut1 = Mutation(x[0]["chain"], x[0]["position"], x[0]["wtAA"], x[0]["prAA"])
            pose_position1 = pose.pdb_info().pdb2pose(mut1.chain, mut1.position)

            mut2 = Mutation(x[1]["chain"], x[1]["position"], x[1]["wtAA"], x[1]["prAA"])
            pose_position2 = pose.pdb_info().pdb2pose(mut2.chain, mut2.position)

            name = f"{count}_{mut1.old}_{mut1.position}_{mut1.new}+{mut2.old}_{mut2.position}_{mut2.new}"

            pose_m1 = predict_ddG.mutate_residue(pose, pose_position1, mut1.new, 8.0, sfxn)
            pose_m2 = predict_ddG.mutate_residue(pose_m1, pose_position2, mut2.new, 8.0, sfxn)
            score_2 = sfxn.score(pose_m2)
            ddG = score_2 - score_1
            writer.writerow(
                {
                    "sr": count,
                    "pdb": pdb_file,
                    "combination": name,
                    "mut1_wtAA": mut1.old,
                    "chain1": mut1.chain,
                    "mut1_position": mut1.position,
                    "mut1_prAA": mut1.new,
                    "mut2_wtAA": mut2.old,
                    "chain2": mut2.chain,
                    "mut2_position": mut2.position,
                    "mut2_prAA": mut2.new,
                    "double_ddG_value": ddG,
                }
            )
            count += 1
        return out_csv


def get_mut_csv(file, total=0):
    with open(file) as f:
        reader = csv.DictReader(f)
        result = []
        count = 0
        if total == 0:
            total = 100000
        for row in reader:
            result.append(row)
            count += 1
            if count >= total:
                break
    return result


def get_inputs():
    parser = argparse.ArgumentParser(
        description="Takes input a PDB file and CSV of mutations and then calculated the Double ddG values for all the mutations.", epilog="Written by Naisarg Patel (https://github.com/naisarg14)"
    )
    parser.add_argument("-p", "--pdb", help="PDB File", metavar="PDB", required=True)
    parser.add_argument("-i", "--input", help="Input Double CSV File", metavar="CSV", required=True)
    parser.add_argument("-n", "--num", help="Number of top Single ddG values to take", metavar="N", type=int, default=50)
    parser.add_argument("-o", "--output", help="Output CSV with the ddG", metavar="CSV")

    args = parser.parse_args()
    pdb_file = args.pdb
    single_csv = args.input
    out_csv = args.output
    total = args.num


    if not os.path.isabs(pdb_file):
        print(f"Assuming current directory for {pdb_file} as root since path not specified.")
        dir = os.getcwd()
        file = pdb_file
    else:
        dir, file = os.path.split(os.path.abspath(pdb_file))

    full_pdb_path = os.path.join(dir, file)

    if not os.path.isfile(full_pdb_path):
        sys.exit("No such file found in current directory, enter full path for other directories.")

    if not file.endswith(".pdb"):
        sys.exit("Given file is not a PDB file, input should be a PDB file.")

    
    if not os.path.isabs(single_csv):
        print(f"Assuming current directory for {single_csv} as root since path not specified.")
        dir = os.getcwd()
        file = single_csv
    else:
        dir, file = os.path.split(os.path.abspath(single_csv))

    full_csv_path = os.path.join(dir, file)

    if not os.path.isfile(full_csv_path):
        sys.exit("No such file found in current directory, enter full path for other directories.")

    if not full_csv_path.endswith(".csv"):
        sys.exit("Given file is not a CSV file, input should be a CSV file.")

    if out_csv and not out_csv.endswith(".csv"):
        out_csv += ".csv"

    return full_pdb_path, full_csv_path, total, out_csv


if __name__ == "__main__":
    main()
