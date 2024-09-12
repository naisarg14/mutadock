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


import predict_ddG
from itertools import combinations
from Amino import get_1
from helpers import backup
import os, sys, argparse, csv

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
    calc_triple_ddg(pdb_file=full_pdb_path, double_csv=full_csv_path, total=total, out_csv=out_csv)


def calc_triple_ddg(pdb_file, double_csv, out_csv=None, total=0):
    if not out_csv:
        out_csv = f"{pdb_file.removesuffix('.pdb')}_triple_ddg.csv"
    backup(out_csv)

    result = get_mut_csv(double_csv, total)

    sys.stdout = open(os.devnull, 'w')
    init('-mute all')
    sys.stdout = sys.__stdout__
    
    with open(out_csv, "w+") as out:
        writer = csv.DictWriter(out, fieldnames=["sr", "pdb", "combination", "mut1_wtAA", "chain1", "mut1_position", "mut1_prAA", "mut2_wtAA", "chain2", "mut2_position", "mut2_prAA", "mut3_wtAA", "chain3", "mut3_position", "mut3_prAA", "triple_ddG_value"])
        writer.writeheader()
        count = 1
        pose = pose_from_pdb(pdb_file)
        sfxn = get_fa_scorefxn()
        score_1 = sfxn.score(pose)
        for x in tqdm(combinations(result, 3)):
            mut1 = x[0]
            mut2 = x[1]
            mut3 = x[2]

            chain1 = mut1["chain"]
            position1 = int(mut1["position"])
            old1 = mut1["wtAA"]
            mutation1 = mut1["prAA"]
            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)

            chain2 = mut2["chain"]
            position2 = int(mut2["position"])
            old2 = mut2["wtAA"]
            mutation2 = mut2["prAA"]
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)

            chain3 = mut3["chain"]
            position3 = int(mut3["position"])
            old3 = mut3["wtAA"]
            mutation3 = mut3["prAA"]
            pose_position3 = pose.pdb_info().pdb2pose(chain3, position3)

            name = f"{count}_{old1}_{position1}_{mutation1}+{old2}_{position2}_{mutation2}+{old3}_{position3}_{mutation3}"

            pose_m1 = predict_ddG.mutate_residue(pose, pose_position1, get_1(mutation1), 8.0, sfxn)
            pose_m2 = predict_ddG.mutate_residue(pose_m1, pose_position2, get_1(mutation2), 8.0, sfxn)
            pose_m3 = predict_ddG.mutate_residue(pose_m2, pose_position3, get_1(mutation3), 8.0, sfxn)
            score_2 = sfxn.score(pose_m3)
            ddG = score_2 - score_1
            writer.writerow(
                {
                    "sr": count,
                    "pdb": pdb_file,
                    "combination": name,
                    "mut1_wtAA": old1,
                    "chain1": chain1,
                    "mut1_position": position1,
                    "mut1_prAA": mutation1,
                    "mut2_wtAA": old2,
                    "chain2": chain2,
                    "mut2_position": position2,
                    "mut2_prAA": mutation2,
                    "mut3_wtAA": old3,
                    "chain3": chain3,
                    "mut3_position": position3,
                    "mut3_prAA": mutation3,
                    "triple_ddG_value": ddG,
                }
            )
            count += 1
        return out_csv


def get_mut_csv(file, total=-1):
    with open(file) as f:
        reader = csv.DictReader(f)
        results = []
        count = 0
        if total == -1:
            total = 100000
        for row in reader:
            mutation = {"chain": row["chain1"], "position": row["mut1_position"], "wtAA": row["mut1_wtAA"], "prAA": row["mut1_prAA"]}
            if mutation not in results:
                results.append(mutation.copy())
                count += 1

            mutation = {"chain": row["chain2"], "position": row["mut2_position"], "wtAA": row["mut2_wtAA"], "prAA": row["mut2_prAA"]}
            if mutation not in results:
                results.append(mutation.copy())
                count += 1

            if count >= total:
                break

    return results


def get_inputs():
    parser = argparse.ArgumentParser(
        description="Takes input a PDB file and CSV of mutations and then calculated the Triple ddG values for all the mutations.", epilog="Written by Naisarg Patel (https://github.com/naisarg14)"
    )
    parser.add_argument("-p", "--pdb", help="PDB File", metavar="PDB", required=True)
    parser.add_argument("-i", "--input", help="Input Double CSV File", metavar="CSV", required=True)
    parser.add_argument("-n", "--num", help="Number of top Double ddG values to take", metavar="N", type=int, default=50)
    parser.add_argument("-o", "--output", help="Output CSV with the ddG", metavar="CSV")

    args = parser.parse_args()
    pdb_file = args.pdb
    double_csv = args.input
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

    if not os.path.isabs(double_csv):
        print(f"Assuming current directory for {double_csv} as root since path not specified.")
        dir = os.getcwd()
        file = double_csv
    else:
        dir, file = os.path.split(os.path.abspath(double_csv))

    full_csv_path = os.path.join(dir, file)

    if not os.path.isfile(full_csv_path):
        sys.exit("No such file found in current directory, enter full path for other directories.")

    if not full_csv_path.endswith(".csv"):
        print("Given file is not a CSV file, input should be a CSV file.")

    if out_csv and not out_csv.endswith(".csv"):
        out_csv += ".csv"

    return full_pdb_path, full_csv_path, total, out_csv


if __name__ == "__main__":
    main()
