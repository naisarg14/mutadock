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


import csv, os, sys
import mutation.helpers
from mutation.Amino import get_1

try:
    from pyrosetta import *
    from pyrosetta.toolbox import mutate_residue
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
    print("This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's docking module.")

def generate_single_mutation(pdb_file, csv_file, total=0, folder=None, text_file=None):
    if int(total) == 0:
        idk = 0
    else:
        idk = 1
    init('-mute all')
    count = 1
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            sr = row["sr"]
            chain = row["chain"]
            position = int(row["position"])
            aa_old = row["wtAA"]
            aa_new = row["prAA"]

            name = f"{os.path.dirname(pdb_file)}/{sr}_{aa_old}-{chain}{position}-{aa_new}.pdb"

            pose = pose_from_pdb(pdb_file)
            pose_position = pose.pdb_info().pdb2pose(chain, position)
            mutate_residue(pose, pose_position, get_1(aa_new))
            
            pose.dump_pdb(name)

            if folder:
                master_folder, file = os.path.split(os.path.abspath(csv_file))
                out_folder, out_file = helpers.move_file(name, folder)
            if text_file:
                with open(text_file, "a+") as file:
                    file.write(f"{out_file}\n")
            if count >= total:
                break
            count += idk
    return out_folder


def generate_double_mutation(pdb_file, csv_file, total=10, folder=None, text_file=None):
    if int(total) == 0:
        idk = 0
    else:
        idk = 1
    init('-mute all')
    count = 1
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            name = f"{os.path.dirname(pdb_file)}/{row['sr']}_{row['combination']}.pdb"

            chain1 = row["chain1"]
            position1 = int(row["mut1_position"])
            aa_new1 = row["mut1_prAA"]

            chain2 = row["chain2"]
            position2 = int(row["mut2_position"])
            aa_new2 = row["mut2_prAA"]

            pose = pose_from_pdb(pdb_file)

            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)

            mutate_residue(pose, pose_position1, get_1(aa_new1))
            mutate_residue(pose, pose_position2, get_1(aa_new2))

            pose.dump_pdb(name)

            if folder:
                master_folder, file = os.path.split(os.path.abspath(csv_file))
                out_folder, out_file = helpers.move_file(name, folder)
            if text_file:
                with open(text_file, "a+") as file:
                    file.write(f"{out_file}\n")

            if count >= total:
                break
            count += idk
    return out_folder


def generate_triple_mutation(pdb_file, csv_file, total=-1, folder=None, text_file=None):
    if int(total) == -1:
        idk = 0
    else:
        idk = 1
    init('-mute all')
    count = 1
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            name = f"{os.path.dirname(pdb_file)}/{row['sr']}_{row['combination']}.pdb"

            chain1 = row["chain1"]
            position1 = int(row["mut1_position"])
            aa_new1 = row["mut1_prAA"]

            chain2 = row["chain2"]
            position2 = int(row["mut2_position"])
            aa_new2 = row["mut2_prAA"]

            chain3 = row["chain3"]
            position3 = int(row["mut3_position"])
            aa_new3 = row["mut3_prAA"]

            pose = pose_from_pdb(pdb_file)

            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)
            pose_position3 = pose.pdb_info().pdb2pose(chain3, position3)

            mutate_residue(pose, pose_position1, get_1(aa_new1))
            mutate_residue(pose, pose_position2, get_1(aa_new2))
            mutate_residue(pose, pose_position3, get_1(aa_new3))

            pose.dump_pdb(name)

            if folder:
                master_folder, file = os.path.split(os.path.abspath(csv_file))
                out_folder, out_file = helpers.move_file(name, folder)
            if text_file:
                with open(text_file, "a+") as file:
                    file.write(f"{out_file}\n")

            if count >= total:
                break
            count += idk
    return out_folder


if __name__ == "__main__":
    main()
