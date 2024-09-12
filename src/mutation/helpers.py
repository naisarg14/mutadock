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


import os
from datetime import datetime
import shutil
from Amino import check_3, get_1

class Mutation():
    def __init__(self, chain, position, aa_old, aa_new):
        self.chain = chain
        self.position = int(position)
        if check_3(aa_old):
            self.old = get_1(aa_old)
        else:
            self.old = aa_old
        if check_3(aa_new):
            self.new = get_1(aa_new)
        else:
            self.new = aa_new


    def __str__(self):
        return f"The aa at {self.chain}{self.position} is mutated from {self.aa_old} to {self.aa_new}."


def backup(file_path):
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


def move_file(file_path, folder):
    master_folder, file = os.path.split(os.path.abspath(file_path))
    folder_path = os.path.join(master_folder, folder)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    shutil.move(file_path, os.path.join(folder_path, file))
    return (folder_path, (os.path.join(folder_path, file)))


def process_vina(file):
    master_folder, out_file = os.path.split(os.path.abspath(file))
    ligand1 = os.path.join(master_folder, f"{out_file}_ligand_1.pdbqt")
    for i in range(2, 10):
        try:
            file_path = os.path.join(master_folder, f"{out_file}_ligand_{i}.pdbqt")
            os.remove(file_path)
        except FileNotFoundError:
            pass
    return ligand1

def file_info(file_path):
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            row_count = sum(1 for _ in file)
        return (True, row_count)
    else:
        return (False, 0)

def permutations(n, r):
    from math import factorial
    if n < r:
        return 0
    else:
        return int(factorial(n) / factorial(n - r))
    

def clean_pdb(pdb_file, pdb_clean=None):
    from datetime import datetime
    if not pdb_clean:
        pdb_clean = pdb_file.replace(".pdb", "_clean.pdb")
        
    try:
        with open(pdb_file, 'r') as infile, open(pdb_clean, 'w+') as outfile:
            outfile.write("REMARK This file was cleaned using the clean_pdb function of mutadock library.\n")
            outfile.write("REMARK All non-ATOM lines were removed.\n")
            outfile.write(f"REMARK Created on {datetime.now().strftime('%b-%d-%Y %H:%M:%S')}\n")
            for line in infile:
                if line.startswith("ATOM") or line.startswith("TER"):
                    outfile.write(line)
                elif line.strip() == "END":
                    outfile.write(line)
                    break
        return (True, pdb_clean)
    except Exception as e:
        return (False, e)


if __name__ == "__main__":
    print("This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's mutation module.")