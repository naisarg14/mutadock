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
    

if __name__ == "__main__":
    print("This file only contains functions used in other scripts.")