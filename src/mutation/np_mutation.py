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


import os, sys, time, argparse
from csv_generator import generate_csv
from ddg_calc import calc_ddg
from csv_sort import sort_csv
from ddg_calc_double import calc_double_ddg
from ddg_calc_triple import calc_triple_ddg
from generate_mutants import generate_single_mutation, generate_double_mutation, generate_triple_mutation
from helpers import backup, file_info, permutations, clean_pdb
from contextlib import contextmanager


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
    start_time = time.time()
    full_pdb_path, num_double_ddg, num_single_mut, num_double_mut, num_triple_ddg, num_triple_mut, append, quiet = get_inputs()

    #Clean the PDB file
    if not quiet: print(f"Cleaning the PDB file {full_pdb_path}")
    cleaned_pdb = f"{full_pdb_path.removesuffix('.pdb')}_clean.pdb"
    with suppress_stdout():
        clean_pdb(full_pdb_path, cleaned_pdb)
    if not quiet: print(f"PDB file cleaned successfully and saved as {cleaned_pdb}.")

    full_pdb_path = cleaned_pdb

    # Generates possible mutations
    if not quiet: print(f"Generating mutations list for {full_pdb_path}")
    out_all = f"{full_pdb_path.removesuffix('.pdb')}_mutations_all.csv"
    out_op = f"{full_pdb_path.removesuffix('.pdb')}_mutations.csv"
    if file_info(out_all)[0] and file_info(out_op)[0] and append:
        mut_length = file_info(out_op)[1]
        mut_out = out_op
        if not quiet: print(f"Files already exist, skipping generation of mutations list.")
    else:
        with suppress_stdout():
            mut_out, _ = generate_csv(full_pdb_path, out_all=out_all, out_op=out_op)
        if not quiet: print(f"CSV for mutations generated successfully and saved as {mut_out}.")


    # Calculate ddG for the generated mutants
    if not quiet: print(f"Calculating ddG values for for {mut_out}")
    out_file = f"{full_pdb_path.removesuffix('.csv')}_ddG.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == mut_length:
        ddg_out = out_file
        if not quiet: print(f"File already exists, skipping calculation of ddG values.")
    else:
        with suppress_stdout():
            ddg_out = calc_ddg(full_pdb_path, mut_out, out_file=out_file)
        if not quiet: print(f"ddG values for mutations calculated successfully and saved as {ddg_out}.")


    # Sorting the ddG values in descending order
    if not quiet: print(f"Sorting the Single ddG values.")
    out_file = f"{ddg_out.removesuffix('.csv')}_sorted.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == mut_length:
        ddg_out_sort = out_file
        if not quiet: print(f"File already exists, skipping sorting of ddG values.")
    else:
        with suppress_stdout():
            ddg_out_sort = sort_csv(ddg_out, out_file=out_file, col_num=8)
        if not quiet: print(f"ddG values are sorted and stored as {ddg_out_sort}.")


    # Calculate Double ddG for the generated mutants
    if not quiet: print(f"Calculating Double ddG values for {ddg_out_sort}.")
    out_csv = f"{full_pdb_path.removesuffix('.pdb')}_double_ddg.csv"
    if file_info(out_csv)[0] and append and file_info(out_csv)[1] == permutations(num_double_ddg, 2):
        double_ddg_out = out_csv
        if not quiet: print(f"File already exists, skipping calculation of double ddG values.")
    else:
        with suppress_stdout():
            double_ddg_out = calc_double_ddg(full_pdb_path, out_csv=out_csv, single_csv=ddg_out_sort, total=num_double_ddg)
        if not quiet: print(f"Double ddG values for mutations calculated successfully and saved as {double_ddg_out}.")


    # Sorting the ddG values in descending order
    if not quiet: print(f"Sorting the Double ddG values.")
    out_file = f"{double_ddg_out.removesuffix('.csv')}_sorted.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == permutations(num_double_ddg, 2):
        double_ddg_out_sort = out_file
        if not quiet: print(f"File already exists, skipping sorting of double ddG values.")
    else:
        with suppress_stdout():
            double_ddg_out_sort = sort_csv(double_ddg_out, out_file=out_file, col_num=11)
        if not quiet: print(f"Double ddG values are sorted and stored as {double_ddg_out_sort}.")


    # Calculate Triple ddG for the generated mutants
    if not quiet: print(f"Calculating Triple ddG values for {double_ddg_out_sort}")
    out_csv = f"{double_ddg_out_sort.removesuffix('.pdb')}_triple_ddg.csv"
    if file_info(out_csv)[0] and append and file_info(out_csv)[1] == permutations(num_triple_ddg, 3):
        triple_ddg_out = out_csv
        if not quiet: print(f"File already exists, skipping calculation of triple ddG values.")
    else:
        with suppress_stdout():
            triple_ddg_out = calc_triple_ddg(pdb_file=full_pdb_path, double_csv=double_ddg_out_sort, out_csv=None, total=num_triple_ddg)
        if not quiet: print(f"Triple ddG values for mutations calculated successfully and saved as {triple_ddg_out}.")


    # Sorting the Triple ddG values in descending order
    if not quiet: print(f"Sorting the Triple ddG values.")
    out_file = f"{triple_ddg_out.removesuffix('.csv')}_sorted.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == permutations(num_double_ddg, 2):
        double_ddg_out_sort = out_file
        if not quiet: print(f"File already exists, skipping sorting of triple ddG values.")
    else:
        with suppress_stdout():
            triple_ddg_out_sort = sort_csv(triple_ddg_out, col_num=11)
        if not quiet: print(f"Triple ddG values are sorted and stored as {triple_ddg_out_sort}")


    # Generate mutants and put them in a folder and csv with names for easy docking
    if not quiet: print(f"Generating mutants for {ddg_out_sort}")
    backup(f"{full_pdb_path.removesuffix('.pdb')}_mutants.txt")
    _, pdb_name = os.path.split(os.path.abspath(full_pdb_path))
    mut_folder = f"mutation_{pdb_name.removesuffix('.pdb')}"
    with suppress_stdout():
        out_folder = generate_single_mutation(full_pdb_path, ddg_out_sort, total=num_single_mut, folder=mut_folder, text_file=f"{full_pdb_path.removesuffix('.pdb')}_mutants.txt")
    if not quiet: print(f"Generated single mutations PDB in {out_folder}")

    if not quiet: print(f"Generating mutants for {double_ddg_out_sort}")
    with suppress_stdout():
        out_folder = generate_double_mutation(full_pdb_path, double_ddg_out_sort, total=num_double_mut, folder=mut_folder, text_file=f"{full_pdb_path.removesuffix('.pdb')}_mutants.txt")
    if not quiet: print(f"Generated double mutations PDB in {out_folder}")

    if not quiet: print(f"Generating mutants for {triple_ddg_out_sort}")
    with suppress_stdout():
        out_folder = generate_triple_mutation(full_pdb_path, triple_ddg_out_sort, total=num_triple_mut, folder=mut_folder, text_file=f"{full_pdb_path.removesuffix('.pdb')}_mutants.txt")
    if not quiet: print(f"Generated triple mutations PDB in {out_folder}")


    end_time = time.time()
    elapsed_time = (end_time - start_time)/60
    print(f"Completed in {elapsed_time:.2f} minutes!")



def get_inputs():
    parser = argparse.ArgumentParser(prog="np_mutation", description=None, epilog="Part of mutadock library. Written by Naisarg Patel (https://github.com/naisarg14)")
    parser.add_argument("-i",'--input', help="PDB File for predicting and generating mutations", metavar="PDB", required=True)
    parser.add_argument("-s", "--double", help="Number of Compounds for Double ddG", metavar="D", type=int, default=100)
    parser.add_argument("-t", "--triple", help="Number of Compounds for Triple ddG", metavar="T", type=int, default=40)
    parser.add_argument("-S", "--spm", help="Number of Compounds for Single Point Mutations", metavar="SM", type=int, default=15)
    parser.add_argument("-D", "--dpm", help="Number of Compounds for Double Point Mutations", metavar="DM", type=int, default=15)
    parser.add_argument("-T", "--tpm", help="Number of Compounds for Triple Point Mutations", metavar="TM", type=int, default=15)
    parser.add_argument("--noappend", help="Generate all files again", action="store_false")
    parser.add_argument("--quiet", help="Run the Docking in quiet mode (default: False).", action="store_true", default=False)
    
    args = parser.parse_args()

    pdb_file = args.input

    if not os.path.isabs(pdb_file):
        print("Assuming current directory as root since path not specified.")
        dir = os.getcwd()
        file = pdb_file
    else:
        dir, file = os.path.split(os.path.abspath(pdb_file))

    full_pdb_path = os.path.join(dir, file)

    if not os.path.isfile(full_pdb_path):
        sys.exit("No such file found in current directory, enter full path for other directories.")

    if not file.endswith(".pdb"):
        sys.exit("Given file is not a PDB file, input should be a PDB file.")


    return full_pdb_path, args.double, args.spm, args.dpm, args.triple, args.tpm, args.noappend, args.quiet




if __name__ == "__main__":
    naisarg()
