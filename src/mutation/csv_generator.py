import csv, os, sys
from Amino import get_dict, get_scfn_250
from Bio.PDB import PDBParser
from helpers import backup
import argparse


def main():
    pdb_file, out_op, out_all = get_inputs()
    generate_csv(pdb_file, out_all, out_op)


def generate_csv(pdb_file, out_all=None, out_op=None):
    if not out_all:
        out_all = f"{pdb_file.removesuffix('.pdb')}_mutations_all.csv"
    if not out_op:
        out_op = f"{pdb_file.removesuffix('.pdb')}_mutations.csv"
    
    if out_all == out_op:
        out_all = out_all.removesuffix(".csv") + "_all.csv"
    residues = get_residues(pdb_file)
    if not residues:
        print("Check File Name")
        return None

    backup(out_op)
    backup(out_all)
    f = open(out_op, "w")
    f_all = open(out_all, "w")

    writer_all = csv.writer(f_all)
    writer_op = csv.writer(f)
    header = [
            "sr",
            "pdb",
            "chain",
            "position",
            "wtAA",
            "prAA",
            "wtProb",
            "prProb",
        ]
    writer_all.writerow(header)
    writer_op.writerow(header)
    count1 = 1
    count2 = 1
    aa_dict = get_dict()
    score_dict = get_scfn_250()
    for residue in residues:
        residue = residues[residue]
        for aa in aa_dict:
            if aa == residue[2]:
                continue
            prProb = float(score_dict[residue[2]][aa] / 100)
            row = [
                count1,
                pdb_file.removesuffix(".pdb"),
                residue[0],
                residue[1],
                (residue[2]),
                (aa),
                float(score_dict[residue[2]][residue[2]] / 100),
                prProb,
            ]
            writer_all.writerow(row)
            if prProb > 0.0:
                row[0] = count2
                writer_op.writerow(row)
                count2 += 1
            count1 += 1
    f.close()
    return out_op, out_all


def get_residues(file):
    residues = {}
    count = 0
    parser = PDBParser(PERMISSIVE=1)
    try:
        structure = parser.get_structure(file, file)
    except FileNotFoundError:
        return None
    for model in structure:
        for chain in model:
            for residue in chain:
                if count == 0:
                    count += 1
                    continue
                full_id = residue.get_full_id()
                chain = full_id[2]
                position = full_id[3][1]
                name = residue.get_resname()
                residues[count] = (chain, position, name)
                count += 1
    return residues

def get_inputs():
    parser = argparse.ArgumentParser(description="This program takes as input a PDB file and generates all possible mutations and also gives the PAM250 score.", epilog="Written by Naisarg Patel (https://github.com/naisarg14)")
    parser.add_argument("-i",'--input', help="PDB File for predicting mutations", metavar="PDB", required=True)
    parser.add_argument("-o", "--positive", help="Output CSV for only non-zero values", metavar="FILE")
    parser.add_argument("-O", "--all", help="Output CSV for all posible mutations", metavar="FILE")

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

    return full_pdb_path, args.positive, args.all


if __name__ == "__main__":
    main()
