################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2026 Naisarg Patel (https://github.com/naisarg14)        #
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


import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Optional

from .Amino import get_dict, get_scfn_250
from .exceptions import CSVGenerationError, MutationError, PDBFileError
from .helpers import backup, convert_cif_pdb

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

try:
    from Bio.PDB import PDBParser
except ImportError:
    msg = "Error with importing biopython module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install biopython using the following command:\n\n"
    msg += "python -m pip install biopython\n"
    msg += "If you already have biopython installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    logger.error(msg)
    sys.exit(2)


def main() -> None:
    """CLI entry point for ``md_csv_generator``."""
    pdb_file, out_op, out_all = get_inputs()
    generate_csv(pdb_file, out_all, out_op)


def generate_csv(
    pdb_file: str, out_all: Optional[str] = None, out_op: Optional[str] = None
) -> tuple[str, str]:
    """Generate mutation CSV files for all possible substitutions in a PDB.

    Args:
        pdb_file: Path to the PDB file.
        out_all: Output path for the full mutations CSV (all prProb values).
        out_op: Output path for the filtered CSV (prProb > 0).

    Returns:
        ``(out_op, out_all)`` paths.

    Raises:
        PDBFileError: If the PDB file cannot be read.
        CSVGenerationError: If the PDB file contains no residues.
    """
    if not out_all:
        out_all = f"{pdb_file.removesuffix('.pdb')}_mutations_all.csv"
    if not out_op:
        out_op = f"{pdb_file.removesuffix('.pdb')}_mutations.csv"

    if out_all == out_op:
        out_all = out_all.removesuffix(".csv") + "_all.csv"
    residues = get_residues(pdb_file)
    if not residues:
        raise CSVGenerationError(f"No residues found in PDB file '{pdb_file}'")

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
    for key in residues:
        residue = residues[key]
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


def get_residues(file: str) -> dict[int, tuple[str, int, str]]:
    """Parse residue information from a PDB file.

    Args:
        file: Path to the PDB file.

    Returns:
        Mapping of sequential index to ``(chain_id, position, residue_name)``.

    Raises:
        PDBFileError: If the file is not found.
    """
    residues: dict[int, tuple[str, int, str]] = {}
    count = 0
    parser = PDBParser(PERMISSIVE=1)
    try:
        structure = parser.get_structure(file, file)
    except FileNotFoundError as e:
        raise PDBFileError(f"PDB file not found: '{file}'") from e
    for model in structure:
        for chain in model:
            for residue in chain:
                if count == 0:
                    count += 1
                    continue
                full_id = residue.get_full_id()
                chain_id = full_id[2]
                position = full_id[3][1]
                name = residue.get_resname()
                residues[count] = (chain_id, position, name)
                count += 1
    return residues


def get_inputs() -> tuple[str, Optional[str], Optional[str]]:
    parser = argparse.ArgumentParser(
        description="This program takes as input a PDB file and generates all possible mutations and also gives the PAM250 score.",
        epilog="Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="PDB File for predicting mutations",
        metavar="PDB",
        required=True,
    )
    parser.add_argument(
        "-o", "--positive", help="Output CSV for only non-zero values", metavar="FILE"
    )
    parser.add_argument(
        "-O", "--all", help="Output CSV for all posible mutations", metavar="FILE"
    )

    args = parser.parse_args()

    pdb_file = args.input

    pdb_path = Path(pdb_file)
    if not pdb_path.is_absolute():
        logger.info("Assuming current directory as root since path not specified.")
        full_pdb_path = Path.cwd() / pdb_file
    else:
        full_pdb_path = pdb_path

    if full_pdb_path.suffix == ".cif":
        try:
            full_pdb_path = Path(convert_cif_pdb(str(full_pdb_path)))
        except MutationError as e:
            logger.error("Error converting CIF to PDB: %s", e)
            sys.exit(2)

    if not full_pdb_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_pdb_path.suffix != ".pdb":
        sys.exit("Given file is not a PDB file, input should be a PDB file.")

    return str(full_pdb_path), args.positive, args.all


if __name__ == "__main__":
    main()
