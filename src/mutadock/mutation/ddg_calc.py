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
import os
import sys
from pathlib import Path
from typing import Optional

from . import predict_ddG
from .Amino import get_1
from .helpers import backup, in_directory

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

try:
    sys.stdout = open(os.devnull, "w")
    from pyrosetta import *

    sys.stdout = sys.__stdout__
except ImportError:
    sys.stdout = sys.__stdout__
    msg = "Error with importing pyrosetta module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install pyrosetta using the following command:\n\n"
    msg += "python -m pip install pyrosetta_installer && python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'\n"
    msg += "Alternative is to install from PyRosetta's official website.\n"
    msg += "If you already have pyrosetta installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    logger.error(msg)
    sys.exit(2)

try:
    from tqdm import tqdm
except ImportError:
    msg = "Error with importing tqdm module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install tqdm using the following command:\n\n"
    msg += "python -m pip install tqdm\n"
    msg += "If you already have tqdm installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    logger.error(msg)
    sys.exit(2)


def main() -> None:
    full_pdb_path, full_csv_path, out_csv = get_inputs()
    calc_ddg(pdb_file=full_pdb_path, in_file=full_csv_path, out_file=out_csv)


def calc_ddg(pdb_file: str, in_file: str, out_file: Optional[str] = None) -> str:

    sys.stdout = open(os.devnull, "w")
    init("-mute all")
    sys.stdout = sys.__stdout__

    if not out_file:
        out_file = f"{in_file.removesuffix('.csv')}_ddG.csv"
    backup(out_file)

    with open(out_file, "w+") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "sr",
                "pdb",
                "chain",
                "position",
                "wtAA",
                "prAA",
                "wtProb",
                "prProb",
                "ddG_value",
            ],
        )
        writer.writeheader()
        with open(in_file) as csvfile:
            reader = csv.DictReader(csvfile)
            pdb_path = Path(pdb_file).resolve()
            with in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)
            sfxn = get_fa_scorefxn()
            score_1 = sfxn.score(pose)
            for row in tqdm(reader):
                chain = row["chain"]
                position = int(row["position"])
                mutation = row["prAA"]
                pose_position = pose.pdb_info().pdb2pose(chain, position)
                aa_str = get_1(mutation)
                assert aa_str is not None, f"Unknown amino acid code: {mutation}"
                mut_pose = predict_ddG.mutate_residue(
                    pose, pose_position, aa_str, 8.0, sfxn
                )
                score_2 = sfxn.score(mut_pose)
                ddG = score_2 - score_1
                row["ddG_value"] = ddG
                writer.writerow(row)
    return out_file


def get_inputs() -> tuple[str, str, Optional[str]]:
    parser = argparse.ArgumentParser(
        description="Takes input a PDB file and CSV of mutations and then calculated the ddG values for all the mutations. Input CSV can be easily generated by using the 'csv_generator.py' in the package.",
        epilog="Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument("-p", "--pdb", help="PDB File", metavar="PDB", required=True)
    parser.add_argument(
        "-i", "--input", help="Input CSV File", metavar="CSV", required=True
    )
    parser.add_argument("-o", "--output", help="Output CSV with the ddG", metavar="CSV")

    args = parser.parse_args()

    pdb_file = args.pdb
    mut_csv = args.input
    out_csv = args.output

    p = Path(pdb_file)
    if not p.is_absolute():
        logger.info(
            f"Assuming current directory for {pdb_file} as root since path not specified."
        )
        full_pdb_path = Path.cwd() / pdb_file
    else:
        full_pdb_path = p.resolve()

    if not full_pdb_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_pdb_path.suffix != ".pdb":
        logger.error("Given file is not a PDB file, input should be a PDB file.")
        sys.exit(
            "Usage: python ddg_calc_double.py <PDB file> <single_csv> <total> (optional) <out_csv>"
        )

    p = Path(mut_csv)
    if not p.is_absolute():
        logger.info(
            f"Assuming current directory for {mut_csv} as root since path not specified."
        )
        full_csv_path = Path.cwd() / mut_csv
    else:
        full_csv_path = p.resolve()

    if not full_csv_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_csv_path.suffix != ".csv":
        logger.error("Given file is not a CSV file, input should be a CSV file.")
        sys.exit(
            "Usage: python ddg_calc_double.py <PDB file> <single_csv> (optional) <out_csv>"
        )

    return str(full_pdb_path), str(full_csv_path), out_csv


if __name__ == "__main__":
    main()
