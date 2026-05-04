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


import csv
import logging
import os
import sys
from pathlib import Path
from typing import Optional

from . import helpers
from .Amino import get_1

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

try:
    sys.stdout = open(os.devnull, "w")
    from pyrosetta import *
    from pyrosetta.toolbox import mutate_residue

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
    """Stub entry point — this module is used as a library, not run directly."""
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's docking module."
    )


def generate_single_mutation(
    pdb_file: str,
    csv_file: str,
    total: int = 0,
    folder: Optional[str] = None,
    text_file: Optional[str] = None,
) -> str:
    """Generate PDB files for the top single-point mutations.

    Reads each row from the sorted ddG CSV, applies the substitution to the
    wild-type structure using PyRosetta's ``mutate_residue``, dumps the mutant
    PDB, and optionally moves it into *folder* and records its path in
    *text_file*.

    Args:
        pdb_file: Path to the cleaned wild-type PDB file.
        csv_file: Path to the sorted single-ddG CSV (output of ``sort_csv``).
        total: Maximum number of mutants to generate.  0 means generate all.
        folder: Sub-folder name (relative to the PDB directory) to collect the
            mutant PDB files.  Skipped if ``None``.
        text_file: Path to a text file where generated mutant paths are
            appended, one per line.  Skipped if ``None``.

    Returns:
        Path to the output folder containing the mutant PDB files.
    """
    if int(total) == 0:
        idk = 0
    else:
        idk = 1
    init("-mute all")
    count = 1
    out_folder: str = ""
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            sr = row["sr"]
            chain = row["chain"]
            position = int(row["position"])
            aa_old = row["wtAA"]
            aa_new = row["prAA"]

            name = str(
                Path(pdb_file).parent / f"{sr}_{aa_old}-{chain}{position}-{aa_new}.pdb"
            )

            pdb_path = Path(pdb_file).resolve()
            with helpers.in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)
            pose_position = pose.pdb_info().pdb2pose(chain, position)
            aa_str = get_1(aa_new)
            assert aa_str is not None, f"Unknown amino acid code: {aa_new}"
            mutate_residue(pose, pose_position, aa_str)

            pose.dump_pdb(name)

            out_file = name
            if folder:
                out_folder, out_file = helpers.move_file(name, folder)
            if text_file:
                with open(text_file, "a+") as file:
                    file.write(f"{out_file}\n")
            if count >= total:
                break
            count += idk
    return out_folder


def generate_double_mutation(
    pdb_file: str,
    csv_file: str,
    total: int = 10,
    folder: Optional[str] = None,
    text_file: Optional[str] = None,
) -> str:
    """Generate PDB files for the top double-point mutations.

    Reads each row from the sorted double-ddG CSV, applies both substitutions
    sequentially to the wild-type structure using PyRosetta's
    ``mutate_residue``, dumps the mutant PDB, and optionally moves it into
    *folder* and records its path in *text_file*.

    Args:
        pdb_file: Path to the cleaned wild-type PDB file.
        csv_file: Path to the sorted double-ddG CSV.
        total: Maximum number of mutants to generate.  0 means generate all.
        folder: Sub-folder name (relative to the PDB directory) to collect the
            mutant PDB files.  Skipped if ``None``.
        text_file: Path to a text file where generated mutant paths are
            appended, one per line.  Skipped if ``None``.

    Returns:
        Path to the output folder containing the mutant PDB files.
    """
    if int(total) == 0:
        idk = 0
    else:
        idk = 1
    init("-mute all")
    count = 1
    out_folder: str = ""
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            name = str(Path(pdb_file).parent / f"{row['sr']}_{row['combination']}.pdb")

            chain1 = row["chain1"]
            position1 = int(row["mut1_position"])
            aa_new1 = row["mut1_prAA"]

            chain2 = row["chain2"]
            position2 = int(row["mut2_position"])
            aa_new2 = row["mut2_prAA"]

            pdb_path = Path(pdb_file).resolve()
            with helpers.in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)

            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)

            aa_str1 = get_1(aa_new1)
            aa_str2 = get_1(aa_new2)
            assert aa_str1 is not None, f"Unknown amino acid code: {aa_new1}"
            assert aa_str2 is not None, f"Unknown amino acid code: {aa_new2}"
            mutate_residue(pose, pose_position1, aa_str1)
            mutate_residue(pose, pose_position2, aa_str2)

            pose.dump_pdb(name)

            out_file = name
            if folder:
                out_folder, out_file = helpers.move_file(name, folder)
            if text_file:
                with open(text_file, "a+") as file:
                    file.write(f"{out_file}\n")

            if count >= total:
                break
            count += idk
    return out_folder


def generate_triple_mutation(
    pdb_file: str,
    csv_file: str,
    total: int = -1,
    folder: Optional[str] = None,
    text_file: Optional[str] = None,
) -> str:
    """Generate PDB files for the top triple-point mutations.

    Reads each row from the sorted triple-ddG CSV, applies all three
    substitutions sequentially to the wild-type structure using PyRosetta's
    ``mutate_residue``, dumps the mutant PDB, and optionally moves it into
    *folder* and records its path in *text_file*.

    Args:
        pdb_file: Path to the cleaned wild-type PDB file.
        csv_file: Path to the sorted triple-ddG CSV.
        total: Maximum number of mutants to generate.  -1 means generate all.
        folder: Sub-folder name (relative to the PDB directory) to collect the
            mutant PDB files.  Skipped if ``None``.
        text_file: Path to a text file where generated mutant paths are
            appended, one per line.  Skipped if ``None``.

    Returns:
        Path to the output folder containing the mutant PDB files.
    """
    if int(total) == -1:
        idk = 0
    else:
        idk = 1
    init("-mute all")
    count = 1
    out_folder: str = ""
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in tqdm(reader):
            name = str(Path(pdb_file).parent / f"{row['sr']}_{row['combination']}.pdb")

            chain1 = row["chain1"]
            position1 = int(row["mut1_position"])
            aa_new1 = row["mut1_prAA"]

            chain2 = row["chain2"]
            position2 = int(row["mut2_position"])
            aa_new2 = row["mut2_prAA"]

            chain3 = row["chain3"]
            position3 = int(row["mut3_position"])
            aa_new3 = row["mut3_prAA"]

            pdb_path = Path(pdb_file).resolve()
            with helpers.in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)

            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)
            pose_position3 = pose.pdb_info().pdb2pose(chain3, position3)

            aa_str1 = get_1(aa_new1)
            aa_str2 = get_1(aa_new2)
            aa_str3 = get_1(aa_new3)
            assert aa_str1 is not None, f"Unknown amino acid code: {aa_new1}"
            assert aa_str2 is not None, f"Unknown amino acid code: {aa_new2}"
            assert aa_str3 is not None, f"Unknown amino acid code: {aa_new3}"
            mutate_residue(pose, pose_position1, aa_str1)
            mutate_residue(pose, pose_position2, aa_str2)
            mutate_residue(pose, pose_position3, aa_str3)

            pose.dump_pdb(name)

            out_file = name
            if folder:
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
