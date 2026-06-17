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
import logging
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

try:
    from .csv_generator import generate_csv
    from .csv_sort import sort_csv
    from .ddg_calc import calc_ddg
    from .ddg_calc_double import calc_double_ddg
    from .ddg_calc_triple import calc_triple_ddg
    from .generate_mutants import (
        generate_double_mutation,
        generate_single_mutation,
        generate_triple_mutation,
    )
    from .helpers import backup, clean_pdb, file_info, permutations
except ImportError:
    src_dir = Path(__file__).resolve().parents[2]
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    from mutadock.mutation.csv_generator import generate_csv
    from mutadock.mutation.csv_sort import sort_csv
    from mutadock.mutation.ddg_calc import calc_ddg
    from mutadock.mutation.ddg_calc_double import calc_double_ddg
    from mutadock.mutation.ddg_calc_triple import calc_triple_ddg
    from mutadock.mutation.generate_mutants import (
        generate_double_mutation,
        generate_single_mutation,
        generate_triple_mutation,
    )
    from mutadock.mutation.helpers import backup, clean_pdb, file_info, permutations

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


@contextmanager
def suppress_stdout() -> Any:
    """Context manager that redirects stdout to /dev/null for the duration."""
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def np_mutation() -> None:
    """Run the full single/double/triple mutation workflow for a PDB file.

    Orchestrates: PDB cleaning → mutation CSV → ddG calculation → sorting →
    double/triple ddG → mutant PDB generation.  All intermediate files are
    written to the same directory as the input PDB.  CLI entry point for
    ``md_mutate``.
    """
    start_time = time.time()
    (
        full_pdb_path,
        num_double_ddg,
        num_single_mut,
        num_double_mut,
        num_triple_ddg,
        num_triple_mut,
        append,
        quiet,
    ) = get_inputs()

    # Clean the PDB file
    if not quiet:
        logger.info(f"Cleaning the input file {full_pdb_path}")
    cleaned_pdb = f"{Path(full_pdb_path).with_suffix('')}_clean.pdb"
    with suppress_stdout():
        clean_pdb(full_pdb_path, cleaned_pdb)
    if not quiet:
        logger.info(f"PDB file cleaned successfully and saved as {cleaned_pdb}.")

    full_pdb_path = cleaned_pdb
    base_name = str(
        Path(full_pdb_path).parent / Path(full_pdb_path).stem.removesuffix("_clean")
    )

    # Generates possible mutations
    if not quiet:
        logger.info(f"Generating mutations list for {full_pdb_path}")
    out_all = f"{base_name}_mutations_all.csv"
    out_op = f"{base_name}_mutations.csv"
    mut_length = 0
    if file_info(out_all)[0] and file_info(out_op)[0] and append:
        mut_length = file_info(out_op)[1]
        mut_out = out_op
        if not quiet:
            logger.info("Files already exist, skipping generation of mutations list.")
    else:
        with suppress_stdout():
            mut_out, _ = generate_csv(full_pdb_path, out_all=out_all, out_op=out_op)
        if not quiet:
            logger.info(
                f"CSV for mutations generated successfully and saved as {mut_out}."
            )

    # Calculate ddG for the generated mutants
    if not quiet:
        logger.info(f"Calculating ddG values for for {mut_out}")
    out_file = f"{base_name}_ddG.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == mut_length:
        ddg_out = out_file
        if not quiet:
            logger.info("File already exists, skipping calculation of ddG values.")
    else:
        with suppress_stdout():
            ddg_out = calc_ddg(full_pdb_path, mut_out, out_file=out_file)
        if not quiet:
            logger.info(
                f"ddG values for mutations calculated successfully and saved as {ddg_out}."
            )

    # Sorting the ddG values in descending order
    if not quiet:
        logger.info("Sorting the Single ddG values.")
    out_file = f"{ddg_out.removesuffix('.csv')}_sorted.csv"
    if file_info(out_file)[0] and append and file_info(out_file)[1] == mut_length:
        ddg_out_sort = out_file
        if not quiet:
            logger.info("File already exists, skipping sorting of ddG values.")
    else:
        with suppress_stdout():
            ddg_out_sort = sort_csv(ddg_out, out_file=out_file, col_num=8)
        if not quiet:
            logger.info(f"ddG values are sorted and stored as {ddg_out_sort}.")

    # Calculate Double ddG for the generated mutants
    if not quiet:
        logger.info(f"Calculating Double ddG values for {ddg_out_sort}.")
    out_csv = f"{base_name}_double_ddg.csv"
    if (
        file_info(out_csv)[0]
        and append
        and file_info(out_csv)[1] == permutations(num_double_ddg, 2)
    ):
        double_ddg_out = out_csv
        if not quiet:
            logger.info(
                "File already exists, skipping calculation of double ddG values."
            )
    else:
        with suppress_stdout():
            double_ddg_out = calc_double_ddg(
                full_pdb_path,
                out_csv=out_csv,
                single_csv=ddg_out_sort,
                total=num_double_ddg,
            )
        if not quiet:
            logger.info(
                f"Double ddG values for mutations calculated successfully and saved as {double_ddg_out}."
            )

    # Sorting the ddG values in descending order
    if not quiet:
        logger.info("Sorting the Double ddG values.")
    out_file = f"{double_ddg_out.removesuffix('.csv')}_sorted.csv"
    if (
        file_info(out_file)[0]
        and append
        and file_info(out_file)[1] == permutations(num_double_ddg, 2)
    ):
        double_ddg_out_sort = out_file
        if not quiet:
            logger.info("File already exists, skipping sorting of double ddG values.")
    else:
        with suppress_stdout():
            double_ddg_out_sort = sort_csv(
                double_ddg_out, out_file=out_file, col_num=11
            )
        if not quiet:
            logger.info(
                f"Double ddG values are sorted and stored as {double_ddg_out_sort}."
            )

    # Calculate Triple ddG for the generated mutants
    if not quiet:
        logger.info(f"Calculating Triple ddG values for {double_ddg_out_sort}")
    out_csv = f"{base_name}_triple_ddg.csv"
    if (
        file_info(out_csv)[0]
        and append
        and file_info(out_csv)[1] == permutations(num_triple_ddg, 3)
    ):
        triple_ddg_out = out_csv
        if not quiet:
            logger.info(
                "File already exists, skipping calculation of triple ddG values."
            )
    else:
        with suppress_stdout():
            triple_ddg_out = calc_triple_ddg(
                pdb_file=full_pdb_path,
                double_csv=double_ddg_out_sort,
                out_csv=None,
                total=num_triple_ddg,
            )
        if not quiet:
            logger.info(
                f"Triple ddG values for mutations calculated successfully and saved as {triple_ddg_out}."
            )

    # Sorting the Triple ddG values in descending order
    if not quiet:
        logger.info("Sorting the Triple ddG values.")
    out_file = f"{triple_ddg_out.removesuffix('.csv')}_sorted.csv"
    if (
        file_info(out_file)[0]
        and append
        and file_info(out_file)[1] == permutations(num_double_ddg, 2)
    ):
        triple_ddg_out_sort = out_file
        if not quiet:
            logger.info("File already exists, skipping sorting of triple ddG values.")
    else:
        with suppress_stdout():
            triple_ddg_out_sort = sort_csv(triple_ddg_out, col_num=11)
        if not quiet:
            logger.info(
                f"Triple ddG values are sorted and stored as {triple_ddg_out_sort}"
            )

    # Generate mutants and put them in a folder and csv with names for easy docking
    if not quiet:
        logger.info(f"Generating mutants for {ddg_out_sort}")
    backup(f"{base_name}_mutants.txt")
    mut_folder = f"mutation_{Path(base_name).stem}"
    with suppress_stdout():
        out_folder = generate_single_mutation(
            full_pdb_path,
            ddg_out_sort,
            total=num_single_mut,
            folder=mut_folder,
            text_file=f"{base_name}_mutants.txt",
        )
    if not quiet:
        logger.info(f"Generated single mutations PDB in {out_folder}")

    if not quiet:
        logger.info(f"Generating mutants for {double_ddg_out_sort}")
    with suppress_stdout():
        out_folder = generate_double_mutation(
            full_pdb_path,
            double_ddg_out_sort,
            total=num_double_mut,
            folder=mut_folder,
            text_file=f"{base_name}_mutants.txt",
        )
    if not quiet:
        logger.info(f"Generated double mutations PDB in {out_folder}")

    if not quiet:
        logger.info(f"Generating mutants for {triple_ddg_out_sort}")
    with suppress_stdout():
        out_folder = generate_triple_mutation(
            full_pdb_path,
            triple_ddg_out_sort,
            total=num_triple_mut,
            folder=mut_folder,
            text_file=f"{base_name}_mutants.txt",
        )
    if not quiet:
        logger.info(f"Generated triple mutations PDB in {out_folder}")

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    logger.info(f"Completed in {elapsed_time:.2f} minutes!")


def get_inputs() -> tuple[str, int, int, int, int, int, bool, bool]:
    """Parse CLI arguments for ``md_mutate``.

    Returns:
        ``(pdb_path, num_double_ddg, num_single_mut, num_double_mut,
        num_triple_ddg, num_triple_mut, append, quiet)``
    """
    parser = argparse.ArgumentParser(
        prog="np_mutation",
        description=None,
        epilog="Part of mutadock library. Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="PDB File for predicting and generating mutations",
        metavar="PDB",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--double",
        help="Number of Compounds for Double ddG (Default: 100)",
        metavar="D",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-t",
        "--triple",
        help="Number of Compounds for Triple ddG (Default: 40)",
        metavar="T",
        type=int,
        default=40,
    )
    parser.add_argument(
        "-S",
        "--spm",
        help="Number of Compounds for Single Point Mutations (Default: 15)",
        metavar="SM",
        type=int,
        default=15,
    )
    parser.add_argument(
        "-D",
        "--dpm",
        help="Number of Compounds for Double Point Mutations (Default: 15)",
        metavar="DM",
        type=int,
        default=15,
    )
    parser.add_argument(
        "-T",
        "--tpm",
        help="Number of Compounds for Triple Point Mutations (Default: 15)",
        metavar="TM",
        type=int,
        default=15,
    )
    parser.add_argument(
        "--noappend",
        help="Generate all files again (Default: False)",
        action="store_false",
    )
    parser.add_argument(
        "--quiet",
        help="Run the Docking in quiet mode (default: False).",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    pdb_file = args.input

    pdb_path = Path(pdb_file)
    if not pdb_path.is_absolute():
        logger.info("Assuming current directory as root since path not specified.")
        full_pdb_path = Path.cwd() / pdb_file
    else:
        full_pdb_path = pdb_path

    if not full_pdb_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_pdb_path.suffix not in (".pdb", ".cif"):
        sys.exit(
            "Given file is not a PDB/CIF file, input should be a .pdb or .cif file."
        )

    return (
        str(full_pdb_path),
        args.double,
        args.spm,
        args.dpm,
        args.triple,
        args.tpm,
        args.noappend,
        args.quiet,
    )


if __name__ == "__main__":
    np_mutation()
