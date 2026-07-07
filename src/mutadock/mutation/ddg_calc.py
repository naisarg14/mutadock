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
import statistics
import sys
from pathlib import Path
from typing import Optional

from . import predict_ddG
from .Amino import get_1
from .helpers import (
    add_ddg_protocol_args,
    backup,
    in_directory,
    read_partial_ddg,
    resolve_ddg_params,
)
from .predict_ddG import REU_TO_KCAL_SCALE

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
    full_pdb_path, full_csv_path, out_csv, params = get_inputs()
    calc_ddg(pdb_file=full_pdb_path, in_file=full_csv_path, out_file=out_csv, **params)


def calc_ddg(
    pdb_file: str,
    in_file: str,
    out_file: Optional[str] = None,
    *,
    replicates: int = 1,
    pack_radius: float = 8.0,
    backbone_minimization: bool = False,
    cartesian: bool = False,
    reu_to_kcal: float = REU_TO_KCAL_SCALE,
    resume: bool = False,
) -> str:
    """Compute single-mutation ΔΔG values for every row of *in_file*.

    ΔΔG = score(mutant) − score(wild-type self-mutation reference), where both
    sides undergo the *same* repack (and optional minimization) protocol so the
    difference is unbiased (a null WT→WT mutation scores ~0, not strongly
    "stabilizing").  Values are in Rosetta Energy Units (REU); ``ddG_kcal`` is
    ``ddG_value * reu_to_kcal`` (see ``predict_ddG.REU_TO_KCAL_SCALE``).

    Args:
        pdb_file: Wild-type PDB.
        in_file: Mutations CSV (needs ``chain``/``position``/``prAA``).
        out_file: Output CSV path (defaults to ``<in_file>_ddG.csv``).
        replicates: Independent estimates per mutation → mean ± SD (Rosetta
            packing is stochastic).
        pack_radius: Å neighbourhood repacked (and minimized) around the site.
        backbone_minimization: Minimize backbone + side chains after packing.
        cartesian: Use the Cartesian ``ref2015_cart`` protocol (implies/enables
            Cartesian minimization — the ``cartesian_ddg`` style).
        reu_to_kcal: REU → kcal/mol scaling factor for the ``ddG_kcal`` column.
        resume: When ``True`` and *out_file* already holds partial results from
            an interrupted run with identical settings, keep those rows and
            compute only the mutations still missing (per-item checkpoint).
            Resume matches on the recorded ``ddG_protocol``/``n_replicates`` but
            *assumes* the same ``pack_radius``/``reu_to_kcal`` (these are not
            stored, so a change is echoed to the log rather than guarded).

    Returns:
        Path to the written CSV.
    """
    if not out_file:
        out_file = f"{in_file.removesuffix('.csv')}_ddG.csv"

    replicates = max(1, int(replicates))
    protocol = predict_ddG.protocol_label(backbone_minimization, cartesian)

    fieldnames = [
        "sr",
        "pdb",
        "chain",
        "position",
        "wtAA",
        "prAA",
        "wtProb",
        "prProb",
        "ddG_value",
        "ddG_sd",
        "ddG_kcal",
        "n_replicates",
        "ddG_protocol",
    ]

    # Per-item checkpoint: reuse rows already written by an interrupted run with
    # the same protocol/replicates, and skip those mutations below.  Only back up
    # (and thus overwrite) the output when starting fresh.
    resume_rows = (
        read_partial_ddg(
            out_file,
            must_match={"ddG_protocol": protocol, "n_replicates": replicates},
        )
        if resume
        else None
    )
    if resume_rows is None:
        backup(out_file)
        done_keys: set = set()
    else:
        done_keys = {(r["chain"], r["position"], r["prAA"]) for r in resume_rows}
        # pack_radius/reu_to_kcal are not stored per row, so resume cannot verify
        # they match the interrupted run — echo them so a parameter change is
        # visible rather than silently mixed into the kept rows.
        logger.info(
            "Resuming single ΔΔG: %d mutation(s) already computed in '%s'. "
            "Assuming unchanged pack_radius=%s, reu_to_kcal=%s.",
            len(done_keys),
            out_file,
            pack_radius,
            reu_to_kcal,
        )

    sys.stdout = open(os.devnull, "w")
    init("-mute all")
    sys.stdout = sys.__stdout__

    with open(out_file, "w+") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        if resume_rows:
            for r in resume_rows:
                writer.writerow({k: r.get(k, "") for k in fieldnames})
            out.flush()
        with open(in_file) as csvfile:
            reader = csv.DictReader(csvfile)
            pdb_path = Path(pdb_file).resolve()
            with in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)
            sfxn = predict_ddG.get_scorefxn(cartesian)
            # Wild-type self-mutation reference, cached per site: identical for
            # all 19 substitutions at a position, so compute it once per site.
            wt_ref_cache: dict[int, float] = {}
            for row in tqdm(reader):
                if (row["chain"], row["position"], row["prAA"]) in done_keys:
                    continue
                chain = row["chain"]
                position = int(row["position"])
                mutation = row["prAA"]
                pose_position = pose.pdb_info().pdb2pose(chain, position)
                aa_str = get_1(mutation)
                assert aa_str is not None, f"Unknown amino acid code: {mutation}"

                if pose_position not in wt_ref_cache:
                    wt_ref_cache[pose_position] = statistics.fmean(
                        predict_ddG.wt_reference_score(
                            pose,
                            [pose_position],
                            sfxn,
                            pack_radius,
                            backbone_minimization=backbone_minimization,
                            cartesian=cartesian,
                        )
                        for _ in range(replicates)
                    )
                ref = wt_ref_cache[pose_position]

                ddg_values = [
                    sfxn.score(
                        predict_ddG.apply_mutations(
                            pose,
                            [(pose_position, aa_str)],
                            sfxn,
                            pack_radius,
                            backbone_minimization=backbone_minimization,
                            cartesian=cartesian,
                        )
                    )
                    - ref
                    for _ in range(replicates)
                ]
                row.update(predict_ddG.summarize_ddg(ddg_values, reu_to_kcal))
                row["ddG_protocol"] = protocol
                writer.writerow(row)
                # Flush each completed item so an interruption leaves a valid,
                # resumable checkpoint rather than a buffer-full of lost work.
                out.flush()
    return out_file


def get_inputs() -> tuple[str, str, Optional[str], dict]:
    parser = argparse.ArgumentParser(
        description="Takes input a PDB file and CSV of mutations and then calculated the ddG values for all the mutations. Input CSV can be easily generated by using the 'csv_generator.py' in the package.",
        epilog="Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument("-p", "--pdb", help="PDB File", metavar="PDB", required=True)
    parser.add_argument(
        "-i", "--input", help="Input CSV File", metavar="CSV", required=True
    )
    parser.add_argument("-o", "--output", help="Output CSV with the ddG", metavar="CSV")
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume an interrupted run: keep rows already in the output CSV "
        "(same protocol/replicates) and compute only the missing mutations. "
        "Assumes an identical --pack-radius/--reu-to-kcal (not recorded, not "
        "checked).",
    )
    add_ddg_protocol_args(parser)

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

    params = resolve_ddg_params(args)
    params["resume"] = args.resume
    return str(full_pdb_path), str(full_csv_path), out_csv, params


if __name__ == "__main__":
    main()
