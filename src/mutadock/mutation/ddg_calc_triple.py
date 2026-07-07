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
from itertools import combinations
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
    full_pdb_path, full_csv_path, total, out_csv, params = get_inputs()
    calc_triple_ddg(
        pdb_file=full_pdb_path,
        double_csv=full_csv_path,
        total=total,
        out_csv=out_csv,
        **params,
    )


def calc_triple_ddg(
    pdb_file: str,
    double_csv: str,
    out_csv: Optional[str] = None,
    total: int = 0,
    *,
    replicates: int = 1,
    pack_radius: float = 8.0,
    backbone_minimization: bool = False,
    cartesian: bool = False,
    reu_to_kcal: float = REU_TO_KCAL_SCALE,
    resume: bool = False,
) -> str:
    """Compute triple-mutation ΔΔG values for combinations of *double_csv*.

    ΔΔG = score(triple mutant) − score(wild-type self-mutation reference), where
    both sides undergo the *same* repack (and optional minimization) protocol so
    the difference is unbiased (a null WT→WT triple mutation scores ~0, not
    strongly "stabilizing").  The reference depends only on the *triple of
    positions*, so it is cached per position-tuple.  Values are in Rosetta
    Energy Units (REU); ``triple_ddG_kcal`` is ``triple_ddG_value * reu_to_kcal``
    (see ``predict_ddG.REU_TO_KCAL_SCALE``).

    Args:
        pdb_file: Wild-type PDB.
        double_csv: Double-mutation CSV whose single-mutation rows are combined
            via ``combinations(..., 3)`` to form triple mutants.
        out_csv: Output CSV path (defaults to ``<pdb_file>_triple_ddg.csv``).
        total: Number of unique single-mutation rows to combine (0 = all).
        replicates: Independent estimates per combination → mean ± SD (Rosetta
            packing is stochastic).
        pack_radius: Å neighbourhood repacked (and minimized) around each site.
        backbone_minimization: Minimize backbone + side chains after packing.
        cartesian: Use the Cartesian ``ref2015_cart`` protocol (implies/enables
            Cartesian minimization — the ``cartesian_ddg`` style).
        reu_to_kcal: REU → kcal/mol scaling factor for the ``triple_ddG_kcal``
            column.
        resume: When ``True`` and *out_csv* already holds partial results from an
            interrupted run with identical settings, keep those rows and compute
            only the combinations still missing (per-item checkpoint).  Matches
            on the recorded ``ddG_protocol``/``n_replicates`` but *assumes* the
            same ``pack_radius``/``reu_to_kcal`` (not stored, so a change is
            echoed to the log rather than guarded).

    Returns:
        Path to the written CSV.
    """
    if not out_csv:
        out_csv = f"{pdb_file.removesuffix('.pdb')}_triple_ddg.csv"

    result = get_mut_csv(double_csv, total)

    replicates = max(1, int(replicates))
    protocol = predict_ddG.protocol_label(backbone_minimization, cartesian)

    fieldnames = [
        "sr",
        "pdb",
        "combination",
        "mut1_wtAA",
        "chain1",
        "mut1_position",
        "mut1_prAA",
        "mut2_wtAA",
        "chain2",
        "mut2_position",
        "mut2_prAA",
        "mut3_wtAA",
        "chain3",
        "mut3_position",
        "mut3_prAA",
        "triple_ddG_value",
        "triple_ddG_sd",
        "triple_ddG_kcal",
        "n_replicates",
        "ddG_protocol",
    ]

    # Per-item checkpoint: reuse combinations already written by an interrupted
    # run with the same protocol/replicates, keyed by the triple's mutation
    # identities.  Only back up (overwrite) the output when starting fresh.
    resume_rows = (
        read_partial_ddg(
            out_csv,
            must_match={"ddG_protocol": protocol, "n_replicates": replicates},
        )
        if resume
        else None
    )
    if resume_rows is None:
        backup(out_csv)
        done_keys: set = set()
    else:
        done_keys = {
            (
                r["chain1"],
                str(r["mut1_position"]),
                r["mut1_prAA"],
                r["chain2"],
                str(r["mut2_position"]),
                r["mut2_prAA"],
                r["chain3"],
                str(r["mut3_position"]),
                r["mut3_prAA"],
            )
            for r in resume_rows
        }
        # pack_radius/reu_to_kcal are not stored per row, so resume cannot verify
        # they match the interrupted run — echo them so a parameter change is
        # visible rather than silently mixed into the kept rows.
        logger.info(
            "Resuming triple ΔΔG: %d combination(s) already computed in '%s'. "
            "Assuming unchanged pack_radius=%s, reu_to_kcal=%s.",
            len(done_keys),
            out_csv,
            pack_radius,
            reu_to_kcal,
        )

    sys.stdout = open(os.devnull, "w")
    init("-mute all")
    sys.stdout = sys.__stdout__

    with open(out_csv, "w+") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        if resume_rows:
            for r in resume_rows:
                writer.writerow({k: r.get(k, "") for k in fieldnames})
            out.flush()
        count = 1
        pdb_path = Path(pdb_file).resolve()
        with in_directory(pdb_path.parent):
            pose = pose_from_pdb(pdb_path.name)
        sfxn = predict_ddG.get_scorefxn(cartesian)
        # Wild-type self-mutation reference, cached per position triple:
        # identical for all substitutions at a triple of sites, so compute it
        # once per triple.
        wt_ref_cache: dict[tuple[int, int, int], float] = {}
        for x in tqdm(combinations(result, 3)):
            # Deterministic combinations order → count stays aligned whether or
            # not this item was already done, so advance it before skipping.
            key = (
                x[0]["chain"],
                str(x[0]["position"]),
                x[0]["prAA"],
                x[1]["chain"],
                str(x[1]["position"]),
                x[1]["prAA"],
                x[2]["chain"],
                str(x[2]["position"]),
                x[2]["prAA"],
            )
            if key in done_keys:
                count += 1
                continue

            mut1 = x[0]
            mut2 = x[1]
            mut3 = x[2]

            chain1 = mut1["chain"]
            position1 = int(mut1["position"])
            old1 = mut1["wtAA"]
            mutation1 = mut1["prAA"]
            pose_position1 = pose.pdb_info().pdb2pose(chain1, position1)

            chain2 = mut2["chain"]
            position2 = int(mut2["position"])
            old2 = mut2["wtAA"]
            mutation2 = mut2["prAA"]
            pose_position2 = pose.pdb_info().pdb2pose(chain2, position2)

            chain3 = mut3["chain"]
            position3 = int(mut3["position"])
            old3 = mut3["wtAA"]
            mutation3 = mut3["prAA"]
            pose_position3 = pose.pdb_info().pdb2pose(chain3, position3)

            name = f"{count}_{old1}_{position1}_{mutation1}+{old2}_{position2}_{mutation2}+{old3}_{position3}_{mutation3}"

            aa1 = get_1(mutation1)
            aa2 = get_1(mutation2)
            aa3 = get_1(mutation3)
            assert aa1 is not None, f"Unknown amino acid code: {mutation1}"
            assert aa2 is not None, f"Unknown amino acid code: {mutation2}"
            assert aa3 is not None, f"Unknown amino acid code: {mutation3}"

            triple = (pose_position1, pose_position2, pose_position3)
            if triple not in wt_ref_cache:
                wt_ref_cache[triple] = statistics.fmean(
                    predict_ddG.wt_reference_score(
                        pose,
                        [pose_position1, pose_position2, pose_position3],
                        sfxn,
                        pack_radius,
                        backbone_minimization=backbone_minimization,
                        cartesian=cartesian,
                    )
                    for _ in range(replicates)
                )
            ref = wt_ref_cache[triple]

            ddg_values = [
                sfxn.score(
                    predict_ddG.apply_mutations(
                        pose,
                        [
                            (pose_position1, aa1),
                            (pose_position2, aa2),
                            (pose_position3, aa3),
                        ],
                        sfxn,
                        pack_radius,
                        backbone_minimization=backbone_minimization,
                        cartesian=cartesian,
                    )
                )
                - ref
                for _ in range(replicates)
            ]
            summary = predict_ddG.summarize_ddg(ddg_values, reu_to_kcal)
            row = {
                "sr": count,
                "pdb": pdb_file,
                "combination": name,
                "mut1_wtAA": old1,
                "chain1": chain1,
                "mut1_position": position1,
                "mut1_prAA": mutation1,
                "mut2_wtAA": old2,
                "chain2": chain2,
                "mut2_position": position2,
                "mut2_prAA": mutation2,
                "mut3_wtAA": old3,
                "chain3": chain3,
                "mut3_position": position3,
                "mut3_prAA": mutation3,
                "triple_ddG_value": summary["ddG_value"],
                "triple_ddG_sd": summary["ddG_sd"],
                "triple_ddG_kcal": summary["ddG_kcal"],
                "n_replicates": summary["n_replicates"],
                "ddG_protocol": protocol,
            }
            writer.writerow(row)
            # Flush each completed combination so an interruption leaves a valid,
            # resumable checkpoint rather than a buffer-full of lost work.
            out.flush()
            count += 1
        return out_csv


def get_mut_csv(file: str, total: int = -1) -> list[dict[str, str]]:
    with open(file) as f:
        reader = csv.DictReader(f)
        results = []
        count = 0
        if total == -1:
            total = 100000
        for row in reader:
            mutation = {
                "chain": row["chain1"],
                "position": row["mut1_position"],
                "wtAA": row["mut1_wtAA"],
                "prAA": row["mut1_prAA"],
            }
            if mutation not in results:
                results.append(mutation.copy())
                count += 1

            mutation = {
                "chain": row["chain2"],
                "position": row["mut2_position"],
                "wtAA": row["mut2_wtAA"],
                "prAA": row["mut2_prAA"],
            }
            if mutation not in results:
                results.append(mutation.copy())
                count += 1

            if count >= total:
                break

    return results


def get_inputs() -> tuple[str, str, int, Optional[str], dict]:
    parser = argparse.ArgumentParser(
        description="Takes input a PDB file and CSV of mutations and then calculated the Triple ddG values for all the mutations.",
        epilog="Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument("-p", "--pdb", help="PDB File", metavar="PDB", required=True)
    parser.add_argument(
        "-i", "--input", help="Input Double CSV File", metavar="CSV", required=True
    )
    parser.add_argument(
        "-n",
        "--num",
        help="Number of top Double ddG values to take",
        metavar="N",
        type=int,
        default=50,
    )
    parser.add_argument("-o", "--output", help="Output CSV with the ddG", metavar="CSV")
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume an interrupted run: keep combinations already in the output "
        "CSV (same protocol/replicates) and compute only the missing ones. "
        "Assumes an identical --pack-radius/--reu-to-kcal (not recorded, not "
        "checked).",
    )
    add_ddg_protocol_args(parser)

    args = parser.parse_args()
    pdb_file = args.pdb
    double_csv = args.input
    out_csv = args.output
    total = args.num

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
        sys.exit("Given file is not a PDB file, input should be a PDB file.")

    p = Path(double_csv)
    if not p.is_absolute():
        logger.info(
            f"Assuming current directory for {double_csv} as root since path not specified."
        )
        full_csv_path = Path.cwd() / double_csv
    else:
        full_csv_path = p.resolve()

    if not full_csv_path.is_file():
        sys.exit(
            "No such file found in current directory, enter full path for other directories."
        )

    if full_csv_path.suffix != ".csv":
        logger.error("Given file is not a CSV file, input should be a CSV file.")

    if out_csv and not out_csv.endswith(".csv"):
        out_csv += ".csv"

    params = resolve_ddg_params(args)
    params["resume"] = args.resume
    return (
        str(full_pdb_path),
        str(full_csv_path),
        total,
        out_csv,
        params,
    )


if __name__ == "__main__":
    main()
