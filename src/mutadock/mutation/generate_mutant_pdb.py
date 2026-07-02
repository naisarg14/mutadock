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

from . import helpers
from .Amino import check_3, get_1
from .exceptions import MutationError

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

try:
    sys.stdout = open(os.devnull, "w")
    from pyrosetta import *
    from pyrosetta import init, pose_from_pdb
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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_aa1(aa: str) -> str:
    """Return 1-letter code for *aa* (accepts 1- or 3-letter input)."""
    if check_3(aa):
        result = get_1(aa)
        if result is None:
            raise MutationError(f"Unknown amino acid: '{aa}'")
        return result
    return aa


def _mutation_label(wt_aa: str, chain: str, position: int, new_aa1: str) -> str:
    """Build a compact label like ASN-A386-ALA for use in file names."""
    return f"{wt_aa}-{chain}{position}-{new_aa1}"


def _read_mutations_csv(csv_file: str) -> list[dict]:
    """Read mutations from a CSV file.

    Accepted column names (case-insensitive):
        chain     — PDB chain identifier
        position  — residue sequence number
        wtAA      — wild-type amino acid (1- or 3-letter; used for naming only)
        prAA      — mutant amino acid (1- or 3-letter)

    Aliases: ``old_aa``/``oldaa``/``from`` for wtAA;
             ``new_aa``/``newaa``/``to``/``mutaa`` for prAA.
    """
    WT_ALIASES = {"wtaa", "old_aa", "oldaa", "from", "wt"}
    NEW_ALIASES = {"praa", "new_aa", "newaa", "to", "mutaa", "mut_aa"}

    mutations: list[dict] = []
    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise MutationError(f"CSV file '{csv_file}' has no header row.")

        # map canonical names to actual header strings (case-insensitive)
        lower_fields = {h.lower().strip(): h for h in reader.fieldnames}

        def _col(aliases: set[str]) -> Optional[str]:
            for alias in aliases:
                if alias in lower_fields:
                    return lower_fields[alias]
            return None

        chain_col = _col({"chain"})
        pos_col = _col({"position", "pos", "residue_number", "resnum"})
        wt_col = _col(WT_ALIASES)
        new_col = _col(NEW_ALIASES)

        if chain_col is None:
            raise MutationError(f"CSV '{csv_file}' is missing a 'chain' column.")
        if pos_col is None:
            raise MutationError(f"CSV '{csv_file}' is missing a 'position' column.")
        if new_col is None:
            raise MutationError(
                f"CSV '{csv_file}' is missing a mutant AA column (prAA / new_aa / to)."
            )

        for row in reader:
            chain = row[chain_col].strip()
            try:
                position = int(row[pos_col])
            except ValueError as e:
                raise MutationError(
                    f"Non-integer position '{row[pos_col]}' in CSV row: {row}"
                ) from e
            wt_aa = row[wt_col].strip() if wt_col else "X"
            new_aa = row[new_col].strip()
            if not chain or not new_aa:
                continue
            mutations.append(
                {"chain": chain, "position": position, "wtAA": wt_aa, "prAA": new_aa}
            )

    return mutations


def _parse_mutation_arg(mutation_str: str) -> dict:
    """Parse ``CHAIN:POSITION:WTAA:NEWAA`` or ``CHAIN:POSITION:NEWAA``.

    Examples::

        A:386:ASN:ALA   →  chain=A  position=386  wtAA=ASN  prAA=ALA
        A:386:ALA       →  chain=A  position=386  wtAA=X    prAA=ALA
    """
    parts = mutation_str.split(":")
    if len(parts) == 4:
        chain, pos_str, wt_aa, new_aa = parts
    elif len(parts) == 3:
        chain, pos_str, new_aa = parts
        wt_aa = "X"
    else:
        raise MutationError(
            f"Invalid mutation format '{mutation_str}'. "
            "Expected CHAIN:POSITION:WTAA:NEWAA or CHAIN:POSITION:NEWAA."
        )
    try:
        position = int(pos_str)
    except ValueError as e:
        raise MutationError(
            f"Non-integer position '{pos_str}' in mutation '{mutation_str}'."
        ) from e
    return {
        "chain": chain.strip(),
        "position": position,
        "wtAA": wt_aa.strip(),
        "prAA": new_aa.strip(),
    }


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------


def generate_pdb(
    pdb_file: str,
    mutations: list[dict],
    output_folder: Optional[str] = None,
    compound: bool = False,
) -> list[str]:
    """Apply mutations to a PDB structure and write the mutant PDB file(s).

    Each dict in *mutations* must contain:
        ``chain``    — PDB chain identifier (e.g. ``"A"``).
        ``position`` — residue sequence number (int).
        ``prAA``     — mutant amino acid (1- or 3-letter code).
        ``wtAA``     — wild-type amino acid (optional; used for file naming only).

    Args:
        pdb_file:       Path to the wild-type PDB file.
        mutations:      List of mutation dicts (see above).
        output_folder:  Folder for output PDB files.  Defaults to the same
                        directory as *pdb_file*.
        compound:       If ``True``, all mutations are applied sequentially to a
                        **single** pose and one combined PDB is written.  If
                        ``False`` (default), each mutation is applied independently
                        to the wild-type structure, producing one PDB per mutation.

    Returns:
        List of absolute paths to the written PDB files.

    Raises:
        MutationError: If a residue is not found in the structure or an
            amino-acid code cannot be resolved.
    """
    if not mutations:
        raise MutationError("No mutations provided.")

    pdb_path = Path(pdb_file).resolve()
    if not pdb_path.is_file():
        raise MutationError(f"PDB file not found: '{pdb_file}'")

    out_dir = Path(output_folder).resolve() if output_folder else pdb_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    init("-mute all")
    outputs: list[str] = []

    if compound:
        with helpers.in_directory(pdb_path.parent):
            pose = pose_from_pdb(pdb_path.name)

        labels: list[str] = []
        for mut in tqdm(mutations, desc="Applying mutations"):
            chain = mut["chain"]
            position = int(mut["position"])
            wt_aa = mut.get("wtAA", "X")
            new_aa1 = _resolve_aa1(mut["prAA"])

            pose_position = pose.pdb_info().pdb2pose(chain, position)
            if pose_position == 0:
                raise MutationError(
                    helpers.format_missing_residue(
                        helpers.pose_residue_map(pose), chain, position, pdb_file
                    )
                )

            mutate_residue(pose, pose_position, new_aa1)
            labels.append(_mutation_label(wt_aa, chain, position, new_aa1))
            logger.info("Applied %s", labels[-1])

        out_name = out_dir / f"{pdb_path.stem}_{'_'.join(labels)}.pdb"
        pose.dump_pdb(str(out_name))
        logger.info("Written: %s", out_name)
        outputs.append(str(out_name))

    else:
        for i, mut in enumerate(tqdm(mutations, desc="Generating PDBs"), 1):
            chain = mut["chain"]
            position = int(mut["position"])
            wt_aa = mut.get("wtAA", "X")
            new_aa1 = _resolve_aa1(mut["prAA"])

            with helpers.in_directory(pdb_path.parent):
                pose = pose_from_pdb(pdb_path.name)

            pose_position = pose.pdb_info().pdb2pose(chain, position)
            if pose_position == 0:
                raise MutationError(
                    helpers.format_missing_residue(
                        helpers.pose_residue_map(pose), chain, position, pdb_file
                    )
                )

            mutate_residue(pose, pose_position, new_aa1)

            label = _mutation_label(wt_aa, chain, position, new_aa1)
            out_name = out_dir / f"{pdb_path.stem}_{label}.pdb"
            pose.dump_pdb(str(out_name))
            logger.info("[%d/%d] Written: %s", i, len(mutations), out_name)
            outputs.append(str(out_name))

    return outputs


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for ``md_generate_pdb``."""
    parser = argparse.ArgumentParser(
        prog="md_generate_pdb",
        description=(
            "Generate mutant PDB file(s) from a single mutation or a CSV list.\n\n"
            "Single mutation:  md_generate_pdb -i wt.pdb -m A:386:ASN:ALA\n"
            "From CSV:         md_generate_pdb -i wt.pdb -c mutations.csv\n"
            "Compound mutant:  md_generate_pdb -i wt.pdb -c mutations.csv --compound"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="PDB",
        help="Wild-type PDB file.",
    )
    parser.add_argument(
        "-m",
        "--mutation",
        metavar="CHAIN:POS:WTAA:NEWAA",
        help=(
            "Single mutation in CHAIN:POSITION:WTAA:NEWAA format "
            "(e.g. A:386:ASN:ALA).  WTAA may be omitted: A:386:ALA."
        ),
    )
    parser.add_argument(
        "-c",
        "--csv",
        metavar="FILE",
        help=(
            "CSV file listing mutations.  Required columns: chain, position, prAA. "
            "Optional: wtAA (used for output naming).  Column names are case-insensitive."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        metavar="FOLDER",
        help="Folder for output PDB files (default: same directory as input PDB).",
    )
    parser.add_argument(
        "--compound",
        action="store_true",
        help=(
            "Apply all mutations from --csv to a single pose and write one combined PDB "
            "(default: one independent PDB per mutation)."
        ),
    )

    args = parser.parse_args()

    if not args.mutation and not args.csv:
        parser.error("Provide either --mutation (-m) or --csv (-c).")
    if args.mutation and args.csv:
        parser.error("--mutation and --csv are mutually exclusive.")

    if not Path(args.input).is_file():
        parser.error(f"PDB file not found: '{args.input}'")

    try:
        if args.mutation:
            mutations = [_parse_mutation_arg(args.mutation)]
        else:
            mutations = _read_mutations_csv(args.csv)
            if not mutations:
                logger.error("No mutations found in '%s'.", args.csv)
                sys.exit(1)

        outputs = generate_pdb(
            args.input,
            mutations,
            output_folder=args.output_folder,
            compound=args.compound,
        )

    except MutationError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    print(f"\nGenerated {len(outputs)} PDB file(s):")
    for path in outputs:
        print(f"  {path}")


if __name__ == "__main__":
    main()
