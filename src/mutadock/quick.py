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

"""``md_quick`` — the 30-second demo pipeline.

Give a PDB ID (or local file), a single mutation, and a ligand code, and get a
ΔΔG for the mutation plus a docked pose of the ligand against the mutant — all
wrapped up in the standard MUTADOCK report.

    md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Optional

from mutadock.docking.exceptions import (
    ConfigError,
    DockingError,
    DockingRunError,
    LigandPreparationError,
    ReceptorPreparationError,
)
from mutadock.docking.vina_helper import (
    add_score_to_csv,
    calculate_geometric_center,
    calculate_radius,
    dock_vina,
    fetch_ligand,
    prepare_ligand,
    prepare_receptor,
    read_config,
    run_autosite,
    vina_split,
)
from mutadock.mutation.Amino import get_3
from mutadock.mutation.csv_generator import get_residues
from mutadock.mutation.ddg_calc import calc_ddg
from mutadock.mutation.exceptions import MutationError
from mutadock.mutation.generate_mutant_pdb import _parse_mutation_arg, generate_pdb
from mutadock.mutation.helpers import (
    add_ddg_protocol_args,
    clean_pdb,
    convert_cif_pdb,
    fetch_pdb,
    resolve_ddg_params,
)
from mutadock.report.report import generate_report

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Extra padding (Å) added to the AutoSite cluster diameter when auto-sizing the
# docking box; matches np_docking.DEFAULT_BOX_MARGIN.
DEFAULT_BOX_MARGIN = 8.0
# A modest default exhaustiveness so the "quick" demo stays quick.
DEFAULT_EXHAUSTIVENESS = 8


@contextmanager
def suppress_stdout() -> Any:
    """Redirect stdout to /dev/null for the duration (quiets PyRosetta/meeko)."""
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def _resolve_input_pdb(
    pdb_id: Optional[str], input_path: Optional[str], out_dir: Path, quiet: bool
) -> Path:
    """Fetch by PDB ID or use a local file; convert CIF → PDB if needed."""
    if pdb_id:
        try:
            raw = Path(fetch_pdb(pdb_id, dest_dir=out_dir))
        except MutationError as e:
            sys.exit(str(e))
    else:
        raw = Path(input_path)  # type: ignore[arg-type]
        if not raw.is_absolute():
            raw = Path.cwd() / raw
        if not raw.is_file():
            sys.exit(f"Input file not found: {raw}")
        if raw.suffix.lower() not in (".pdb", ".cif"):
            sys.exit("Input must be a .pdb or .cif file.")

    if raw.suffix.lower() == ".cif":
        if not quiet:
            logger.info("Converting CIF %s to PDB ...", raw.name)
        converted = out_dir / f"{raw.stem}.pdb"
        with suppress_stdout():
            convert_cif_pdb(str(raw), str(converted))
        raw = converted
    return raw


def _validate_and_fill_mutation(cleaned_pdb: str, mut: dict) -> dict:
    """Check the residue exists, normalize AA codes, and fill wtAA from the PDB.

    Raises:
        MutationError: If the chain/position is absent or the mutant AA is not a
            standard amino acid.
    """
    residues = get_residues(cleaned_pdb)
    by_site = {(c, p): name for c, p, name in residues.values()}

    chain, pos = mut["chain"], mut["position"]
    if (chain, pos) not in by_site:
        chains = sorted({c for c, _ in by_site})
        positions = sorted(p for c, p in by_site if c == chain)
        if positions:
            hint = f"chain {chain} spans residues {positions[0]}–{positions[-1]}"
        else:
            hint = f"available chains: {', '.join(chains) or '(none)'}"
        raise MutationError(
            f"Residue {chain}{pos} not found in the structure ({hint})."
        )

    pr3 = get_3(mut["prAA"])
    if pr3 is None:
        raise MutationError(
            f"Unknown mutant amino acid '{mut['prAA']}' — expected a standard "
            "1- or 3-letter code (e.g. HIS or H)."
        )
    mut["prAA"] = pr3

    # Fill wtAA from the actual structure so report labels are correct.
    wt_actual = by_site[(chain, pos)]
    if mut.get("wtAA") in (None, "", "X"):
        mut["wtAA"] = wt_actual
    return mut


def _compute_ddg(
    cleaned_pdb: str, mut: dict, base_name: str, quiet: bool, **ddg_params: Any
) -> tuple[float, str]:
    """Write a one-row mutations CSV, run calc_ddg, and return (ddG, csv_path)."""
    mut_csv = f"{base_name}_mutation.csv"
    with open(mut_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["sr", "pdb", "chain", "position", "wtAA", "prAA", "wtProb", "prProb"]
        )
        writer.writerow(
            [
                1,
                cleaned_pdb.removesuffix(".pdb"),
                mut["chain"],
                mut["position"],
                mut["wtAA"],
                mut["prAA"],
                0,
                0,
            ]
        )

    # discover_run prefers *_ddG_sorted.csv; with one row sorted == unsorted.
    ddg_csv = f"{base_name}_ddG_sorted.csv"
    with suppress_stdout():
        calc_ddg(cleaned_pdb, mut_csv, out_file=ddg_csv, **ddg_params)

    with open(ddg_csv) as f:
        row = next(csv.DictReader(f))
        ddg = float(row["ddG_value"])
    return ddg, ddg_csv


def _dock(
    mutant_pdb: str,
    ligand_file: str,
    out_dir: Path,
    config: Optional[str],
    exhaustiveness: int,
    quiet: bool,
) -> tuple[float, str]:
    """Dock *ligand_file* into *mutant_pdb*; return (affinity, pose_sdf_path).

    The box comes from *config* when given, otherwise from AutoSite run on the
    prepared mutant receptor.
    """
    mutant_path = Path(mutant_pdb)
    prepared_receptor = str(mutant_path.with_suffix(".pdbqt"))
    if not quiet:
        logger.info("Preparing receptor %s ...", mutant_path.name)
    prepare_receptor(input_pdb=mutant_pdb, output_pdbqt=prepared_receptor)

    n_poses, n_poses_write, overwrite = 20, 5, True
    if config is not None:
        center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite = (
            read_config(config)
        )
    else:
        if not quiet:
            logger.info("Running AutoSite to locate the binding pocket ...")
        cluster_pdb = run_autosite(prepared_receptor)
        center = list(calculate_geometric_center(cluster_pdb))
        dim = calculate_radius(cluster_pdb) * 2 + DEFAULT_BOX_MARGIN
        box_size = [dim, dim, dim]
        if not quiet:
            logger.info(
                "Box center %s, size %.1f Å", [round(c, 1) for c in center], dim
            )

    lig_path = Path(ligand_file)
    prepared_ligand = str(out_dir / f"{lig_path.stem}.pdbqt")
    if not quiet:
        logger.info("Preparing ligand %s ...", lig_path.name)
    with suppress_stdout():
        prepare_ligand(in_file=ligand_file, out_file=prepared_ligand)

    stem = f"{mutant_path.stem}_{lig_path.stem}_out"
    out_pdbqt = str(out_dir / f"{stem}.pdbqt")
    out_sdf = str(out_dir / f"{stem}.sdf")
    log_file = str(out_dir / f"{mutant_path.stem}_{lig_path.stem}_log.txt")

    if not quiet:
        logger.info("Docking (exhaustiveness=%d) ...", exhaustiveness)
    with suppress_stdout():
        dock_vina(
            prepared_receptor,
            prepared_ligand,
            out_pdbqt,
            log_file,
            center=center,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
            n_poses=n_poses,
            n_poses_write=n_poses_write,
            overwrite=overwrite,
        )
        score, _ = vina_split(input_file=out_pdbqt, output_file=out_sdf)

    add_score_to_csv(out_sdf, str(out_dir / "docking_results.csv"), score)
    return score, out_sdf


def run_quick(
    *,
    pdb_id: Optional[str],
    input_path: Optional[str],
    mutation: str,
    ligand_code: Optional[str],
    ligand_file: Optional[str],
    out_dir: Path,
    config: Optional[str] = None,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    make_report: bool = True,
    quiet: bool = False,
    replicates: int = 1,
    pack_radius: float = 8.0,
    backbone_minimization: bool = False,
    cartesian: bool = False,
    reu_to_kcal: Optional[float] = None,
) -> dict:
    """Run the full fetch → mutate → ΔΔG → fetch-ligand → dock → report pipeline.

    Returns a dict with keys ``ddg``, ``affinity``, ``mutant_pdb``, ``pose_sdf``,
    and (when a report is generated) ``report``.
    """
    start = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Structure ----------------------------------------------------------
    raw_pdb = _resolve_input_pdb(pdb_id, input_path, out_dir, quiet)
    stem = raw_pdb.stem
    base_name = str(out_dir / stem)
    cleaned_pdb = f"{base_name}_clean.pdb"
    if not quiet:
        logger.info("Cleaning %s ...", raw_pdb.name)
    with suppress_stdout():
        clean_pdb(str(raw_pdb), cleaned_pdb)

    # 2. Mutation + ΔΔG -----------------------------------------------------
    try:
        mut = _parse_mutation_arg(mutation)
        mut = _validate_and_fill_mutation(cleaned_pdb, mut)
    except MutationError as e:
        sys.exit(str(e))
    label = f"{mut['wtAA']}-{mut['chain']}{mut['position']}-{mut['prAA']}"
    if not quiet:
        logger.info("Computing ΔΔG for %s ...", label)
    ddg_params: dict = {
        "replicates": replicates,
        "pack_radius": pack_radius,
        "backbone_minimization": backbone_minimization,
        "cartesian": cartesian,
    }
    if reu_to_kcal is not None:
        ddg_params["reu_to_kcal"] = reu_to_kcal
    ddg, _ = _compute_ddg(cleaned_pdb, mut, base_name, quiet, **ddg_params)
    if not quiet:
        logger.info("ΔΔG(%s) = %.2f REU", label, ddg)

    # 3. Mutant structure ---------------------------------------------------
    mut_folder = str(out_dir / f"mutation_{stem}")
    if not quiet:
        logger.info("Building mutant structure ...")
    with suppress_stdout():
        mutant_pdbs = generate_pdb(cleaned_pdb, [mut], output_folder=mut_folder)
    mutant_pdb = mutant_pdbs[0]

    # 4. Ligand -------------------------------------------------------------
    if ligand_file:
        lig = ligand_file
        if not Path(lig).is_file():
            sys.exit(f"Ligand file not found: {lig}")
    else:
        try:
            lig = fetch_ligand(ligand_code, dest_dir=out_dir)  # type: ignore[arg-type]
        except LigandPreparationError as e:
            sys.exit(str(e))

    # 5. Dock ---------------------------------------------------------------
    affinity: Optional[float] = None
    pose_sdf: Optional[str] = None
    try:
        affinity, pose_sdf = _dock(
            mutant_pdb, lig, out_dir, config, exhaustiveness, quiet
        )
        if not quiet:
            logger.info("Docking affinity = %.3f kcal/mol", affinity)
    except (
        ReceptorPreparationError,
        LigandPreparationError,
        DockingRunError,
        DockingError,
        ConfigError,
    ) as e:
        logger.error("Docking failed: %s", e)
        logger.error("Continuing — the ΔΔG result and report are still produced.")

    # 6. Report -------------------------------------------------------------
    result: dict = {
        "ddg": ddg,
        "affinity": affinity,
        "mutant_pdb": mutant_pdb,
        "pose_sdf": pose_sdf,
        "label": label,
    }
    if make_report:
        try:
            written = generate_report(str(out_dir), base_name=stem, quiet=quiet)
            if "html" in written:
                result["report"] = written["html"]
        except Exception as e:  # pragma: no cover - report must never break the run
            logger.warning("Report generation failed: %s", e)

    elapsed = (time.time() - start) / 60
    logger.info("Completed in %.2f minutes!", elapsed)
    return result


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="md_quick",
        description="Fetch a structure, apply one mutation, compute its ΔΔG, "
        "fetch a ligand, dock it against the mutant, and write a report.",
        epilog="Part of mutadock library. Written by Naisarg Patel "
        "(https://github.com/naisarg14)",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--pdb-id",
        dest="pdb_id",
        metavar="ID",
        help="4-character RCSB PDB ID (e.g. 4QJR).",
    )
    source.add_argument(
        "-i", "--input", metavar="PDB", help="Local .pdb/.cif structure file."
    )
    parser.add_argument(
        "-m",
        "--mutation",
        required=True,
        metavar="SPEC",
        help="Mutation as CHAIN:POSITION:WTAA:NEWAA or CHAIN:POSITION:NEWAA "
        "(e.g. A:386:ASN:HIS or A:386:H).",
    )
    lig = parser.add_mutually_exclusive_group(required=True)
    lig.add_argument(
        "--ligand-code",
        dest="ligand_code",
        metavar="CODE",
        help="PubChem CID or name (e.g. 2244, imatinib, cid:5291, name:aspirin).",
    )
    lig.add_argument(
        "--ligand-file",
        dest="ligand_file",
        metavar="FILE",
        help="Local ligand file (.sdf/.mol2/.pdbqt) instead of fetching.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        metavar="DIR",
        default=None,
        help="Output directory (default: ./mdquick_<id>).",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        metavar="CONFIG",
        help="Vina config file with center/box (default: AutoSite auto-box).",
    )
    parser.add_argument(
        "-e",
        "--exhaustiveness",
        type=int,
        default=DEFAULT_EXHAUSTIVENESS,
        metavar="N",
        help=f"Vina exhaustiveness (default: {DEFAULT_EXHAUSTIVENESS}).",
    )
    parser.add_argument(
        "--no-report", dest="no_report", action="store_true", help="Skip report."
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress informational logging."
    )
    add_ddg_protocol_args(parser)
    return parser.parse_args(argv)


def md_quick(argv: Optional[list[str]] = None) -> None:
    """CLI entry point for ``md_quick``."""
    args = _parse_args(argv)

    if args.output_dir:
        out_dir = Path(args.output_dir)
    elif args.pdb_id:
        out_dir = Path.cwd() / f"mdquick_{args.pdb_id.strip().upper()}"
    else:
        out_dir = Path.cwd() / f"mdquick_{Path(args.input).stem}"

    result = run_quick(
        pdb_id=args.pdb_id,
        input_path=args.input,
        mutation=args.mutation,
        ligand_code=args.ligand_code,
        ligand_file=args.ligand_file,
        out_dir=out_dir,
        config=args.config,
        exhaustiveness=args.exhaustiveness,
        make_report=not args.no_report,
        quiet=args.quiet,
        **resolve_ddg_params(args),
    )

    affinity = result.get("affinity")
    affinity_str = (
        f"{affinity:.3f} kcal/mol" if affinity is not None else "(docking failed)"
    )
    logger.info("=" * 60)
    logger.info("md_quick summary")
    logger.info("  Mutation : %s", result["label"])
    logger.info("  ΔΔG      : %.2f REU (negative = stabilizing)", result["ddg"])
    logger.info("  Affinity : %s", affinity_str)
    if result.get("pose_sdf"):
        logger.info("  Pose     : %s", result["pose_sdf"])
    if result.get("report"):
        logger.info("  Report   : %s", result["report"])
    logger.info("=" * 60)


if __name__ == "__main__":
    md_quick()
