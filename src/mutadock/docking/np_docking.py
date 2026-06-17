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
import shutil
import sys
import time
from contextlib import contextmanager
from itertools import product
from pathlib import Path
from typing import Any, Optional

from tqdm import tqdm

try:
    from .exceptions import (
        ConfigError,
        DockingError,
        DockingRunError,
        LigandPreparationError,
        ReceptorPreparationError,
    )
    from .vina_helper import (
        add_score_to_csv,
        backup,
        calculate_geometric_center,
        calculate_radius,
        dock_vina,
        prepare_ligand,
        prepare_receptor,
        read_config,
        run_autosite,
        vina_split,
    )
except ImportError:
    from exceptions import (  # type: ignore[no-redef]
        ConfigError,
        DockingError,
        DockingRunError,
        LigandPreparationError,
        ReceptorPreparationError,
    )
    from vina_helper import (  # type: ignore[no-redef]
        add_score_to_csv,
        backup,
        calculate_geometric_center,
        calculate_radius,
        dock_vina,
        prepare_ligand,
        prepare_receptor,
        read_config,
        run_autosite,
        vina_split,
    )

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


def np_docking() -> None:
    """Run batch receptor–ligand docking over all combinations from text files.

    Reads receptor and ligand paths from text files, prepares PDBQT inputs,
    runs AutoDock Vina for each combination, and appends scores to a CSV.
    Skips already-completed pairs unless ``--ignore_existing`` is set.
    CLI entry point for ``md_dock``.
    """
    start_time = time.time()
    receptors, ligands, config, autosite, quiet, completed_name, ignore_existing = (
        prepare_inputs()
    )

    center: Optional[list[float]] = None
    box_size: Optional[list[float]] = None
    exhaustiveness: int = 32
    n_poses: int = 20
    n_poses_write: int = 5
    overwrite: bool = True
    use_autosite_binary: bool = False

    if config is not None:
        try:
            center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite = (
                read_config(config)
            )
        except ConfigError as e:
            sys.exit(f"Error while reading the config file {config}\nError: {e}")

    if autosite is not None:
        try:
            center = list(calculate_geometric_center(autosite))
            radius = calculate_radius(autosite)
            box_dim = max(50.0, radius) * 2
            box_size = [box_dim, box_dim, box_dim]
        except Exception as e:
            sys.exit(
                f"Error while calculating the geometric center of the autosite file {autosite}\nError: {e}"
            )

    if config is None and autosite is None and center is None:
        if shutil.which("autosite") is not None:
            use_autosite_binary = True
            if not quiet:
                logger.info(
                    "autosite binary found. Will run autosite per receptor to determine the binding site."
                )
        else:
            sys.exit(
                "Error: No config file, autosite PDB, or center coordinates provided, "
                "and the autosite binary was not found on PATH.\n"
                "Please provide one of:\n"
                "  -c CONFIG     Vina configuration file with center_x/y/z\n"
                "  -a AUTOSITE   AutoSite-generated cluster PDB file\n"
                "  Or ensure the autosite binary is available on PATH."
            )

    try:
        with open(completed_name) as f:
            completed = f.readlines()
            completed = [x.replace("\n", "") for x in completed]
    except FileNotFoundError:
        completed = []

    if ignore_existing:
        completed = []

    combinations = list(product(receptors, ligands))

    combinations = [x for x in combinations if str(x) not in completed]

    if not quiet and not ignore_existing:
        logger.info(
            f"Found {len(completed)} completed receptor-ligand combinations. {len(combinations)} combinations to be docked."
        )

    # Maps prepared_receptor path -> (center, box_size) on success, None on failure.
    autosite_cache: dict[str, Optional[tuple[list[float], list[float]]]] = {}

    output_dir = None
    for combination in tqdm(combinations):
        receptor = combination[0]
        ligand = combination[1]
        prepared_receptor = f"{receptor}qt"

        if ligand.endswith(".sdf"):
            prepared_ligand = f"{ligand.removesuffix('.sdf')}.pdbqt"
        elif ligand.endswith(".mol2"):
            prepared_ligand = f"{ligand.removesuffix('.mol2')}.pdbqt"
        elif ligand.endswith(".pdbqt"):
            prepared_ligand = ligand
        else:
            logger.error(f"Unsupported ligand format for {ligand}, skipping.")
            continue

        if not quiet:
            logger.info(f"Docking the receptor {receptor} to the ligand {ligand}")
        receptor_path = Path(receptor).resolve()
        ligand_path = Path(ligand).resolve()

        output_dir = receptor_path.parent / "out"
        output_dir.mkdir(parents=True, exist_ok=True)

        out_pdb = str(output_dir / f"{receptor_path.stem}_{ligand_path.stem}_out.pdb")
        log_file = str(output_dir / f"{receptor_path.stem}_{ligand_path.stem}_log.txt")
        csv_file = str(output_dir / "docking_results.csv")

        backup(out_pdb)
        backup(log_file)

        loop_center = center
        loop_box_size = box_size

        try:
            if not quiet:
                logger.info(f"Docking for {ligand} with {receptor}")
            if not quiet:
                logger.info(
                    "Press Ctrl+D (EOFE Error) to skip this receptor-ligand combination."
                )
            if not Path(prepared_receptor).exists() or ignore_existing:
                if not quiet:
                    logger.info(f"Preparing receptor {receptor}")
                try:
                    with suppress_stdout():
                        prepare_receptor(
                            input_pdb=receptor,
                            output_pdbqt=prepared_receptor,
                        )
                except ReceptorPreparationError as e:
                    logger.error(
                        f"Error while preparing receptor {receptor}\nError: {e}\nSkipping this receptor-ligand combination."
                    )
                    continue

            if use_autosite_binary:
                if prepared_receptor in autosite_cache:
                    cached = autosite_cache[prepared_receptor]
                    if cached is None:
                        continue
                    loop_center, loop_box_size = cached
                else:
                    receptor_pdbqt_path = Path(prepared_receptor)
                    autosite_out_dir = (
                        receptor_pdbqt_path.parent
                        / f"{receptor_pdbqt_path.stem}_autosite_out"
                    )
                    cluster_pdb = (
                        autosite_out_dir / f"{receptor_pdbqt_path.stem}_cl_001.pdb"
                    )

                    if not cluster_pdb.exists() or ignore_existing:
                        if autosite_out_dir.exists():
                            shutil.rmtree(autosite_out_dir)
                        if not quiet:
                            logger.info(f"Running autosite on {prepared_receptor}")
                        try:
                            cluster_pdb = Path(run_autosite(prepared_receptor))
                        except DockingRunError as e:
                            logger.error(
                                f"autosite failed for {receptor}\nError: {e}\nSkipping this receptor-ligand combination."
                            )
                            autosite_cache[prepared_receptor] = None
                            continue
                    elif not quiet:
                        logger.info(f"Reusing existing AutoSite output: {cluster_pdb}")

                    try:
                        loop_center = list(calculate_geometric_center(str(cluster_pdb)))
                        print(f"Calculated box center: {loop_center}")
                        radius = calculate_radius(str(cluster_pdb))
                        box_dim = max(50.0, radius) * 2
                        print(f"Calculated box size: {box_dim} (radius: {radius})")
                        loop_box_size = [box_dim, box_dim, box_dim]
                        autosite_cache[prepared_receptor] = (loop_center, loop_box_size)
                    except Exception as e:
                        logger.error(
                            f"Failed to calculate binding site from {cluster_pdb}\nError: {e}\nSkipping this receptor-ligand combination."
                        )
                        autosite_cache[prepared_receptor] = None
                        continue

                    config_path = (
                        receptor_pdbqt_path.parent
                        / f"{receptor_pdbqt_path.stem}_autosite.conf"
                    )
                    with open(config_path, "w") as cfg:
                        cfg.write(f"center_x = {loop_center[0]}\n")
                        cfg.write(f"center_y = {loop_center[1]}\n")
                        cfg.write(f"center_z = {loop_center[2]}\n")
                        cfg.write(f"size_x = {loop_box_size[0]}\n")
                        cfg.write(f"size_y = {loop_box_size[1]}\n")
                        cfg.write(f"size_z = {loop_box_size[2]}\n")
                        cfg.write(f"exhaustiveness = {exhaustiveness}\n")
                        cfg.write(f"n_poses = {n_poses}\n")
                        cfg.write(f"n_poses_write = {n_poses_write}\n")
                        cfg.write(f"overwrite = {overwrite}\n")
                    if not quiet:
                        logger.info(f"Wrote AutoSite config: {config_path}")

            if not Path(prepared_ligand).exists() or ignore_existing:
                if not quiet:
                    logger.info(f"Preparing ligand {ligand}")
                try:
                    with suppress_stdout():
                        prepare_ligand(in_file=ligand, out_file=prepared_ligand)
                except LigandPreparationError as e:
                    logger.error(
                        f"Error while preparing ligand {ligand}\nError: {e}\nSkipping this receptor-ligand combination."
                    )
                    continue

            if not quiet:
                logger.info("Starting docking")
            try:
                with suppress_stdout():
                    dock_vina(
                        prepared_receptor,
                        prepared_ligand,
                        out_pdb,
                        log_file,
                        center=loop_center,
                        box_size=loop_box_size,
                        exhaustiveness=exhaustiveness,
                        n_poses=n_poses,
                        n_poses_write=n_poses_write,
                        overwrite=overwrite,
                    )
            except DockingRunError as e:
                logger.error(
                    f"Error while docking {ligand} to {receptor}\nError: {e}\nSkipping this receptor-ligand combination."
                )
                continue
            if not quiet:
                logger.info("Docking Completed, writing log file")

            if not quiet:
                logger.info("Getting the ligand 1 after docking")
            ligand_1 = out_pdb.replace(".pdbqt", "_ligand_1.sdf")
            with suppress_stdout():
                score, _ = vina_split(input_file=out_pdb, output_file=out_pdb)
            if not quiet:
                logger.info("Adding affinity to CSV")
            try:
                add_score_to_csv(out_pdb, csv_file, score)
            except DockingError as e:
                logger.error(
                    f"Error while adding affinity to CSV file {csv_file}\nError: {e}"
                )

            with open(completed_name, "a+") as file:
                file.write(f"{combination}\n")

            if not quiet:
                logger.info(
                    f"Docking completed, log file is {log_file}, ligand_1 is {ligand_1.replace('.pdbqt', '.sdf')}, docking affinity is {score}. \n"
                )

        except EOFError:
            continue

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60

    if output_dir is not None:
        logger.info(f"All Outputs are saved in the folder: {output_dir}")

    logger.info(f"Completed in {elapsed_time:.2f} minutes!")


def prepare_inputs() -> (
    tuple[list[str], list[str], Optional[str], Optional[str], bool, str, bool]
):
    """Parse CLI arguments for ``md_dock``.

    Returns:
        ``(receptors, ligands, config, autosite, quiet, completed_name,
        ignore_existing)`` where *receptors* and *ligands* are lists of
        absolute file paths read from the provided text files.
    """
    parser = argparse.ArgumentParser(
        prog="np_dock",
        description=None,
        epilog="Part of mutadock library. Written by Naisarg Patel (https://github.com/naisarg14)",
    )
    parser.add_argument(
        "-r", "--receptor_txt", help="Text File with all receptors", metavar="RECEPTOR"
    )
    parser.add_argument(
        "-l", "--ligand_txt", help="Text File with all ligands", metavar="LIGAND"
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Text File with all Vina Configuration Settings",
        metavar="CONFIG",
    )
    parser.add_argument(
        "-a",
        "--autosite",
        default=None,
        help="PDB generated by autosite for binding site of protein.",
        metavar="AUTOSITE",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Run the Docking in quiet mode (default: False).",
    )
    parser.add_argument(
        "-i",
        "--ignore_existing",
        action="store_true",
        help="Run the Docking while ingoring existing files. All dockings will be performed again. (default: False).",
    )

    args = parser.parse_args()

    receptor_txt = args.receptor_txt
    ligand_txt = args.ligand_txt
    config = args.config

    logger.info(f"Receptor file: {receptor_txt}")
    logger.info(f"Ligands file: {ligand_txt}")
    logger.info(f"Config file: {config}")
    for f in [receptor_txt, ligand_txt]:
        if f is None:
            sys.exit("Error: Required Input Files not provided.")
        if not Path(f).exists():
            logger.error(f"{f} not found or cannot be opened.")
    try:
        with open(receptor_txt) as rec:
            receptors = rec.readlines()
    except OSError as err:
        sys.exit(f"Error reading the file {receptor_txt}: ".format(receptor_txt, err))
    try:
        with open(ligand_txt) as lig:
            ligands = lig.readlines()
    except OSError as err:
        sys.exit(f"Error reading the file {ligand_txt}: ".format(ligand_txt, err))

    for i in range(len(receptors)):
        receptors[i] = receptors[i].strip()
        if not Path(receptors[i]).is_absolute():
            receptors[i] = str(Path.cwd() / receptors[i])
    for i in range(len(ligands)):
        ligands[i] = ligands[i].strip()
        if not Path(ligands[i]).is_absolute():
            ligands[i] = str(Path.cwd() / ligands[i])

    rec_path = Path(receptor_txt).resolve()
    lig_path = Path(ligand_txt).resolve()
    completed_file_name = f"{rec_path.stem}_{lig_path.stem}_completed.txt"
    completed_name = str(rec_path.parent / completed_file_name)

    return (
        receptors,
        ligands,
        config,
        args.autosite,
        args.quiet,
        completed_name,
        args.ignore_existing,
    )


naisarg = np_docking


if __name__ == "__main__":
    np_docking()
