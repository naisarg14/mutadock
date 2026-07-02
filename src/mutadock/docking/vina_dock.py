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
import sys
from pathlib import Path
from typing import Optional

try:
    from .exceptions import DockingRunError
except ImportError:
    from exceptions import DockingRunError  # type: ignore[no-redef]

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def vina_dock(
    receptor: str,
    ligand: str,
    output: str,
    center: Optional[list[float]] = None,
    box_size: Optional[list[int]] = None,
    exhaustiveness: int = 32,
    n_poses: int = 20,
    n_poses_write: int = 5,
    overwrite: bool = True,
) -> None:
    """Perform molecular docking using the AutoDock Vina Python bindings.

    Sets up a Vina scoring function, loads the receptor and ligand, computes
    the search map over the specified box, runs local optimisation followed by
    global docking, and writes the requested number of poses to *output*.

    Args:
        receptor: Path to the prepared receptor PDBQT file.
        ligand: Path to the prepared ligand PDBQT file.
        output: Destination path for the docking poses (PDBQT format).
        center: Search box center ``[x, y, z]`` in Ångströms (default
            ``[0, 0, 0]``).
        box_size: Search box dimensions ``[x, y, z]`` in Ångströms (default
            ``[30, 30, 30]``).
        exhaustiveness: Exhaustiveness of the global search (default 32).
            Higher values are slower but more thorough.
        n_poses: Number of binding poses to generate internally (default 20).
        n_poses_write: Number of top poses to write to *output* (default 5).
        overwrite: If ``True``, overwrite an existing *output* file.

    Raises:
        DockingRunError: If inputs are missing or the Vina run fails.
    """
    try:
        from vina import Vina
    except ModuleNotFoundError:
        msg = "Error importing Vina module.\n"
        msg += "Install via: python -m pip install vina\n"
        msg += "See: https://autodock-vina.readthedocs.io/en/latest/installation.html\n"
        msg += "Contact: naisarg.patel14@hotmail.com"
        logger.error(msg)
        sys.exit(2)

    if center is None:
        center = [0, 0, 0]
    if box_size is None:
        box_size = [30, 30, 30]

    if not Path(receptor).exists():
        raise DockingRunError(f"Receptor file not found: {receptor}")
    if not Path(ligand).exists():
        raise DockingRunError(f"Ligand file not found: {ligand}")

    try:
        logger.info("Docking using mutadock library by Naisarg Patel (@naisarg14)")
        v = Vina(sf_name="vina")

        v.set_receptor(receptor)
        v.set_ligand_from_file(ligand)
        logger.info(str(v))

        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(output, n_poses=n_poses_write, overwrite=overwrite)

        logger.info(f"Docking completed successfully. Output: {output}")

    except DockingRunError:
        raise
    except Exception as e:
        raise DockingRunError(f"Docking failed: {e}") from e


def main() -> None:
    """CLI entry point for ``md_vina_dock``."""
    parser = argparse.ArgumentParser(
        description="Run molecular docking using AutoDock Vina Python bindings."
    )

    parser.add_argument(
        "--receptor",
        type=str,
        required=True,
        help="Path to the receptor file (PDBQT format)",
    )
    parser.add_argument(
        "--ligand",
        type=str,
        required=True,
        help="Path to the ligand file (PDBQT format)",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path for output file"
    )
    parser.add_argument(
        "--center",
        nargs=3,
        type=float,
        default=[0, 0, 0],
        help="Center coordinates of search box (default: 0 0 0)",
    )
    parser.add_argument(
        "--box_size",
        nargs=3,
        type=float,
        default=[30, 30, 30],
        help="Size of search box in Angstroms (default: 30 30 30)",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=32,
        help="Exhaustiveness of search (default: 32)",
    )
    parser.add_argument(
        "--n_poses",
        type=int,
        default=20,
        help="Number of poses to generate (default: 20)",
    )
    parser.add_argument(
        "--n_poses_write",
        type=int,
        default=5,
        help="Number of poses to write (default: 5)",
    )
    parser.add_argument(
        "--overwrite",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Overwrite existing output file (use --no-overwrite to disable)",
    )

    args = parser.parse_args()

    try:
        vina_dock(
            receptor=args.receptor,
            ligand=args.ligand,
            output=args.output,
            center=args.center,
            box_size=args.box_size,
            exhaustiveness=args.exhaustiveness,
            n_poses=args.n_poses,
            n_poses_write=args.n_poses_write,
            overwrite=args.overwrite,
        )
    except DockingRunError as e:
        logger.error(str(e))
        sys.exit(1)
    else:
        logger.info("Docking completed successfully!")
        sys.exit(0)


if __name__ == "__main__":
    main()
