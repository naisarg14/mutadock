################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2024 Naisarg Patel (https://github.com/naisarg14)        #
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

import sys, argparse
import logging
from typing import List, Tuple
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def vina_dock(
    receptor: str, 
    ligand: str, 
    output: str, 
    center: List[float] = [0, 0, 0], 
    box_size: List[int] = [30, 30, 30], 
    exhaustiveness: int = 32, 
    n_poses: int = 20, 
    n_poses_write: int = 5, 
    overwrite: bool = True
) -> Tuple[bool, str]:
    """
    Perform molecular docking using AutoDock Vina.
    
    Returns:
        Tuple of (success: bool, message: str)
    """
    try:
        from vina import Vina
    except ModuleNotFoundError:
        msg = "Error importing Vina module.\n"
        msg += "Install via: python -m pip install vina\n"
        msg += "See: https://autodock-vina.readthedocs.io/en/latest/installation.html\n"
        msg += "Contact: naisarg.patel14@hotmail.com"
        logger.error(msg)
        return (False, msg)
    
    # Validate inputs
    if not Path(receptor).exists():
        return (False, f"Receptor file not found: {receptor}")
    if not Path(ligand).exists():
        return (False, f"Ligand file not found: {ligand}")
    
    try:
        logger.info("Docking using mutadock library by Naisarg Patel (@naisarg14)")
        v = Vina(sf_name='vina')
        
        v.set_receptor(receptor)
        v.set_ligand_from_file(ligand)
        logger.info(str(v))
        
        v.compute_vina_maps(center=center, box_size=box_size)
        v.optimize()
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(output, n_poses=n_poses_write, overwrite=overwrite)
        
        logger.info(f"Docking completed successfully. Output: {output}")
        return (True, "Docking successful")
        
    except Exception as e:
        error_msg = f"Docking failed: {str(e)}"
        logger.error(error_msg)
        return (False, error_msg)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run molecular docking using AutoDock Vina Python bindings.")
    
    # Required arguments
    parser.add_argument("--receptor", type=str, required=True, help="Path to the receptor file (PDBQT format)")
    parser.add_argument("--ligand", type=str, required=True, help="Path to the ligand file (PDBQT format)")
    parser.add_argument("--output", type=str, required=True, help="Path for output file")
    
    # Optional arguments
    parser.add_argument("--center", nargs=3, type=float, default=[0, 0, 0], help="Center coordinates of search box (default: 0 0 0)")
    parser.add_argument("--box_size", nargs=3, type=float, default=[30, 30, 30], help="Size of search box in Angstroms (default: 30 30 30)")
    parser.add_argument("--exhaustiveness", type=int, default=32, help="Exhaustiveness of search (default: 32)")
    parser.add_argument("--n_poses", type=int, default=20, help="Number of poses to generate (default: 20)")
    parser.add_argument("--n_poses_write", type=int, default=5, help="Number of poses to write (default: 5)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output file")
    
    args = parser.parse_args()
    
    # Run docking
    success, message = vina_dock(
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
    
    if not success:
        logger.error(f"Docking failed: {message}")
        sys.exit(1)
    else:
        logger.info("Docking completed successfully!")
        sys.exit(0)


if __name__ == "__main__":
    main()