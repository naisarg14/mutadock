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

def vina_dock(receptor, ligand, output, center=[0, 0, 0], box_size=[30, 30, 30], exhaustiveness=32, n_poses=20, n_poses_write=5, overwrite=True):
    try:
        from vina import Vina
    except ModuleNotFoundError:
        msg = "Error with importing modules for docking using Vina.\n"
        msg += "Easiest way to fix this is to install vina using the following command:\n\n"
        msg += "python -m pip install vina\n"
        msg += "If you already have vina installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        print(msg)
        sys.exit(2)
    try:
        print("Docking done using mutadock library developed by Naisarg Patel (Github:@naisarg14)")
        v = Vina(sf_name='vina')

        v.set_receptor(receptor)
        v.set_ligand_from_file(ligand)
        print(v)
        v.compute_vina_maps(center=center, box_size=box_size)

        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(output, n_poses=n_poses_write, overwrite=overwrite)

    except Exception as e:
        return (False, e)
    return (True, "")


def main():
    parser = argparse.ArgumentParser(description="Run docking using AutoDock vina Python bindings.")

    parser.add_argument("--receptor", type=str, help="Path to the receptor file.")
    parser.add_argument("--ligand", type=str, help="Path to the ligand file.")
    parser.add_argument("--output", type=str, help="Path for saving the output file.")
    parser.add_argument("--log_file", type=str, help="Path for saving the log file.")

    parser.add_argument("--center", nargs=3, type=float, default=[0, 0, 0], help="X-dimension of the center of search box (default: [0, 0, 0]).")

    parser.add_argument("--box_size", nargs=3, type=float, default=[30, 30, 30], help="Size of the search box (default: [30, 30, 30]).")
    parser.add_argument("--exhaustiveness", type=int, default=32, help="Exhaustiveness of the search (default: 32).")
    parser.add_argument("--n_poses", type=int, default=20, help="Number of poses to generate (default: 20).")
    parser.add_argument("--n_poses_write", type=int, default=5, help="Number of poses to write to the output (default: 5).")
    parser.add_argument("--nooverwrite", action="store_false", default=True, help="Do not overwrite existing files.")

    args = parser.parse_args()

    vina_dock(
        receptor=args.receptor,
        ligand=args.ligand,
        output=args.output,
        center=args.center,
        box_size=args.box_size,
        exhaustiveness=args.exhaustiveness,
        n_poses=args.n_poses,
        n_poses_write=args.n_poses_write,
        overwrite=args.nooverwrite,
        )


if __name__ == "__main__":
    main()