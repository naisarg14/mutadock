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

"""MUTADOCK — mutation studies and receptor-ligand docking.

Top-level package exposing two subpackages:

* :mod:`mutadock.mutation` — amino-acid substitution analysis (mutation CSVs,
  ddG calculations, mutant PDB generation via PyRosetta).
* :mod:`mutadock.docking` — AutoDock Vina integration for receptor/ligand
  preparation and batch docking workflows.
"""

__version__ = "2.0.0"
