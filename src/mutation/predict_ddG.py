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


from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory


def mutate_residue(pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn):
    if pose.is_fullatom() == False:
        IOError("mutate_residue only works with fullatom poses")

    test_pose = Pose()
    test_pose.assign(pose)

    task = TaskFactory.create_packer_task(test_pose)

    aa_bool = vector1_bool()

    mutant_aa = aa_from_oneletter_code(mutant_aa)

    for i in range(1, 21):
        aa_bool.append(i == mutant_aa)

    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        dist = center.distance_squared(test_pose.residue(i).nbr_atom_xyz())

        if i != mutant_position and dist > pow(float(pack_radius), 2):
            task.nonconst_residue_task(i).prevent_repacking()
        elif i != mutant_position and dist <= pow(float(pack_radius), 2):
            task.nonconst_residue_task(i).restrict_to_repacking()

    packer = PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose


if __name__ == "__main__":
    import sys
    sys.exit("This file only contains functions used in other libraries.")
