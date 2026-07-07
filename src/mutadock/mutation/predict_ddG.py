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

import logging
import statistics
from typing import Any

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# REU -> kcal/mol scaling factor.
#
# Rosetta ΔΔG values are in *Rosetta Energy Units (REU)*, NOT kcal/mol. To get
# numbers that are roughly comparable to experiment, REU are commonly scaled by
# a constant. For the ref2015 `cartesian_ddg` protocol the community convention
# is to divide REU by ~2.94 (i.e. multiply by ~0.34); reported values in the
# literature span roughly 0.29–0.40 depending on protocol and dataset.
#
#   reference: Park, H. et al. (2016) "Simultaneous Optimization of Biomolecular
#   Energy Functions on Features from Small Molecules and Macromolecules."
#   J. Chem. Theory Comput. 12(12):6201-6212. doi:10.1021/acs.jctc.6b00819
#
# >>> CHANGE THIS VALUE to use a different scaling, or override per-run with the
# >>> `--reu-to-kcal` CLI flag (md_mutate / md_quick / md_ddg_*). Reported
# >>> ddG_kcal = ddG_value(REU) * REU_TO_KCAL_SCALE.
# ---------------------------------------------------------------------------
REU_TO_KCAL_SCALE = 0.34  # ≈ 1 / 2.94

try:
    from pyrosetta import *
    from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.protocols.minimization_packing import (
        MinMover,
        PackRotamersMover,
    )
    from pyrosetta.rosetta.utility import vector1_bool
    from pyrosetta.teaching import *
except ImportError:
    import sys

    msg = "Error with importing pyrosetta module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install pyrosetta using the following command:\n\n"
    msg += "python -m pip install pyrosetta_installer && python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'\n"
    msg += "Alternative is to install from PyRosetta's official website.\n"
    msg += "If you already have pyrosetta installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    logger.error(msg)
    sys.exit(2)


def get_scorefxn(cartesian: bool = False) -> Any:
    """Return the score function for the requested ΔΔG protocol.

    Args:
        cartesian: If ``True`` return the Cartesian-space ``ref2015_cart``
            score function (required for Cartesian backbone minimization, as in
            Rosetta's ``cartesian_ddg`` protocol); otherwise the default
            full-atom ``ref2015`` score function.
    """
    if cartesian:
        return create_score_function("ref2015_cart")
    return get_fa_scorefxn()


def mutate_residue(
    pose: Any,
    mutant_position: int,
    mutant_aa: str,
    pack_radius: float,
    pack_scorefxn: Any,
    backbone_minimization: bool = False,
    cartesian: bool = False,
) -> Any:
    """Apply a single point mutation and repack (and optionally minimize).

    Creates a copy of *pose*, substitutes residue *mutant_position* with
    *mutant_aa*, and runs ``PackRotamersMover`` on all residues within
    *pack_radius* Ångströms of the mutation site.  When *backbone_minimization*
    is set, the neighbourhood's backbone and side chains are additionally
    energy-minimized (Cartesian-space when *cartesian* is ``True``), approaching
    Rosetta's ``cartesian_ddg`` protocol.

    Args:
        pose: PyRosetta ``Pose`` object (must be full-atom).
        mutant_position: Rosetta pose numbering of the residue to mutate.
        mutant_aa: One-letter code of the target amino acid.
        pack_radius: Distance threshold (Å) for side-chain repacking (and, when
            enabled, backbone minimization).
        pack_scorefxn: PyRosetta score function used by ``PackRotamersMover``
            and the optional ``MinMover``.
        backbone_minimization: If ``True``, minimize backbone + side chains of
            the neighbourhood after packing.  Slower but less noisy / more
            accurate.
        cartesian: If ``True`` the minimization is done in Cartesian space
            (requires *pack_scorefxn* to be a ``*_cart`` function, e.g. from
            ``get_scorefxn(cartesian=True)``).

    Returns:
        A new ``Pose`` with the mutation applied and neighbours repacked
        (and minimized when requested).

    Raises:
        OSError: If *pose* is not a full-atom pose.
    """
    if not pose.is_fullatom():
        raise OSError("mutate_residue only works with fullatom poses")

    test_pose = Pose()
    test_pose.assign(pose)

    task = TaskFactory.create_packer_task(test_pose)

    aa_bool = vector1_bool()

    mutant_aa = aa_from_oneletter_code(mutant_aa)

    for i in range(1, 21):
        aa_bool.append(i == mutant_aa)

    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)
    center = pose.residue(mutant_position).nbr_atom_xyz()
    radius_sq = pow(float(pack_radius), 2)
    # Track the neighbourhood (mutation site + residues within pack_radius) so
    # the optional minimization operates on the same shell that was repacked.
    in_shell = [False] * (pose.total_residue() + 1)
    in_shell[mutant_position] = True
    for i in range(1, pose.total_residue() + 1):
        dist = center.distance_squared(test_pose.residue(i).nbr_atom_xyz())

        if i != mutant_position and dist > radius_sq:
            task.nonconst_residue_task(i).prevent_repacking()
        elif i != mutant_position and dist <= radius_sq:
            task.nonconst_residue_task(i).restrict_to_repacking()
            in_shell[i] = True

    packer = PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    if backbone_minimization:
        movemap = MoveMap()
        for i in range(1, test_pose.total_residue() + 1):
            if in_shell[i]:
                movemap.set_bb(i, True)
                movemap.set_chi(i, True)
        min_mover = MinMover(
            movemap, pack_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, True
        )
        if cartesian:
            min_mover.cartesian(True)
        min_mover.apply(test_pose)

    return test_pose


def apply_mutations(
    pose: Any,
    mutations: list[tuple[int, str]],
    pack_scorefxn: Any,
    pack_radius: float,
    backbone_minimization: bool = False,
    cartesian: bool = False,
) -> Any:
    """Apply a list of ``(pose_position, aa1)`` mutations sequentially.

    Each mutation repacks (and optionally minimizes) its own neighbourhood via
    :func:`mutate_residue`.  Passing the wild-type residues here produces the
    *self-mutation reference* used to remove packing bias (see
    :func:`wt_reference_score`).
    """
    work = pose
    for position, aa1 in mutations:
        work = mutate_residue(
            work,
            position,
            aa1,
            pack_radius,
            pack_scorefxn,
            backbone_minimization=backbone_minimization,
            cartesian=cartesian,
        )
    return work


def wt_reference_score(
    pose: Any,
    positions: list[int],
    pack_scorefxn: Any,
    pack_radius: float,
    backbone_minimization: bool = False,
    cartesian: bool = False,
) -> float:
    """Score of the *wild-type self-mutation* reference for *positions*.

    Applies the identical repack (+minimization) protocol to the wild-type
    residues at *positions*, so the ΔΔG ``score(mutant) - score(reference)`` is a
    fair difference.  Without this, repacking only the mutant lowers its energy
    relative to the untouched wild type, biasing every ΔΔG toward "stabilizing"
    (a null WT→WT mutation would otherwise score strongly negative instead of
    ~0).  The reference depends only on the *set of positions*, so callers
    should cache it per site across a saturation scan.
    """
    wt = [(p, pose.residue(p).name1()) for p in positions]
    ref_pose = apply_mutations(
        pose,
        wt,
        pack_scorefxn,
        pack_radius,
        backbone_minimization=backbone_minimization,
        cartesian=cartesian,
    )
    return pack_scorefxn.score(ref_pose)


def protocol_label(backbone_minimization: bool, cartesian: bool) -> str:
    """Name the ΔΔG protocol for provenance: ``fast`` / ``min`` / ``cartesian``.

    ``fast`` (single repack, no minimization) is screening-only — absolute ΔΔG
    values are unreliable because a failed wild-type repack can invert a site's
    ranking; ``min`` and ``cartesian`` add minimization and are trustworthy.
    """
    if cartesian:
        return "cartesian"
    if backbone_minimization:
        return "min"
    return "fast"


def summarize_ddg(values: list[float], scale: float = REU_TO_KCAL_SCALE) -> dict:
    """Summarize a list of replicate ΔΔG estimates (in REU).

    Args:
        values: Per-replicate ΔΔG values in Rosetta Energy Units.
        scale: REU → kcal/mol scaling factor applied to the mean.

    Returns:
        Dict with ``ddG_value`` (mean REU), ``ddG_sd`` (population-sample SD,
        0.0 for a single replicate), ``ddG_kcal`` (mean × *scale*), and
        ``n_replicates``.
    """
    mean = statistics.fmean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "ddG_value": mean,
        "ddG_sd": sd,
        "ddG_kcal": mean * scale,
        "n_replicates": len(values),
    }


if __name__ == "__main__":
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's docking module."
    )
