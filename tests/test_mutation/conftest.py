"""
conftest.py — test_mutation package
=====================================
This file's *module-level* code runs during pytest startup, BEFORE any test
module in this package is imported. That guarantees the heavyweight optional
dependencies (BioPython, PyRosetta, tqdm) are already stubbed in sys.modules
when each test file does its top-level imports — without requiring the real
packages to be installed.
"""

import sys
import types
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _stub(name: str, **attrs) -> types.ModuleType:
    """Insert a lightweight module stub into sys.modules (no-op if present)."""
    if name not in sys.modules:
        mod = types.ModuleType(name)
        for k, v in attrs.items():
            setattr(mod, k, v)
        sys.modules[name] = mod
    return sys.modules[name]


# ---------------------------------------------------------------------------
# BioPython  (needed by csv_generator.py and helpers.py → convert_cif_pdb)
# ---------------------------------------------------------------------------
# Only stub BioPython when the real package is unavailable.  These tests patch
# ``mutadock.mutation.csv_generator.PDBParser`` locally, so they pass either way;
# stubbing unconditionally used to clobber ``sys.modules["Bio.PDB"]`` with a
# MagicMock that leaked into sibling test packages (e.g. test_report, which does
# real Cα parsing) and made their PDBParser return an empty structure.
try:
    import Bio.PDB  # noqa: F401
except ImportError:
    _pdb_parser_cls = MagicMock(name="PDBParser")
    _mmcif_parser_cls = MagicMock(name="MMCIFParser")
    _pdbio_cls = MagicMock(name="PDBIO")

    _bio = _stub("Bio")
    _bio_pdb = _stub(
        "Bio.PDB",
        PDBParser=_pdb_parser_cls,
        MMCIFParser=_mmcif_parser_cls,
        PDBIO=_pdbio_cls,
    )
    _bio.PDB = _bio_pdb

# ---------------------------------------------------------------------------
# PyRosetta core  (needed by generate_mutants.py, ddg_calc.py, predict_ddG.py)
# ---------------------------------------------------------------------------
_init_mock = MagicMock(name="init")
_pfp_mock = MagicMock(name="pose_from_pdb")
_gfa_mock = MagicMock(name="get_fa_scorefxn")
_pose_cls_mock = MagicMock(name="Pose")

_csf_mock = MagicMock(name="create_score_function")

_pyrosetta = _stub(
    "pyrosetta",
    init=_init_mock,
    pose_from_pdb=_pfp_mock,
    get_fa_scorefxn=_gfa_mock,
    create_score_function=_csf_mock,
    Pose=_pose_cls_mock,
    __all__=[
        "init",
        "pose_from_pdb",
        "get_fa_scorefxn",
        "create_score_function",
        "Pose",
    ],
)

# ---------------------------------------------------------------------------
# PyRosetta sub-packages  (needed by predict_ddG.py)
# ---------------------------------------------------------------------------
_mutate_residue_mock = MagicMock(name="toolbox_mutate_residue")

_stub("pyrosetta.toolbox", mutate_residue=_mutate_residue_mock)
_stub("pyrosetta.teaching")

_rosetta = _stub("pyrosetta.rosetta")
_ros_util = _stub(
    "pyrosetta.rosetta.utility", vector1_bool=MagicMock(name="vector1_bool")
)
_ros_core = _stub("pyrosetta.rosetta.core")
_ros_chem = _stub(
    "pyrosetta.rosetta.core.chemical",
    aa_from_oneletter_code=MagicMock(name="aa_from_oneletter_code"),
)
_ros_pack = _stub("pyrosetta.rosetta.core.pack")
_ros_task = _stub(
    "pyrosetta.rosetta.core.pack.task", TaskFactory=MagicMock(name="TaskFactory")
)
_ros_kinematics = _stub(
    "pyrosetta.rosetta.core.kinematics", MoveMap=MagicMock(name="MoveMap")
)
_ros_proto = _stub("pyrosetta.rosetta.protocols")
_ros_minpack = _stub(
    "pyrosetta.rosetta.protocols.minimization_packing",
    PackRotamersMover=MagicMock(name="PackRotamersMover"),
    MinMover=MagicMock(name="MinMover"),
)

# Expose sub-packages as attributes on their parents so attribute access works
_pyrosetta.toolbox = sys.modules["pyrosetta.toolbox"]
_pyrosetta.teaching = sys.modules["pyrosetta.teaching"]
_pyrosetta.rosetta = _rosetta
_rosetta.utility = _ros_util
_rosetta.core = _ros_core
_ros_core.chemical = _ros_chem
_ros_core.pack = _ros_pack
_ros_core.kinematics = _ros_kinematics
_ros_pack.task = _ros_task
_rosetta.protocols = _ros_proto
_ros_proto.minimization_packing = _ros_minpack

# ---------------------------------------------------------------------------
# tqdm  (needed by generate_mutants.py and ddg_calc.py)
# ---------------------------------------------------------------------------
_stub("tqdm", tqdm=lambda x, **kw: x)
