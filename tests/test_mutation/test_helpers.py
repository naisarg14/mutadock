"""
Tests for mutadock.mutation.helpers.clean_pdb
---------------------------------------------
Bio / PyRosetta / tqdm stubs are injected by conftest.py before this module is
imported, so ``from mutadock.mutation.helpers import ...`` does not pull in any
heavyweight optional dependency (helpers.py performs those imports lazily inside
functions, not at module level).

Run from the project root:
    pytest tests/test_mutation/test_helpers.py -q
"""

from pathlib import Path

from mutadock.mutation.helpers import clean_pdb

# ---------------------------------------------------------------------------
# Synthetic PDB fixtures
# ---------------------------------------------------------------------------

# A multi-model style file with a bare ``END`` marker *between* the two ATOM
# blocks.  The pre-fix implementation broke out of its loop at that first
# ``END``, silently discarding everything after it (model 2).  A correct
# implementation keeps all ATOM records.
_MULTI_MODEL_PDB = (
    "HEADER    TEST MULTI-MODEL\n"
    "MODEL        1\n"
    "ATOM      1  N   MET A   1      11.104  13.207  10.567  1.00  0.00           N\n"
    "ATOM      2  CA  MET A   1      12.560  13.207  10.567  1.00  0.00           C\n"
    "TER       3      MET A   1\n"
    "ENDMDL\n"
    "END\n"  # bare END before the second model -> old code truncated here
    "MODEL        2\n"
    "ATOM      4  N   GLY A   1      21.104  23.207  20.567  1.00  0.00           N\n"
    "ATOM      5  CA  GLY A   1      22.560  23.207  20.567  1.00  0.00           C\n"
    "TER       6      GLY A   1\n"
    "ENDMDL\n"
    "END\n"
)

_PDB_WITH_HETATM = (
    "HEADER    TEST HETATM\n"
    "ATOM      1  N   MET A   1      11.104  13.207  10.567  1.00  0.00           N\n"
    "ATOM      2  CA  MET A   1      12.560  13.207  10.567  1.00  0.00           C\n"
    "HETATM    3  O   HOH A 101      31.104  33.207  30.567  1.00  0.00           O\n"
    "HETATM    4 ZN    ZN A 102      41.104  43.207  40.567  1.00  0.00          ZN\n"
    "TER       5      MET A   1\n"
    "END\n"
)


def _atom_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.startswith("ATOM")]


# ---------------------------------------------------------------------------
# Bug 2.15 (1): multi-model / stray-END files must not be truncated
# ---------------------------------------------------------------------------


def test_multi_model_not_truncated_at_first_end(tmp_path: Path):
    src = tmp_path / "multi.pdb"
    src.write_text(_MULTI_MODEL_PDB)

    out = clean_pdb(str(src), str(tmp_path / "multi_clean.pdb"))
    cleaned = Path(out).read_text()

    atoms = _atom_lines(cleaned)
    # Both models' ATOM records must survive (2 from model 1, 2 from model 2).
    assert len(atoms) == 4, f"expected 4 ATOM records, got {len(atoms)}:\n{cleaned}"
    assert any("MET A   1" in ln for ln in atoms)  # model 1 content
    assert any("GLY A   1" in ln for ln in atoms)  # model 2 content (dropped by bug)


def test_single_terminal_end_marker(tmp_path: Path):
    # There must be exactly one END, and it must be the final record, so that
    # downstream PDB parsers never truncate the structure at an early marker.
    src = tmp_path / "multi.pdb"
    src.write_text(_MULTI_MODEL_PDB)

    out = clean_pdb(str(src), str(tmp_path / "multi_clean.pdb"))
    lines = [ln for ln in Path(out).read_text().splitlines() if ln.strip()]

    end_lines = [ln for ln in lines if ln.strip() == "END"]
    assert len(end_lines) == 1, f"expected exactly one END, got {len(end_lines)}"
    assert lines[-1].strip() == "END", "END must be the last record"


# ---------------------------------------------------------------------------
# Bug 2.15 (2): HETATM stripping is the deliberate default
# ---------------------------------------------------------------------------


def test_hetatm_stripped_by_default(tmp_path: Path):
    src = tmp_path / "het.pdb"
    src.write_text(_PDB_WITH_HETATM)

    out = clean_pdb(str(src), str(tmp_path / "het_clean.pdb"))
    cleaned = Path(out).read_text()

    assert "HETATM" not in cleaned, "HETATM records must be stripped by default"
    assert len(_atom_lines(cleaned)) == 2  # the two protein atoms remain


def test_hetatm_kept_when_requested(tmp_path: Path):
    src = tmp_path / "het.pdb"
    src.write_text(_PDB_WITH_HETATM)

    out = clean_pdb(str(src), str(tmp_path / "het_keep_clean.pdb"), keep_hetatm=True)
    cleaned = Path(out).read_text()

    het_lines = [ln for ln in cleaned.splitlines() if ln.startswith("HETATM")]
    assert len(het_lines) == 2, "both HETATM records must be preserved"
    assert len(_atom_lines(cleaned)) == 2


# ---------------------------------------------------------------------------
# Regression guard: normal single-model output shape is preserved
# ---------------------------------------------------------------------------


def test_default_output_path_and_header(tmp_path: Path):
    src = tmp_path / "plain.pdb"
    src.write_text(
        "ATOM      1  N   MET A   1      11.104  13.207  10.567  1.00  0.00           N\n"
        "TER       2      MET A   1\n"
        "END\n"
    )

    out = clean_pdb(str(src))  # default destination: <stem>_clean.pdb
    assert out == str(src).replace(".pdb", "_clean.pdb")

    cleaned = Path(out).read_text()
    assert cleaned.startswith("REMARK This file was cleaned")
    assert len(_atom_lines(cleaned)) == 1
    assert cleaned.rstrip().endswith("END")
