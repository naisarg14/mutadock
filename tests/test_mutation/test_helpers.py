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
from unittest.mock import patch

import pytest

from mutadock.mutation.exceptions import MutationError
from mutadock.mutation.helpers import (
    clean_pdb,
    fetch_pdb,
    format_missing_residue,
    structure_warnings,
)

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


# ---------------------------------------------------------------------------
# 3.4: fetch_pdb — fetch a structure from RCSB by PDB ID
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad_id", ["AB", "toolong", "12-4", ""])
def test_fetch_pdb_rejects_invalid_id(bad_id):
    # Malformed IDs must be rejected *before* any network access.
    with pytest.raises(MutationError, match="Invalid PDB ID"):
        fetch_pdb(bad_id)


def test_fetch_pdb_downloads_when_absent(tmp_path: Path):
    def fake_urlretrieve(url, dest):
        Path(dest).write_text("ATOM  ...\nEND\n")

    with patch(
        "mutadock.mutation.helpers.urllib.request.urlretrieve",
        side_effect=fake_urlretrieve,
    ) as mock_get:
        out = fetch_pdb("4qjr", dest_dir=tmp_path)  # lower-case on purpose

    # ID is upper-cased for the filename and the download happened once.
    assert Path(out) == tmp_path / "4QJR.pdb"
    assert Path(out).is_file()
    mock_get.assert_called_once()
    called_url = mock_get.call_args[0][0]
    assert called_url == "https://files.rcsb.org/download/4QJR.pdb"


def test_fetch_pdb_reuses_existing_file(tmp_path: Path):
    existing = tmp_path / "4QJR.pdb"
    existing.write_text("ATOM  ...\nEND\n")

    with patch("mutadock.mutation.helpers.urllib.request.urlretrieve") as mock_get:
        out = fetch_pdb("4QJR", dest_dir=tmp_path)

    assert Path(out) == existing
    mock_get.assert_not_called()  # cached file must not be re-downloaded


# ---------------------------------------------------------------------------
# 3.4: structure_warnings — flag features that affect residue numbering
# ---------------------------------------------------------------------------


def _atom(serial=1, altloc=" ", icode=" ", chain="A", resseq=1, resname="MET"):
    """Build a column-correct PDB ATOM line (altLoc col 17, iCode col 27)."""
    return (
        f"ATOM  {serial:>5} N   {altloc}{resname:>3} {chain}{resseq:>4}{icode}"
        "   11.104  13.207  10.567  1.00  0.00           N\n"
    )


def _write(tmp_path: Path, name: str, text: str) -> str:
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_structure_warnings_multi_model(tmp_path: Path):
    pdb = _write(
        tmp_path,
        "nmr.pdb",
        "MODEL        1\n" + _atom(1) + "ENDMDL\n"
        "MODEL        2\n" + _atom(2) + "ENDMDL\nEND\n",
    )
    warnings = structure_warnings(pdb)
    assert any("2 models" in w for w in warnings)


def test_structure_warnings_altloc(tmp_path: Path):
    pdb = _write(
        tmp_path, "alt.pdb", _atom(1, altloc="A") + _atom(2, altloc="B") + "END\n"
    )
    warnings = structure_warnings(pdb)
    assert any("alternate location" in w for w in warnings)


def test_structure_warnings_insertion_code_counts_residues(tmp_path: Path):
    # Two atoms of the SAME insertion-coded residue -> counted as one residue.
    pdb = _write(
        tmp_path,
        "ins.pdb",
        _atom(1, icode="A", resseq=52) + _atom(2, icode="A", resseq=52) + "END\n",
    )
    warnings = structure_warnings(pdb)
    assert any("1 residue(s) with insertion codes" in w for w in warnings)


def test_structure_warnings_clean_structure(tmp_path: Path):
    pdb = _write(tmp_path, "clean.pdb", _atom(1) + _atom(2, resseq=2) + "END\n")
    assert structure_warnings(pdb) == []


def test_structure_warnings_missing_file_is_silent(tmp_path: Path):
    assert structure_warnings(str(tmp_path / "does_not_exist.pdb")) == []


# ---------------------------------------------------------------------------
# 3.4: format_missing_residue — actionable "not found" messages
# ---------------------------------------------------------------------------


def test_format_missing_residue_lists_chains_when_chain_absent():
    available = {"A": [1, 2, 3], "B": [10, 11]}
    msg = format_missing_residue(available, "Z", 5, "prot.pdb")
    assert "Chain 'Z' not found" in msg
    assert "A, B" in msg  # available chains listed, sorted


def test_format_missing_residue_shows_span_when_position_absent():
    available = {"A": [218, 219, 220, 461]}
    msg = format_missing_residue(available, "A", 999, "prot.pdb")
    assert "Residue 999 not found in chain 'A'" in msg
    assert "218-461" in msg  # numbering span reported
