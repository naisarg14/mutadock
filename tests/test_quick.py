"""Tests for ``mutadock.quick`` (the md_quick pipeline).

The heavy steps (RCSB/PubChem fetch, PyRosetta ΔΔG, receptor prep, AutoSite,
Vina) are mocked; these tests exercise the orchestration, mutation validation,
and error handling.
"""

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mutadock import quick
from mutadock.docking.exceptions import DockingRunError
from mutadock.mutation.exceptions import MutationError

RESIDUES = {1: ("A", 386, "ASN"), 2: ("A", 387, "HIS")}


# --------------------------------------------------------------------------- #
# _validate_and_fill_mutation                                                 #
# --------------------------------------------------------------------------- #
class TestValidateMutation:
    def _pdb(self, tmp_path):
        p = tmp_path / "p.pdb"
        p.write_text("ATOM\n")
        return str(p)

    def test_fills_wt_and_normalizes_praa(self, tmp_path):
        with patch("mutadock.quick.get_residues", return_value=RESIDUES):
            mut = quick._validate_and_fill_mutation(
                self._pdb(tmp_path),
                {"chain": "A", "position": 386, "wtAA": "X", "prAA": "H"},
            )
        assert mut["wtAA"] == "ASN"  # filled from structure
        assert mut["prAA"] == "HIS"  # 1-letter normalized to 3-letter

    def test_missing_residue_raises(self, tmp_path):
        with patch("mutadock.quick.get_residues", return_value=RESIDUES):
            with pytest.raises(MutationError):
                quick._validate_and_fill_mutation(
                    self._pdb(tmp_path),
                    {"chain": "A", "position": 999, "wtAA": "X", "prAA": "ALA"},
                )

    def test_unknown_mutant_aa_raises(self, tmp_path):
        with patch("mutadock.quick.get_residues", return_value=RESIDUES):
            with pytest.raises(MutationError):
                quick._validate_and_fill_mutation(
                    self._pdb(tmp_path),
                    {"chain": "A", "position": 386, "wtAA": "ASN", "prAA": "ZZZ"},
                )


# --------------------------------------------------------------------------- #
# run_quick orchestration                                                     #
# --------------------------------------------------------------------------- #
def _common_patches(dock):
    """patch.multiple kwargs for mutadock.quick with the heavy steps stubbed."""

    def fake_fetch_pdb(pdb_id, dest_dir=None):
        p = Path(dest_dir) / f"{pdb_id}.pdb"
        p.write_text("HEADER\nEND\n")
        return str(p)

    def fake_clean(src, dst):
        Path(dst).write_text("ATOM\n")

    def fake_calc(pdb, mut_csv, out_file=None, **kwargs):
        with open(out_file, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "sr",
                    "pdb",
                    "chain",
                    "position",
                    "wtAA",
                    "prAA",
                    "wtProb",
                    "prProb",
                    "ddG_value",
                ]
            )
            w.writerow([1, pdb, "A", 386, "ASN", "HIS", 0, 0, -27.3])
        return out_file

    def fake_gen(pdb, muts, output_folder=None):
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        m = Path(output_folder) / "1_ASN-A386-H.pdb"
        m.write_text("ATOM\n")
        return [str(m)]

    def fake_lig(code, dest_dir=None):
        p = Path(dest_dir) / "aspirin.sdf"
        p.write_text("sdf")
        return str(p)

    def fake_report(od, base_name=None, quiet=False):
        return {"html": Path(od) / "reports" / "report.html"}

    return {
        "fetch_pdb": fake_fetch_pdb,
        "clean_pdb": fake_clean,
        "get_residues": MagicMock(return_value=RESIDUES),
        "calc_ddg": fake_calc,
        "generate_pdb": fake_gen,
        "fetch_ligand": fake_lig,
        "_dock": dock,
        "generate_report": MagicMock(side_effect=fake_report),
    }


def test_run_quick_happy_path(tmp_path):
    dock = MagicMock(return_value=(-7.5, str(tmp_path / "pose.sdf")))
    with patch.multiple("mutadock.quick", **_common_patches(dock)):
        res = quick.run_quick(
            pdb_id="4QJR",
            input_path=None,
            mutation="A:386:ASN:HIS",
            ligand_code="aspirin",
            ligand_file=None,
            out_dir=tmp_path / "out",
            quiet=True,
        )
    assert res["ddg"] == pytest.approx(-27.3)
    assert res["affinity"] == pytest.approx(-7.5)
    assert res["label"] == "ASN-A386-HIS"
    assert res["report"].name == "report.html"
    dock.assert_called_once()


def test_run_quick_docking_failure_still_reports(tmp_path):
    dock = MagicMock(side_effect=DockingRunError("vina blew up"))
    with patch.multiple("mutadock.quick", **_common_patches(dock)):
        res = quick.run_quick(
            pdb_id="4QJR",
            input_path=None,
            mutation="A:386:ASN:HIS",
            ligand_code="aspirin",
            ligand_file=None,
            out_dir=tmp_path / "out",
            quiet=True,
        )
    # ΔΔG and report survive a docking failure; affinity is absent.
    assert res["ddg"] == pytest.approx(-27.3)
    assert res["affinity"] is None
    assert "report" in res


def test_run_quick_no_report(tmp_path):
    dock = MagicMock(return_value=(-7.5, str(tmp_path / "pose.sdf")))
    patches = _common_patches(dock)
    with patch.multiple("mutadock.quick", **patches):
        res = quick.run_quick(
            pdb_id="4QJR",
            input_path=None,
            mutation="A:386:ASN:HIS",
            ligand_code="aspirin",
            ligand_file=None,
            out_dir=tmp_path / "out",
            make_report=False,
            quiet=True,
        )
    assert "report" not in res
    patches["generate_report"].assert_not_called()
