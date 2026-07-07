"""
Tests for mutadock.mutation.ddg_calc
-----------------------------
PyRosetta, tqdm, and the predict_ddG sub-module are mocked before import.

Run from the project root:
    pytest tests/
"""

import csv
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Path setup.
# Bio / PyRosetta / tqdm stubs are already injected by conftest.py before
# this module is imported.
# ---------------------------------------------------------------------------
from mutadock.mutation import ddg_calc

# ---------------------------------------------------------------------------
# CSV fixture helpers
# ---------------------------------------------------------------------------

_SINGLE_FIELDS = ["sr", "pdb", "chain", "position", "wtAA", "prAA", "wtProb", "prProb"]
_DDG_FIELDS = _SINGLE_FIELDS + ["ddG_value"]


def _write_input_csv(path: str, rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_SINGLE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_output_csv(path: str) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _row(
    sr=1, chain="A", position=10, wt="ALA", pr="GLY", wt_prob=0.02, pr_prob=0.01
) -> dict:
    return {
        "sr": sr,
        "pdb": "test",
        "chain": chain,
        "position": position,
        "wtAA": wt,
        "prAA": pr,
        "wtProb": wt_prob,
        "prProb": pr_prob,
    }


# ---------------------------------------------------------------------------
# Patch context manager
# ---------------------------------------------------------------------------


@contextmanager
def _patch_calc_ddg(wt_score: float = 10.0, mut_score: float = 7.5):
    """Patch all pyrosetta deps *in ddg_calc's namespace* with real-valued mocks.

    Yielded dict keys: pose, sfxn, mut_pose, get_scorefxn.
    Patches init, pose_from_pdb, predict_ddG.get_scorefxn, backup, in_directory.
    The caller is responsible for patching predict_ddG.mutate_residue.

    ΔΔG = score(mutant) − score(wild-type self-mutation reference).  Within each
    site the reference is scored first, then the mutant, so the first score
    returns *wt_score* and subsequent scores return *mut_score* → ddG for a
    single-row call is ``mut_score − wt_score``.
    """
    mock_pose = MagicMock(name="pose")
    pdb_info = MagicMock()
    pdb_info.pdb2pose.return_value = 5
    mock_pose.pdb_info.return_value = pdb_info

    mock_mut_pose = MagicMock(name="mut_pose")

    mock_sfxn = MagicMock(name="sfxn")
    call_state = {"n": 0}

    def _score(p):
        call_state["n"] += 1
        return wt_score if call_state["n"] == 1 else mut_score

    mock_sfxn.score.side_effect = _score
    # calc_ddg now obtains its score function via predict_ddG.get_scorefxn(cartesian)
    # (so the Cartesian protocol can swap in ref2015_cart), not get_fa_scorefxn.
    mock_get_scorefxn = MagicMock(name="get_scorefxn", return_value=mock_sfxn)

    with (
        patch("mutadock.mutation.ddg_calc.init"),
        patch("mutadock.mutation.ddg_calc.pose_from_pdb", return_value=mock_pose),
        patch("mutadock.mutation.predict_ddG.get_scorefxn", mock_get_scorefxn),
        patch("mutadock.mutation.ddg_calc.backup"),
        patch("mutadock.mutation.ddg_calc.in_directory"),
    ):
        yield {
            "pose": mock_pose,
            "sfxn": mock_sfxn,
            "mut_pose": mock_mut_pose,
            "get_scorefxn": mock_get_scorefxn,
        }


# ---------------------------------------------------------------------------
# Tests: calc_ddg()
# ---------------------------------------------------------------------------


class TestCalcDdg(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb_file = str(self.tmpdir / "wt.pdb")
        Path(self.pdb_file).touch()
        self.in_csv = str(self.tmpdir / "mutations.csv")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------ #
    # Output-file naming                                                   #
    # ------------------------------------------------------------------ #

    def test_default_output_filename_derived_from_input(self):
        """Without out_file, the output path appends '_ddG' before .csv."""
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                result = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        expected = self.in_csv.removesuffix(".csv") + "_ddG.csv"
        self.assertEqual(result, expected)

    def test_explicit_out_file_is_used(self):
        """When out_file is given, that exact path is used."""
        _write_input_csv(self.in_csv, [_row()])
        out = str(self.tmpdir / "custom_out.csv")
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                result = ddg_calc.calc_ddg(self.pdb_file, self.in_csv, out_file=out)

        self.assertEqual(result, out)

    # ------------------------------------------------------------------ #
    # Output CSV structure                                                 #
    # ------------------------------------------------------------------ #

    def test_output_csv_is_created(self):
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        self.assertTrue(Path(out).is_file())

    def test_output_has_ddg_value_column(self):
        """Output CSV must contain the 'ddG_value' column."""
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        self.assertGreater(len(rows), 0)
        self.assertIn("ddG_value", rows[0])

    def test_output_preserves_input_columns(self):
        """All original columns must survive in the output."""
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        for col in _SINGLE_FIELDS:
            self.assertIn(col, rows[0], msg=f"Column '{col}' missing from output")

    def test_all_input_rows_appear_in_output(self):
        """Row count in output must equal row count in input."""
        input_rows = [_row(sr=i, position=i) for i in range(1, 6)]
        _write_input_csv(self.in_csv, input_rows)
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        self.assertEqual(len(rows), len(input_rows))

    # ------------------------------------------------------------------ #
    # ddG arithmetic                                                       #
    # ------------------------------------------------------------------ #

    def test_ddg_equals_mutant_score_minus_wt_score(self):
        """ddG_value must be score(mutant) − score(wildtype)."""
        wt_score = 12.0
        mut_score = 8.5
        expected_ddg = mut_score - wt_score  # -3.5

        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg(wt_score, mut_score) as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(float(rows[0]["ddG_value"]), expected_ddg, places=5)

    def test_positive_ddg_for_destabilising_mutation(self):
        """A mutation that raises the score (destabilising) → positive ddG."""
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg(wt_score=5.0, mut_score=9.0) as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        self.assertGreater(float(rows[0]["ddG_value"]), 0)

    def test_negative_ddg_for_stabilising_mutation(self):
        """A mutation that lowers the score (stabilising) → negative ddG."""
        _write_input_csv(self.in_csv, [_row()])
        with _patch_calc_ddg(wt_score=10.0, mut_score=3.0) as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                out = ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        rows = _read_output_csv(out)
        self.assertLess(float(rows[0]["ddG_value"]), 0)

    # ------------------------------------------------------------------ #
    # Interactions with dependencies                                       #
    # ------------------------------------------------------------------ #

    def test_scorefxn_built_once(self):
        """Score function is built once (via get_scorefxn) and reused."""
        input_rows = [_row(sr=i, position=i) for i in range(1, 4)]
        _write_input_csv(self.in_csv, input_rows)
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        mocks["get_scorefxn"].assert_called_once()

    def test_mutate_residue_called_for_reference_and_mutant(self):
        """WT self-reference is computed once per site (cached), plus one mutant per row.

        ``mutate_residue`` runs once for the wild-type self-mutation reference
        and once per mutant.  The mocked ``pdb2pose`` maps every row to the same
        pose position (5), so the reference is computed once and reused — giving
        ``1 + N`` calls for N rows and exercising the per-site reference cache.
        """
        input_rows = [_row(sr=i, position=i) for i in range(1, 4)]
        _write_input_csv(self.in_csv, input_rows)
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ) as mock_mut:
                ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        self.assertEqual(mock_mut.call_count, 1 + len(input_rows))

    def test_backup_called_for_output_file(self):
        """backup() must be called for the output CSV path."""
        _write_input_csv(self.in_csv, [_row()])
        out = str(self.tmpdir / "out_ddG.csv")
        with _patch_calc_ddg() as mocks:
            with (
                patch.object(
                    ddg_calc.predict_ddG,
                    "mutate_residue",
                    return_value=mocks["mut_pose"],
                ),
                patch("mutadock.mutation.ddg_calc.backup") as mock_backup,
            ):
                ddg_calc.calc_ddg(self.pdb_file, self.in_csv, out_file=out)

        backed_up_paths = [c.args[0] for c in mock_backup.call_args_list]
        self.assertIn(out, backed_up_paths)

    def test_chain_and_position_passed_to_pdb2pose(self):
        """pdb2pose must receive the chain and position from each CSV row."""
        row_data = _row(chain="B", position=25)
        _write_input_csv(self.in_csv, [row_data])
        with _patch_calc_ddg() as mocks:
            with patch.object(
                ddg_calc.predict_ddG, "mutate_residue", return_value=mocks["mut_pose"]
            ):
                ddg_calc.calc_ddg(self.pdb_file, self.in_csv)

        pdb2pose_calls = mocks["pose"].pdb_info.return_value.pdb2pose.call_args_list
        self.assertTrue(
            any(c.args == ("B", 25) for c in pdb2pose_calls),
            msg=f"Expected pdb2pose('B', 25) in {pdb2pose_calls}",
        )


if __name__ == "__main__":
    unittest.main()
