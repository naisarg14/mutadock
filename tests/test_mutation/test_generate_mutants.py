"""
Tests for mutation.generate_mutants
-------------------------------------
PyRosetta and tqdm are mocked at the sys.modules level *before* the module is
imported so the tests run without the real (heavyweight) packages.

Run from the project root:
    pytest tests/
"""

import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Path setup.
# Bio / PyRosetta / tqdm stubs are already injected by conftest.py before
# this module is imported.
# ---------------------------------------------------------------------------
from mutation import generate_mutants

# Grab references to the toolbox mock so tests can reset / inspect it.
_pyrosetta_mock = sys.modules["pyrosetta"]
_toolbox_mock = sys.modules["pyrosetta.toolbox"]


# ---------------------------------------------------------------------------
# CSV fixture helpers
# ---------------------------------------------------------------------------


def _write_single_csv(path: str, rows: list[dict]) -> None:
    """Write a minimal single-mutation CSV."""
    fieldnames = [
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
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_double_csv(path: str, rows: list[dict]) -> None:
    """Write a minimal double-mutation CSV."""
    fieldnames = [
        "sr",
        "combination",
        "chain1",
        "mut1_position",
        "mut1_wtAA",
        "mut1_prAA",
        "chain2",
        "mut2_position",
        "mut2_wtAA",
        "mut2_prAA",
        "ddG_value",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_triple_csv(path: str, rows: list[dict]) -> None:
    """Write a minimal triple-mutation CSV."""
    fieldnames = [
        "sr",
        "combination",
        "chain1",
        "mut1_position",
        "mut1_wtAA",
        "mut1_prAA",
        "chain2",
        "mut2_position",
        "mut2_wtAA",
        "mut2_prAA",
        "chain3",
        "mut3_position",
        "mut3_wtAA",
        "mut3_prAA",
        "ddG_value",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _single_row(sr=1, chain="A", position=10, wt="ALA", pr="GLY") -> dict:
    return {
        "sr": sr,
        "pdb": "test",
        "chain": chain,
        "position": position,
        "wtAA": wt,
        "prAA": pr,
        "wtProb": 0.02,
        "prProb": 0.01,
        "ddG_value": -1.5,
    }


def _double_row(sr=1) -> dict:
    return {
        "sr": sr,
        "combination": "ALA-A10-GLY_SER-B20-THR",
        "chain1": "A",
        "mut1_position": 10,
        "mut1_wtAA": "ALA",
        "mut1_prAA": "GLY",
        "chain2": "B",
        "mut2_position": 20,
        "mut2_wtAA": "SER",
        "mut2_prAA": "THR",
        "ddG_value": -2.0,
    }


def _triple_row(sr=1) -> dict:
    return {
        "sr": sr,
        "combination": "ALA-A10-GLY_SER-B20-THR_VAL-C30-ILE",
        "chain1": "A",
        "mut1_position": 10,
        "mut1_wtAA": "ALA",
        "mut1_prAA": "GLY",
        "chain2": "B",
        "mut2_position": 20,
        "mut2_wtAA": "SER",
        "mut2_prAA": "THR",
        "chain3": "C",
        "mut3_position": 30,
        "mut3_wtAA": "VAL",
        "mut3_prAA": "ILE",
        "ddG_value": -3.0,
    }


# ---------------------------------------------------------------------------
# Shared mock-pose factory
# ---------------------------------------------------------------------------


def _make_pose(chain="A", position=10, pose_position=5):
    pose = MagicMock(name="pose")
    pdb_info = MagicMock()
    pdb_info.pdb2pose.return_value = pose_position
    pose.pdb_info.return_value = pdb_info
    return pose


# ---------------------------------------------------------------------------
# Tests: generate_single_mutation()
# ---------------------------------------------------------------------------


class TestGenerateSingleMutation(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb_file = str(self.tmpdir / "wt.pdb")
        Path(self.pdb_file).touch()

        self.csv_file = str(self.tmpdir / "mutations.csv")
        self.mock_pose = _make_pose()
        _pyrosetta_mock.pose_from_pdb.return_value = self.mock_pose

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)
        _pyrosetta_mock.init.reset_mock()
        _pyrosetta_mock.pose_from_pdb.reset_mock()
        _toolbox_mock.mutate_residue.reset_mock()

    def test_init_is_called_with_mute_all(self):
        """PyRosetta init() must be called with '-mute all'."""
        _write_single_csv(self.csv_file, [_single_row()])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )
        _pyrosetta_mock.init.assert_called_once_with("-mute all")

    def test_pose_from_pdb_called_per_row(self):
        """pose_from_pdb should be called once for each CSV row processed."""
        rows = [_single_row(sr=i, position=i + 10) for i in range(1, 4)]
        _write_single_csv(self.csv_file, rows)
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=3, folder="mutants"
            )
        self.assertEqual(_pyrosetta_mock.pose_from_pdb.call_count, 3)

    def test_mutate_residue_called_per_row(self):
        """toolbox.mutate_residue must be called once per processed row."""
        rows = [_single_row(sr=i) for i in range(1, 3)]
        _write_single_csv(self.csv_file, rows)
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=2, folder="mutants"
            )
        self.assertEqual(_toolbox_mock.mutate_residue.call_count, 2)

    def test_dump_pdb_called_per_row(self):
        """pose.dump_pdb must be called once per processed row."""
        rows = [_single_row(sr=i) for i in range(1, 3)]
        _write_single_csv(self.csv_file, rows)
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=2, folder="mutants"
            )
        self.assertEqual(self.mock_pose.dump_pdb.call_count, 2)

    def test_output_filename_uses_naming_convention(self):
        """Output PDB name must follow {sr}_{wtAA}-{chain}{position}-{prAA}.pdb."""
        row = _single_row(sr=3, chain="A", position=10, wt="ALA", pr="GLY")
        _write_single_csv(self.csv_file, [row])
        dump_calls = []
        self.mock_pose.dump_pdb.side_effect = lambda path: dump_calls.append(path)

        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )

        self.assertEqual(len(dump_calls), 1)
        filename = Path(dump_calls[0]).name
        self.assertEqual(filename, "3_ALA-A10-GLY.pdb")

    def test_move_file_called_when_folder_given(self):
        """helpers.move_file must be called when folder is provided."""
        _write_single_csv(self.csv_file, [_single_row()])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ) as mock_move:
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )
        mock_move.assert_called_once()

    def test_text_file_written_when_given(self):
        """When text_file is provided, the mutant path is appended to it."""
        _write_single_csv(self.csv_file, [_single_row()])
        text_file = str(self.tmpdir / "paths.txt")
        out_path = str(self.tmpdir / "mutants" / "out.pdb")

        with patch.object(
            generate_mutants.helpers, "move_file", return_value=(self.tmpdir, out_path)
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file,
                self.csv_file,
                total=1,
                folder="mutants",
                text_file=text_file,
            )

        self.assertTrue(Path(text_file).is_file())
        content = Path(text_file).read_text()
        self.assertIn(out_path, content)

    def test_total_n_processes_exactly_n_rows(self):
        """total=N must stop processing after N rows."""
        rows = [_single_row(sr=i, position=i) for i in range(1, 6)]
        _write_single_csv(self.csv_file, rows)

        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_single_mutation(
                self.pdb_file, self.csv_file, total=3, folder="mutants"
            )
        self.assertEqual(_pyrosetta_mock.pose_from_pdb.call_count, 3)


# ---------------------------------------------------------------------------
# Tests: generate_double_mutation()
# ---------------------------------------------------------------------------


class TestGenerateDoubleMutation(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb_file = str(self.tmpdir / "wt.pdb")
        Path(self.pdb_file).touch()
        self.csv_file = str(self.tmpdir / "double.csv")

        self.mock_pose = _make_pose()
        _pyrosetta_mock.pose_from_pdb.return_value = self.mock_pose

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)
        _pyrosetta_mock.init.reset_mock()
        _pyrosetta_mock.pose_from_pdb.reset_mock()
        _toolbox_mock.mutate_residue.reset_mock()

    def test_mutate_residue_called_twice_per_row(self):
        """Each double mutation must apply mutate_residue exactly twice."""
        _write_double_csv(self.csv_file, [_double_row()])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_double_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )
        self.assertEqual(_toolbox_mock.mutate_residue.call_count, 2)

    def test_output_filename_uses_combination_column(self):
        """Output PDB name must be {sr}_{combination}.pdb."""
        row = _double_row(sr=7)
        _write_double_csv(self.csv_file, [row])
        dump_calls = []
        self.mock_pose.dump_pdb.side_effect = lambda p: dump_calls.append(p)

        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_double_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )

        filename = Path(dump_calls[0]).name
        self.assertEqual(filename, f"7_{row['combination']}.pdb")

    def test_total_limits_processing(self):
        rows = [_double_row(sr=i) for i in range(1, 5)]
        _write_double_csv(self.csv_file, rows)
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_double_mutation(
                self.pdb_file, self.csv_file, total=2, folder="mutants"
            )
        # 2 rows × 2 mutations each
        self.assertEqual(_toolbox_mock.mutate_residue.call_count, 4)

    def test_init_is_called(self):
        _write_double_csv(self.csv_file, [_double_row()])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_double_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )
        _pyrosetta_mock.init.assert_called_once_with("-mute all")


# ---------------------------------------------------------------------------
# Tests: generate_triple_mutation()
# ---------------------------------------------------------------------------


class TestGenerateTripleMutation(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb_file = str(self.tmpdir / "wt.pdb")
        Path(self.pdb_file).touch()
        self.csv_file = str(self.tmpdir / "triple.csv")

        self.mock_pose = _make_pose()
        _pyrosetta_mock.pose_from_pdb.return_value = self.mock_pose

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)
        _pyrosetta_mock.init.reset_mock()
        _pyrosetta_mock.pose_from_pdb.reset_mock()
        _toolbox_mock.mutate_residue.reset_mock()

    def test_mutate_residue_called_three_times_per_row(self):
        """Each triple mutation must apply mutate_residue exactly three times."""
        _write_triple_csv(self.csv_file, [_triple_row()])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_triple_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )
        self.assertEqual(_toolbox_mock.mutate_residue.call_count, 3)

    def test_output_filename_uses_combination_column(self):
        row = _triple_row(sr=5)
        _write_triple_csv(self.csv_file, [row])
        dump_calls = []
        self.mock_pose.dump_pdb.side_effect = lambda p: dump_calls.append(p)

        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_triple_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )

        filename = Path(dump_calls[0]).name
        self.assertEqual(filename, f"5_{row['combination']}.pdb")

    def test_total_n_limits_rows(self):
        rows = [_triple_row(sr=i) for i in range(1, 5)]
        _write_triple_csv(self.csv_file, rows)
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_triple_mutation(
                self.pdb_file, self.csv_file, total=2, folder="mutants"
            )
        # 2 rows × 3 mutations each
        self.assertEqual(_toolbox_mock.mutate_residue.call_count, 6)

    def test_all_three_chains_are_used(self):
        """pdb2pose must be queried for all three mutation sites."""
        row = _triple_row()
        _write_triple_csv(self.csv_file, [row])
        with patch.object(
            generate_mutants.helpers,
            "move_file",
            return_value=(self.tmpdir, "/tmp/out.pdb"),
        ):
            generate_mutants.generate_triple_mutation(
                self.pdb_file, self.csv_file, total=1, folder="mutants"
            )

        pdb2pose_calls = self.mock_pose.pdb_info.return_value.pdb2pose.call_args_list
        queried_chains = {c.args[0] for c in pdb2pose_calls}
        self.assertIn("A", queried_chains)
        self.assertIn("B", queried_chains)
        self.assertIn("C", queried_chains)


if __name__ == "__main__":
    unittest.main()
