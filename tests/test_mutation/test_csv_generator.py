"""
Tests for mutadock.mutation.csv_generator
---------------------------------
Covers: get_residues(), generate_csv(), and PAM250 filtering logic.

Run from the project root:
    pytest tests/
"""

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Path setup — add src/ so "mutadock.mutation.*" imports resolve without installation.
# Bio / PyRosetta / tqdm stubs are already injected by conftest.py before
# this module is imported, so no additional stub setup is needed here.
# ---------------------------------------------------------------------------
from mutadock.mutation.Amino import get_dict
from mutadock.mutation.csv_generator import generate_csv, get_residues
from mutadock.mutation.exceptions import CSVGenerationError, PDBFileError
from mutadock.mutation.helpers import DATA_DIR, load_matrix

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_mock_residue(chain: str, position: int, resname: str) -> MagicMock:
    """Return a BioPython Residue-like mock."""
    r = MagicMock()
    r.get_full_id.return_value = ("struct", 0, chain, (" ", position, " "))
    r.get_resname.return_value = resname
    return r


def _make_mock_structure(residues: list) -> MagicMock:
    """Wrap a flat list of residue mocks into a structure → model → chain tree."""
    mock_chain = MagicMock()
    mock_chain.__iter__ = MagicMock(return_value=iter(residues))
    mock_model = MagicMock()
    mock_model.__iter__ = MagicMock(return_value=iter([mock_chain]))
    mock_structure = MagicMock()
    mock_structure.__iter__ = MagicMock(return_value=iter([mock_model]))
    return mock_structure


# ---------------------------------------------------------------------------
# Tests: get_residues()
# ---------------------------------------------------------------------------


class TestGetResidues(unittest.TestCase):

    def test_raises_for_nonexistent_file(self):
        """FileNotFoundError from PDBParser.get_structure → PDBFileError."""
        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.side_effect = FileNotFoundError
            with self.assertRaises(PDBFileError):
                get_residues("/nonexistent/path/fake_protein.pdb")

    def test_returns_dict_for_valid_structure(self):
        """Two residues in PDB → both stored (no residue is dropped)."""
        residues_in = [
            _make_mock_residue("A", 1, "ALA"),
            _make_mock_residue("A", 2, "GLY"),
        ]
        mock_structure = _make_mock_structure(residues_in)

        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.return_value = mock_structure
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                tmp_path = f.name
            try:
                result = get_residues(tmp_path)
            finally:
                Path(tmp_path).unlink()

        self.assertIsNotNone(result)
        self.assertIsInstance(result, dict)
        # Both residues end up in the dict — the first is no longer dropped.
        self.assertEqual(len(result), 2)

    def test_first_residue_is_included(self):
        """The N-terminal residue of the first chain must be present (regression: 2.3)."""
        residues_in = [
            _make_mock_residue("A", 1, "ALA"),  # first residue — must NOT be dropped
            _make_mock_residue("A", 2, "GLY"),
            _make_mock_residue("A", 3, "VAL"),
        ]
        mock_structure = _make_mock_structure(residues_in)

        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.return_value = mock_structure
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                tmp_path = f.name
            try:
                result = get_residues(tmp_path)
            finally:
                Path(tmp_path).unlink()

        self.assertEqual(len(result), 3)
        # First stored residue is the ALA at position 1.
        first_key = min(result.keys())
        self.assertEqual(result[first_key], ("A", 1, "ALA"))
        # Numbering still starts at 1.
        self.assertEqual(first_key, 1)
        # Every residue is represented.
        names = {v[2] for v in result.values()}
        self.assertEqual(names, {"ALA", "GLY", "VAL"})

    def test_residue_tuple_is_chain_position_name(self):
        """Each dict value is (chain, position, resname)."""
        residues_in = [
            _make_mock_residue("B", 10, "LYS"),
            _make_mock_residue("B", 11, "TRP"),
        ]
        mock_structure = _make_mock_structure(residues_in)

        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.return_value = mock_structure
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                tmp_path = f.name
            try:
                result = get_residues(tmp_path)
            finally:
                Path(tmp_path).unlink()

        self.assertEqual(len(result), 2)
        first_key = min(result.keys())
        self.assertEqual(result[first_key], ("B", 10, "LYS"))
        last_key = max(result.keys())
        chain, position, name = result[last_key]
        self.assertEqual(chain, "B")
        self.assertEqual(position, 11)
        self.assertEqual(name, "TRP")

    def test_empty_structure_returns_empty_dict(self):
        """A structure with no residues → empty dict (not None)."""
        mock_structure = _make_mock_structure([])

        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.return_value = mock_structure
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                tmp_path = f.name
            try:
                result = get_residues(tmp_path)
            finally:
                Path(tmp_path).unlink()

        self.assertIsNotNone(result)
        self.assertEqual(len(result), 0)

    def test_multiple_chains(self):
        """Residues across two chains are all collected."""
        chain_a_res = [
            _make_mock_residue("A", 1, "ALA"),
            _make_mock_residue("A", 2, "SER"),
        ]
        chain_b_res = [
            _make_mock_residue("B", 1, "PHE"),
            _make_mock_residue("B", 2, "ILE"),
        ]

        mock_chain_a = MagicMock()
        mock_chain_a.__iter__ = MagicMock(return_value=iter(chain_a_res))
        mock_chain_b = MagicMock()
        mock_chain_b.__iter__ = MagicMock(return_value=iter(chain_b_res))
        mock_model = MagicMock()
        mock_model.__iter__ = MagicMock(return_value=iter([mock_chain_a, mock_chain_b]))
        mock_structure = MagicMock()
        mock_structure.__iter__ = MagicMock(return_value=iter([mock_model]))

        with patch("mutadock.mutation.csv_generator.PDBParser") as MockParser:
            MockParser.return_value.get_structure.return_value = mock_structure
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                tmp_path = f.name
            try:
                result = get_residues(tmp_path)
            finally:
                Path(tmp_path).unlink()

        # 4 residues across both chains → all 4 stored
        self.assertEqual(len(result), 4)


# ---------------------------------------------------------------------------
# Tests: generate_csv()
# ---------------------------------------------------------------------------


class TestGenerateCsv(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb_file = str(self.tmpdir / "test.pdb")
        # Create a dummy PDB so backup() doesn't crash on a missing file
        Path(self.pdb_file).touch()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # --- helpers ---

    def _one_residue_dict(self, chain="A", position=5, resname="ALA"):
        """Minimal residues dict: one entry (as if first was already skipped)."""
        return {1: (chain, position, resname)}

    def _read_csv(self, path):
        with open(path, newline="") as f:
            return list(csv.DictReader(f))

    # --- tests ---

    def test_raises_when_get_residues_returns_none(self):
        with (
            patch("mutadock.mutation.csv_generator.get_residues", return_value=None),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            with self.assertRaises(CSVGenerationError):
                generate_csv(self.pdb_file)

    def test_raises_when_residues_empty(self):
        with (
            patch("mutadock.mutation.csv_generator.get_residues", return_value={}),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            with self.assertRaises(CSVGenerationError):
                generate_csv(self.pdb_file)

    def test_returns_tuple_of_two_paths(self):
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            result = generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)

    def test_creates_both_output_files(self):
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)
        self.assertTrue(Path(out_op).is_file())
        self.assertTrue(Path(out_all).is_file())

    def test_output_files_have_correct_headers(self):
        expected_headers = [
            "sr",
            "pdb",
            "chain",
            "position",
            "wtAA",
            "prAA",
            "wtProb",
            "prProb",
        ]
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        for path in (out_op, out_all):
            with open(path, newline="") as f:
                reader = csv.reader(f)
                headers = next(reader)
            self.assertEqual(headers, expected_headers, msg=f"Wrong headers in {path}")

    def test_all_file_contains_every_mutation(self):
        """_all.csv must contain 19 rows for one residue (20 aa minus itself)."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(resname="ALA"),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)
        rows = self._read_csv(out_all)
        self.assertEqual(len(rows), 19)

    def test_op_file_only_contains_positive_pam250(self):
        """_mutations.csv keeps only rows where PAM250 score > 0."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(resname="ALA"),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        rows = self._read_csv(out_op)
        score_dict = load_matrix(DATA_DIR / "PAM250")
        for row in rows:
            wt = row["wtAA"]
            pr = row["prAA"]
            self.assertGreater(
                score_dict[wt][pr],
                0,
                msg=f"{wt}→{pr} has score ≤ 0 but appeared in op file",
            )

    def test_no_self_mutation_in_output(self):
        """A residue must never appear as its own mutation target."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(resname="GLY"),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        for path in (out_op, out_all):
            for row in self._read_csv(path):
                self.assertNotEqual(
                    row["wtAA"], row["prAA"], msg=f"Self-mutation found in {path}"
                )

    def test_default_output_naming(self):
        """Without explicit out_* args, filenames are derived from pdb_file."""
        pdb = str(self.tmpdir / "myprotein.pdb")
        Path(pdb).touch()
        expected_op = str(self.tmpdir / "myprotein_mutations.csv")
        expected_all = str(self.tmpdir / "myprotein_mutations_all.csv")

        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            result = generate_csv(pdb)

        self.assertIsNotNone(result)
        op_path, all_path = result
        self.assertEqual(op_path, expected_op)
        self.assertEqual(all_path, expected_all)

    def test_same_name_collision_gets_all_suffix(self):
        """If out_all == out_op, out_all is renamed to end with _all.csv."""
        same_path = str(self.tmpdir / "same.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            result = generate_csv(self.pdb_file, out_all=same_path, out_op=same_path)

        self.assertIsNotNone(result)
        op_path, all_path = result
        self.assertNotEqual(op_path, all_path)
        self.assertTrue(all_path.endswith("_all.csv"))

    def test_wtprob_matches_pam250_diagonal(self):
        """wtProb must equal PAM250[wt][wt] / 100."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(resname="LEU"),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        score_dict = load_matrix(DATA_DIR / "PAM250")
        expected_wt_prob = score_dict["LEU"]["LEU"] / 100.0
        for row in self._read_csv(out_all):
            self.assertAlmostEqual(float(row["wtProb"]), expected_wt_prob)

    def test_chain_and_position_propagate_correctly(self):
        """chain/position in CSV rows must match the residue dict."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        residues = {1: ("C", 42, "SER")}
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues", return_value=residues
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        for row in self._read_csv(out_all):
            self.assertEqual(row["chain"], "C")
            self.assertEqual(row["position"], "42")

    def test_serial_numbers_are_sequential(self):
        """sr column in _all.csv must be 1, 2, 3, … without gaps."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        rows = self._read_csv(out_all)
        sr_values = [int(r["sr"]) for r in rows]
        self.assertEqual(sr_values, list(range(1, len(sr_values) + 1)))

    def test_backup_is_called_for_both_outputs(self):
        """backup() must be invoked for both output paths."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues",
                return_value=self._one_residue_dict(),
            ),
            patch("mutadock.mutation.csv_generator.backup") as mock_backup,
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        called_paths = {call.args[0] for call in mock_backup.call_args_list}
        self.assertIn(out_op, called_paths)
        self.assertIn(out_all, called_paths)

    def test_nonstandard_residue_does_not_raise(self):
        """A residue absent from the matrix (e.g. MSE/HOH) is skipped, not fatal (2.9)."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        residues = {1: ("A", 1, "MSE"), 2: ("A", 2, "HOH")}
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues", return_value=residues
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            # Must not raise KeyError despite non-standard residue names,
            # and a warning must be logged for the skipped residue.
            with self.assertLogs("mutadock.mutation.csv_generator", level="WARNING"):
                result = generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        self.assertIsNotNone(result)
        # No standard mutations could be enumerated → only the header exists.
        self.assertEqual(self._read_csv(out_all), [])
        self.assertEqual(self._read_csv(out_op), [])

    def test_nonstandard_residue_does_not_corrupt_serial_numbers(self):
        """A skipped non-standard residue must not create gaps in the sr column."""
        out_op = str(self.tmpdir / "op.csv")
        out_all = str(self.tmpdir / "all.csv")
        # A valid residue, a non-standard one, then another valid residue.
        residues = {
            1: ("A", 1, "ALA"),
            2: ("A", 2, "MSE"),  # skipped
            3: ("A", 3, "GLY"),
        }
        with (
            patch(
                "mutadock.mutation.csv_generator.get_residues", return_value=residues
            ),
            patch("mutadock.mutation.csv_generator.backup"),
        ):
            generate_csv(self.pdb_file, out_all=out_all, out_op=out_op)

        rows = self._read_csv(out_all)
        # Two valid residues × 19 substitutions each = 38 rows, no MSE rows.
        self.assertEqual(len(rows), 38)
        self.assertNotIn("MSE", {r["wtAA"] for r in rows})
        # Serial numbers remain contiguous 1..38 with no gap from the skip.
        sr_values = [int(r["sr"]) for r in rows]
        self.assertEqual(sr_values, list(range(1, len(sr_values) + 1)))


# ---------------------------------------------------------------------------
# Tests: PAM250 matrix sanity (load_matrix)
# ---------------------------------------------------------------------------


class TestPAM250Matrix(unittest.TestCase):

    def test_matrix_covers_all_20_amino_acids(self):
        score_dict = load_matrix(DATA_DIR / "PAM250")
        aa_dict = get_dict()
        self.assertEqual(set(score_dict.keys()), set(aa_dict.keys()))

    def test_diagonal_is_positive(self):
        """Each amino acid scores highest against itself."""
        score_dict = load_matrix(DATA_DIR / "PAM250")
        for aa, row in score_dict.items():
            self.assertGreater(row[aa], 0, msg=f"Diagonal for {aa} is not positive")

    def test_matrix_is_symmetric(self):
        """PAM250 should be symmetric: score[a][b] == score[b][a]."""
        score_dict = load_matrix(DATA_DIR / "PAM250")
        for aa1, row in score_dict.items():
            for aa2, val in row.items():
                self.assertEqual(
                    val, score_dict[aa2][aa1], msg=f"Asymmetry: {aa1}↔{aa2}"
                )


if __name__ == "__main__":
    unittest.main()
