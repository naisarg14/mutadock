"""
Tests for mutadock.mutation.np_mutation
---------------------------------------
The heavy optional dependencies (BioPython, PyRosetta, tqdm) are stubbed by
conft.py before this module is imported, so np_mutation can be imported without
the real packages installed.

Run from the project root:
    pytest tests/test_mutation/test_np_mutation.py
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mutadock.mutation import np_mutation

# ---------------------------------------------------------------------------
# Flow tests for the double/triple-ddG "skip if exists" predicate (bug 2.7)
# ---------------------------------------------------------------------------
#
# The skip check for the *sorted triple* file must compare the file's line
# count against the number of rows a COMPLETE triple-ddG CSV would contain
# (``combinations(unique_singles, 3)`` + a header) — not a permutation count,
# and not the double count.  These tests drive the full ``np_mutation()``
# orchestrator with every pipeline stage mocked and use the number of
# ``sort_csv`` calls as the discriminator:
#
#   * single-sort + double-sort always run       -> 2 calls
#   * triple-sort runs only when NOT skipped      -> +1 call (3 total)
#
# ``_expected_double_lines`` / ``_expected_triple_lines`` are patched to known,
# distinct values so the predicate's wiring (skips only on an exact match to
# the TRIPLE count) is observable without touching the filesystem.

EXPECTED_DOUBLE_LINES = 12
EXPECTED_TRIPLE_LINES = 6


class TestTripleSkipPredicate(unittest.TestCase):
    # Chosen so that P(double,2) != P(triple,3).
    NUM_DOUBLE_DDG = 4  # P(4, 2) = 12
    NUM_TRIPLE_DDG = 3  # P(3, 3) = 6

    def setUp(self):
        # Pure string path; no filesystem access happens because file_info is
        # patched.  base_name derived by np_mutation() will be "/work/prot".
        self.pdb = "/work/prot.pdb"
        self.base_name = "/work/prot"
        self.triple_sorted = f"{self.base_name}_triple_ddg_sorted.csv"

    def _run(self, triple_sorted_count: int) -> int:
        """Drive np_mutation() and return the number of sort_csv() calls.

        Only the sorted-triple file is reported as already existing (with
        ``triple_sorted_count`` lines); every other file is reported absent, so
        all earlier stages execute normally.
        """
        inputs = (
            self.pdb,  # full_pdb_path
            self.NUM_DOUBLE_DDG,  # num_double_ddg
            15,  # num_single_mut
            15,  # num_double_mut
            self.NUM_TRIPLE_DDG,  # num_triple_ddg
            15,  # num_triple_mut
            True,  # append
            True,  # quiet
            None,  # output_dir
            True,  # no_report
            {},  # ddg_params
        )

        def fake_file_info(path):
            if path == self.triple_sorted:
                return (True, triple_sorted_count)
            return (False, 0)

        with (
            patch.object(np_mutation, "get_inputs", return_value=inputs),
            patch.object(np_mutation, "file_info", side_effect=fake_file_info),
            patch.object(np_mutation, "clean_pdb"),
            patch.object(
                np_mutation,
                "generate_csv",
                return_value=(f"{self.base_name}_mutations.csv", None),
            ),
            patch.object(
                np_mutation, "calc_ddg", return_value=f"{self.base_name}_ddG.csv"
            ),
            patch.object(
                np_mutation,
                "calc_double_ddg",
                return_value=f"{self.base_name}_double_ddg.csv",
            ),
            patch.object(
                np_mutation,
                "calc_triple_ddg",
                return_value=f"{self.base_name}_triple_ddg.csv",
            ),
            patch.object(
                np_mutation,
                "sort_csv",
                side_effect=lambda in_file, *a, **k: (
                    k.get("out_file")
                    or f"{str(in_file).removesuffix('.csv')}_sorted.csv"
                ),
            ) as mock_sort,
            patch.object(np_mutation, "generate_single_mutation"),
            patch.object(np_mutation, "generate_double_mutation"),
            patch.object(np_mutation, "generate_triple_mutation"),
            patch.object(np_mutation, "backup"),
            patch.object(
                np_mutation,
                "_expected_double_lines",
                return_value=EXPECTED_DOUBLE_LINES,
            ),
            patch.object(
                np_mutation,
                "_expected_triple_lines",
                return_value=EXPECTED_TRIPLE_LINES,
            ),
        ):
            np_mutation.np_mutation()

        return mock_sort.call_count

    def test_skip_fires_when_count_equals_triple_combination(self):
        """Sorted-triple file with the exact triple-combination line count ->
        triple sort is skipped."""
        count = self._run(EXPECTED_TRIPLE_LINES)
        # single + double sort ran; triple sort was skipped.
        self.assertEqual(count, 2)

    def test_skip_does_not_fire_when_count_equals_double_count(self):
        """Sorted-triple file whose line count matches the DOUBLE count ->
        triple sort still runs.

        This is the case the original bug got wrong: it compared the triple
        file against the double count and would have wrongly skipped.
        """
        count = self._run(EXPECTED_DOUBLE_LINES)
        # single + double + triple sort all ran.
        self.assertEqual(count, 3)


class TestResumeWiring(unittest.TestCase):
    """np_mutation must thread ``resume=append`` into all three ΔΔG calcs, so an
    interrupted run resumes per-item (append=True) rather than recomputing."""

    def _run_and_capture(self, append: bool):
        base = "/work/prot"
        inputs = (f"{base}.pdb", 4, 15, 15, 3, 15, append, True, None, True, {})
        # Every output absent -> every calc stage executes (nothing skipped).
        with (
            patch.object(np_mutation, "get_inputs", return_value=inputs),
            patch.object(np_mutation, "file_info", return_value=(False, 0)),
            patch.object(np_mutation, "clean_pdb"),
            patch.object(
                np_mutation,
                "generate_csv",
                return_value=(f"{base}_mutations.csv", None),
            ),
            patch.object(
                np_mutation, "calc_ddg", return_value=f"{base}_ddG.csv"
            ) as m_single,
            patch.object(
                np_mutation, "calc_double_ddg", return_value=f"{base}_double_ddg.csv"
            ) as m_double,
            patch.object(
                np_mutation, "calc_triple_ddg", return_value=f"{base}_triple_ddg.csv"
            ) as m_triple,
            patch.object(
                np_mutation,
                "sort_csv",
                side_effect=lambda in_file, *a, **k: (
                    k.get("out_file")
                    or f"{str(in_file).removesuffix('.csv')}_sorted.csv"
                ),
            ),
            patch.object(np_mutation, "generate_single_mutation"),
            patch.object(np_mutation, "generate_double_mutation"),
            patch.object(np_mutation, "generate_triple_mutation"),
            patch.object(np_mutation, "backup"),
            patch.object(np_mutation, "_expected_double_lines", return_value=-1),
            patch.object(np_mutation, "_expected_triple_lines", return_value=-1),
        ):
            np_mutation.np_mutation()
        return m_single, m_double, m_triple

    def test_resume_true_passed_when_append(self):
        for m in self._run_and_capture(append=True):
            self.assertIs(m.call_args.kwargs.get("resume"), True, m._mock_name)

    def test_resume_false_passed_when_no_append(self):
        for m in self._run_and_capture(append=False):
            self.assertIs(m.call_args.kwargs.get("resume"), False, m._mock_name)


# ---------------------------------------------------------------------------
# Argparse flag test for --no-append -> append (bug 2.16)
# ---------------------------------------------------------------------------


class TestNoAppendFlag(unittest.TestCase):

    _APPEND_INDEX = 6  # position of `append` in the get_inputs() return tuple

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.pdb = self.tmpdir / "prot.pdb"
        self.pdb.touch()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _get_append(self, argv: list[str]) -> bool:
        with patch.object(sys, "argv", ["np_mutation", *argv]):
            result = np_mutation.get_inputs()
        return result[self._APPEND_INDEX]

    def test_append_true_by_default(self):
        """Without the flag, append defaults to True (reuse existing files)."""
        self.assertIs(self._get_append(["-i", str(self.pdb)]), True)

    def test_no_append_flag_sets_append_false(self):
        """Passing --no-append sets append to False (regenerate everything)."""
        self.assertIs(self._get_append(["-i", str(self.pdb), "--no-append"]), False)


# ---------------------------------------------------------------------------
# 3.4: --pdb-id fetches from RCSB instead of -i
# ---------------------------------------------------------------------------


class TestPdbIdArg(unittest.TestCase):

    _PATH_INDEX = 0  # position of the resolved PDB path in get_inputs()'s tuple

    def _get_inputs(self, argv: list[str]):
        with patch.object(sys, "argv", ["np_mutation", *argv]):
            return np_mutation.get_inputs()

    def test_pdb_id_resolves_to_fetched_path(self):
        """--pdb-id 4QJR fetches from RCSB and returns the downloaded path."""
        with patch.object(
            np_mutation, "fetch_pdb", return_value="/downloads/4QJR.pdb"
        ) as mock_fetch:
            result = self._get_inputs(["--pdb-id", "4QJR"])
        mock_fetch.assert_called_once_with("4QJR")
        self.assertEqual(result[self._PATH_INDEX], "/downloads/4QJR.pdb")

    def test_input_and_pdb_id_are_mutually_exclusive(self):
        """Supplying both -i and --pdb-id is an argparse error (exit)."""
        with self.assertRaises(SystemExit):
            self._get_inputs(["-i", "prot.pdb", "--pdb-id", "4QJR"])

    def test_one_of_input_or_pdb_id_is_required(self):
        """Supplying neither -i nor --pdb-id is an argparse error (exit)."""
        with self.assertRaises(SystemExit):
            self._get_inputs([])


# ---------------------------------------------------------------------------
# CIF input is converted to PDB before the rest of the pipeline runs
# ---------------------------------------------------------------------------


class TestCifConversion(unittest.TestCase):
    """The pipeline only understands PDB; a .cif must be converted up front
    rather than text-cleaned as if it were already PDB."""

    def _drive(self, input_path: str, convert_return: str = "/work/prot.pdb"):
        """Run np_mutation() with all stages mocked.

        Returns ``(convert_mock, clean_mock)`` for assertions.
        """
        inputs = (input_path, 4, 15, 15, 3, 15, True, True, None, True, {})
        with (
            patch.object(np_mutation, "get_inputs", return_value=inputs),
            patch.object(np_mutation, "file_info", return_value=(False, 0)),
            patch.object(
                np_mutation, "convert_cif_pdb", return_value=convert_return
            ) as mock_conv,
            patch.object(np_mutation, "structure_warnings", return_value=[]),
            patch.object(np_mutation, "clean_pdb") as mock_clean,
            patch.object(
                np_mutation,
                "generate_csv",
                return_value=("/work/prot_mutations.csv", None),
            ),
            patch.object(np_mutation, "calc_ddg", return_value="/work/prot_ddG.csv"),
            patch.object(
                np_mutation,
                "calc_double_ddg",
                return_value="/work/prot_double_ddg.csv",
            ),
            patch.object(
                np_mutation,
                "calc_triple_ddg",
                return_value="/work/prot_triple_ddg.csv",
            ),
            patch.object(
                np_mutation,
                "sort_csv",
                side_effect=lambda in_file, *a, **k: (
                    k.get("out_file")
                    or f"{str(in_file).removesuffix('.csv')}_sorted.csv"
                ),
            ),
            patch.object(np_mutation, "generate_single_mutation"),
            patch.object(np_mutation, "generate_double_mutation"),
            patch.object(np_mutation, "generate_triple_mutation"),
            patch.object(np_mutation, "backup"),
        ):
            np_mutation.np_mutation()
        return mock_conv, mock_clean

    def test_cif_input_is_converted_before_cleaning(self):
        mock_conv, mock_clean = self._drive("/work/prot.cif")
        # convert_cif_pdb was called with the .cif input
        mock_conv.assert_called_once()
        self.assertEqual(mock_conv.call_args[0][0], "/work/prot.cif")
        # clean_pdb operated on the converted .pdb, never the raw .cif
        clean_input = str(mock_clean.call_args[0][0])
        self.assertTrue(clean_input.endswith(".pdb"))
        self.assertNotIn(".cif", clean_input)

    def test_pdb_input_is_not_converted(self):
        mock_conv, mock_clean = self._drive("/work/prot.pdb")
        mock_conv.assert_not_called()
        self.assertEqual(str(mock_clean.call_args[0][0]), "/work/prot.pdb")


if __name__ == "__main__":
    unittest.main()
