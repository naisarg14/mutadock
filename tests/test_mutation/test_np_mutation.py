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
from mutadock.mutation.helpers import permutations

# ---------------------------------------------------------------------------
# Flow tests for the triple-ddG "skip if exists" predicate (bug 2.7)
# ---------------------------------------------------------------------------
#
# The buggy skip check for the *sorted triple* file compared the file's line
# count against ``permutations(num_double_ddg, 2)`` instead of
# ``permutations(num_triple_ddg, 3)``.  These tests drive the full
# ``np_mutation()`` orchestrator with every pipeline stage mocked, and use the
# number of ``sort_csv`` calls as the discriminator:
#
#   * single-sort + double-sort always run       -> 2 calls
#   * triple-sort runs only when NOT skipped      -> +1 call (3 total)
#
# With num_double_ddg=4 (P(4,2)=12) and num_triple_ddg=3 (P(3,3)=6) the two
# expected counts differ, so the double- vs triple-count comparison is
# observable.


# Sanity: the two permutation counts must differ for the discriminator to work.
assert permutations(4, 2) != permutations(3, 3)  # 12 != 6


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
        ):
            np_mutation.np_mutation()

        return mock_sort.call_count

    def test_skip_fires_when_count_equals_triple_permutation(self):
        """Sorted-triple file with P(triple, 3) lines -> triple sort skipped."""
        count = self._run(permutations(self.NUM_TRIPLE_DDG, 3))  # 6
        # single + double sort ran; triple sort was skipped.
        self.assertEqual(count, 2)

    def test_skip_does_not_fire_when_count_equals_double_permutation(self):
        """Sorted-triple file with P(double, 2) lines -> triple sort still runs.

        This is the case the original bug got wrong: it would have (incorrectly)
        skipped the triple sort because it compared against the double count.
        """
        count = self._run(permutations(self.NUM_DOUBLE_DDG, 2))  # 12
        # single + double + triple sort all ran.
        self.assertEqual(count, 3)


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


if __name__ == "__main__":
    unittest.main()
