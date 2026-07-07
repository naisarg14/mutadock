"""
Tests for per-item checkpoint/resume of the ΔΔG calculators
-----------------------------------------------------------
Bio / PyRosetta / tqdm stubs are injected by conftest.py before this module is
imported, so the heavyweight optional dependencies need not be installed.

Covers:
  * ``helpers.read_partial_ddg`` guards (absent / header-only / torn trailing
    line / protocol mismatch / replicates mismatch).
  * ``calc_ddg`` / ``calc_double_ddg`` / ``calc_triple_ddg`` resume: after an
    interrupted run is truncated to a prefix, resuming reproduces the exact same
    row *set*, ``sr`` numbering, and ``combination`` names as a from-scratch run
    (the stochastic ΔΔG floats are not compared — only identity/ordering).

Run from the project root:
    pytest tests/test_mutation/test_ddg_resume.py -q
"""

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from mutadock.mutation import ddg_calc, ddg_calc_double, ddg_calc_triple
from mutadock.mutation.helpers import read_partial_ddg

_HDR = (
    "sr,pdb,chain,position,wtAA,prAA,wtProb,prProb,"
    "ddG_value,ddG_sd,ddG_kcal,n_replicates,ddG_protocol"
)


# ---------------------------------------------------------------------------
# read_partial_ddg guards (pure — no pyrosetta)
# ---------------------------------------------------------------------------


class TestReadPartialDdg(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, lines):
        p = self.tmp / name
        p.write_text("\n".join(lines) + "\n")
        return str(p)

    def test_absent_file_returns_none(self):
        self.assertIsNone(read_partial_ddg(str(self.tmp / "nope.csv")))

    def test_header_only_returns_none(self):
        p = self._write("h.csv", [_HDR])
        self.assertIsNone(read_partial_ddg(p))

    def test_torn_trailing_line_dropped(self):
        """A complete run leaves N full rows; a mid-write crash may leave a
        short final line.  That torn row is skipped (recomputed), not returned."""
        p = self._write(
            "t.csv",
            [
                _HDR,
                "1,x,A,10,ALA,GLY,0.02,0.01,-1.0,0.0,-0.34,1,min",
                "2,x,A,11,ALA,SER,0.02,0.01,-2.0,0.0,-0.68,1,min",
                "3,x,A,12,ALA,THR,0.02,0.01",  # torn: missing final columns
            ],
        )
        rows = read_partial_ddg(
            p, must_match={"ddG_protocol": "min", "n_replicates": 1}
        )
        self.assertIsNotNone(rows)
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["position"] for r in rows}, {"10", "11"})

    def test_protocol_mismatch_returns_none(self):
        p = self._write(
            "p.csv",
            [_HDR, "1,x,A,10,ALA,GLY,0.02,0.01,-1.0,0.0,-0.34,1,min"],
        )
        self.assertIsNone(read_partial_ddg(p, must_match={"ddG_protocol": "cartesian"}))

    def test_replicates_mismatch_returns_none(self):
        p = self._write(
            "r.csv",
            [_HDR, "1,x,A,10,ALA,GLY,0.02,0.01,-1.0,0.0,-0.34,1,min"],
        )
        self.assertIsNone(read_partial_ddg(p, must_match={"n_replicates": 3}))


# ---------------------------------------------------------------------------
# Resume identity for the three calculators
# ---------------------------------------------------------------------------


def _mock_pose():
    pose = MagicMock(name="pose")
    info = MagicMock()
    info.pdb2pose.side_effect = lambda c, p: (hash((c, int(p))) % 997) + 1
    pose.pdb_info.return_value = info
    return pose


def _patches(module):
    sfxn = MagicMock(name="sfxn")
    state = {"n": 0}

    def _score(_p):
        state["n"] += 1
        return float(state["n"])

    sfxn.score.side_effect = _score
    return [
        patch.object(module, "init"),
        patch.object(module, "pose_from_pdb", return_value=_mock_pose()),
        patch.object(module.predict_ddG, "get_scorefxn", return_value=sfxn),
        patch.object(module.predict_ddG, "wt_reference_score", return_value=0.0),
        patch.object(module.predict_ddG, "apply_mutations", return_value=MagicMock()),
        patch.object(module, "in_directory"),
    ]


def _run(module, fn, *args, **kwargs):
    ps = _patches(module)
    for p in ps:
        p.start()
    try:
        return fn(*args, **kwargs)
    finally:
        for p in ps:
            p.stop()


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _truncate(src, dst, keep):
    rows = _rows(src)
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows[:keep]:
            w.writerow(r)


class TestResumeIdentity(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.pdb = self.tmp / "wt.pdb"
        self.pdb.touch()
        self.single_in = self.tmp / "muts.csv"
        aa = ["GLY", "SER", "THR", "VAL", "LEU", "ILE"]
        with open(self.single_in, "w", newline="") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
                    "sr",
                    "pdb",
                    "chain",
                    "position",
                    "wtAA",
                    "prAA",
                    "wtProb",
                    "prProb",
                ],
            )
            w.writeheader()
            for i in range(6):
                w.writerow(
                    {
                        "sr": i + 1,
                        "pdb": "wt",
                        "chain": "A",
                        "position": 10 + i,
                        "wtAA": "ALA",
                        "prAA": aa[i],
                        "wtProb": 0.02,
                        "prProb": 0.01,
                    }
                )

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _assert_same(self, scratch, resumed, id_cols):
        a, b = _rows(scratch), _rows(resumed)
        self.assertEqual(len(a), len(b))
        for ra, rb in zip(a, b, strict=True):
            self.assertEqual(ra["sr"], rb["sr"])
            self.assertEqual(
                tuple(ra[c] for c in id_cols), tuple(rb[c] for c in id_cols)
            )
            if "combination" in ra:
                self.assertEqual(ra["combination"], rb["combination"])

    def test_single_resume_matches_scratch(self):
        scratch = str(self.tmp / "s_scratch.csv")
        resumed = str(self.tmp / "s_resumed.csv")
        _run(
            ddg_calc,
            ddg_calc.calc_ddg,
            str(self.pdb),
            str(self.single_in),
            out_file=scratch,
        )
        _truncate(scratch, resumed, keep=4)
        _run(
            ddg_calc,
            ddg_calc.calc_ddg,
            str(self.pdb),
            str(self.single_in),
            out_file=resumed,
            resume=True,
        )
        self._assert_same(scratch, resumed, ["chain", "position", "prAA"])

    def test_double_resume_matches_scratch(self):
        scratch = str(self.tmp / "d_scratch.csv")
        resumed = str(self.tmp / "d_resumed.csv")
        _run(
            ddg_calc_double,
            ddg_calc_double.calc_double_ddg,
            str(self.pdb),
            single_csv=str(self.single_in),
            out_csv=scratch,
            total=6,
        )
        _truncate(scratch, resumed, keep=5)
        _run(
            ddg_calc_double,
            ddg_calc_double.calc_double_ddg,
            str(self.pdb),
            single_csv=str(self.single_in),
            out_csv=resumed,
            total=6,
            resume=True,
        )
        self._assert_same(
            scratch,
            resumed,
            [
                "chain1",
                "mut1_position",
                "mut1_prAA",
                "chain2",
                "mut2_position",
                "mut2_prAA",
            ],
        )

    def test_triple_resume_matches_scratch(self):
        double = str(self.tmp / "d_for_triple.csv")
        _run(
            ddg_calc_double,
            ddg_calc_double.calc_double_ddg,
            str(self.pdb),
            single_csv=str(self.single_in),
            out_csv=double,
            total=6,
        )
        scratch = str(self.tmp / "t_scratch.csv")
        resumed = str(self.tmp / "t_resumed.csv")
        _run(
            ddg_calc_triple,
            ddg_calc_triple.calc_triple_ddg,
            str(self.pdb),
            double_csv=double,
            out_csv=scratch,
            total=6,
        )
        _truncate(scratch, resumed, keep=max(1, len(_rows(scratch)) // 2))
        _run(
            ddg_calc_triple,
            ddg_calc_triple.calc_triple_ddg,
            str(self.pdb),
            double_csv=double,
            out_csv=resumed,
            total=6,
            resume=True,
        )
        self._assert_same(
            scratch,
            resumed,
            [
                "chain1",
                "mut1_position",
                "mut1_prAA",
                "chain2",
                "mut2_position",
                "mut2_prAA",
                "chain3",
                "mut3_position",
                "mut3_prAA",
            ],
        )

    def test_resume_without_flag_recomputes_from_scratch(self):
        """resume=False (default) must ignore any existing file and overwrite
        it via backup — preserving the pre-feature behaviour."""
        out = str(self.tmp / "s.csv")
        _run(
            ddg_calc,
            ddg_calc.calc_ddg,
            str(self.pdb),
            str(self.single_in),
            out_file=out,
        )
        with patch("mutadock.mutation.ddg_calc.backup") as mock_backup:
            _run(
                ddg_calc,
                ddg_calc.calc_ddg,
                str(self.pdb),
                str(self.single_in),
                out_file=out,
            )  # resume defaults to False
        mock_backup.assert_called()  # fresh run still backs up the old file


if __name__ == "__main__":
    unittest.main()
