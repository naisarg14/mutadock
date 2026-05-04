"""
Tests for docking.np_docking
------------------------------
All external docking helpers (prepare_receptor, prepare_ligand, dock_vina,
vina_split, add_score_to_csv, backup, read_config) are patched in
np_docking's own namespace so the real vina/meeko/AutoDockTools packages
need not be installed.

Run from the project root:
    pytest tests/
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from docking import np_docking

# ---------------------------------------------------------------------------
# Helpers for prepare_inputs() — patch sys.argv
# ---------------------------------------------------------------------------


def _argv(*extra):
    """Build a sys.argv list for np_docking's argparse."""
    return ["np_dock"] + list(extra)


# ===========================================================================
# TestPrepareInputs
# ===========================================================================


class TestPrepareInputs(unittest.TestCase):
    """Tests for prepare_inputs() via sys.argv patching."""

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        # Create receptor and ligand text files
        self.rec_txt = str(self.tmpdir / "recs.txt")
        self.lig_txt = str(self.tmpdir / "ligs.txt")
        self.rec_pdb = str(self.tmpdir / "receptor.pdb")
        self.lig_sdf = str(self.tmpdir / "ligand.sdf")
        Path(self.rec_pdb).touch()
        Path(self.lig_sdf).touch()
        with open(self.rec_txt, "w") as f:
            f.write(self.rec_pdb + "\n")
        with open(self.lig_txt, "w") as f:
            f.write(self.lig_sdf + "\n")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_receptors_loaded_from_file(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt)):
            receptors, _, _, _, _, _, _ = np_docking.prepare_inputs()
        self.assertEqual(len(receptors), 1)
        self.assertIn("receptor.pdb", receptors[0])

    def test_ligands_loaded_from_file(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt)):
            _, ligands, _, _, _, _, _ = np_docking.prepare_inputs()
        self.assertEqual(len(ligands), 1)
        self.assertIn("ligand.sdf", ligands[0])

    def test_relative_paths_converted_to_absolute(self):
        # Write relative path into txt file
        rel_rec_txt = str(self.tmpdir / "rel_recs.txt")
        with open(rel_rec_txt, "w") as f:
            f.write("receptor.pdb\n")
        with patch("sys.argv", _argv("-r", rel_rec_txt, "-l", self.lig_txt)):
            receptors, _, _, _, _, _, _ = np_docking.prepare_inputs()
        self.assertTrue(Path(receptors[0]).is_absolute())

    def test_completed_filename_derived_from_stems(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt)):
            _, _, _, _, _, completed_name, _ = np_docking.prepare_inputs()
        self.assertIn("recs", completed_name)
        self.assertIn("ligs", completed_name)
        self.assertIn("completed", completed_name)

    def test_quiet_flag_defaults_to_false(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt)):
            _, _, _, _, quiet, _, _ = np_docking.prepare_inputs()
        self.assertFalse(quiet)

    def test_quiet_flag_set_when_passed(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt, "-q")):
            _, _, _, _, quiet, _, _ = np_docking.prepare_inputs()
        self.assertTrue(quiet)

    def test_ignore_existing_flag_set_when_passed(self):
        with patch("sys.argv", _argv("-r", self.rec_txt, "-l", self.lig_txt, "-i")):
            _, _, _, _, _, _, ignore = np_docking.prepare_inputs()
        self.assertTrue(ignore)


# ===========================================================================
# TestNaisarg
# ===========================================================================


class TestNaisarg(unittest.TestCase):
    """Integration tests for the main naisarg() orchestration function."""

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())
        self.rec = str(self.tmpdir / "receptor.pdb")
        self.lig = str(self.tmpdir / "ligand.sdf")
        self.completed = str(self.tmpdir / "recs_ligs_completed.txt")
        # Minimal config file so naisarg() has center/box_size defined
        self.config = str(self.tmpdir / "config.txt")
        Path(self.rec).touch()
        Path(self.lig).touch()
        with open(self.config, "w") as f:
            f.write("center_x = 0.0\ncenter_y = 0.0\ncenter_z = 0.0\n")
            f.write("size_x = 30\nsize_y = 30\nsize_z = 30\n")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _prepare_inputs_return(self, ignore_existing=False):
        """A prepare_inputs() return value using this test's files."""
        # read_config returns strings; naisarg unpacks them for dock_vina
        return (
            [self.rec],
            [self.lig],
            self.config,  # config file path
            None,  # autosite
            True,  # quiet
            self.completed,
            ignore_existing,
        )

    def _all_patches(self):
        """Context manager that patches every vina_helper function in np_docking."""
        return (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(),
            ),
            patch("docking.np_docking.prepare_receptor", return_value=(True, "")),
            patch("docking.np_docking.prepare_ligand", return_value=(True, "")),
            patch("docking.np_docking.dock_vina", return_value=(True, "")),
            patch(
                "docking.np_docking.vina_split", return_value=(-8.5, "out_ligand_1.sdf")
            ),
            patch("docking.np_docking.add_score_to_csv", return_value=(True, "name")),
            patch("docking.np_docking.backup", return_value=False),
        )

    def _run(self, **extra_patches):
        """Run naisarg() under all patches, optionally overriding some."""
        patches = list(self._all_patches())
        with (
            patches[0] as mock_pi,
            patches[1] as mock_rec,
            patches[2] as mock_lig,
            patches[3] as mock_dock,
            patches[4] as mock_split,
            patches[5] as mock_csv,
            patches[6] as mock_backup,
        ):
            np_docking.naisarg()
        return mock_pi, mock_rec, mock_lig, mock_dock, mock_split, mock_csv, mock_backup

    # ------------------------------------------------------------------ #
    # Workflow coverage                                                    #
    # ------------------------------------------------------------------ #

    def test_single_combination_docked(self):
        _, mock_rec, mock_lig, mock_dock, _, _, _ = self._run()
        mock_rec.assert_called_once()
        mock_lig.assert_called_once()
        mock_dock.assert_called_once()

    def test_vina_split_called_after_dock(self):
        _, _, _, _, mock_split, _, _ = self._run()
        mock_split.assert_called_once()

    def test_add_score_to_csv_called_after_split(self):
        _, _, _, _, _, mock_csv, _ = self._run()
        mock_csv.assert_called_once()

    def test_combination_written_to_completed_file(self):
        self._run()
        self.assertTrue(Path(self.completed).is_file())
        content = Path(self.completed).read_text()
        self.assertIn("receptor.pdb", content)
        self.assertIn("ligand.sdf", content)

    def test_completed_combination_not_docked_again(self):
        # Run once to populate the completed file
        self._run()
        # Second run: dock should NOT be called again
        with (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(),
            ),
            patch("docking.np_docking.prepare_receptor", return_value=(True, "")),
            patch("docking.np_docking.prepare_ligand", return_value=(True, "")),
            patch(
                "docking.np_docking.dock_vina", return_value=(True, "")
            ) as mock_dock2,
            patch("docking.np_docking.vina_split", return_value=(-8.5, "f.sdf")),
            patch("docking.np_docking.add_score_to_csv", return_value=(True, "n")),
            patch("docking.np_docking.backup", return_value=False),
        ):
            np_docking.naisarg()
        mock_dock2.assert_not_called()

    def test_ignore_existing_redocks_all(self):
        # Write a completed entry first
        Path(self.completed).write_text(f"('{self.rec}', '{self.lig}')\n")
        with (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(ignore_existing=True),
            ),
            patch("docking.np_docking.prepare_receptor", return_value=(True, "")),
            patch("docking.np_docking.prepare_ligand", return_value=(True, "")),
            patch("docking.np_docking.dock_vina", return_value=(True, "")) as mock_dock,
            patch("docking.np_docking.vina_split", return_value=(-8.5, "f.sdf")),
            patch("docking.np_docking.add_score_to_csv", return_value=(True, "n")),
            patch("docking.np_docking.backup", return_value=False),
        ):
            np_docking.naisarg()
        mock_dock.assert_called_once()

    # ------------------------------------------------------------------ #
    # Error / skip behaviour                                               #
    # ------------------------------------------------------------------ #

    def test_receptor_prep_failure_skips_combination(self):
        with (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(),
            ),
            patch(
                "docking.np_docking.prepare_receptor",
                return_value=(False, "bad receptor"),
            ),
            patch("docking.np_docking.prepare_ligand", return_value=(True, "")),
            patch("docking.np_docking.dock_vina", return_value=(True, "")) as mock_dock,
            patch("docking.np_docking.vina_split", return_value=(-8.5, "f.sdf")),
            patch("docking.np_docking.add_score_to_csv", return_value=(True, "n")),
            patch("docking.np_docking.backup", return_value=False),
        ):
            np_docking.naisarg()
        mock_dock.assert_not_called()

    def test_ligand_prep_failure_skips_combination(self):
        with (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(),
            ),
            patch("docking.np_docking.prepare_receptor", return_value=(True, "")),
            patch(
                "docking.np_docking.prepare_ligand", return_value=(False, "bad ligand")
            ),
            patch("docking.np_docking.dock_vina", return_value=(True, "")) as mock_dock,
            patch("docking.np_docking.vina_split", return_value=(-8.5, "f.sdf")),
            patch("docking.np_docking.add_score_to_csv", return_value=(True, "n")),
            patch("docking.np_docking.backup", return_value=False),
        ):
            np_docking.naisarg()
        mock_dock.assert_not_called()

    def test_docking_failure_skips_split_and_csv(self):
        with (
            patch(
                "docking.np_docking.prepare_inputs",
                return_value=self._prepare_inputs_return(),
            ),
            patch("docking.np_docking.prepare_receptor", return_value=(True, "")),
            patch("docking.np_docking.prepare_ligand", return_value=(True, "")),
            patch("docking.np_docking.dock_vina", return_value=(False, "vina error")),
            patch(
                "docking.np_docking.vina_split", return_value=(-8.5, "f.sdf")
            ) as mock_split,
            patch(
                "docking.np_docking.add_score_to_csv", return_value=(True, "n")
            ) as mock_csv,
            patch("docking.np_docking.backup", return_value=False),
        ):
            np_docking.naisarg()
        mock_split.assert_not_called()
        mock_csv.assert_not_called()


if __name__ == "__main__":
    unittest.main()
