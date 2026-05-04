"""
Tests for docking.vina_dock
-----------------------------
The AutoDock Vina Python bindings ('vina' package) are stubbed via
patch.dict(sys.modules) so the real package need not be installed.

Run from the project root:
    pytest tests/
"""

import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from docking import vina_dock

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_vina_stub():
    """Return (stub_module, Vina_class_mock, vina_instance_mock)."""
    mock_instance = MagicMock(name="vina_instance")
    mock_cls = MagicMock(name="Vina", return_value=mock_instance)
    stub = types.ModuleType("vina")
    stub.Vina = mock_cls
    return stub, mock_cls, mock_instance


# ===========================================================================
# TestVinaDock
# ===========================================================================


class TestVinaDock(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.receptor = str(Path(self.tmpdir) / "rec.pdbqt")
        self.ligand = str(Path(self.tmpdir) / "lig.pdbqt")
        self.output = str(Path(self.tmpdir) / "out.pdbqt")
        Path(self.receptor).touch()
        Path(self.ligand).touch()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------ #
    # Input validation                                                     #
    # ------------------------------------------------------------------ #

    def test_missing_receptor_returns_false(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock("ghost_rec.pdbqt", self.ligand, self.output)
        self.assertIs(result[0], False)
        self.assertIn("Receptor file not found", result[1])

    def test_missing_ligand_returns_false(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock(self.receptor, "ghost_lig.pdbqt", self.output)
        self.assertIs(result[0], False)
        self.assertIn("Ligand file not found", result[1])

    # ------------------------------------------------------------------ #
    # Success path                                                         #
    # ------------------------------------------------------------------ #

    def test_success_returns_true_and_message(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertTrue(result[0])
        self.assertEqual(result[1], "Docking successful")

    def test_set_receptor_called_with_path(self):
        stub, _, mock_v = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        mock_v.set_receptor.assert_called_once_with(self.receptor)

    def test_set_ligand_called_with_path(self):
        stub, _, mock_v = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        mock_v.set_ligand_from_file.assert_called_once_with(self.ligand)

    def test_compute_vina_maps_called_with_center_and_box(self):
        stub, _, mock_v = _make_vina_stub()
        center = [1.0, 2.0, 3.0]
        box = [20, 20, 20]
        with patch.dict(sys.modules, {"vina": stub}):
            vina_dock.vina_dock(
                self.receptor, self.ligand, self.output, center=center, box_size=box
            )
        mock_v.compute_vina_maps.assert_called_once_with(center=center, box_size=box)

    def test_dock_called_with_exhaustiveness_and_n_poses(self):
        stub, _, mock_v = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            vina_dock.vina_dock(
                self.receptor, self.ligand, self.output, exhaustiveness=16, n_poses=10
            )
        mock_v.dock.assert_called_once_with(exhaustiveness=16, n_poses=10)

    def test_write_poses_called_with_n_poses_write_and_overwrite(self):
        stub, _, mock_v = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            vina_dock.vina_dock(
                self.receptor,
                self.ligand,
                self.output,
                n_poses_write=3,
                overwrite=False,
            )
        mock_v.write_poses.assert_called_once_with(
            self.output, n_poses=3, overwrite=False
        )

    # ------------------------------------------------------------------ #
    # Failure paths                                                        #
    # ------------------------------------------------------------------ #

    def test_vina_not_installed_returns_false(self):
        """Simulate missing vina package by putting None sentinel in sys.modules."""
        with patch.dict(sys.modules, {"vina": None}):
            result = vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertIs(result[0], False)

    def test_exception_during_docking_returns_false_with_message(self):
        stub, _, mock_v = _make_vina_stub()
        mock_v.dock.side_effect = RuntimeError("search failed")
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertIs(result[0], False)
        self.assertIn("Docking failed", result[1])

    def test_exception_during_write_poses_returns_false(self):
        stub, _, mock_v = _make_vina_stub()
        mock_v.write_poses.side_effect = OSError("disk full")
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertIs(result[0], False)


if __name__ == "__main__":
    unittest.main()
