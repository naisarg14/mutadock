"""
Tests for mutadock.docking.vina_dock
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

from mutadock.docking import vina_dock, vina_helper
from mutadock.docking.exceptions import DockingRunError

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

    def test_missing_receptor_raises(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            with self.assertRaises(DockingRunError) as ctx:
                vina_dock.vina_dock("ghost_rec.pdbqt", self.ligand, self.output)
        self.assertIn("Receptor file not found", str(ctx.exception))

    def test_missing_ligand_raises(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            with self.assertRaises(DockingRunError) as ctx:
                vina_dock.vina_dock(self.receptor, "ghost_lig.pdbqt", self.output)
        self.assertIn("Ligand file not found", str(ctx.exception))

    # ------------------------------------------------------------------ #
    # Success path                                                         #
    # ------------------------------------------------------------------ #

    def test_success_does_not_raise(self):
        stub, _, _ = _make_vina_stub()
        with patch.dict(sys.modules, {"vina": stub}):
            result = vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertIsNone(result)

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

    def test_vina_not_installed_exits(self):
        """Simulate missing vina package by putting None sentinel in sys.modules."""
        with patch.dict(sys.modules, {"vina": None}):
            with self.assertRaises(SystemExit):
                vina_dock.vina_dock(self.receptor, self.ligand, self.output)

    def test_exception_during_docking_raises_docking_run_error(self):
        stub, _, mock_v = _make_vina_stub()
        mock_v.dock.side_effect = RuntimeError("search failed")
        with patch.dict(sys.modules, {"vina": stub}):
            with self.assertRaises(DockingRunError) as ctx:
                vina_dock.vina_dock(self.receptor, self.ligand, self.output)
        self.assertIn("Docking failed", str(ctx.exception))

    def test_exception_during_write_poses_raises_docking_run_error(self):
        stub, _, mock_v = _make_vina_stub()
        mock_v.write_poses.side_effect = OSError("disk full")
        with patch.dict(sys.modules, {"vina": stub}):
            with self.assertRaises(DockingRunError):
                vina_dock.vina_dock(self.receptor, self.ligand, self.output)


# ===========================================================================
# TestOverwriteRoundTrip
# ===========================================================================
#
# Bug 2.6: dock_vina (vina_helper) emits an overwrite flag that vina_dock.py's
# argparse must accept, and the parsed value must reach write_poses with the
# intended boolean. This exercises BOTH sides so a one-sided fix is caught.


class TestOverwriteRoundTrip(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.receptor = str(Path(self.tmpdir) / "rec.pdbqt")
        self.ligand = str(Path(self.tmpdir) / "lig.pdbqt")
        self.output = str(Path(self.tmpdir) / "out.pdbqt")
        self.log = str(Path(self.tmpdir) / "log.txt")
        Path(self.receptor).touch()
        Path(self.ligand).touch()

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _emitted_commands(self, overwrite):
        """Return the argv list dock_vina would run vina_dock.py with."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        with patch("subprocess.run", return_value=mock_result) as mock_run:
            vina_helper.dock_vina(
                self.receptor,
                self.ligand,
                self.output,
                self.log,
                overwrite=overwrite,
            )
        return mock_run.call_args[0][0]

    def _parse_overwrite(self, cmd):
        """Feed the emitted args (minus python + script path) into vina_dock's
        parser via main() and capture the overwrite value handed to vina_dock."""
        argv = ["vina_dock.py"] + list(cmd[2:])
        with patch.object(vina_dock, "vina_dock") as mock_dock:
            with patch.object(sys, "argv", argv):
                with self.assertRaises(SystemExit) as ctx:
                    vina_dock.main()
        self.assertEqual(ctx.exception.code, 0)
        return mock_dock.call_args.kwargs["overwrite"]

    def test_overwrite_true_roundtrips(self):
        cmd = self._emitted_commands(True)
        self.assertIn("--overwrite", cmd)
        self.assertTrue(self._parse_overwrite(cmd))

    def test_overwrite_false_roundtrips(self):
        cmd = self._emitted_commands(False)
        self.assertIn("--no-overwrite", cmd)
        self.assertFalse(self._parse_overwrite(cmd))


if __name__ == "__main__":
    unittest.main()
