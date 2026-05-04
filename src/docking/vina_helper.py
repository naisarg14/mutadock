################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2026 Naisarg Patel (https://github.com/naisarg14)        #
#                                                                              #
#          Project: https://github.com/naisarg14/mutadock                      #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
################################################################################


import logging
from pathlib import Path
from typing import Any, Optional

from .exceptions import (
    ConfigError,
    DockingError,
    DockingRunError,
    LigandPreparationError,
    PDBFileError,
    ReceptorPreparationError,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def backup(file_path: str) -> bool:
    """Move *file_path* into a timestamped backup in a ``backups/`` sub-folder.

    The backup folder is created next to the file if it does not exist.  The
    file is renamed to ``<stem>_<modified-timestamp><suffix>`` so multiple
    backups of the same file never collide.

    Args:
        file_path: Path to the file to back up.

    Returns:
        ``True`` if the file was moved, ``False`` if it did not exist.
    """
    from datetime import datetime

    path = Path(file_path)
    if not path.exists():
        return False
    target_directory = path.parent / "backups"
    target_directory.mkdir(exist_ok=True)
    modified_time = path.stat().st_mtime
    timestamp = datetime.fromtimestamp(modified_time).strftime("%b-%d-%Y_%H.%M")
    target_file = target_directory / f"{path.stem}_{timestamp}{path.suffix}"
    path.rename(target_file)
    return True


def read_pdb_file(file_path: str) -> list[dict[str, Any]]:
    """Parse ATOM/HETATM records from a PDB file into a list of dicts.

    Each dict contains: ``atom_serial_number``, ``atom_name``,
    ``residue_name``, ``chain_id``, ``residue_sequence_number``, ``x``,
    ``y``, ``z``, ``occupancy``, ``temp_factor``, and ``extra_factor``.

    Args:
        file_path: Path to the PDB file.

    Returns:
        List of atom dicts.

    Raises:
        PDBFileError: If the file cannot be read or parsed.
    """
    import re

    try:
        atoms = []
        pattern = r"^ATOM\s+(\d+)\s+([A-Z]+)\s+([A-Z]{2,3})\s+([A-Z]?)\s*(\d+)\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(\d+\.\d{1,2})\s+(\d+\.\d{1,3})(?:\s+(\d+\.\d{1,3}))?\s+[A-Z]$"
        with open(file_path) as file:
            for line in file:
                if line.startswith(("ATOM", "HETATM")):
                    match = re.match(pattern, line)
                if match:
                    parts = {
                        "atom_serial_number": int(match.group(1)),
                        "atom_name": match.group(2).strip(),
                        "residue_name": match.group(3).strip(),
                        "chain_id": match.group(4).strip() or None,
                        "residue_sequence_number": int(match.group(5)),
                        "x": float(match.group(6)),
                        "y": float(match.group(7)),
                        "z": float(match.group(8)),
                        "occupancy": float(match.group(9)),
                        "temp_factor": float(match.group(10)),
                        "extra_factor": (
                            float(match.group(11)) if match.group(11) else None
                        ),
                    }
                    atoms.append(parts)
        return atoms
    except Exception as e:
        raise PDBFileError(f"Failed to read PDB file '{file_path}': {e}") from e


def calculate_geometric_center(pdb_file: str) -> tuple[float, float, float]:
    """Return the geometric center (x, y, z) of all atoms in a PDB file.

    Args:
        pdb_file: Path to the PDB file.

    Returns:
        A 3-tuple ``(cx, cy, cz)`` of mean coordinates.

    Raises:
        PDBFileError: If the PDB file cannot be parsed.
    """
    atoms = read_pdb_file(pdb_file)
    num_atoms = len(atoms)
    x_sum = sum(atom["x"] for atom in atoms)
    y_sum = sum(atom["y"] for atom in atoms)
    z_sum = sum(atom["z"] for atom in atoms)

    return (x_sum / num_atoms, y_sum / num_atoms, z_sum / num_atoms)


def calculate_radius(pdb_file: str) -> float:
    """Return the maximum distance from the geometric center to any atom.

    Args:
        pdb_file: Path to the PDB file.

    Returns:
        Max distance in Ångströms.

    Raises:
        PDBFileError: If the PDB file cannot be parsed.
    """
    import math

    atoms = read_pdb_file(pdb_file)

    num_atoms = len(atoms)
    cx = sum(a["x"] for a in atoms) / num_atoms
    cy = sum(a["y"] for a in atoms) / num_atoms
    cz = sum(a["z"] for a in atoms) / num_atoms

    max_distance: float = 0.0
    for atom in atoms:
        distance = math.sqrt(
            (atom["x"] - cx) ** 2 + (atom["y"] - cy) ** 2 + (atom["z"] - cz) ** 2
        )
        if distance > max_distance:
            max_distance = distance

    return max_distance


def vina_split(input_file: str, output_file: Optional[str] = None) -> tuple[float, str]:
    """Extract the best docking pose from a Vina PDBQT output and write it as SDF.

    Reads the first MODEL block from *input_file*, converts it to an SDF
    molecule via meeko, appends the docking score as an SD tag, and writes
    the result to *output_file*.

    Args:
        input_file: Path to the Vina PDBQT output file.
        output_file: Destination SDF path.  Defaults to replacing ``.pdbqt``
            with ``_ligand_1.sdf``.

    Returns:
        A 2-tuple ``(score, output_file)`` where *score* is the binding
        affinity in kcal/mol.
    """
    try:
        import sys

        from meeko import PDBQTMolecule, RDKitMolCreate
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing ligand files for Docking.\n"
        msg += "Easaies way to fix this is to install meeko using the following command:\n\n"
        msg += "python -m pip install meeko\n"
        msg += "If you already have meeko installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        logger.error(msg)
        sys.exit(2)

    if output_file is None:
        output_file = input_file.replace(".pdbqt", "_ligand_1.sdf")

    pdbqt_string = ""
    with open(input_file) as infile:
        for line in infile:
            if "vina result" in line.lower():
                score = [
                    float(x)
                    for x in line.split()
                    if x.replace(".", "", 1).replace("-", "", 1).isdigit()
                ][0]
            pdbqt_string += line
            if line.startswith("ENDMDL"):
                break
    molecule = PDBQTMolecule(pdbqt_string)
    sdf_string, failures = RDKitMolCreate.write_sd_string(molecule)

    if len(failures) > 0:
        msg = "\nCould not convert to RDKit. Maybe this library was not used for preparing\n"
        msg += "the input PDBQT for docking, and the SMILES string is missing?\n"
        msg += "Except for standard protein sidechains, all ligands and flexible residues\n"
        msg += "require a REMARK SMILES line in the PDBQT, which is added automatically by meeko."
        raise RuntimeError(msg)

    footer_string = f"> <Docking Score>\n{score}\n> <Credits>\nCreated using a script in mutadock library written by Naisarg Patel (https://github.com/naisarg14/mutadock).\n$$$$\n"
    with open(output_file, "w") as outfile:
        outfile.write(sdf_string.replace("$$$$", footer_string))

    return (score, output_file)


def prepare_ligand(in_file: str, out_file: Optional[str] = None) -> str:
    """Convert an SDF or MOL2 ligand file to PDBQT format using meeko.

    Args:
        in_file: Path to the input ligand file (``.sdf`` or ``.mol2``).
        out_file: Destination PDBQT path.  Defaults to replacing the input
            extension with ``.pdbqt``.

    Returns:
        The PDBQT string written to *out_file*.

    Raises:
        LigandPreparationError: If preparation fails for any reason.
    """
    try:
        import sys

        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing ligand files for Docking.\n"
        msg += "Easaies way to fix this is to install meeko and rdkit using the following command:\n\n"
        msg += "python -m pip install meeko rdkit\n"
        msg += "If you already have meeko and rdkit installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        logger.error(msg)
        sys.exit(2)

    if out_file is None:
        if in_file.endswith(".sdf"):
            out_file = f"{in_file.removesuffix('.sdf')}.pdbqt"
        elif in_file.endswith(".mol2"):
            out_file = f"{in_file.removesuffix('.mol2')}.pdbqt"
        else:
            raise LigandPreparationError("Input file is not in SDF or MOL2 format.")

    try:
        if in_file.endswith(".sdf"):
            mol = Chem.SDMolSupplier(in_file)[0]
        if in_file.endswith(".mol2"):
            mol = Chem.MolFromMol2File(in_file)

        mol = Chem.AddHs(mol)
        mp = MoleculePreparation()
        molecule_setups = mp.prepare(mol)
        pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(
            molecule_setups[0]
        )

        if not success:
            raise RuntimeError(f"Could not convert to PDBQT: {error_msg}")
        with open(out_file, "w") as output_file:
            output_file.write(pdbqt_string)
    except LigandPreparationError:
        raise
    except Exception as e:
        raise LigandPreparationError(
            f"Failed to prepare ligand '{in_file}': {e}"
        ) from e

    return pdbqt_string


def prepare_receptor(
    receptor_filename: str,
    outputfilename: Optional[str] = None,
    ph: float = 7.4,
) -> None:
    """Convert a PDB or CIF receptor file to PDBQT format.

    Fixes the structure with PDBFixer (missing residues/atoms, nonstandard
    residues, hydrogens at *ph*), writes a temporary PDB, then converts it
    to PDBQT via meeko's ``mk_receptor`` CLI.  The temporary file is always
    removed when the function exits.

    Args:
        receptor_filename: Path to the input receptor file (``.pdb`` or ``.cif``).
        outputfilename: Destination PDBQT path.  Defaults to replacing the
            input extension with ``.pdbqt``.
        ph: pH used when adding missing hydrogens (default 7.4).

    Raises:
        ReceptorPreparationError: If preparation or conversion fails.
    """
    import subprocess
    import sys
    import tempfile

    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing receptor files for Docking.\n"
        msg += "Install pdbfixer and openmm:\n\n"
        msg += "  conda install -c conda-forge pdbfixer openmm\n"
        msg += "  or: python -m pip install pdbfixer openmm\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        logger.error(msg)
        sys.exit(2)

    if outputfilename is None:
        outputfilename = str(Path(receptor_filename).with_suffix(".pdbqt"))

    tmp_fd, tmp_path = tempfile.mkstemp(suffix="_fixed.pdb")
    try:
        fixer = PDBFixer(filename=receptor_filename)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)

        with open(tmp_fd, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        result = subprocess.run(
            ["mk_receptor", "-i", tmp_path, "-o", outputfilename],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"mk_receptor failed:\n{result.stderr.strip()}")
    except ReceptorPreparationError:
        raise
    except Exception as e:
        raise ReceptorPreparationError(
            f"Failed to prepare receptor '{receptor_filename}': {e}"
        ) from e
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def add_score_to_csv(out_pdb: str, csv_file: str, score: float) -> str:
    """Append a docking result row to the aggregate CSV file.

    Reads the last row to determine the next serial number, then appends a
    row with columns ``sr``, ``name``, and ``affinity``.

    Args:
        out_pdb: Path to the docking output file (used to derive the run name).
        csv_file: Path to the CSV file to append to (created if absent).
        score: Binding affinity in kcal/mol.

    Returns:
        The run name derived from *out_pdb*.

    Raises:
        DockingError: If writing to the CSV file fails.
    """
    import csv

    try:
        with open(csv_file) as lc:
            final_line = lc.readlines()[-1]
            count = int(final_line.split(",")[0]) + 1
    except (FileNotFoundError, ValueError):
        count = 1

    name = Path(out_pdb).name.removesuffix("_out.pdb")
    try:
        with open(csv_file, "a+", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=["sr", "name", "affinity"])
            writer.writerow({"sr": count, "name": name, "affinity": score})
    except Exception as e:
        raise DockingError(f"Failed to write to CSV '{csv_file}': {e}") from e

    return name


def read_config(
    config_file: str,
) -> tuple[list[float], list[float], int, int, int, bool]:
    """Parse a Vina configuration file and return docking parameters.

    Expects ``key = value`` lines; lines starting with ``#`` are ignored.
    Falls back to safe defaults for any missing keys.

    Args:
        config_file: Path to the Vina ``.conf`` / ``.txt`` config file.

    Returns:
        ``(center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite)``

    Raises:
        ConfigError: If the file cannot be read or parsed.
    """
    try:
        config = {}
        with open(config_file) as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, value = line.split("=")
                    config[key.strip()] = value.strip()

        center = [
            float(config.get("center_x", "0.0")),
            float(config.get("center_y", "0.0")),
            float(config.get("center_z", "0.0")),
        ]
        box_size = [
            float(config.get("size_x", "30.0")),
            float(config.get("size_y", "30.0")),
            float(config.get("size_z", "30.0")),
        ]
        exhaustiveness = int(config.get("exhaustiveness", "32"))
        n_poses = int(config.get("n_poses", "20"))
        n_poses_write = int(config.get("n_poses_write", "5"))
        overwrite = config.get("overwrite", "True").lower() in ("true", "1", "yes")

        return (center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite)

    except Exception as e:
        raise ConfigError(f"Failed to read config file '{config_file}': {e}") from e


def dock_vina(
    receptor: str,
    ligand: str,
    output: str,
    log_file: str,
    config: Optional[str] = None,
    autosite: Optional[str] = None,
    center: Optional[list[float]] = None,
    box_size: Optional[list[float]] = None,
    exhaustiveness: int = 32,
    n_poses: int = 20,
    n_poses_write: int = 5,
    overwrite: bool = True,
) -> None:
    """Run AutoDock Vina by spawning ``vina_dock.py`` as a subprocess.

    Docking parameters can be supplied directly or read from a config file.
    If *autosite* is provided its geometric center overrides *center*.

    Args:
        receptor: Path to the prepared receptor PDBQT file.
        ligand: Path to the prepared ligand PDBQT file.
        output: Path for the docking output PDBQT file.
        log_file: Path where stdout/stderr from Vina will be written.
        config: Optional Vina config file; overrides *center*, *box_size*,
            and other parameters when provided.
        autosite: Optional PDB file from AutoSite; its geometric center is
            used as the docking box center when provided.
        center: Search box center ``[x, y, z]`` in Ångströms.
        box_size: Search box dimensions ``[x, y, z]`` in Ångströms.
        exhaustiveness: Exhaustiveness of the global search (default 32).
        n_poses: Number of poses to generate (default 20).
        n_poses_write: Number of poses to write to *output* (default 5).
        overwrite: Whether to overwrite an existing output file.

    Raises:
        ConfigError: If the config file cannot be read.
        PDBFileError: If the autosite PDB cannot be parsed.
        DockingRunError: If the Vina subprocess exits with a non-zero code.
    """
    if config is not None:
        center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite = (
            read_config(config)
        )

    if autosite is not None:
        center = list(calculate_geometric_center(autosite))

    if center is None:
        center = [0.0, 0.0, 0.0]
    if box_size is None:
        box_size = [30.0, 30.0, 30.0]

    import subprocess

    vina_dock_script = Path(__file__).parent / "vina_dock.py"

    commands = [
        "python3",
        str(vina_dock_script),
        "--receptor",
        receptor,
        "--ligand",
        ligand,
        "--output",
        output,
        "--center",
        str(center[0]),
        str(center[1]),
        str(center[2]),
        "--box_size",
        str(box_size[0]),
        str(box_size[1]),
        str(box_size[2]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--n_poses",
        str(n_poses),
        "--n_poses_write",
        str(n_poses_write),
    ]
    if not overwrite:
        commands.append("--nooverwrite")
    with open(log_file, "w+") as lfile:
        result = subprocess.run(commands, stdout=lfile, stderr=lfile, text=True)

    if result.returncode != 0:
        raise DockingRunError(f"Check the error in {log_file}")


if __name__ == "__main__":
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's docking module."
    )
