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
import os
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

try:
    from .exceptions import (
        ConfigError,
        DockingError,
        DockingRunError,
        LigandPreparationError,
        PDBFileError,
        ReceptorPreparationError,
    )
except ImportError:
    from exceptions import (  # type: ignore[no-redef]
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

# PubChem PUG REST endpoint for a 3D SDF conformer, by CID or name.
PUBCHEM_SDF_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{route}/SDF?record_type=3d"
)


def _timeout_from_env(var: str, default: float) -> Optional[float]:
    """Resolve a subprocess timeout (in seconds) from *var*, else *default*.

    A non-positive or unparseable value disables the timeout (returns ``None``),
    giving users an escape hatch for legitimately long-running jobs.
    """
    raw = os.environ.get(var)
    if raw is None:
        return default
    try:
        value = float(raw)
    except ValueError:
        logger.warning(
            "Ignoring invalid %s=%r; falling back to %ss.", var, raw, default
        )
        return default
    return value if value > 0 else None


# Per-subprocess wall-clock ceilings (seconds) so a hung external tool
# (mk_prepare_receptor, AutoSite, Vina) can't stall an entire batch.  Each is
# overridable via the matching ``MUTADOCK_*`` environment variable; set it to 0
# (or a negative value) to disable that timeout entirely.
RECEPTOR_PREP_TIMEOUT = _timeout_from_env("MUTADOCK_RECEPTOR_PREP_TIMEOUT", 900.0)
AUTOSITE_TIMEOUT = _timeout_from_env("MUTADOCK_AUTOSITE_TIMEOUT", 1800.0)
VINA_TIMEOUT = _timeout_from_env("MUTADOCK_VINA_TIMEOUT", 3600.0)

# Default Vina RNG seed for every docking mutadock performs.
#
# Vina's own default is 0, which means "choose a random seed", so before this
# constant existed two identical mutadock runs could return different
# affinities with nothing in the output to explain it. That is not acceptable
# for a tool whose results get published: the seed materially changes the
# answer, not just its last decimal. Measured on BCR-ABL/imatinib at
# exhaustiveness 32, changing only the seed moves a single Delta-binding by up
# to 2.0 kcal/mol, because the seed can change WHICH binding mode is found
# rather than merely how precisely one mode is scored.
#
# So a fixed seed is the default, it is recorded on every output row, and it is
# overridable (``--seed`` on md_dock / md_vina_dock / md_quick, or ``seed =``
# in a Vina config file) for anyone who wants replicates. A fixed default makes
# a run reproducible; it does NOT make one run sufficient -- report replicates
# across several seeds before quoting a Delta-binding as a result.
DEFAULT_VINA_SEED = 19


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
        mutadock.docking.exceptions.PDBFileError: If the file cannot be read or parsed.
    """
    import re

    try:
        atoms = []
        pattern = r"^(?:ATOM|HETATM)\s+(\d+)\s+([A-Z]+)\s+([A-Z]{2,3})\s+([A-Z]?)\s*(\d+)\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(\d+\.\d{1,2})\s+(\d+\.\d{1,3})(?:\s+(\d+\.\d{1,3}))?\s+[A-Z]$"
        with open(file_path) as file:
            for line in file:
                match = None
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
        mutadock.docking.exceptions.PDBFileError: If the PDB file cannot be parsed.
    """
    atoms = read_pdb_file(pdb_file)
    num_atoms = len(atoms)
    if num_atoms == 0:
        raise PDBFileError(
            f"Cannot compute geometric center: no atoms found in '{pdb_file}'."
        )
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
        mutadock.docking.exceptions.PDBFileError: If the PDB file cannot be parsed.
    """
    import math

    atoms = read_pdb_file(pdb_file)

    num_atoms = len(atoms)
    if num_atoms == 0:
        raise PDBFileError(f"Cannot compute radius: no atoms found in '{pdb_file}'.")
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
        import re

        from meeko import PDBQTMolecule, RDKitMolCreate
    except ModuleNotFoundError:
        msg = "Error with importing modules for preparing ligand files for Docking.\n"
        msg += "Easiest way to fix this is to install meeko using the following command:\n\n"
        msg += "python -m pip install meeko\n"
        msg += "If you already have meeko installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        logger.error(msg)
        sys.exit(2)

    if output_file is None:
        output_file = input_file.replace(".pdbqt", "_ligand_1.sdf")

    pdbqt_string = ""
    score = None
    with open(input_file) as infile:
        for line in infile:
            if "vina result" in line.lower():
                number_match = re.search(
                    r"vina result:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
                    line,
                    re.IGNORECASE,
                )
                if number_match:
                    score = float(number_match.group(1))
            pdbqt_string += line
            if line.startswith("ENDMDL"):
                break

    if score is None:
        raise DockingError(f"No 'VINA RESULT' affinity found in '{input_file}'.")

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


def _validate_sdf(path: Path, label: str) -> None:
    """Raise ``LigandPreparationError`` unless *path* is a usable 3D SDF.

    PubChem returns an HTTP error body (not an SDF) for unknown compounds; a
    plain ``urlretrieve`` writes that body to disk without raising, so the
    download must be validated rather than trusted.
    """
    text = path.read_text(errors="ignore")
    if "$$$$" not in text or ("V2000" not in text and "V3000" not in text):
        path.unlink(missing_ok=True)
        raise LigandPreparationError(
            f"Downloaded file for '{label}' is not a valid SDF "
            "(PubChem likely returned an error page)."
        )
    try:
        from rdkit import Chem
    except ImportError:
        return  # RDKit unavailable — the textual check above is the best we can do
    try:
        mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)), None)
        if mol is None or mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
            raise ValueError("no parseable molecule with 3D coordinates")
    except LigandPreparationError:
        raise
    except Exception as e:
        path.unlink(missing_ok=True)
        raise LigandPreparationError(
            f"Downloaded SDF for '{label}' is unusable: {e}"
        ) from e


def fetch_ligand(
    code: str,
    dest_dir: Optional[str | Path] = None,
    name: Optional[str] = None,
) -> str:
    """Download a 3D ligand SDF from PubChem by CID or compound name.

    Args:
        code: A PubChem identifier. Digits are treated as a CID; anything else
            is treated as a compound name. Explicit ``cid:`` / ``name:``
            prefixes override the auto-detection (e.g. ``"cid:2244"``,
            ``"name:aspirin"``).
        dest_dir: Directory to save the SDF into. Defaults to the current
            working directory. An existing non-empty file is reused.
        name: Base filename (without extension). Defaults to ``CID_<n>`` for a
            CID or a filesystem-safe form of the name.

    Returns:
        Path to the downloaded (or already-present) 3D ``.sdf`` file.

    Raises:
        LigandPreparationError: If *code* is malformed, the download fails, or
            the downloaded file is not a usable 3D SDF.
    """
    import re
    import urllib.error
    import urllib.parse
    import urllib.request

    raw = code.strip()
    if not raw:
        raise LigandPreparationError("Empty ligand code.")
    lowered = raw.lower()
    if lowered.startswith("cid:"):
        kind, value = "cid", raw[4:].strip()
    elif lowered.startswith("name:"):
        kind, value = "name", raw[5:].strip()
    elif raw.isdigit():
        kind, value = "cid", raw
    else:
        kind, value = "name", raw
    if not value:
        raise LigandPreparationError(f"Invalid ligand code '{code}'.")

    if kind == "cid":
        if not value.isdigit():
            raise LigandPreparationError(
                f"Invalid PubChem CID '{value}': expected digits (e.g. 2244)."
            )
        route = f"cid/{value}"
        stem = name or f"CID_{value}"
    else:
        route = f"name/{urllib.parse.quote(value)}"
        stem = name or (re.sub(r"[^0-9A-Za-z._-]+", "_", value).strip("_") or "ligand")

    directory = Path(dest_dir) if dest_dir else Path.cwd()
    directory.mkdir(parents=True, exist_ok=True)
    dest = directory / f"{stem}.sdf"
    if dest.is_file() and dest.stat().st_size > 0:
        logger.info(
            "Ligand '%s' already present at '%s'; skipping download.", value, dest
        )
        return str(dest)

    url = PUBCHEM_SDF_URL.format(route=route)
    logger.info("Fetching ligand %s '%s' from PubChem ...", kind, value)
    try:
        urllib.request.urlretrieve(url, dest)
    except urllib.error.HTTPError as e:
        raise LigandPreparationError(
            f"PubChem has no 3D record for {kind} '{value}' (HTTP {e.code}). "
            "Check the CID/name, or supply your own SDF/MOL2 ligand file."
        ) from e
    except Exception as e:
        raise LigandPreparationError(
            f"Failed to fetch ligand '{value}' from PubChem ({url}): {e}"
        ) from e

    _validate_sdf(dest, value)
    logger.info("Saved ligand to '%s'.", dest)
    return str(dest)


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
        msg += "Easiest way to fix this is to install meeko and rdkit using the following command:\n\n"
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

        # addCoords=True is REQUIRED. Without it RDKit adds the hydrogens as
        # topology only, with no 3D coordinates, and meeko then writes every one
        # of them at the origin (0, 0, 0). For a ligand docked in a box centred
        # anywhere else that puts the polar/donor hydrogens tens of Angstroms
        # outside the search volume, and Vina infers donor typing from whether a
        # polar H sits within bonding distance of its N/O -- so those atoms stop
        # being treated as H-bond donors and the score changes. Measured on
        # imatinib/2HYY: WT -12.22 -> -12.55 kcal/mol, 5x the seed-noise floor.
        mol = Chem.AddHs(mol, addCoords=True)
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
    input_pdb: str,
    output_pdbqt: Optional[str] = None,
) -> None:
    """Prepare a receptor PDB file for docking using PDBFixer and meeko.

    Fixes missing residues/atoms and adds hydrogens with PDBFixer, then converts
    the result to a rigid-receptor PDBQT with meeko's ``mk_prepare_receptor.py``
    command-line tool (which, unlike meeko's ligand ``MoleculePreparation`` API,
    correctly handles multi-chain / multi-fragment receptors).  A temporary
    ``_fixed.pdb`` file is created and removed regardless of success or failure.

    Args:
        input_pdb: Path to the input PDB file.
        output_pdbqt: Destination PDBQT path.  Defaults to *input_pdb* with
            the extension replaced by ``.pdbqt``.

    Raises:
        ReceptorPreparationError: If any preparation step fails or the meeko
            ``mk_prepare_receptor.py`` tool is not on PATH.
    """
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ModuleNotFoundError as e:
        raise ReceptorPreparationError(f"Missing required module: {e}") from e

    mk_receptor = shutil.which("mk_prepare_receptor.py") or shutil.which(
        "mk_prepare_receptor"
    )
    if mk_receptor is None:
        raise ReceptorPreparationError(
            "meeko's 'mk_prepare_receptor.py' was not found on PATH. It ships "
            "with meeko (pip install meeko); ensure the environment's scripts "
            "directory is on PATH."
        )

    if output_pdbqt is None:
        output_pdbqt = str(Path(input_pdb).with_suffix(".pdbqt"))

    input_path = Path(input_pdb)
    tmp_pdb = str(input_path.parent / f"{input_path.stem}_fixed.pdb")

    try:
        fixer = PDBFixer(filename=input_pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        with open(tmp_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        result = subprocess.run(
            [sys.executable, mk_receptor, "--read_pdb", tmp_pdb, "-p", output_pdbqt],
            capture_output=True,
            text=True,
            timeout=RECEPTOR_PREP_TIMEOUT,
        )
        if result.returncode != 0 or not Path(output_pdbqt).is_file():
            raise ValueError(
                "mk_prepare_receptor failed: "
                f"{result.stderr.strip() or result.stdout.strip()}"
            )

    except ReceptorPreparationError:
        raise
    except Exception as e:
        raise ReceptorPreparationError(
            f"Failed to prepare receptor '{input_pdb}': {e}"
        ) from e
    finally:
        if Path(tmp_pdb).exists():
            Path(tmp_pdb).unlink()


def add_score_to_csv(
    pose_file: str,
    csv_file: str,
    score: float,
    center: Optional[Sequence[float]] = None,
    box_size: Optional[Sequence[float]] = None,
    exhaustiveness: Optional[int] = None,
    seed: Optional[int] = None,
) -> str:
    """Append a docking result row to the aggregate CSV file.

    Reads the last row to determine the next serial number, then appends a
    row with columns ``sr``, ``name``, ``affinity`` and — when supplied — the
    search parameters the score was computed with.

    Recording these per row is not cosmetic. A Vina affinity is only
    interpretable relative to the box it was computed in, so two rows scored in
    different boxes are not comparable and no ``Δaffinity`` between them is
    meaningful. The same argument applies to the RNG seed: at fixed settings a
    seed change can move a single affinity by ~2 kcal/mol by finding a different
    binding mode, so a row without its seed cannot be reproduced or replicated.
    Without these columns a results CSV cannot be audited after the fact, and
    silently inconsistent datasets are indistinguishable from valid ones.

    Args:
        pose_file: Path to the docking pose file (used to derive the run
            name); the ``_out`` suffix and any extension are stripped, so
            ``rec_lig_out.sdf``/``.pdbqt``/``.pdb`` all yield ``rec_lig``.
        csv_file: Path to the CSV file to append to (created if absent).
        score: Binding affinity in kcal/mol.
        center: Optional ``(x, y, z)`` search-box center actually used.
        box_size: Optional ``(x, y, z)`` search-box dimensions actually used.
        exhaustiveness: Optional Vina exhaustiveness actually used.
        seed: Optional Vina RNG seed actually used.

    Returns:
        The run name derived from *pose_file*.

    Raises:
        DockingError: If writing to the CSV file fails.
    """
    import csv

    path = Path(csv_file)
    file_exists = path.exists() and path.stat().st_size > 0

    count = 1
    if file_exists:
        try:
            with open(csv_file) as lc:
                lines = [ln.strip() for ln in lc.readlines() if ln.strip()]
            # Walk backwards to the most recent data row, skipping the header
            # row (whose first field is not an integer).
            for line in reversed(lines):
                first_field = line.split(",")[0]
                try:
                    count = int(first_field) + 1
                    break
                except ValueError:
                    continue
        except Exception:
            count = 1

    # Strip the trailing "_out" and any extension so the run name is the same
    # regardless of which pose artifact (SDF / PDBQT) is passed in.
    name = Path(pose_file).stem.removesuffix("_out")

    def _fmt(v: object) -> object:
        return round(float(v), 3) if isinstance(v, int | float) else ""

    # Column order is append-only: new run parameters go on the END so that a
    # CSV written by an older version stays readable by column name.
    param_cols = [
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "exhaustiveness",
        "seed",
    ]
    c: list[Optional[float]] = list(center) if center is not None else [None] * 3
    b: list[Optional[float]] = list(box_size) if box_size is not None else [None] * 3
    param_vals = [
        _fmt(c[0]) if len(c) > 0 and c[0] is not None else "",
        _fmt(c[1]) if len(c) > 1 and c[1] is not None else "",
        _fmt(c[2]) if len(c) > 2 and c[2] is not None else "",
        _fmt(b[0]) if len(b) > 0 and b[0] is not None else "",
        _fmt(b[1]) if len(b) > 1 and b[1] is not None else "",
        _fmt(b[2]) if len(b) > 2 and b[2] is not None else "",
        exhaustiveness if exhaustiveness is not None else "",
        seed if seed is not None else "",
    ]
    full_row = [count, name, score, *param_vals]

    # Match the existing file's column count instead of producing a ragged file.
    # This covers the legacy 3-column CSVs and, now, files written before the
    # seed column existed: the row is truncated to the width already on disk.
    # Truncation loses a parameter rather than corrupting the file, so it is
    # logged -- silently dropping the seed is exactly the kind of missing
    # provenance this column was added to prevent.
    existing_cols: Optional[int] = None
    if file_exists:
        try:
            with open(csv_file) as hc:
                for ln in hc:
                    if ln.strip():
                        existing_cols = len(ln.split(","))
                        break
        except Exception:
            existing_cols = None

    row = full_row
    if existing_cols is not None and existing_cols < len(full_row):
        dropped = param_cols[existing_cols - 3 :]
        row = full_row[:existing_cols]
        logger.warning(
            "%s has %d columns; appending a %d-column row would make it ragged, "
            "so %s %s not recorded for this run. Start a new CSV to capture them.",
            Path(csv_file).name,
            existing_cols,
            len(full_row),
            ", ".join(dropped),
            "is" if len(dropped) == 1 else "are",
        )

    try:
        with open(csv_file, "a", newline="") as out:
            writer = csv.writer(out)
            if not file_exists:
                writer.writerow(["sr", "name", "affinity", *param_cols])
            writer.writerow(row)
    except Exception as e:
        raise DockingError(f"Failed to write to CSV '{csv_file}': {e}") from e

    return name


def read_config(
    config_file: str,
) -> tuple[list[float], list[float], int, int, int, bool, int]:
    """Parse a Vina configuration file and return docking parameters.

    Expects ``key = value`` lines; lines starting with ``#`` are ignored.
    Falls back to safe defaults for any missing keys.

    Args:
        config_file: Path to the Vina ``.conf`` / ``.txt`` config file.

    Returns:
        ``(center, box_size, exhaustiveness, n_poses, n_poses_write, overwrite,
        seed)``.  *seed* falls back to :data:`DEFAULT_VINA_SEED` when the config
        does not name one, so a config written before seeds were recorded still
        yields a reproducible run rather than a random one.

    Raises:
        ConfigError: If the file cannot be read or parsed.
    """
    try:
        config = {}
        with open(config_file) as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, value = line.split("=", 1)
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
        seed = int(config.get("seed", str(DEFAULT_VINA_SEED)))

        return (
            center,
            box_size,
            exhaustiveness,
            n_poses,
            n_poses_write,
            overwrite,
            seed,
        )

    except Exception as e:
        raise ConfigError(f"Failed to read config file '{config_file}': {e}") from e


def run_autosite(receptor_pdbqt: str) -> str:
    """Run the AutoSite binary on a prepared receptor PDBQT file.

    Creates an output directory named ``{stem}_autosite_out`` next to the
    receptor file and runs ``autosite -r receptor.pdbqt -o out_dir``.

    The directory is created here because AutoSite requires it to ALREADY EXIST
    and asserts otherwise (``MakeGrids.setWorkingFolder``: "The directory
    specified doesn't exist or is not accessible"). It is created with
    ``exist_ok=False``, so a pre-existing directory raises ``FileExistsError``
    rather than being silently reused with stale maps in it — callers must remove
    it first (see ``_prepare_autosite`` in benchmark2/benchmark_platinum.py and
    the AutoSite branch of np_docking).

    NOTE ON FAILURE MODES: AutoSite can exit 0 having found NOTHING. On a
    receptor with no detectable pocket it prints "analysing 0 clusters" and
    writes no ``*_cl_001.pdb``, which this function reports as a
    ``DockingRunError`` about the missing cluster file rather than a crash. That
    is a legitimate result for some structures, not necessarily a bug — verified
    on PLATINUM systems 2CC7, 2CCB, 2CCC and 1MAR, which yield 0 clusters.

    Args:
        receptor_pdbqt: Path to the prepared receptor PDBQT file.

    Returns:
        Path to the AutoSite cluster PDB file (``{stem}_cl_001.pdb``).

    Raises:
        FileExistsError: If ``{stem}_autosite_out`` already exists.
        DockingRunError: If the autosite binary is not found, the subprocess
            fails or times out, or the expected cluster PDB is absent after the
            run (including the "0 clusters found" case above).
    """
    import shutil
    import subprocess

    autosite_cmd = shutil.which("autosite")
    if autosite_cmd is None:
        raise DockingRunError("autosite binary not found on PATH.")

    receptor_path = Path(receptor_pdbqt)
    out_dir = receptor_path.parent / f"{receptor_path.stem}_autosite_out"

    out_dir.mkdir(exist_ok=False)

    try:
        result = subprocess.run(
            [autosite_cmd, "-r", receptor_path.name, "-o", out_dir.name],
            capture_output=True,
            text=True,
            cwd=str(receptor_path.parent),
            timeout=AUTOSITE_TIMEOUT,
        )
    except subprocess.TimeoutExpired as e:
        raise DockingRunError(
            f"autosite timed out after {AUTOSITE_TIMEOUT}s on "
            f"'{receptor_path.name}' and was terminated."
        ) from e
    if result.returncode != 0:
        raise DockingRunError(f"autosite failed:\n{result.stderr.strip()}")

    cluster_pdb = out_dir / f"{receptor_path.stem}_cl_001.pdb"
    if not cluster_pdb.exists():
        raise DockingRunError(f"Expected AutoSite output not found: {cluster_pdb}")

    return str(cluster_pdb)


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
    seed: int = DEFAULT_VINA_SEED,
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
        seed: Vina RNG seed (default :data:`DEFAULT_VINA_SEED`).  A *config*
            file naming ``seed`` overrides this, like every other parameter.

    Raises:
        ConfigError: If the config file cannot be read.
        mutadock.docking.exceptions.PDBFileError: If the autosite PDB cannot be parsed.
        DockingRunError: If the Vina subprocess exits with a non-zero code.
    """
    if config is not None:
        (
            center,
            box_size,
            exhaustiveness,
            n_poses,
            n_poses_write,
            overwrite,
            seed,
        ) = read_config(config)

    if autosite is not None:
        center = list(calculate_geometric_center(autosite))

    if center is None:
        center = [0.0, 0.0, 0.0]
    if box_size is None:
        box_size = [30.0, 30.0, 30.0]

    import subprocess

    vina_dock_script = Path(__file__).parent / "vina_dock.py"

    commands = [
        sys.executable,
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
        "--seed",
        str(seed),
    ]
    commands.append("--overwrite" if overwrite else "--no-overwrite")
    with open(log_file, "w+") as lfile:
        try:
            result = subprocess.run(
                commands,
                stdout=lfile,
                stderr=lfile,
                text=True,
                timeout=VINA_TIMEOUT,
            )
        except subprocess.TimeoutExpired as e:
            lfile.write(f"\nVina timed out after {VINA_TIMEOUT}s and was terminated.\n")
            raise DockingRunError(
                f"Vina timed out after {VINA_TIMEOUT}s for "
                f"'{Path(receptor).name}' + '{Path(ligand).name}'."
            ) from e

    if result.returncode != 0:
        raise DockingRunError(f"Check the error in {log_file}")


if __name__ == "__main__":
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's docking module."
    )
