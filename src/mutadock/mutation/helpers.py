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


import argparse
import csv
import logging
import os
import re
import shutil
import urllib.request
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Optional

from .Amino import check_3, get_1, get_3
from .exceptions import MutationError, ResidueMismatchError

DATA_DIR = Path(__file__).parent.parent.parent / "data"
NCBI_MATRIX_URL = "https://ftp.ncbi.nih.gov/blast/matrices/"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def add_ddg_protocol_args(parser: argparse.ArgumentParser) -> None:
    """Add the shared ΔΔG protocol flags to *parser*.

    Used by every ΔΔG CLI (``md_ddg_single``, the double/triple calculators,
    ``md_mutate``) so the protocol knobs are identical everywhere.  Pair with
    :func:`resolve_ddg_params` to turn the parsed args into ``calc_*`` kwargs.
    """
    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        metavar="N",
        help="Independent ΔΔG estimates per mutation for mean ± SD (Rosetta "
        "packing is stochastic; default 1).",
    )
    parser.add_argument(
        "--pack-radius",
        dest="pack_radius",
        type=float,
        default=8.0,
        metavar="Å",
        help="Å neighbourhood repacked/minimized around the mutation " "(default 8.0).",
    )
    parser.add_argument(
        "--protocol",
        choices=["fast", "min", "cartesian"],
        default="min",
        help="ΔΔG protocol (default: 'min'). 'min' = repack + backbone/side-chain "
        "minimization (reliable; recommended). 'cartesian' = ref2015_cart + "
        "Cartesian minimization (most accurate, cartesian_ddg style — Park et al. "
        "2016; slowest). 'fast' = single repack, NO minimization — SCREENING "
        "ONLY: fast but absolute ΔΔG values are unreliable (a failed wild-type "
        "repack can invert a site's ranking).",
    )
    parser.add_argument(
        "--reu-to-kcal",
        dest="reu_to_kcal",
        type=float,
        default=None,
        metavar="F",
        help="Override the REU→kcal/mol scaling factor "
        "(default: predict_ddG.REU_TO_KCAL_SCALE ≈ 0.34).",
    )


def resolve_ddg_params(args: argparse.Namespace) -> dict:
    """Resolve the shared ΔΔG protocol flags into ``calc_*`` keyword arguments.

    Protocol → (cartesian, backbone_minimization):

    * ``fast``      → ``(False, False)`` — single repack, no minimization
      (screening only; absolute ΔΔG unreliable).
    * ``min``       → ``(False, True)``  — repack + minimization (default).
    * ``cartesian`` → ``(True, True)``   — ref2015_cart + Cartesian minimization.

    ``reu_to_kcal`` is only included when the user supplied ``--reu-to-kcal`` so
    the ``calc_*`` default is used otherwise.
    """
    cartesian = args.protocol == "cartesian"
    backbone_min = args.protocol in ("min", "cartesian")
    params: dict = {
        "replicates": args.replicates,
        "pack_radius": args.pack_radius,
        "backbone_minimization": backbone_min,
        "cartesian": cartesian,
    }
    if args.reu_to_kcal is not None:
        params["reu_to_kcal"] = args.reu_to_kcal
    return params


class Mutation:
    """Represents a single amino-acid substitution (chain, position, old→new)."""

    def __init__(self, chain: str, position: int, aa_old: str, aa_new: str) -> None:
        """
        Args:
            chain: PDB chain identifier (e.g. ``"A"``).
            position: Residue sequence number.
            aa_old: Wild-type amino acid (1- or 3-letter code).
            aa_new: Mutant amino acid (1- or 3-letter code).
        """
        self.chain = chain
        self.position = int(position)
        if check_3(aa_old):
            _old = get_1(aa_old)
            self.old = _old if _old is not None else aa_old
        else:
            self.old = aa_old
        if check_3(aa_new):
            _new = get_1(aa_new)
            self.new = _new if _new is not None else aa_new
        else:
            self.new = aa_new

    def __str__(self) -> str:
        return f"The aa at {self.chain}{self.position} is mutated from {self.old} to {self.new}."


@contextmanager
def in_directory(path: str | Path):
    """Temporarily change the working directory.

    Rosetta tokenizes file paths on whitespace, so paths containing spaces
    cause it to fail even when the file exists.  This context manager
    changes to the target directory, yields, then restores the original
    directory — allowing callers to pass only the basename to Rosetta.
    """
    old_dir = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(old_dir)


def backup(file_path: str) -> bool:
    """Move *file_path* to a timestamped copy inside a ``backups/`` sub-folder.

    Args:
        file_path: Path to the file to back up.

    Returns:
        ``True`` if the file was moved, ``False`` if it did not exist.
    """
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


def move_file(file_path: str, folder: str) -> tuple[str, str]:
    """Move *file_path* into *folder* (created next to the file if absent).

    Args:
        file_path: Path to the file to move.
        folder: Sub-folder name relative to the file's parent directory.

    Returns:
        ``(folder_path, destination_path)`` as absolute strings.
    """
    path = Path(file_path).resolve()
    folder_path = path.parent / folder
    folder_path.mkdir(exist_ok=True)
    dest = folder_path / path.name
    shutil.move(file_path, dest)
    return (str(folder_path), str(dest))


def process_vina(file: str) -> str:
    """Delete extra Vina pose files and return the path to pose 1.

    Removes ``<file>_ligand_2.pdbqt`` … ``<file>_ligand_9.pdbqt`` (if they
    exist) and returns the expected path of the first pose.

    Args:
        file: Base path used by Vina when naming pose files.

    Returns:
        Path string for ``<file>_ligand_1.pdbqt``.
    """
    path = Path(file).resolve()
    ligand1 = path.parent / f"{path.name}_ligand_1.pdbqt"
    for i in range(2, 10):
        candidate = path.parent / f"{path.name}_ligand_{i}.pdbqt"
        try:
            candidate.unlink()
        except FileNotFoundError:
            pass
    return str(ligand1)


def file_info(file_path: str) -> tuple[bool, int]:
    """Return existence and line count for a file.

    Args:
        file_path: Path to inspect.

    Returns:
        ``(exists, line_count)`` — ``(False, 0)`` if the file does not exist.
    """
    path = Path(file_path)
    if path.is_file():
        row_count = sum(1 for _ in path.open("r"))
        return (True, row_count)
    return (False, 0)


def read_partial_ddg(
    out_file: str,
    *,
    must_match: Optional[dict[str, object]] = None,
) -> Optional[list[dict[str, str]]]:
    """Read the completed rows of a partially-written ΔΔG CSV for per-item resume.

    The ΔΔG output CSV acts as its own checkpoint: every fully-written row is a
    completed item.  This returns the list of complete data rows already present
    so a caller can rewrite them verbatim and skip recomputing them, computing
    only the items that are missing.

    A row is treated as *complete* only when its final column (``ddG_protocol``)
    is non-empty, which drops a torn trailing line left by a mid-write
    interruption; the caller recomputes that one item cleanly.

    Args:
        out_file: Path to the (possibly partial) ΔΔG CSV.
        must_match: Column→value pairs every completed row must equal (compared
            as strings) for the file to be resumable, e.g.
            ``{"ddG_protocol": protocol, "n_replicates": replicates}``.  This
            prevents silently mixing rows computed under different settings.

    Returns:
        The list of complete row dicts, or ``None`` when the file is absent,
        unreadable, empty/header-only, or disagrees with *must_match* — in which
        case the caller should start a fresh calculation.
    """
    path = Path(out_file)
    if not path.is_file():
        return None
    try:
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "ddG_protocol" not in reader.fieldnames:
                return None
            rows: list[dict[str, str]] = []
            for row in reader:
                if not row.get("ddG_protocol"):
                    # Torn/incomplete trailing line — recompute this item.
                    continue
                if must_match:
                    for col, val in must_match.items():
                        if str(row.get(col, "")) != str(val):
                            # Settings changed since the partial run: mixing rows
                            # would be invalid, so start over from scratch.
                            return None
                rows.append(row)
    except OSError:
        return None
    return rows or None


def permutations(n: int, r: int) -> int:
    """Return P(n, r) — the number of ordered r-permutations of n items.

    Args:
        n: Pool size.
        r: Selection size.

    Returns:
        ``n! / (n-r)!``, or ``0`` when ``n < r``.
    """
    from math import factorial

    if n < r:
        return 0
    else:
        return int(factorial(n) / factorial(n - r))


def clean_pdb(
    pdb_file: str, pdb_clean: Optional[str] = None, keep_hetatm: bool = False
) -> str:
    """Strip a PDB file down to ATOM/TER records, terminated by a single END.

    Every ``ATOM`` (and ``TER``) record from the input is retained, regardless
    of how many models the file contains.  Any ``END`` / ``ENDMDL`` markers in
    the source are dropped and a single ``END`` line is appended at the very
    end.  This matters because ``END`` is the PDB end-of-file marker: re-emitting
    a stray or per-model ``END`` mid-file would cause any downstream parser
    (e.g. PyRosetta's ``pose_from_pdb``) to silently truncate the structure at
    the first one — the exact bug this function must avoid.

    ``HETATM`` records (waters, ions, ligands, cofactors) are stripped by
    default.  This is deliberate: mutation scoring runs on the apo protein, so
    heteroatoms are intentionally excluded.  Pass ``keep_hetatm=True`` only if
    you explicitly need them retained.

    Args:
        pdb_file: Path to the source PDB file.
        pdb_clean: Destination path.  Defaults to ``<stem>_clean.pdb``.
        keep_hetatm: If ``True``, ``HETATM`` records are preserved.  Defaults
            to ``False`` (strip them — the correct behaviour for apo scoring).

    Returns:
        Path to the cleaned PDB file.

    Raises:
        MutationError: If the file cannot be read or written.
    """
    if not pdb_clean:
        pdb_clean = pdb_file.replace(".pdb", "_clean.pdb")

    try:
        with open(pdb_file) as infile, open(pdb_clean, "w+") as outfile:
            outfile.write(
                "REMARK This file was cleaned using the clean_pdb function of mutadock library.\n"
            )
            outfile.write("REMARK All non-ATOM lines were removed.\n")
            outfile.write(
                f"REMARK Created on {datetime.now().strftime('%b-%d-%Y %H:%M:%S')}\n"
            )
            for line in infile:
                if line.startswith("ATOM") or line.startswith("TER"):
                    outfile.write(line)
                elif keep_hetatm and line.startswith("HETATM"):
                    outfile.write(line)
            # A single terminal END; intermediate END/ENDMDL markers are dropped
            # so downstream parsers never truncate at an early end-of-file marker.
            outfile.write("END\n")
        return pdb_clean
    except Exception as e:
        raise MutationError(f"Failed to clean PDB file '{pdb_file}': {e}") from e


def convert_cif_pdb(cif_file: str, pdb_file: Optional[str] = None) -> str:
    """Convert a mmCIF file to PDB format using BioPython.

    Args:
        cif_file: Path to the source ``.cif`` file.
        pdb_file: Destination path.  Defaults to replacing ``.cif`` with ``.pdb``.

    Returns:
        Path to the written PDB file.

    Raises:
        MutationError: If conversion fails.
    """
    from Bio.PDB import PDBIO, MMCIFParser

    if not pdb_file:
        pdb_file = cif_file.replace(".cif", ".pdb")
    parser = MMCIFParser()
    try:
        structure = parser.get_structure("structure", cif_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        return pdb_file
    except Exception as e:
        raise MutationError(f"Failed to convert CIF to PDB '{cif_file}': {e}") from e


def fetch_pdb(pdb_id: str, dest_dir: Optional[str | Path] = None) -> str:
    """Download a structure from RCSB by its 4-character PDB ID.

    Args:
        pdb_id: A 4-character PDB accession (e.g. ``"4QJR"``); case-insensitive.
        dest_dir: Directory to save ``<PDB_ID>.pdb`` into.  Defaults to the
            current working directory.  An existing file is reused (not
            re-downloaded).

    Returns:
        Path to the downloaded (or already-present) ``.pdb`` file.

    Raises:
        MutationError: If *pdb_id* is malformed or the download fails.
    """
    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9A-Z]{4}", pdb_id):
        raise MutationError(
            f"Invalid PDB ID '{pdb_id}': expected a 4-character code such as 4QJR."
        )
    directory = Path(dest_dir) if dest_dir else Path.cwd()
    directory.mkdir(parents=True, exist_ok=True)
    dest = directory / f"{pdb_id}.pdb"
    if dest.is_file():
        logger.info("PDB %s already present at '%s'; skipping download.", pdb_id, dest)
        return str(dest)

    url = RCSB_PDB_URL.format(pdb_id=pdb_id)
    logger.info("Fetching PDB %s from %s ...", pdb_id, url)
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        raise MutationError(
            f"Failed to fetch PDB '{pdb_id}' from RCSB ({url}): {e}"
        ) from e
    logger.info("Saved to '%s'.", dest)
    return str(dest)


def structure_warnings(pdb_file: str) -> list[str]:
    """Detect structural features that silently affect residue numbering.

    Scans a PDB file for multi-model ensembles (e.g. NMR), alternate location
    indicators (altlocs), and residue insertion codes — each of which can make
    the integer residue numbering used for mutations ambiguous.  Must be run on
    the *original* structure: ``clean_pdb`` drops ``MODEL`` records, so a cleaned
    file would never report a multi-model warning.

    Args:
        pdb_file: Path to the PDB file to inspect.

    Returns:
        A (possibly empty) list of human-readable warning strings.
    """
    model_count = 0
    altlocs: set[str] = set()
    insertion_residues: set[tuple[str, str, str]] = set()
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("MODEL"):
                    model_count += 1
                elif line.startswith(("ATOM", "HETATM")):
                    altloc = line[16:17].strip()
                    if altloc:
                        altlocs.add(altloc)
                    icode = line[26:27].strip()
                    if icode:
                        chain = line[21:22]
                        resseq = line[22:26].strip()
                        insertion_residues.add((chain, resseq, icode))
    except OSError:
        return []

    warnings: list[str] = []
    if model_count > 1:
        warnings.append(
            f"Structure contains {model_count} models (e.g. an NMR ensemble); "
            "residue numbering can differ between models and only one is scored."
        )
    if altlocs:
        warnings.append(
            f"Structure contains alternate location indicators "
            f"({', '.join(sorted(altlocs))}); alternate conformations may affect "
            "mutation scoring."
        )
    if insertion_residues:
        warnings.append(
            f"Structure contains {len(insertion_residues)} residue(s) with "
            "insertion codes; insertion codes are not reflected in the integer "
            "residue numbering used for mutations."
        )
    return warnings


def pose_residue_map(pose: object) -> dict[str, list[int]]:
    """Return ``{chain: [residue numbers]}`` for a PyRosetta pose.

    Uses the pose's PDB info so the numbers match the PDB numbering that
    ``pdb2pose`` expects.  Used to build actionable "residue not found" errors.
    """
    info = pose.pdb_info()  # type: ignore[attr-defined]
    residues: dict[str, list[int]] = {}
    for i in range(1, pose.total_residue() + 1):  # type: ignore[attr-defined]
        residues.setdefault(info.chain(i), []).append(info.number(i))
    return residues


def format_missing_residue(
    available: dict[str, list[int]], chain: str, position: int, pdb_file: str
) -> str:
    """Build an actionable message for a residue that isn't in the structure.

    Args:
        available: ``{chain: [residue numbers]}`` present in the structure.
        chain: Requested chain ID.
        position: Requested residue number.
        pdb_file: Structure path, for context in the message.

    Returns:
        A message listing available chains (if the chain is absent) or the
        numbering span of the requested chain (if only the position is absent).
    """
    if chain not in available:
        chains = ", ".join(sorted(available)) or "none"
        return (
            f"Chain '{chain}' not found in '{pdb_file}'. "
            f"Available chains: {chains}."
        )
    positions = available[chain]
    lo, hi = min(positions), max(positions)
    return (
        f"Residue {position} not found in chain '{chain}' of '{pdb_file}'. "
        f"Chain '{chain}' spans positions {lo}-{hi} ({len(positions)} residues)."
    )


def load_matrix(
    file_path: str | Path = DATA_DIR / "PAM250",
) -> dict[str, dict[str, int]]:
    """Parse an NCBI-format substitution matrix file into 3-letter-coded dicts.

    Lines starting with ``#`` are ignored.  The first non-comment line is
    treated as the column header (one-letter codes).  Rows whose leading
    letter does not map to a standard amino acid (B, Z, X, ``*``) are skipped.
    """
    path = Path(file_path)
    columns: list[str] = []
    result: dict[str, dict[str, int]] = {}

    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if not columns:
                columns = parts
                continue
            three = get_3(parts[0])
            if three is None:
                continue
            result[three] = {}
            for col, val in zip(columns, parts[1:], strict=False):
                col_three = get_3(col)
                if col_three is None:
                    continue
                result[three][col_three] = int(val)

    return result


def download_matrix(name: str) -> Path:
    """Download *name* from the NCBI BLAST matrices FTP into ``data/``."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    dest = DATA_DIR / name
    url = NCBI_MATRIX_URL + name
    logger.info("Downloading matrix '%s' from %s ...", name, url)
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        raise MutationError(f"Failed to download matrix '{name}': {e}") from e
    logger.info("Saved to '%s'", dest)
    return dest


def resolve_matrix(
    name: str = "PAM250", custom_file: Optional[str] = None
) -> dict[str, dict[str, int]]:
    """Return a scoring matrix, downloading from NCBI if not already in ``data/``.

    If *custom_file* is given it is loaded directly, ignoring *name*.
    Otherwise the matrix is looked up in ``data/``; if absent it is downloaded
    from ``https://ftp.ncbi.nih.gov/blast/matrices/``.
    """
    if custom_file:
        path = Path(custom_file)
        if not path.is_file():
            raise MutationError(f"Custom matrix file not found: '{custom_file}'")
        return load_matrix(path)

    name = name.upper()
    path = DATA_DIR / name
    if not path.is_file():
        path = download_matrix(name)
    return load_matrix(path)


_rosetta_init_done: bool = False


def mutate_and_score(
    pdb_file: str,
    chain: str,
    aa_number: int,
    orig_aa: str,
    mutated_aa: str,
    output_file: str,
    pack_radius: float = 8.0,
    backbone_minimization: bool = False,
    cartesian: bool = False,
) -> tuple[float, float, float]:
    """Apply a single-point mutation, score both WT and mutant, and write the mutant PDB.

    Initialises PyRosetta once per process (subsequent calls reuse the existing
    session), loads *pdb_file* as a full-atom pose, and computes ΔΔG the same
    way :func:`ddg_calc.calc_ddg` does: both the wild-type reference and the
    mutant undergo the *identical* repack (and optional minimization) protocol
    via :func:`predict_ddG.wt_reference_score` / :func:`predict_ddG.apply_mutations`,
    so the difference is not biased by repacking only one side. The mutant
    coordinates are dumped to *output_file*.

    The residue actually present at (*chain*, *aa_number*) is verified against
    *orig_aa* before mutating. If it doesn't match *orig_aa* but does match
    *mutated_aa*, the structure is treated as already being the mutant (e.g. a
    deposited mutant/co-crystal structure rather than wild-type): the mutant
    score/PDB come from a self-mutation (consistent packing, no identity
    change) and the wild-type reference comes from one reverse mutation back
    to *orig_aa* — avoiding a wasted reverse-then-forward round trip through
    the wild-type identity. If it matches neither, :class:`ResidueMismatchError`
    is raised.

    Args:
        pdb_file: Path to the cleaned wild-type PDB (ATOM records only).
        chain: PDB chain identifier, e.g. ``"A"``.
        aa_number: Residue sequence number to mutate.
        orig_aa: Expected wild-type amino acid (1- or 3-letter code); verified
            against the structure at (*chain*, *aa_number*) before mutating.
        mutated_aa: Target amino acid (1- or 3-letter code).
        output_file: Destination path for the mutant PDB file.
            Parent directories are created if absent.
        pack_radius: Å neighbourhood repacked (and minimized) around the site.
        backbone_minimization: Minimize backbone + side chains after packing.
        cartesian: Use the Cartesian ``ref2015_cart`` protocol (implies/enables
            Cartesian minimization — the ``cartesian_ddg`` style).

    Returns:
        ``(wt_score, mut_score, ddG)`` where ``ddG = mut_score - wt_score``
        in Rosetta energy units.

    Raises:
        MutationError: If PyRosetta is unavailable, the residue is not found
            in the pose, or the amino-acid code cannot be resolved.
        ResidueMismatchError: If the residue at (*chain*, *aa_number*) matches
            neither *orig_aa* nor *mutated_aa*.
    """
    global _rosetta_init_done
    import os
    import sys as _sys

    try:
        _sys.stdout = open(os.devnull, "w")
        from pyrosetta import get_fa_scorefxn, init, pose_from_pdb  # noqa: F401

        _sys.stdout = _sys.__stdout__
    except ImportError as err:
        _sys.stdout = _sys.__stdout__
        raise MutationError(
            "PyRosetta not found. Install it with:\n"
            "  python -m pip install pyrosetta_installer && "
            "python3 -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'"
        ) from err

    from . import predict_ddG  # local import avoids circular dependency at module load

    if not _rosetta_init_done:
        _sys.stdout = open(os.devnull, "w")
        init("-mute all")
        _sys.stdout = _sys.__stdout__
        _rosetta_init_done = True

    pdb_path = Path(pdb_file).resolve()
    with in_directory(pdb_path.parent):
        pose = pose_from_pdb(pdb_path.name)

    sfxn = predict_ddG.get_scorefxn(cartesian)

    pose_position: int = pose.pdb_info().pdb2pose(chain, aa_number)
    if pose_position == 0:
        raise MutationError(
            format_missing_residue(pose_residue_map(pose), chain, aa_number, pdb_file)
        )

    mut_1: Optional[str] = get_1(mutated_aa) if check_3(mutated_aa) else mutated_aa
    if mut_1 is None:
        raise MutationError(f"Unknown amino acid code: '{mutated_aa}'")

    orig_1: Optional[str] = get_1(orig_aa) if check_3(orig_aa) else orig_aa
    native_aa = pose.residue(pose_position).name1()

    if orig_1 is not None and native_aa != orig_1 and native_aa != mut_1:
        raise ResidueMismatchError(
            f"Expected wild-type '{orig_1}' at {chain}{aa_number} in '{pdb_file}', "
            f"but the structure has '{native_aa}'. The (chain, position) numbering "
            "likely doesn't match this structure's own residue numbering."
        )

    if orig_1 is not None and native_aa == mut_1 and native_aa != orig_1:
        # The structure already has the target mutant residue at this position
        # (e.g. a deposited mutant/co-crystal structure rather than wild-type).
        # Treat it as the native MT endpoint: self-mutate (consistent packing,
        # no-op identity-wise) for the mutant score/PDB, and reverse-mutate to
        # the wild-type letter for the reference — avoids a wasted
        # reverse-then-forward round trip through the wild-type identity.
        logger.info(
            "%s chain=%s pos=%d: structure already has target residue '%s' "
            "(expected wild-type '%s'); treating as native mutant and "
            "reverse-mutating for the wild-type reference.",
            pdb_path.name,
            chain,
            aa_number,
            native_aa,
            orig_1,
        )
        mut_pose = predict_ddG.apply_mutations(
            pose,
            [(pose_position, native_aa)],
            sfxn,
            pack_radius,
            backbone_minimization=backbone_minimization,
            cartesian=cartesian,
        )
        mut_score: float = sfxn.score(mut_pose)

        wt_pose = predict_ddG.apply_mutations(
            pose,
            [(pose_position, orig_1)],
            sfxn,
            pack_radius,
            backbone_minimization=backbone_minimization,
            cartesian=cartesian,
        )
        wt_score: float = sfxn.score(wt_pose)
    else:
        wt_score = predict_ddG.wt_reference_score(
            pose,
            [pose_position],
            sfxn,
            pack_radius,
            backbone_minimization=backbone_minimization,
            cartesian=cartesian,
        )

        mut_pose = predict_ddG.apply_mutations(
            pose,
            [(pose_position, mut_1)],
            sfxn,
            pack_radius,
            backbone_minimization=backbone_minimization,
            cartesian=cartesian,
        )
        mut_score = sfxn.score(mut_pose)

    ddg: float = mut_score - wt_score

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mut_pose.dump_pdb(str(out_path))

    logger.info(
        "%s chain=%s pos=%d %s->%s | WT=%.3f Mut=%.3f ddG=%.3f",
        pdb_path.name,
        chain,
        aa_number,
        orig_aa,
        mut_1,
        wt_score,
        mut_score,
        ddg,
    )
    return (wt_score, mut_score, ddg)


if __name__ == "__main__":
    logger.info(
        "This is a dependency file for mutadock (https://github.com/naisarg14/mutadock) library's mutation module."
    )
