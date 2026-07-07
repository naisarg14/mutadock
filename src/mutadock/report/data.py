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
"""Data discovery and loading for MUTADOCK reports.

This module is the shared data contract for the ``mutadock.report`` subpackage.
It turns a completed MUTADOCK run directory into a single, well-typed
:class:`RunData` object that the figure, HTML and PPTX generators build on.

The discovery is deliberately tolerant: files are located by glob/substring
patterns (never by assuming an exact stem), pre-sorted files are *re-sorted*
ascending by their own ΔΔG column (the pipeline has a known quirk that sorts the
triple CSV by the wrong column), and empty-but-present CSVs are treated as *no
data* (``None``) rather than empty frames that downstream code might misread.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime
from importlib.metadata import version as _pkg_version
from pathlib import Path
from typing import Any

import pandas as pd

# ---------------------------------------------------------------------------
# ΔΔG column names (matched by name, not position)
# ---------------------------------------------------------------------------
SINGLE_DDG_COL = "ddG_value"
DOUBLE_DDG_COL = "double_ddG_value"
TRIPLE_DDG_COL = "triple_ddG_value"
AFFINITY_COL = "affinity"

# Optional companion columns (present on runs produced after the replicate /
# kcal-scaling feature; absent on older runs — the report degrades gracefully).
SINGLE_DDG_SD_COL = "ddG_sd"
SINGLE_DDG_KCAL_COL = "ddG_kcal"


@dataclass
class RunData:
    """A fully-loaded MUTADOCK run, ready for reporting.

    Every ``*_df`` field is ``None`` when the corresponding file is missing or
    contains zero data rows. All ΔΔG frames are re-sorted ascending by their own
    value column (most-stabilizing first); ``docking_df`` is sorted ascending by
    ``affinity`` (best binder first).
    """

    output_dir: Path
    base_name: str
    input_pdb: Path | None
    mutations_df: pd.DataFrame | None
    single_df: pd.DataFrame | None
    double_df: pd.DataFrame | None
    triple_df: pd.DataFrame | None
    docking_df: pd.DataFrame | None
    mutant_pdbs: list[Path]
    top_mutant_pdb: Path | None
    top_single_row: dict | None
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
_SINGLE_MUTANT_RE = re.compile(r"^\d+_[A-Z]{3}-.*\.pdb$")

# Suffix tokens stripped (repeatedly, from the end) to recover the run stem.
_STEM_SUFFIXES = (
    "_sorted",
    "_double_ddg",
    "_triple_ddg",
    "_ddg",
    "_mutations_all",
    "_mutations",
    "_clean",
    "_all",
)


def _to_py(value: Any) -> Any:
    """Coerce a numpy/pandas scalar to a plain Python object."""
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return item()
        except (ValueError, TypeError):
            return value
    return value


def _row_to_dict(row: pd.Series) -> dict:
    """Convert a DataFrame row to a plain dict with native Python scalars."""
    return {k: _to_py(v) for k, v in row.to_dict().items()}


def _infer_stem(filename: str) -> str:
    """Recover the run stem (e.g. ``4QJR``) from a discovered filename."""
    stem = Path(filename).stem
    changed = True
    while changed:
        changed = False
        low = stem.lower()
        for suf in _STEM_SUFFIXES:
            if low.endswith(suf):
                stem = stem[: len(stem) - len(suf)]
                changed = True
                low = stem.lower()
    return stem


def _prefer_sorted(candidates: list[Path]) -> Path | None:
    """Return the ``*_sorted`` variant if present, else the first candidate."""
    if not candidates:
        return None
    for path in candidates:
        if path.stem.lower().endswith("_sorted"):
            return path
    return candidates[0]


def _classify_csvs(output_dir: Path) -> dict[str, Path | None]:
    """Group the CSV files in *output_dir* by role (single/double/triple/...).

    Matching is done on lower-cased filenames so it is robust on both
    case-sensitive and case-insensitive filesystems.
    """
    single: list[Path] = []
    double: list[Path] = []
    triple: list[Path] = []
    mutations: list[Path] = []
    mutations_all: list[Path] = []
    docking: list[Path] = []

    for path in sorted(output_dir.glob("*.csv")):
        low = path.name.lower()
        if "double_ddg" in low:
            double.append(path)
        elif "triple_ddg" in low:
            triple.append(path)
        elif "_ddg" in low:
            single.append(path)
        elif "mutations_all" in low:
            mutations_all.append(path)
        elif "mutations" in low:
            mutations.append(path)
        elif "docking" in low:
            docking.append(path)

    # Prefer the plain mutations file over the (larger) *_all variant.
    mutations_pick = _prefer_sorted(mutations) or _prefer_sorted(mutations_all)

    return {
        "single": _prefer_sorted(single),
        "double": _prefer_sorted(double),
        "triple": _prefer_sorted(triple),
        "mutations": mutations_pick,
        "docking": _prefer_sorted(docking),
    }


def _load_csv(path: Path | None) -> pd.DataFrame | None:
    """Load a CSV, returning ``None`` for missing/empty (zero data-row) files."""
    if path is None or not path.exists():
        return None
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None
    except (OSError, ValueError):
        return None
    if df.empty:
        return None
    return df


def _sort_ascending(df: pd.DataFrame | None, column: str) -> pd.DataFrame | None:
    """Re-sort *df* ascending by *column* (by name). No-op if column absent."""
    if df is None or column not in df.columns:
        return df
    return df.sort_values(column, ascending=True, kind="mergesort").reset_index(
        drop=True
    )


def _collect_mutant_pdbs(output_dir: Path) -> list[Path]:
    """Union mutant PDBs from ``*_mutants.txt`` and any ``mutation_*`` folder."""
    ordered: list[Path] = []
    seen: set[str] = set()

    def _add(path: Path) -> None:
        key = str(path.resolve())
        if key not in seen and path.exists():
            seen.add(key)
            ordered.append(path)

    for txt in sorted(output_dir.glob("*_mutants.txt")):
        try:
            lines = txt.read_text().splitlines()
        except OSError:
            continue
        for line in lines:
            line = line.strip()
            if line:
                _add(Path(line))

    for folder in sorted(output_dir.glob("mutation_*")):
        if folder.is_dir():
            for pdb in sorted(folder.glob("*.pdb")):
                _add(pdb)

    return ordered


def _resolve_top_mutant(
    top_row: dict | None,
    output_dir: Path,
    mutant_pdbs: list[Path],
) -> Path | None:
    """Find the mutant PDB corresponding to the most-stabilizing single row."""
    folders = [f for f in sorted(output_dir.glob("mutation_*")) if f.is_dir()]

    # 1) Construct the expected single-mutant filename and match it exactly.
    if top_row is not None:
        try:
            expected = (
                f"{int(top_row['sr'])}_{top_row['wtAA']}-"
                f"{top_row['chain']}{int(top_row['position'])}-"
                f"{top_row['prAA']}.pdb"
            )
        except (KeyError, ValueError, TypeError):
            expected = None
        if expected:
            for folder in folders:
                candidate = folder / expected
                if candidate.exists():
                    return candidate

    # 2) Fall back to the first file that looks like a single-mutant PDB.
    for folder in folders:
        for pdb in sorted(folder.glob("*.pdb")):
            if _SINGLE_MUTANT_RE.match(pdb.name):
                return pdb

    # 3) Fall back to the first known mutant PDB, else None.
    return mutant_pdbs[0] if mutant_pdbs else None


def _pkg_ver(name: str) -> str:
    """Best-effort package version lookup (no heavy import / init)."""
    try:
        return _pkg_version(name)
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def discover_run(
    output_dir: str | Path,
    base_name: str | None = None,
    timestamp: str | None = None,
) -> RunData:
    """Discover and load a completed MUTADOCK run from *output_dir*.

    Args:
        output_dir: Directory containing the run outputs (str or ``Path``).
        base_name: Run stem (e.g. ``"4QJR"``). Inferred from filenames if
            omitted.
        timestamp: Value for ``metadata["generated_at"]``. Defaults to
            ``datetime.now().isoformat(timespec="seconds")`` when ``None``.

    Returns:
        A populated :class:`RunData`.
    """
    output_dir = Path(output_dir)

    picks = _classify_csvs(output_dir)

    mutations_df = _load_csv(picks["mutations"])
    single_df = _sort_ascending(_load_csv(picks["single"]), SINGLE_DDG_COL)
    double_df = _sort_ascending(_load_csv(picks["double"]), DOUBLE_DDG_COL)
    triple_df = _sort_ascending(_load_csv(picks["triple"]), TRIPLE_DDG_COL)
    docking_df = _sort_ascending(_load_csv(picks["docking"]), AFFINITY_COL)

    # Infer the run stem from the cleanest available filename.
    if base_name is None:
        for role in ("mutations", "single", "double", "triple"):
            pick = picks[role]
            if pick is not None:
                base_name = _infer_stem(pick.name)
                break
        if base_name is None:
            base_name = output_dir.name

    # Preferred input structure: <stem>_clean.pdb, then <stem>.pdb.
    input_pdb: Path | None = None
    for cand in (
        output_dir / f"{base_name}_clean.pdb",
        output_dir / f"{base_name}.pdb",
    ):
        if cand.exists():
            input_pdb = cand
            break

    mutant_pdbs = _collect_mutant_pdbs(output_dir)

    top_single_row: dict | None = None
    if single_df is not None and not single_df.empty:
        top_single_row = _row_to_dict(single_df.iloc[0])

    top_mutant_pdb = _resolve_top_mutant(top_single_row, output_dir, mutant_pdbs)

    if timestamp is None:
        timestamp = datetime.now().isoformat(timespec="seconds")

    discovered = sorted(p.name for p in output_dir.glob("*") if p.is_file())

    metadata: dict = {
        "generated_at": timestamp,
        "input": str(input_pdb) if input_pdb is not None else None,
        "output_dir": str(output_dir),
        "files": discovered,
        "versions": {
            "pyrosetta": _pkg_ver("pyrosetta"),
            "vina": _pkg_ver("vina"),
        },
    }

    return RunData(
        output_dir=output_dir,
        base_name=base_name,
        input_pdb=input_pdb,
        mutations_df=mutations_df,
        single_df=single_df,
        double_df=double_df,
        triple_df=triple_df,
        docking_df=docking_df,
        mutant_pdbs=mutant_pdbs,
        top_mutant_pdb=top_mutant_pdb,
        top_single_row=top_single_row,
        metadata=metadata,
    )


def mutation_label(row: dict) -> str:
    """Return a human label for a single-mutation row, e.g. ``ASN-A386-HIS``."""
    wt = row.get("wtAA", "")
    chain = row.get("chain", "")
    prot = row.get("prAA", "")
    pos = row.get("position", "")
    try:
        pos = int(pos)
    except (ValueError, TypeError):
        pass
    return f"{wt}-{chain}{pos}-{prot}"


def summary_stats(run: RunData) -> dict:
    """Compute headline numbers shared by the HTML and PPTX reports."""
    mutations_df = run.mutations_df
    single_df = run.single_df
    docking_df = run.docking_df

    # n_scored: prefer mutations_df, fall back to single_df.
    if mutations_df is not None:
        n_scored = int(len(mutations_df))
    elif single_df is not None:
        n_scored = int(len(single_df))
    else:
        n_scored = 0

    n_stabilizing = 0
    best_ddg: float | None = None
    best_ddg_label: str | None = None
    best_ddg_kcal: float | None = None
    reu_to_kcal: float | None = None
    if single_df is not None and SINGLE_DDG_COL in single_df.columns:
        n_stabilizing = int((single_df[SINGLE_DDG_COL] < 0).sum())
        best_ddg = float(single_df[SINGLE_DDG_COL].min())
        best_row = _row_to_dict(single_df.iloc[0])
        best_ddg_label = mutation_label(best_row)
        # Recover the run's actual REU→kcal factor from the data (reflects any
        # --reu-to-kcal override) using any row with a non-zero ΔΔG.
        if SINGLE_DDG_KCAL_COL in single_df.columns:
            reu = pd.to_numeric(single_df[SINGLE_DDG_COL], errors="coerce")
            kcal = pd.to_numeric(single_df[SINGLE_DDG_KCAL_COL], errors="coerce")
            nonzero = reu.abs() > 1e-9
            if bool(nonzero.any()):
                reu_to_kcal = float((kcal[nonzero] / reu[nonzero]).median())
                best_ddg_kcal = best_ddg * reu_to_kcal

    n_replicates: int | None = None
    protocol: str | None = None
    if single_df is not None and "n_replicates" in single_df.columns:
        reps = pd.to_numeric(single_df["n_replicates"], errors="coerce").dropna()
        if not reps.empty:
            n_replicates = int(reps.max())
    if single_df is not None and "ddG_protocol" in single_df.columns:
        vals = single_df["ddG_protocol"].dropna()
        if not vals.empty:
            protocol = str(vals.iloc[0])

    n_docked = 0
    best_affinity: float | None = None
    best_affinity_name: str | None = None
    if docking_df is not None and AFFINITY_COL in docking_df.columns:
        n_docked = int(len(docking_df))
        best_affinity = float(docking_df[AFFINITY_COL].min())
        if "name" in docking_df.columns:
            best_affinity_name = str(docking_df.iloc[0]["name"])

    return {
        "n_scored": n_scored,
        "n_stabilizing": n_stabilizing,
        "best_ddg": best_ddg,
        "best_ddg_label": best_ddg_label,
        "best_ddg_kcal": best_ddg_kcal,
        "reu_to_kcal": reu_to_kcal,
        "n_replicates": n_replicates,
        "protocol": protocol,
        "n_docked": n_docked,
        "best_affinity": best_affinity,
        "best_affinity_name": best_affinity_name,
    }
