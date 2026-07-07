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
"""Lightweight structure helpers for MUTADOCK reports.

Dependency-light on purpose: only ``Bio.PDB`` and the standard library. No
matplotlib here — plotting lives in ``figures.py`` (owned by another module).
"""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBParser

# Ordered CA record: (resseq, x, y, z, resname)
CATrace = dict[str, list[tuple[int, float, float, float, str]]]


def read_pdb_text(pdb_path: str | Path) -> str:
    """Return the full text of *pdb_path* (for inlining into an NGL viewer).

    Raises:
        FileNotFoundError: if the file does not exist.
    """
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")
    return path.read_text()


def parse_ca_trace(pdb_path: str | Path) -> CATrace:
    """Extract the per-chain CA trace from a PDB file.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        A mapping of chain id to an ordered list of
        ``(resseq, x, y, z, resname)`` tuples for the CA atoms of standard
        (non-hetero, non-water) residues. Chains without CA atoms are omitted.

    Raises:
        FileNotFoundError: if the file does not exist.
    """
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")

    parser = PDBParser(QUIET=True)
    trace: CATrace = {}
    try:
        structure = parser.get_structure("s", str(path))
    except Exception:
        # Malformed file → gracefully return whatever (nothing) we have.
        return trace

    try:
        model = next(iter(structure))
    except StopIteration:
        return trace

    for chain in model:
        residues: list[tuple[int, float, float, float, str]] = []
        for residue in chain:
            # residue.id == (hetfield, resseq, icode); skip hetero/water.
            hetfield = residue.id[0]
            if hetfield != " ":
                continue
            if "CA" not in residue:
                continue
            atom = residue["CA"]
            x, y, z = (float(c) for c in atom.get_coord())
            residues.append((int(residue.id[1]), x, y, z, residue.get_resname()))
        if residues:
            trace[chain.id] = residues

    return trace


def mutation_sites(row: dict) -> list[tuple[str, int]]:
    """Extract ``(chain, position)`` pairs from a mutation row.

    Handles single rows (``chain``/``position``) and double/triple rows
    (``chain1``/``mut1_position``, ``chain2``/``mut2_position``,
    ``chain3``/``mut3_position``). Only pairs actually present are returned.
    """
    sites: list[tuple[str, int]] = []

    key_pairs = [
        ("chain", "position"),
        ("chain1", "mut1_position"),
        ("chain2", "mut2_position"),
        ("chain3", "mut3_position"),
    ]
    for chain_key, pos_key in key_pairs:
        if chain_key in row and pos_key in row:
            chain = row[chain_key]
            pos = row[pos_key]
            if chain is None or pos is None:
                continue
            try:
                sites.append((str(chain), int(pos)))
            except (ValueError, TypeError):
                continue

    return sites
