"""Fixtures for the report tests.

Builds a minimal-but-schema-correct synthetic MUTADOCK run under ``tmp_path``,
and exposes the real sample run directory when it happens to be present.

The synthetic stem is ``prot`` (NOT ``4QJR``) on purpose, to prove that
``base_name`` inference is not hard-coded.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

STEM = "prot"

# A real MUTADOCK sample run, generated for these tests. Present in the dev
# scratchpad; absent in CI — tests using it are skipped when missing.
REAL_SAMPLE_DIR = Path(
    "/tmp/claude-1000/-mnt-c-N-Drive-16--Mutadock/"
    "3468880e-17aa-4ebd-b429-baa6208e73b4/scratchpad/sample_run"
)


def _atom(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str = "C",
) -> str:
    """Format a single PDB ATOM record with correct fixed-width columns."""
    return (
        f"ATOM  {serial:>5} {name:<4}{'':1}{resname:>3} {chain:1}{resseq:>4}"
        f"{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"{'':>10}{element:>2}"
    )


def _write_pdb(path: Path, residues: list[tuple[str, str, int, str]]) -> None:
    """Write a tiny valid PDB with one CA atom per given residue.

    ``residues`` items are ``(chain, resname, resseq, element)`` — element is
    always carbon here (CA atoms).
    """
    lines: list[str] = ["REMARK synthetic test structure"]
    serial = 1
    x = 0.0
    for chain, resname, resseq, _elem in residues:
        lines.append(
            _atom(serial, "CA", resname, chain, resseq, x, x + 1.0, x + 2.0, "C")
        )
        serial += 1
        x += 3.0
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def sample_run_dir(tmp_path: Path) -> Path:
    """Create a synthetic run directory and return its path."""
    d = tmp_path

    # --- mutations.csv (no ddG column) --------------------------------------
    mutations = (
        "sr,pdb,chain,position,wtAA,prAA,wtProb,prProb\n"
        f"1,{d}/{STEM}_clean,A,10,ALA,GLY,0.02,0.01\n"
        f"2,{d}/{STEM}_clean,A,20,SER,THR,0.03,0.02\n"
        f"3,{d}/{STEM}_clean,B,30,LEU,ILE,0.02,0.01\n"
    )
    (d / f"{STEM}_mutations.csv").write_text(mutations)

    # --- single ddG (deliberately NOT ascending) ---------------------------
    single = (
        "sr,pdb,chain,position,wtAA,prAA,wtProb,prProb,ddG_value\n"
        f"3,{d}/{STEM}_clean,B,30,LEU,ILE,0.02,0.01,2.0\n"
        f"1,{d}/{STEM}_clean,A,10,ALA,GLY,0.02,0.01,-5.0\n"
        f"2,{d}/{STEM}_clean,A,20,SER,THR,0.03,0.02,-8.5\n"
    )
    (d / f"{STEM}_ddG_sorted.csv").write_text(single)

    # --- double ddG (deliberately NOT ascending) ---------------------------
    double = (
        "sr,pdb,combination,mut1_wtAA,chain1,mut1_position,mut1_prAA,"
        "mut2_wtAA,chain2,mut2_position,mut2_prAA,double_ddG_value\n"
        f"1,{d}/{STEM}_clean.pdb,1_S_20_T+A_10_G,S,A,20,T,A,A,10,G,-3.0\n"
        f"2,{d}/{STEM}_clean.pdb,2_S_20_T+L_30_I,S,A,20,T,L,B,30,I,-6.0\n"
    )
    (d / f"{STEM}_double_ddg_sorted.csv").write_text(double)

    # --- triple ddG (EMPTY: header only) -----------------------------------
    triple_header = (
        "sr,pdb,combination,mut1_wtAA,chain1,mut1_position,mut1_prAA,"
        "mut2_wtAA,chain2,mut2_position,mut2_prAA,"
        "mut3_wtAA,chain3,mut3_position,mut3_prAA,triple_ddG_value\n"
    )
    (d / f"{STEM}_clean_triple_ddg_sorted.csv").write_text(triple_header)

    # --- docking (deliberately NOT ascending) ------------------------------
    docking = (
        "sr,name,affinity\n"
        "1,2_SER-A20-THR_Ligand,-7.0\n"
        "2,1_ALA-A10-GLY_Ligand,-9.5\n"
        "3,prot_clean_Ligand,-5.0\n"
    )
    (d / "docking_results.csv").write_text(docking)

    # --- mutant PDBs (2 chains each so parse_ca_trace returns data) --------
    mut_dir = d / f"mutation_{STEM}"
    mut_dir.mkdir()
    two_chain = [
        ("A", "SER", 20, "C"),
        ("A", "ALA", 21, "C"),
        ("B", "LEU", 30, "C"),
    ]
    # Top single row is sr=2 (ddG=-8.5) → this file must exist by name.
    _write_pdb(mut_dir / "2_SER-A20-THR.pdb", two_chain)
    _write_pdb(mut_dir / "1_ALA-A10-GLY.pdb", two_chain)
    _write_pdb(mut_dir / "1_S_20_T+A_10_G.pdb", two_chain)

    # --- mutants.txt --------------------------------------------------------
    (d / f"{STEM}_mutants.txt").write_text(
        "\n".join(
            str(mut_dir / n)
            for n in ("2_SER-A20-THR.pdb", "1_ALA-A10-GLY.pdb", "1_S_20_T+A_10_G.pdb")
        )
        + "\n"
    )

    # --- cleaned wild-type structure ---------------------------------------
    _write_pdb(d / f"{STEM}_clean.pdb", two_chain)

    return d


@pytest.fixture
def real_sample_dir() -> Path:
    """Path to the real sample run; skips the test if it is not available."""
    if not REAL_SAMPLE_DIR.is_dir() or not os.listdir(REAL_SAMPLE_DIR):
        pytest.skip("real sample run directory not available")
    return REAL_SAMPLE_DIR
