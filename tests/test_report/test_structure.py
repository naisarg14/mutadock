"""Tests for mutadock.report.structure."""

from __future__ import annotations

from pathlib import Path

import pytest

from mutadock.report.structure import mutation_sites, parse_ca_trace, read_pdb_text


def test_parse_ca_trace_synthetic(sample_run_dir: Path) -> None:
    pdb = sample_run_dir / "prot_clean.pdb"
    trace = parse_ca_trace(pdb)

    # Two chains present.
    assert set(trace) == {"A", "B"}

    # Chain A has two CA residues, chain B has one.
    assert len(trace["A"]) == 2
    assert len(trace["B"]) == 1

    resseq, x, y, z, resname = trace["A"][0]
    assert resseq == 20
    assert resname == "SER"
    assert isinstance(x, float) and isinstance(y, float) and isinstance(z, float)


def test_parse_ca_trace_missing() -> None:
    with pytest.raises(FileNotFoundError):
        parse_ca_trace("/no/such/file.pdb")


def test_read_pdb_text(sample_run_dir: Path) -> None:
    pdb = sample_run_dir / "prot_clean.pdb"
    text = read_pdb_text(pdb)
    assert "ATOM" in text
    assert text == pdb.read_text()


def test_read_pdb_text_missing() -> None:
    with pytest.raises(FileNotFoundError):
        read_pdb_text("/no/such/file.pdb")


def test_mutation_sites_single() -> None:
    row = {"chain": "A", "position": 386, "wtAA": "ASN", "prAA": "HIS"}
    assert mutation_sites(row) == [("A", 386)]


def test_mutation_sites_double() -> None:
    row = {
        "chain1": "A",
        "mut1_position": 386,
        "chain2": "B",
        "mut2_position": 30,
    }
    assert mutation_sites(row) == [("A", 386), ("B", 30)]


def test_mutation_sites_triple() -> None:
    row = {
        "chain1": "A",
        "mut1_position": 10,
        "chain2": "A",
        "mut2_position": 20,
        "chain3": "B",
        "mut3_position": 30,
    }
    assert mutation_sites(row) == [("A", 10), ("A", 20), ("B", 30)]
