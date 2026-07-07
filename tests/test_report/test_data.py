"""Tests for mutadock.report.data."""

from __future__ import annotations

from pathlib import Path

from mutadock.report.data import discover_run, mutation_label, summary_stats


def test_discover_synthetic(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)

    # base_name inferred (NOT hard-coded to 4QJR).
    assert run.base_name == "prot"

    # Frames loaded.
    assert run.mutations_df is not None
    assert run.single_df is not None
    assert run.double_df is not None
    assert run.docking_df is not None

    # Empty triple → None.
    assert run.triple_df is None

    # Single re-sorted ascending → first row has the min ddG.
    assert run.single_df.iloc[0]["ddG_value"] == -8.5
    assert run.single_df["ddG_value"].is_monotonic_increasing

    # Double re-sorted ascending.
    assert run.double_df.iloc[0]["double_ddG_value"] == -6.0
    assert run.double_df["double_ddG_value"].is_monotonic_increasing

    # Docking re-sorted ascending by affinity.
    assert run.docking_df.iloc[0]["affinity"] == -9.5
    assert run.docking_df["affinity"].is_monotonic_increasing

    # input_pdb prefers the cleaned structure.
    assert run.input_pdb is not None
    assert run.input_pdb.name == "prot_clean.pdb"

    # Top single row + mutant resolution.
    assert run.top_single_row is not None
    assert run.top_single_row["ddG_value"] == -8.5
    assert isinstance(run.top_single_row["sr"], int)
    assert run.top_mutant_pdb is not None
    assert run.top_mutant_pdb.exists()
    assert run.top_mutant_pdb.name == "2_SER-A20-THR.pdb"

    # mutant_pdbs all exist and are de-duped.
    assert len(run.mutant_pdbs) == 3
    assert all(p.exists() for p in run.mutant_pdbs)


def test_summary_stats_synthetic(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    stats = summary_stats(run)

    assert stats["n_scored"] == 3
    assert stats["n_stabilizing"] == 2
    assert stats["best_ddg"] == -8.5
    assert stats["best_ddg_label"] == "SER-A20-THR"
    assert stats["n_docked"] == 3
    assert stats["best_affinity"] == -9.5
    assert stats["best_affinity_name"] == "1_ALA-A10-GLY_Ligand"


def test_metadata_populated(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir, timestamp="2026-07-02T12:00:00")
    md = run.metadata

    assert md["generated_at"] == "2026-07-02T12:00:00"
    assert md["input"] == str(run.input_pdb)
    assert md["output_dir"] == str(sample_run_dir)
    assert isinstance(md["files"], list) and md["files"]
    assert md["files"] == sorted(md["files"])
    assert "docking_results.csv" in md["files"]
    assert set(md["versions"]) == {"pyrosetta", "vina"}


def test_generated_at_default(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    assert isinstance(run.metadata["generated_at"], str)
    assert run.metadata["generated_at"]


def test_mutation_label() -> None:
    row = {"wtAA": "ASN", "chain": "A", "position": 386, "prAA": "HIS"}
    assert mutation_label(row) == "ASN-A386-HIS"
    # position given as float-like string still renders cleanly.
    assert mutation_label(
        {"wtAA": "S", "chain": "A", "position": "20", "prAA": "T"}
    ) == ("S-A20-T")


def test_degradation_docking_only(tmp_path: Path) -> None:
    (tmp_path / "docking_results.csv").write_text(
        "sr,name,affinity\n1,x_Ligand,-5.0\n2,y_Ligand,-8.0\n"
    )
    run = discover_run(tmp_path)

    assert run.mutations_df is None
    assert run.single_df is None
    assert run.double_df is None
    assert run.triple_df is None
    assert run.docking_df is not None
    assert run.docking_df.iloc[0]["affinity"] == -8.0

    # Should not crash even with no ddG data.
    stats = summary_stats(run)
    assert stats["n_scored"] == 0
    assert stats["n_stabilizing"] == 0
    assert stats["best_ddg"] is None
    assert stats["n_docked"] == 2


def test_real_sample(real_sample_dir: Path) -> None:
    run = discover_run(real_sample_dir)

    assert run.base_name == "4QJR"
    assert run.single_df is not None
    assert 900 < len(run.single_df) < 1050  # ~979 rows

    # Ascending → first is most negative.
    assert run.single_df["ddG_value"].is_monotonic_increasing
    assert run.single_df.iloc[0]["ddG_value"] == run.single_df["ddG_value"].min()

    # Triple file is header-only in the sample → None.
    assert run.triple_df is None
    assert run.double_df is not None

    # Top mutant resolves to a real file.
    assert run.top_mutant_pdb is not None
    assert run.top_mutant_pdb.exists()
