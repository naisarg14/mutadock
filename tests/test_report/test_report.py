"""Tests for ``mutadock.report.report.generate_report`` output layout.

Verifies that reports land in a ``reports/`` sub-folder with a ``figures/``
directory of standalone PNGs alongside the HTML and PPTX files.
"""

from pathlib import Path

from mutadock.report.report import (
    FIGURES_DIRNAME,
    REPORT_HTML_NAME,
    REPORT_PPTX_NAME,
    REPORTS_DIRNAME,
    generate_report,
)


def test_reports_written_into_reports_subfolder(sample_run_dir: Path) -> None:
    written = generate_report(sample_run_dir, quiet=True)

    reports = sample_run_dir / REPORTS_DIRNAME
    assert reports.is_dir()
    assert (reports / REPORT_HTML_NAME).is_file()
    assert (reports / REPORT_PPTX_NAME).is_file()
    # Nothing should be written directly into the run directory.
    assert not (sample_run_dir / REPORT_HTML_NAME).exists()
    assert not (sample_run_dir / REPORT_PPTX_NAME).exists()

    assert written["html"] == reports / REPORT_HTML_NAME
    assert written["pptx"] == reports / REPORT_PPTX_NAME


def test_figures_saved_as_standalone_pngs(sample_run_dir: Path) -> None:
    written = generate_report(sample_run_dir, quiet=True)

    figures_dir = sample_run_dir / REPORTS_DIRNAME / FIGURES_DIRNAME
    assert figures_dir.is_dir()
    assert written["figures"] == figures_dir

    pngs = sorted(p.name for p in figures_dir.glob("*.png"))
    # The synthetic run has single ddG, docking, and a top mutant structure.
    assert "ddg_distribution.png" in pngs
    assert "top_mutations.png" in pngs
    assert "docking_affinity.png" in pngs
    assert "top_mutant_structure.png" in pngs
    # Each PNG is a real, non-empty image (PNG magic bytes).
    for p in figures_dir.glob("*.png"):
        data = p.read_bytes()
        assert data.startswith(b"\x89PNG\r\n\x1a\n")
        assert len(data) > 1000


def test_html_still_self_contained(sample_run_dir: Path) -> None:
    """Saving standalone PNGs must not stop the HTML embedding its own copies."""
    generate_report(sample_run_dir, formats=("html",), quiet=True)
    html = (sample_run_dir / REPORTS_DIRNAME / REPORT_HTML_NAME).read_text(
        encoding="utf-8"
    )
    # Charts are embedded as base64 data URIs, not <img src="figures/...">.
    assert "data:image/png;base64," in html
    assert 'src="figures/' not in html


def test_explicit_report_dir_override(sample_run_dir: Path, tmp_path: Path) -> None:
    dest = tmp_path / "custom_out"
    generate_report(sample_run_dir, report_dir=dest, quiet=True)
    assert (dest / REPORT_HTML_NAME).is_file()
    assert (dest / FIGURES_DIRNAME).is_dir()
    # Default reports/ folder is not created when an explicit dir is given.
    assert not (sample_run_dir / REPORTS_DIRNAME).exists()
