"""Tests for ``mutadock.report.html_report``.

The report is exercised against the conftest synthetic run (which always runs)
and, when present, the real sample run. We verify a full, self-contained HTML
document is produced: summary tiles, base64 images, tables, and — for the
structure case — a genuinely *inlined* NGL library (not merely the literal
``NGL.Stage`` call, which would appear even if the blob were HTML-escaped and
broken).
"""

from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path

from mutadock.report import figures as F
from mutadock.report.data import discover_run
from mutadock.report.html_report import load_ngl_js, render_html


def _figures_for(run) -> dict:
    """Generate the real figure set for a run (any value may be None)."""
    struct = None
    if run.top_mutant_pdb is not None:
        struct = F.structure_backbone(run.top_mutant_pdb)
    return {
        "ddg_distribution": F.ddg_distribution(run.single_df),
        "top_mutations": F.top_mutations_bar(run.single_df),
        "affinity": F.affinity_plot(run.docking_df),
        "structure": struct,
    }


def _assert_parses(html: str) -> None:
    """Feed the HTML through stdlib html.parser; must not raise."""
    HTMLParser().feed(html)


# ---------------------------------------------------------------------------
# load_ngl_js
# ---------------------------------------------------------------------------
def test_load_ngl_js_is_large_and_defines_ngl() -> None:
    js = load_ngl_js()
    assert isinstance(js, str)
    assert len(js) > 1_000_000  # ~1.28 MB vendored build
    assert "loadFile" in js


# ---------------------------------------------------------------------------
# Synthetic run (always available via conftest)
# ---------------------------------------------------------------------------
def test_render_synthetic_full_document(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    html = render_html(run, _figures_for(run))

    assert html.lstrip().lower().startswith("<!doctype")
    assert "<html" in html
    assert run.base_name in html
    # summary tiles
    assert "Mutations scored" in html
    # at least one inlined image and a table
    assert "data:image/png;base64," in html
    assert "<table" in html
    _assert_parses(html)


def test_render_synthetic_has_working_ngl(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    assert run.top_mutant_pdb is not None
    html = render_html(run, _figures_for(run))

    assert 'id="md-viewport"' in html
    assert "NGL.Stage" in html
    assert "atob(" in html
    # Discriminating check: the minified NGL library contains raw ``&&`` and
    # ``<`` operators. If ``|safe`` were missing, autoescape would turn these
    # into ``&amp;&amp;`` / ``&lt;`` and the (still >1MB) file would be broken.
    assert "&amp;&amp;" not in html
    assert "loadFile" in html
    # Proof NGL is actually inlined, not linked.
    assert len(html) > 1_000_000
    # Inlining the whole library must not close the <script> early: only the
    # two authored closing tags may appear (the NGL blob is neutralised).
    assert html.count("</script") == 2


def test_render_default_title(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    html = render_html(run, _figures_for(run))
    assert f"MUTADOCK report — {run.base_name}" in html


def test_render_custom_title(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    html = render_html(run, _figures_for(run), title="My Study")
    assert "My Study" in html


def test_inlined_pdb_decodes_to_structure(sample_run_dir: Path) -> None:
    import base64
    import re

    run = discover_run(sample_run_dir)
    html = render_html(run, _figures_for(run))
    m = re.search(r'atob\("([^"]+)"\)', html)
    assert m is not None
    decoded = base64.b64decode(m.group(1)).decode("utf-8")
    assert decoded.startswith(("ATOM", "HEADER", "REMARK"))


# ---------------------------------------------------------------------------
# Graceful degradation: no structure, no docking
# ---------------------------------------------------------------------------
def test_render_degraded_no_structure_no_docking(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    run.top_mutant_pdb = None
    run.docking_df = None
    figures = {
        "ddg_distribution": F.ddg_distribution(run.single_df),
        "top_mutations": F.top_mutations_bar(run.single_df),
        "affinity": None,
        "structure": None,
    }
    html = render_html(run, figures)

    assert html.lstrip().lower().startswith("<!doctype")
    assert "NGL.Stage" not in html
    assert 'id="md-viewport"' not in html
    assert "Docking affinity" not in html
    # NGL is not inlined → small file.
    assert len(html) < 1_000_000
    _assert_parses(html)


def test_render_no_figures_at_all(sample_run_dir: Path) -> None:
    run = discover_run(sample_run_dir)
    run.top_mutant_pdb = None
    figures = {
        "ddg_distribution": None,
        "top_mutations": None,
        "affinity": None,
        "structure": None,
    }
    html = render_html(run, figures)
    # Tables (from single/double/docking data) still render even with no figs.
    assert "<table" in html
    assert "data:image/png;base64," not in html
    _assert_parses(html)


# ---------------------------------------------------------------------------
# Real sample run (skipped when the sample dir is absent)
# ---------------------------------------------------------------------------
def test_render_real_sample(real_sample_dir: Path) -> None:
    run = discover_run(real_sample_dir)
    html = render_html(run, _figures_for(run))

    assert html.lstrip().lower().startswith("<!doctype")
    assert run.base_name in html
    assert "data:image/png;base64," in html
    assert "<table" in html
    if run.top_mutant_pdb is not None:
        assert "NGL.Stage" in html
        assert "atob(" in html
        assert "&amp;&amp;" not in html
        assert len(html) > 1_000_000
    _assert_parses(html)
