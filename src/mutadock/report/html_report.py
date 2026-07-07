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
"""Self-contained, single-file HTML report generator for MUTADOCK.

This module assembles a completed :class:`~mutadock.report.data.RunData` and a
dict of pre-rendered base64 figures into one standalone HTML document a PI can
open offline or receive by email. Everything is inlined: styles, images (as
``data:`` URIs), the mutant structure (base64-encoded PDB), and — only when a
structure is shown — the entire vendored NGL library, so the interactive 3D
viewer works with no network access.

The public entry point is :func:`render_html`; the orchestrator writes its
returned string to disk. :func:`load_ngl_js` exposes the vendored NGL asset.
"""

from __future__ import annotations

import base64
from importlib.resources import files

from jinja2 import Environment

from .data import RunData, mutation_label, summary_stats
from .structure import mutation_sites, read_pdb_text

__all__ = ["render_html", "load_ngl_js"]

_TEMPLATE_NAME = "templates/report.html.j2"


def load_ngl_js() -> str:
    """Return the text of the vendored NGL UMD build (defines global ``NGL``).

    The text is safe to inline directly inside an inline ``<script>`` element:
    any literal ``</script`` sequence is neutralised to ``<\\/script`` (a JS
    no-op) so the HTML parser cannot terminate the element early.
    """
    text = (
        files("mutadock.report")
        .joinpath("assets/ngl.min.js")
        .read_text(encoding="utf-8")
    )
    return text.replace("</script", "<\\/script")


def _load_template_source() -> str:
    """Read the Jinja2 template source via importlib.resources (wheel-safe)."""
    return files("mutadock.report").joinpath(_TEMPLATE_NAME).read_text(encoding="utf-8")


def _fmt_ddg(value: object) -> str:
    """Format a ΔΔG / affinity value to 2 dp, tolerating None / bad input."""
    try:
        return f"{float(value):.2f}"  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return "—"


def _fmt_sd(value: object) -> str:
    """Format a standard deviation as ``± X.XX``; empty for missing / zero SD."""
    try:
        v = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return ""
    return f"± {v:.2f}" if v > 0 else ""


def _single_table_rows(run: RunData, top_n: int) -> list[dict]:
    """Top-*n* single-mutation rows as template-ready dicts."""
    df = run.single_df
    if df is None:
        return []
    rows: list[dict] = []
    for i, record in enumerate(df.head(top_n).to_dict("records"), start=1):
        rows.append(
            {
                "rank": i,
                "label": mutation_label(record),
                "chain": record.get("chain", ""),
                "position": record.get("position", ""),
                "wtAA": record.get("wtAA", ""),
                "prAA": record.get("prAA", ""),
                "ddg": _fmt_ddg(record.get("ddG_value")),
                "ddg_sd": _fmt_sd(record.get("ddG_sd")),
                "ddg_kcal": _fmt_ddg(record.get("ddG_kcal")),
            }
        )
    return rows


def _double_table_rows(run: RunData, top_n: int) -> list[dict]:
    """Top-*n* double-mutation rows as template-ready dicts."""
    df = run.double_df
    if df is None:
        return []
    rows: list[dict] = []
    for i, record in enumerate(df.head(top_n).to_dict("records"), start=1):
        rows.append(
            {
                "rank": i,
                "combination": record.get("combination", ""),
                "ddg": _fmt_ddg(record.get("double_ddG_value")),
                "ddg_sd": _fmt_sd(record.get("double_ddG_sd")),
                "ddg_kcal": _fmt_ddg(record.get("double_ddG_kcal")),
            }
        )
    return rows


def _docking_table_rows(run: RunData, top_n: int) -> list[dict]:
    """Top-*n* docking rows as template-ready dicts."""
    df = run.docking_df
    if df is None:
        return []
    rows: list[dict] = []
    for record in df.head(top_n).to_dict("records"):
        rows.append(
            {
                "name": record.get("name", ""),
                "affinity": _fmt_ddg(record.get("affinity")),
            }
        )
    return rows


def _stat_tiles(run: RunData) -> list[dict]:
    """Build the summary stat tiles, omitting any whose value is None."""
    stats = summary_stats(run)
    tiles: list[dict] = []

    tiles.append({"label": "Mutations scored", "value": str(stats["n_scored"])})
    tiles.append(
        {"label": "Stabilizing (ΔΔG < 0)", "value": str(stats["n_stabilizing"])}
    )
    if stats["best_ddg"] is not None:
        tiles.append(
            {
                "label": "Best ΔΔG (REU)",
                "value": _fmt_ddg(stats["best_ddg"]),
                "sub": stats["best_ddg_label"],
            }
        )
    if stats.get("best_ddg_kcal") is not None:
        tiles.append(
            {
                "label": "Best ΔΔG (kcal/mol)",
                "value": _fmt_ddg(stats["best_ddg_kcal"]),
                "sub": stats["best_ddg_label"],
            }
        )
    # Only show docking-related tiles when there is an actual affinity.
    if stats["best_affinity"] is not None:
        tiles.append({"label": "Docked", "value": str(stats["n_docked"])})
        tiles.append(
            {
                "label": "Best affinity (kcal/mol)",
                "value": _fmt_ddg(stats["best_affinity"]),
                "sub": stats["best_affinity_name"],
            }
        )
    return tiles


def render_html(
    run: RunData,
    figures: dict,
    *,
    top_n: int = 15,
    title: str | None = None,
) -> str:
    """Render a complete, self-contained HTML report for a MUTADOCK run.

    Args:
        run: The loaded run (see :func:`mutadock.report.data.discover_run`).
        figures: Pre-rendered base64 PNG strings (no ``data:`` prefix) keyed by
            ``ddg_distribution``, ``top_mutations``, ``affinity`` and
            ``structure``. Any ``None`` value omits that block.
        top_n: Maximum rows shown in each table.
        title: Report title; defaults to ``f"MUTADOCK report — {base_name}"``.

    Returns:
        The full HTML document as a single string.
    """
    figures = figures or {}
    if title is None:
        title = f"MUTADOCK report — {run.base_name}"

    metadata = run.metadata or {}
    versions = metadata.get("versions", {}) or {}

    # ΔΔG units / provenance note. ΔΔG is in Rosetta Energy Units (REU), not
    # kcal/mol; state the run's actual scaling factor and replicate count.
    stats = summary_stats(run)
    ddg_note_parts = [
        "ΔΔG is in Rosetta Energy Units (REU), not kcal/mol; "
        "lower (more negative) = more stabilizing."
    ]
    if stats.get("reu_to_kcal") is not None:
        ddg_note_parts.append(
            f"kcal/mol ≈ REU × {stats['reu_to_kcal']:.3g} "
            "(scaling factor; approximate, see docs)."
        )
    n_rep = stats.get("n_replicates")
    if n_rep is not None and n_rep > 1:
        ddg_note_parts.append(f"Values are mean ± SD over {n_rep} replicates.")
    elif n_rep == 1:
        ddg_note_parts.append(
            "Single replicate (no SD); pass --replicates N for mean ± SD."
        )
    protocol = stats.get("protocol")
    if protocol:
        ddg_note_parts.append(f"Protocol: {protocol}.")
    # The fast (no-minimization) protocol is screening-only: a failed wild-type
    # repack can invert a site's ranking, so the top hits may be artifacts.
    ddg_warning = (
        "SCREENING ONLY — the 'fast' protocol skips minimization, so absolute "
        "ΔΔG values (and the top-mutations ranking) can be unreliable. "
        "Re-run with --protocol cartesian (or the default 'min') for "
        "trustworthy values."
        if protocol == "fast"
        else ""
    )
    ddg_note = " ".join(ddg_note_parts)

    # 3D viewer: only build (and inline NGL) when a structure is available.
    has_structure = run.top_mutant_pdb is not None
    pdb_b64 = ""
    ngl_sele = ""
    ngl_js = ""
    mutation_caption = ""
    if has_structure:
        assert run.top_mutant_pdb is not None
        pdb_text = read_pdb_text(run.top_mutant_pdb)
        pdb_b64 = base64.b64encode(pdb_text.encode("utf-8")).decode("ascii")
        sites = mutation_sites(run.top_single_row) if run.top_single_row else []
        ngl_sele = " or ".join(f"{pos} and :{chain}" for chain, pos in sites)
        if run.top_single_row is not None:
            mutation_caption = mutation_label(run.top_single_row)
        ngl_js = load_ngl_js()

    context = {
        "title": title,
        "base_name": run.base_name,
        "generated_at": metadata.get("generated_at", ""),
        "input_name": (
            metadata.get("input", "").rsplit("/", 1)[-1]
            if metadata.get("input")
            else None
        ),
        "tiles": _stat_tiles(run),
        "ddg_note": ddg_note,
        "ddg_warning": ddg_warning,
        "fig_ddg": figures.get("ddg_distribution"),
        "fig_top": figures.get("top_mutations"),
        "fig_affinity": figures.get("affinity"),
        "fig_structure": figures.get("structure"),
        "single_rows": _single_table_rows(run, top_n),
        "double_rows": _double_table_rows(run, top_n),
        "docking_rows": _docking_table_rows(run, top_n),
        "top_n": top_n,
        "has_structure": has_structure,
        "pdb_b64": pdb_b64,
        "ngl_sele": ngl_sele,
        "ngl_js": ngl_js,
        "mutation_caption": mutation_caption,
        "meta": {
            "generated_at": metadata.get("generated_at", ""),
            "input": metadata.get("input"),
            "output_dir": metadata.get("output_dir"),
            "files": metadata.get("files", []) or [],
            "pyrosetta": versions.get("pyrosetta", "unknown"),
            "vina": versions.get("vina", "unknown"),
        },
    }

    env = Environment(autoescape=True, trim_blocks=True, lstrip_blocks=True)
    template = env.from_string(_load_template_source())
    return template.render(**context)
