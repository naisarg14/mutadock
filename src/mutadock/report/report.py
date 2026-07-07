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

"""Orchestrator + ``md_report`` CLI for MUTADOCK run reports.

Discovers a run's output files, builds the shared set of figures once, and
writes them into a ``reports/`` folder inside the run directory: a
self-contained ``report.html``, a ``report.pptx`` deck, and a ``figures/``
sub-folder holding each chart as a standalone PNG.
"""

from __future__ import annotations

import argparse
import base64
import logging
import sys
from pathlib import Path

from mutadock.report.data import RunData, discover_run
from mutadock.report.structure import mutation_sites

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Reports go into a dedicated sub-folder so they never collide with (or get
# mistaken for) the run's raw CSVs / PDBs.
REPORTS_DIRNAME = "reports"
FIGURES_DIRNAME = "figures"
REPORT_HTML_NAME = "report.html"
REPORT_PPTX_NAME = "report.pptx"

# Standalone-PNG file names for each figure key produced by build_figures().
_FIGURE_FILENAMES = {
    "ddg_distribution": "ddg_distribution.png",
    "top_mutations": "top_mutations.png",
    "affinity": "docking_affinity.png",
    "structure": "top_mutant_structure.png",
}


def _write_figures(figures: dict, figures_dir: Path) -> dict[str, Path]:
    """Decode each non-empty base64 figure and write it as a PNG file.

    Args:
        figures: Mapping of figure key → base64-PNG string (values may be
            ``None`` when there is no data for that figure).
        figures_dir: Destination directory (created only if there is at least
            one figure to write).

    Returns:
        Mapping of figure key → written PNG path.
    """
    to_write = {k: v for k, v in figures.items() if v}
    written: dict[str, Path] = {}
    if not to_write:
        return written
    figures_dir.mkdir(parents=True, exist_ok=True)
    for key, b64 in to_write.items():
        path = figures_dir / _FIGURE_FILENAMES.get(key, f"{key}.png")
        path.write_bytes(base64.b64decode(b64))
        written[key] = path
    return written


def build_figures(run: RunData, *, top_n: int = 15, dark: bool = False) -> dict:
    """Build the shared base64-PNG figure set used by both renderers.

    Figures are generated once here and passed to both the HTML and PPTX
    renderers so the (relatively expensive) matplotlib work happens a single
    time per report.

    Args:
        run: The loaded run.
        top_n: How many top mutations the bar chart should show.
        dark: Render dark-mode figure variants.

    Returns:
        Dict keyed by ``ddg_distribution``, ``top_mutations``, ``affinity`` and
        ``structure``; any value may be ``None`` when there is no data.
    """
    # Imported lazily so that merely importing this module (e.g. for the CLI
    # wiring in np_mutation/np_docking) does not pull in matplotlib.
    from mutadock.report import figures as F

    sites = mutation_sites(run.top_single_row) if run.top_single_row else None
    structure_png = None
    if run.top_mutant_pdb is not None:
        structure_png = F.structure_backbone(run.top_mutant_pdb, sites, dark=dark)

    return {
        "ddg_distribution": F.ddg_distribution(run.single_df, dark=dark),
        "top_mutations": F.top_mutations_bar(run.single_df, n=top_n, dark=dark),
        "affinity": F.affinity_plot(run.docking_df, dark=dark),
        "structure": structure_png,
    }


def generate_report(
    output_dir: str | Path,
    base_name: str | None = None,
    formats: tuple[str, ...] = ("html", "pptx"),
    *,
    report_dir: str | Path | None = None,
    top_n: int = 15,
    title: str | None = None,
    timestamp: str | None = None,
    quiet: bool = False,
) -> dict[str, Path]:
    """Discover a run under *output_dir* and write the requested report files.

    Args:
        output_dir: Directory holding a MUTADOCK run's CSVs / mutant PDBs.
        base_name: Optional run stem; inferred from the files when omitted.
        formats: Which reports to emit — any of ``"html"`` and ``"pptx"``.
        report_dir: Where to write the report files (defaults to a
            ``reports/`` sub-folder inside *output_dir*). Figures are written
            to a ``figures/`` sub-folder within this directory.
        top_n: Rows shown in tables / bars.
        title: Report title; a sensible default is used when omitted.
        timestamp: Overrides the generation timestamp (mainly for tests).
        quiet: Suppress informational logging.

    Returns:
        Mapping of format name → written file path (only for formats emitted).
    """
    run = discover_run(output_dir, base_name=base_name, timestamp=timestamp)

    has_any = any(
        df is not None
        for df in (
            run.mutations_df,
            run.single_df,
            run.double_df,
            run.triple_df,
            run.docking_df,
        )
    )
    if not has_any:
        logger.warning(
            "No MUTADOCK result files found in %s — nothing to report.", output_dir
        )
        return {}

    figures = build_figures(run, top_n=top_n)

    # Default: a "reports/" sub-folder inside the run directory. An explicit
    # report_dir overrides this and is used as-is.
    if report_dir is not None:
        rdir = Path(report_dir)
    else:
        rdir = Path(output_dir) / REPORTS_DIRNAME
    rdir.mkdir(parents=True, exist_ok=True)

    written: dict[str, Path] = {}

    # Save each chart as a standalone PNG under reports/figures/ so users can
    # drop them straight into their own slides/papers; the HTML still embeds
    # its own copies so it stays self-contained.
    figure_paths = _write_figures(figures, rdir / FIGURES_DIRNAME)
    if figure_paths:
        written["figures"] = rdir / FIGURES_DIRNAME
        if not quiet:
            logger.info(
                "Figures written: %s (%d files)",
                rdir / FIGURES_DIRNAME,
                len(figure_paths),
            )

    if "html" in formats:
        from mutadock.report.html_report import render_html

        html = render_html(run, figures, top_n=top_n, title=title)
        out_path = rdir / REPORT_HTML_NAME
        # Report files are regenerated every run; overwrite rather than backup
        # so the run directory does not accumulate stale copies.
        out_path.write_text(html, encoding="utf-8")
        written["html"] = out_path
        if not quiet:
            logger.info("HTML report written: %s", out_path)

    if "pptx" in formats:
        from mutadock.report.ppt_report import build_pptx

        prs = build_pptx(run, figures, top_n=top_n, title=title)
        out_path = rdir / REPORT_PPTX_NAME
        prs.save(str(out_path))
        written["pptx"] = out_path
        if not quiet:
            logger.info("PPTX report written: %s", out_path)

    return written


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="md_report",
        description="Generate a self-contained HTML + PPTX report from a "
        "MUTADOCK run directory.",
        epilog="Part of mutadock library. Written by Naisarg Patel "
        "(https://github.com/naisarg14)",
    )
    parser.add_argument(
        "-d",
        "--dir",
        required=True,
        help="Run directory containing MUTADOCK output CSVs / mutant PDBs.",
        metavar="DIR",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Where to write report files (default: a 'reports/' folder "
        "inside the run directory). Figures are saved under '<here>/figures/'.",
        metavar="DIR",
    )
    fmt = parser.add_mutually_exclusive_group()
    fmt.add_argument(
        "--html-only",
        dest="html_only",
        action="store_true",
        help="Only write report.html.",
    )
    fmt.add_argument(
        "--pptx-only",
        dest="pptx_only",
        action="store_true",
        help="Only write report.pptx.",
    )
    parser.add_argument(
        "--top",
        dest="top_n",
        type=int,
        default=15,
        help="Number of top mutations to tabulate/plot (default: 15).",
        metavar="N",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom report title.",
        metavar="TITLE",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress informational logging.",
    )
    return parser.parse_args(argv)


def md_report(argv: list[str] | None = None) -> None:
    """CLI entry point for ``md_report``."""
    args = _parse_args(argv)

    run_dir = Path(args.dir)
    if not run_dir.is_dir():
        sys.exit(f"Run directory not found: {run_dir}")

    if args.html_only:
        formats: tuple[str, ...] = ("html",)
    elif args.pptx_only:
        formats = ("pptx",)
    else:
        formats = ("html", "pptx")

    try:
        written = generate_report(
            run_dir,
            formats=formats,
            report_dir=args.output_dir,
            top_n=args.top_n,
            title=args.title,
            quiet=args.quiet,
        )
    except Exception as e:  # pragma: no cover - defensive CLI guard
        sys.exit(f"Report generation failed: {e}")

    if not written:
        sys.exit(
            f"No MUTADOCK result files found in {run_dir}. "
            "Run md_mutate and/or md_dock first."
        )

    for fmt, path in written.items():
        logger.info("%s -> %s", fmt.upper(), path)


if __name__ == "__main__":
    md_report()
