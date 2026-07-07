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
"""Publication-quality, colorblind-safe charts for MUTADOCK reports.

Every public figure function returns a **base64-encoded PNG string** (no
``data:`` URI prefix) or ``None`` when there is nothing to plot, so both the
HTML and PPTX renderers can share the exact same images and skip empty
sections. Charts follow the project's validated dataviz standard: a small
colorblind-safe palette, recessive chrome, and a light surface tuned to be
excellent by default (with a straightforward dark variant).

The non-interactive Agg backend is selected before ``pyplot`` is imported so
this module never touches a GUI.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .data import mutation_label  # noqa: E402
from .structure import parse_ca_trace  # noqa: E402

# ---------------------------------------------------------------------------
# Palette / theme tokens (from the project's dataviz standard)
# ---------------------------------------------------------------------------
# Categorical chain colors, in this fixed order.
_CHAIN_COLORS_LIGHT = [
    "#2a78d6",  # blue
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
    "#e87ba4",  # magenta
    "#eb6834",  # orange
]
_CHAIN_COLORS_DARK = [
    "#3987e5",  # blue
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#2bb32b",  # green
    "#8a7be0",  # violet
    "#e66767",  # red
    "#e87ba4",  # magenta
    "#eb6834",  # orange
]

_LIGHT: dict[str, Any] = {
    "blue": "#2a78d6",
    "red": "#e34948",
    "blue_seq": "#256abf",
    "highlight": "#eb6834",
    "ink_primary": "#0b0b0b",
    "ink_secondary": "#52514e",
    "ink_muted": "#898781",
    "grid": "#e1e0d9",
    "baseline": "#c3c2b7",
    "surface": "#fcfcfb",
    "chains": _CHAIN_COLORS_LIGHT,
}

_DARK: dict[str, Any] = {
    "blue": "#3987e5",
    "red": "#e66767",
    "blue_seq": "#3987e5",
    "highlight": "#eb6834",
    "ink_primary": "#ffffff",
    "ink_secondary": "#c3c2b7",
    "ink_muted": "#c3c2b7",
    "grid": "#33332f",
    "baseline": "#52514e",
    "surface": "#1a1a19",
    "chains": _CHAIN_COLORS_DARK,
}


def _theme(dark: bool) -> dict[str, Any]:
    """Return the colour token mapping for the requested surface variant."""
    return _DARK if dark else _LIGHT


def _rc(theme: dict[str, Any]) -> dict[str, Any]:
    """Build the ``rcParams`` overrides that keep every figure consistent."""
    return {
        "figure.facecolor": theme["surface"],
        "axes.facecolor": theme["surface"],
        "savefig.facecolor": theme["surface"],
        "font.size": 10,
        "text.color": theme["ink_primary"],
        "axes.edgecolor": theme["baseline"],
        "axes.labelcolor": theme["ink_secondary"],
        "axes.titlecolor": theme["ink_primary"],
        "xtick.color": theme["ink_muted"],
        "ytick.color": theme["ink_muted"],
        "xtick.labelcolor": theme["ink_muted"],
        "ytick.labelcolor": theme["ink_muted"],
        "grid.color": theme["grid"],
        "grid.linewidth": 0.6,
        "axes.grid": False,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.dpi": 150,
    }


def _style_axes(ax: plt.Axes, theme: dict[str, Any]) -> None:
    """Recede the chrome: mute the remaining visible spines."""
    for spine in ax.spines.values():
        spine.set_color(theme["baseline"])
        spine.set_linewidth(0.8)


# ---------------------------------------------------------------------------
# Base64 helper
# ---------------------------------------------------------------------------
def fig_to_base64(fig: plt.Figure) -> str:
    """Render *fig* to a PNG and return its base64 text (figure is closed)."""
    buf = io.BytesIO()
    fig.savefig(
        buf,
        format="png",
        dpi=150,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


# ---------------------------------------------------------------------------
# Internal utilities
# ---------------------------------------------------------------------------
def _has_rows(df: pd.DataFrame | None, column: str | None = None) -> bool:
    """True if *df* has at least one row (and *column* if given)."""
    if df is None or len(df) == 0:
        return False
    if column is not None and column not in df.columns:
        return False
    return True


def _shorten_ligand(name: Any) -> str:
    """Trim a trailing ``_Ligand`` suffix from a docking name for readability."""
    text = str(name)
    if text.endswith("_Ligand"):
        text = text[: -len("_Ligand")]
    return text


def _robust_hist_range(data: np.ndarray) -> tuple[float, float, int]:
    """Return ``(lo, hi, n_beyond)`` for a histogram view of *data*.

    ΔΔG is heavy-tailed (severe clashes can reach thousands of REU), which would
    squash the decision-relevant region around zero. When — and only when — a
    genuine upper tail exists (beyond a Tukey *far* fence, q3 + 3·IQR), the view
    is clipped to that fence and ``n_beyond`` counts the hidden points. On
    well-behaved (bounded) data nothing is clipped and ``n_beyond`` is ``0``.
    """
    lo = float(data.min())
    hi_max = float(data.max())
    q1, q3 = (float(v) for v in np.percentile(data, [25, 75]))
    fence = q3 + 3.0 * (q3 - q1)
    if hi_max > fence > lo:
        return lo, fence, int((data > fence).sum())
    return lo, hi_max, 0


# ---------------------------------------------------------------------------
# 1. ddG distribution
# ---------------------------------------------------------------------------
def ddg_distribution(
    single_df: pd.DataFrame | None, *, dark: bool = False
) -> str | None:
    """Diverging histogram of single-mutation ΔΔG values (blue = stabilizing)."""
    if not _has_rows(single_df, "ddG_value"):
        return None
    assert single_df is not None
    values = pd.to_numeric(single_df["ddG_value"], errors="coerce").dropna()
    if values.empty:
        return None

    data = values.to_numpy(dtype=float)
    n = data.size
    n_stab = int((data < 0).sum())
    mean = float(data.mean())
    lo, hi, n_beyond = _robust_hist_range(data)
    plot_data = data[data <= hi] if n_beyond else data

    theme = _theme(dark)
    with plt.rc_context(_rc(theme)):
        fig, ax = plt.subplots(figsize=(7.2, 4.5))

        counts, edges = np.histogram(plot_data, bins=40, range=(lo, hi))
        centers = (edges[:-1] + edges[1:]) / 2.0
        widths = np.diff(edges)
        colors = [theme["blue"] if c < 0 else theme["red"] for c in centers]
        ax.bar(
            centers,
            counts,
            width=widths,
            color=colors,
            edgecolor=theme["surface"],
            linewidth=0.3,
        )

        ax.axvline(0.0, color=theme["baseline"], linewidth=1.5, zorder=1)
        top = float(counts.max()) if counts.size else 1.0
        if lo <= mean <= hi:
            ax.axvline(
                mean,
                color=theme["ink_muted"],
                linewidth=1.2,
                linestyle="--",
                zorder=2,
            )
            ax.annotate(
                f"mean = {mean:.1f}",
                xy=(mean, top),
                xytext=(4, -2),
                textcoords="offset points",
                ha="left",
                va="top",
                color=theme["ink_secondary"],
                fontsize=9,
            )
        if n_beyond:
            ax.annotate(
                f"+{n_beyond} with ΔΔG > {hi:.0f} (not shown)",
                xy=(0.99, 0.84),
                xycoords="axes fraction",
                ha="right",
                va="top",
                color=theme["ink_secondary"],
                fontsize=8.5,
            )
        # Guard the degenerate single-value case (e.g. md_quick's one mutation),
        # where hi == lo would make the xlim transformation singular.
        pad = (hi - lo) * 0.02 or 1.0
        ax.set_xlim(lo - pad, hi + pad)

        ax.grid(axis="y", color=theme["grid"], linewidth=0.6)
        ax.set_title(
            "Single-mutation ΔΔG distribution", fontsize=13, pad=18, loc="left"
        )
        ax.set_xlabel("ΔΔG (REU)")
        ax.set_ylabel("count")
        ax.annotate(
            f"n = {n}, {n_stab} stabilizing (ΔΔG < 0)",
            xy=(0.0, 1.02),
            xycoords="axes fraction",
            ha="left",
            va="bottom",
            color=theme["ink_secondary"],
            fontsize=9.5,
        )
        _style_axes(ax, theme)
        fig.tight_layout()
        return fig_to_base64(fig)


# ---------------------------------------------------------------------------
# 2. Top stabilizing mutations
# ---------------------------------------------------------------------------
def top_mutations_bar(
    single_df: pd.DataFrame | None, n: int = 15, *, dark: bool = False
) -> str | None:
    """Horizontal bars of the top-*n* most-stabilizing single mutations."""
    if not _has_rows(single_df, "ddG_value"):
        return None
    assert single_df is not None
    df = single_df.head(max(int(n), 1)).copy()
    df = df[pd.to_numeric(df["ddG_value"], errors="coerce").notna()]
    if df.empty:
        return None

    labels = [mutation_label(row) for row in df.to_dict("records")]
    values = pd.to_numeric(df["ddG_value"], errors="coerce").to_numpy(float)
    k = len(values)

    theme = _theme(dark)
    with plt.rc_context(_rc(theme)):
        height = max(3.0, 0.34 * k + 1.4)
        fig, ax = plt.subplots(figsize=(7.2, min(height, 9.0)))

        y = np.arange(k)
        colors = [theme["blue"] if v < 0 else theme["red"] for v in values]
        # Replicate SD as horizontal error bars, when present and non-zero.
        xerr = None
        if "ddG_sd" in df.columns:
            sd = (
                pd.to_numeric(df["ddG_sd"], errors="coerce").fillna(0.0).to_numpy(float)
            )
            if bool((sd > 0).any()):
                xerr = sd
        ax.barh(
            y,
            values,
            color=colors,
            height=0.68,
            zorder=3,
            xerr=xerr,
            error_kw={
                "ecolor": theme["ink_secondary"],
                "elinewidth": 0.8,
                "capsize": 2,
            },
        )
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()  # most-stabilizing (row 0) at TOP

        span = float(np.abs(values).max()) or 1.0
        pad = span * 0.16
        for yi, v in zip(y, values, strict=False):
            if v < 0:
                ax.text(
                    v - span * 0.015,
                    float(yi),
                    f"{v:.1f}",
                    ha="right",
                    va="center",
                    color=theme["ink_secondary"],
                    fontsize=8.5,
                )
            else:
                ax.text(
                    v + span * 0.015,
                    float(yi),
                    f"{v:.1f}",
                    ha="left",
                    va="center",
                    color=theme["ink_secondary"],
                    fontsize=8.5,
                )
        left = min(0.0, float(values.min())) - pad
        right = max(0.0, float(values.max())) + pad
        ax.set_xlim(left, right)

        ax.axvline(0.0, color=theme["baseline"], linewidth=1.2, zorder=2)
        ax.grid(axis="x", color=theme["grid"], linewidth=0.6)
        ax.set_title(
            f"Top {k} stabilizing single mutations",
            fontsize=13,
            pad=12,
            loc="left",
        )
        ax.set_xlabel("ΔΔG (REU)")
        ax.spines["left"].set_visible(False)
        _style_axes(ax, theme)
        ax.spines["left"].set_visible(False)
        fig.tight_layout()
        return fig_to_base64(fig)


# ---------------------------------------------------------------------------
# 3. Docking affinity
# ---------------------------------------------------------------------------
def affinity_plot(
    docking_df: pd.DataFrame | None, n: int = 20, *, dark: bool = False
) -> str | None:
    """Lollipop chart of the top-*n* best (most-negative) docking affinities."""
    if not _has_rows(docking_df, "affinity"):
        return None
    assert docking_df is not None
    df = docking_df.head(max(int(n), 1)).copy()
    df = df[pd.to_numeric(df["affinity"], errors="coerce").notna()]
    if df.empty:
        return None

    values = pd.to_numeric(df["affinity"], errors="coerce").to_numpy(float)
    if "name" in df.columns:
        labels = [_shorten_ligand(v) for v in df["name"].tolist()]
    else:
        labels = [str(i + 1) for i in range(len(values))]
    k = len(values)

    theme = _theme(dark)
    with plt.rc_context(_rc(theme)):
        height = max(3.0, 0.34 * k + 1.4)
        fig, ax = plt.subplots(figsize=(7.2, min(height, 9.0)))

        y = np.arange(k)
        ax.hlines(y, 0, values, color=theme["blue_seq"], linewidth=2.2, zorder=2)
        ax.scatter(values, y, s=46, color=theme["blue_seq"], zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()  # best (row 0) at TOP

        span = float(np.abs(values).max()) or 1.0
        for yi, v in zip(y, values, strict=False):
            ax.text(
                v - span * 0.02,
                float(yi),
                f"{v:.1f}",
                ha="right",
                va="center",
                color=theme["ink_secondary"],
                fontsize=8.5,
            )
        ax.set_xlim(float(values.min()) - span * 0.18, 0.0)

        ax.axvline(0.0, color=theme["baseline"], linewidth=1.2, zorder=1)
        ax.grid(axis="x", color=theme["grid"], linewidth=0.6)
        ax.set_title("Docking affinity (kcal/mol)", fontsize=13, pad=12, loc="left")
        ax.set_xlabel("affinity (kcal/mol, lower = stronger)")
        ax.spines["left"].set_visible(False)
        _style_axes(ax, theme)
        ax.spines["left"].set_visible(False)
        fig.tight_layout()
        return fig_to_base64(fig)


# ---------------------------------------------------------------------------
# 4. Structure backbone (3D Cα trace)
# ---------------------------------------------------------------------------
def structure_backbone(
    pdb_path: str | Path,
    sites: list[tuple[str, int]] | None = None,
    *,
    dark: bool = False,
) -> str | None:
    """3D Cα backbone trace, one line per chain, mutated sites highlighted."""
    trace = parse_ca_trace(pdb_path)
    if not trace:
        return None

    theme = _theme(dark)
    chain_colors = theme["chains"]

    all_x: list[float] = []
    all_y: list[float] = []
    all_z: list[float] = []

    with plt.rc_context(_rc(theme)):
        fig = plt.figure(figsize=(6.5, 6.0))
        ax = fig.add_subplot(projection="3d")
        ax.set_facecolor(theme["surface"])

        for idx, (chain_id, residues) in enumerate(sorted(trace.items())):
            xs = [r[1] for r in residues]
            ys = [r[2] for r in residues]
            zs = [r[3] for r in residues]
            all_x.extend(xs)
            all_y.extend(ys)
            all_z.extend(zs)
            ax.plot(
                xs,
                ys,
                zs,
                color=chain_colors[idx % len(chain_colors)],
                linewidth=1.6,
                label=f"Chain {chain_id}",
            )

        # Highlight requested (chain, position) sites.
        if sites:
            wanted = {(str(c), int(p)) for c, p in sites}
            for chain_id, residues in trace.items():
                for resseq, x, y, z, _resname in residues:
                    if (str(chain_id), int(resseq)) in wanted:
                        ax.scatter(
                            [x],
                            [y],
                            [z],
                            s=120,
                            color=theme["highlight"],
                            edgecolors=theme["ink_primary"],
                            linewidths=0.8,
                            depthshade=False,
                            zorder=6,
                        )
                        ax.text(
                            x,
                            y,
                            z,
                            f"  {chain_id}{resseq}",
                            color=theme["ink_secondary"],
                            fontsize=9,
                        )

        # Undistorted aspect: scale the box by the true per-axis data spans.
        if all_x:
            dx = max(all_x) - min(all_x) or 1.0
            dy = max(all_y) - min(all_y) or 1.0
            dz = max(all_z) - min(all_z) or 1.0
            ax.set_box_aspect((dx, dy, dz))

        # Recessive 3D chrome: transparent panes, no grid, no ticks.
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((0.0, 0.0, 0.0, 0.0))
            axis.line.set_color(theme["grid"])
            axis.line.set_linewidth(0.6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

        site_note = ""
        if sites:
            site_note = "  (" + ", ".join(f"{c}{p}" for c, p in sites) + ")"
        ax.set_title(
            f"Top mutant — Cα backbone{site_note}",
            fontsize=13,
            pad=10,
            loc="left",
            color=theme["ink_primary"],
        )
        if len(trace) > 1:
            ax.legend(loc="upper left", fontsize=9)
        return fig_to_base64(fig)
