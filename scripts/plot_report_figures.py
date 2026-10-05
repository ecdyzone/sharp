#!/usr/bin/env python3
"""Report tables → figures (PNG + SVG), for the wiki, a paper or slides.

Pure presentation: every number is read from the tidy TSVs that
`build_report_tables.py` writes, so this script can be re-run after editing a
table by hand (e.g. typing in S(H)ARP's numbers) without touching the
benchmark outputs. A figure whose table is missing is skipped with a warning,
so the MiBIG and the raw-comparison tables can live in different directories
or be built at different times.

Figures (each written as <name>.png and <name>.svg):

    mibig_metrics          headline rates per tool + how much each tool called
    mibig_by_class         MiBIG clusters per BGC class, and how many each tool found
    mibig_threshold_sweep  the same rates as a scored tool's cutoff rises
    raw_totals             regions called, share of sequence called, region length
    raw_by_class           regions per BGC class, per series
    raw_by_genus           mean regions per genome for the most-sampled genera
    raw_agreement          how much of each series' territory another series also called

**Placeholders.** A series whose rows have `status=placeholder` and nothing but
zeros (S(H)ARP until it has results) keeps its colour and its slot, is marked
"pending" in the plot and "(pending)" in the legend, and is drawn as soon as
any of its numbers is non-zero — so hand-entered numbers show up without
having to edit the status column.

Colours follow the tool, never its position: antiSMASH blue, DeepBGC orange,
SHARP aqua; extra cutoffs of the same tool reuse its hue, darker and hatched.
The palette was checked with a colour-vision-deficiency validator (adjacent
CVD ΔE ≥ 9).

Usage:
    python scripts/plot_report_figures.py \\
        --tables-dir data/processed/report/tables \\
        --output-dir data/processed/report/figures

    # only some figures, PNG only
    python scripts/plot_report_figures.py \\
        --tables-dir data/processed/report/tables \\
        --output-dir data/processed/report/figures \\
        --only mibig_by_class raw_by_class --formats png
"""
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter, PercentFormatter  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

LOG = logging.getLogger("plot_report_figures")

# ═══════════════════════════════ theme ═════════════════════════════════════
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
REFERENCE = "#bdbcb4"            # MiBIG / "known" bars: de-emphasis gray

# Tool identity. Fixed per tool, never by rank, so a tool keeps its colour
# when another is added or removed.
TOOL_COLORS = {
    "antiSMASH": "#2a78d6",      # blue
    "DeepBGC":   "#eb6834",      # orange
    "SHARP":     "#1baf7a",      # aqua
}
# Any other tool takes the next unused slot of the same validated palette.
SPARE_COLORS = ["#4a3aa7", "#e87ba4", "#008300", "#eda100", "#e34948"]

# Sequential ramp (one hue, light → dark) for the agreement heatmap.
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
              "#184f95", "#0d366b"]


def apply_theme() -> None:
    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": TEXT, "axes.labelcolor": TEXT_2,
        "axes.titlecolor": TEXT, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
        "axes.edgecolor": BASELINE, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False,
        "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
        "axes.grid.axis": "y", "axes.axisbelow": True,
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "axes.labelsize": 9.5,
        "xtick.labelsize": 9, "ytick.labelsize": 8.5,
        "legend.fontsize": 9, "legend.frameon": False,
        "hatch.color": SURFACE, "hatch.linewidth": 1.4,
        "svg.fonttype": "none",       # keep text as text in the SVG
    })


# ═══════════════════════════════ series ════════════════════════════════════

@dataclass(frozen=True)
class Style:
    """How one series is drawn."""
    label: str
    color: str
    hatch: str | None
    pending: bool

    @property
    def legend(self) -> str:
        return f"{self.label} (pending)" if self.pending else self.label


def _luminance(rgba) -> float:
    """WCAG relative luminance, to pick ink or white text on a filled cell."""
    def lin(c: float) -> float:
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(c) for c in to_rgb(rgba))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _darken(hex_color: str, amount: float) -> str:
    r, g, b = to_rgb(hex_color)
    return to_hex((r * (1 - amount), g * (1 - amount), b * (1 - amount)))


def is_pending(df: pd.DataFrame, value_cols: list[str]) -> bool:
    """A placeholder that nobody has filled in yet: flagged so, and all zero."""
    if df.empty or "status" not in df or not (df["status"] == "placeholder").all():
        return False
    values = df[value_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return bool((values == 0).all().all())


def styles_for(df: pd.DataFrame, key: str, value_cols: list[str],
               tool_col: str | None = None) -> list[Style]:
    """One Style per series, in table order.

    `key` names the series column ("tool" or "label"); `tool_col` the tool
    family it belongs to, when a tool has several series (score cutoffs). The
    first series of a family takes the family colour; later ones a darker,
    hatched variant of it.
    """
    tool_col = tool_col or key
    spare = iter(SPARE_COLORS)
    family_color: dict[str, str] = {}
    seen: dict[str, int] = {}
    out = []
    for name in dict.fromkeys(df[key]):
        rows = df[df[key] == name]
        tool = str(rows[tool_col].iloc[0])
        if tool not in family_color:
            family_color[tool] = TOOL_COLORS.get(tool) or next(spare, MUTED)
        n = seen.get(tool, 0)
        seen[tool] = n + 1
        color = family_color[tool] if n == 0 else _darken(family_color[tool], 0.28 * n)
        out.append(Style(str(name), color, None if n == 0 else "///",
                         is_pending(rows, value_cols)))
    return out


# ═══════════════════════════════ helpers ═══════════════════════════════════

def pct(v: float, digits: int | None = None) -> str:
    """94% · 6.1% · 0.4% — one decimal below 10%, where it carries the story."""
    if digits is None:
        digits = 1 if 0 < abs(v) < 0.1 else 0
    return f"{100 * v:.{digits}f}%"


def compact(v: float) -> str:
    """1,284 · 12.9k · 219k · 3.47M"""
    v = float(v)
    if abs(v) >= 1e6:
        return f"{v / 1e6:.3g}M"
    if abs(v) >= 1e4:
        return f"{v / 1e3:.0f}k"
    if abs(v) >= 1e3:
        return f"{v / 1e3:.1f}k"
    return f"{v:,.0f}"


def title(fig: plt.Figure, text: str, subtitle: str | None = None) -> None:
    fig.text(0.012, 0.985, text, ha="left", va="top", fontsize=13,
             fontweight="bold", color=TEXT)
    if subtitle:
        fig.text(0.012, 0.935, subtitle, ha="left", va="top", fontsize=9.5,
                 color=TEXT_2)


def grouped_bars(
    ax: plt.Axes, categories: list[str], styles: list[Style],
    values: dict[str, list[float]], *, fmt: Callable[[float], str] | None,
    cap: float | None = None, reference: tuple[str, list[float]] | None = None,
    label_size: float = 7.0,
) -> None:
    """Bars grouped by category, one per series (plus an optional gray
    reference bar first). Values above `cap` are drawn clipped at it with their
    true value and an arrow, so one huge bar cannot flatten every other one.
    Pending series get a 'pending' note in their slot instead of a bar."""
    bars_per_group = len(styles) + (1 if reference else 0)
    width = min(0.8 / bars_per_group, 0.24)
    x = np.arange(len(categories))
    slots: list[tuple[str, list[float], str, str | None, bool]] = []
    if reference:
        slots.append((reference[0], reference[1], REFERENCE, None, False))
    slots += [(s.legend, values[s.label], s.color, s.hatch, s.pending) for s in styles]

    top = cap if cap is not None else max(
        (v for _, vals, *_ in slots for v in vals), default=1) or 1
    for i, (legend, vals, color, hatch, pending) in enumerate(slots):
        pos = x + (i - (bars_per_group - 1) / 2) * width
        raw = np.asarray(vals, dtype=float)
        shown = np.minimum(raw, cap) if cap is not None else raw
        ax.bar(pos, shown, width, color=color, hatch=hatch, label=legend,
               edgecolor=SURFACE, linewidth=1.0, zorder=2)
        for p, v, h in zip(pos, raw, shown):
            if pending:
                ax.text(p, top * 0.015, "pending", rotation=90, ha="center",
                        va="bottom", fontsize=6.5, color=MUTED, zorder=3)
            elif cap is not None and v > cap:
                ax.text(p, h - top * 0.02, f"↑ {fmt(v) if fmt else compact(v)}",
                        rotation=90,
                        ha="center", va="top", fontsize=label_size,
                        color=SURFACE, fontweight="bold", zorder=3)
            elif fmt is not None:
                ax.text(p, h + top * 0.012, fmt(v), ha="center", va="bottom",
                        fontsize=label_size, color=TEXT_2, zorder=3)
    ax.set_xticks(x, categories)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, len(categories) - 0.5)


def legend_below(fig: plt.Figure, ax: plt.Axes, ncol: int | None = None) -> None:
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=ncol or len(labels),
               bbox_to_anchor=(0.5, 0.0), handlelength=1.2, columnspacing=1.6)


def legend_top(fig: plt.Figure, ax: plt.Axes, y: float = 0.895) -> None:
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", ncol=len(labels),
               bbox_to_anchor=(0.005, y), handlelength=1.2, columnspacing=1.6)


# ═══════════════════════════════ MiBIG ═════════════════════════════════════

MIBIG_RATES = [
    ("detection_recall", "Clusters found\n(detection recall)"),
    ("reciprocal_recall", "Found with tight\nboundaries (reciprocal)"),
    ("nucleotide_recall", "Cluster bp\ncovered"),
    ("median_prediction_coverage", "Median share of a\nmatch that is cluster"),
    ("matched_prediction_frac", "Calls on a\nknown cluster"),
    ("nucleotide_precision", "Called bp on a\nknown cluster"),
]


def _mibig_scope(df: pd.DataFrame) -> str:
    row = df.iloc[0]
    return (f"{int(row['n_clusters']):,} MiBIG clusters on "
            f"{int(row['n_contigs']):,} sequences")


def fig_mibig_metrics(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["mibig_metrics"]
    value_cols = [c for c, _ in MIBIG_RATES] + ["n_predictions", "predicted_bp"]
    styles = styles_for(df, "tool", value_cols)
    by_tool = df.set_index("tool")

    fig = plt.figure(figsize=(11, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.62,
                          wspace=0.42, left=0.12, right=0.985, top=0.80,
                          bottom=0.07)
    ax = fig.add_subplot(gs[0, :])
    grouped_bars(ax, [lbl for _, lbl in MIBIG_RATES], styles,
                 {s.label: [float(by_tool.loc[s.label, c]) for c, _ in MIBIG_RATES]
                  for s in styles},
                 fmt=lambda v: pct(v))
    ax.set_ylim(0, 1.1)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_title("a   How well each tool recovers the known clusters", pad=10)
    ax.axvline(3.5, color=BASELINE, linewidth=0.8, zorder=1)
    ax.text(1.5, 1.07, "recall side: of the known clusters…", ha="center",
            fontsize=8.5, color=MUTED)
    ax.text(4.5, 1.07, "call side: of what the tool called…", ha="center",
            fontsize=8.5, color=MUTED)

    known_n = int(df["n_clusters"].iloc[0])
    known_mb = float(df["gt_bp"].iloc[0]) / 1e6
    for col, (field, ref, ref_label, unit, label) in enumerate([
        ("n_predictions", known_n, "MiBIG clusters", "",
         "b   Regions called"),
        ("predicted_bp", known_mb, "MiBIG clusters", " Mb",
         "c   Sequence called (Mb)"),
    ]):
        axv = fig.add_subplot(gs[1, col])
        scale = 1e6 if field == "predicted_bp" else 1
        vals = [ref] + [float(by_tool.loc[s.label, field]) / scale for s in styles]
        names = [ref_label] + [s.legend for s in styles]
        colors = [REFERENCE] + [s.color for s in styles]
        y = np.arange(len(vals))[::-1]
        axv.barh(y, vals, height=0.62, color=colors, edgecolor=SURFACE, zorder=2)
        top = max(vals) or 1
        for yy, v, s_pending in zip(y, vals, [False] + [s.pending for s in styles]):
            txt = "pending" if s_pending else (
                f"{v:,.1f}{unit}" if unit else f"{v:,.0f}")
            axv.text(v + top * 0.015, yy, txt, va="center", fontsize=8,
                     color=MUTED if s_pending else TEXT_2)
        axv.set_yticks(y, names)
        axv.set_xlim(0, top * 1.18)
        axv.grid(axis="x"); axv.grid(axis="y", visible=False)
        axv.tick_params(axis="y", length=0)
        axv.set_title(label, pad=8)
        axv.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    legend_top(fig, ax, y=0.875)
    title(fig, "Benchmark against MiBIG",
          _mibig_scope(df) + " · a cluster is found when one call covers ≥50% of it")
    return fig


def fig_mibig_by_class(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["mibig_by_class"]
    classes = list(dict.fromkeys(df["class"]))
    styles = styles_for(df, "tool", ["n_recovered"])
    known = df.drop_duplicates("class").set_index("class")["n_clusters"]
    values = {s.label: df[df["tool"] == s.label].set_index("class")
              .reindex(classes)["n_recovered"].astype(float).tolist()
              for s in styles}
    labels = [f"{c}\n(n={int(known[c])})" for c in classes]

    fig, ax = plt.subplots(figsize=(11, 5.6))
    fig.subplots_adjust(left=0.06, right=0.985, top=0.84, bottom=0.12)
    grouped_bars(ax, labels, styles, values, fmt=lambda v: f"{v:.0f}",
                 reference=("MiBIG clusters", known.reindex(classes).astype(float).tolist()))
    ax.set_ylabel("Clusters")
    ax.set_ylim(0, float(known.max()) * 1.12)
    legend_top(fig, ax, y=0.90)
    title(fig, "Known clusters found per BGC class",
          f"{int(known.sum()):,} MiBIG clusters · gray = clusters in MiBIG, "
          "colours = clusters each tool found")
    return fig


SWEEP_PANELS = [
    ("detection_recall", "a   Clusters found", True),
    ("matched_prediction_frac", "b   Calls on a known cluster", True),
    ("n_predictions", "c   Regions called", False),
]


def fig_mibig_threshold_sweep(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["mibig_threshold_sweep"]
    styles = styles_for(df, "tool", [c for c, _, _ in SWEEP_PANELS])
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.4))
    fig.subplots_adjust(left=0.06, right=0.985, top=0.72, bottom=0.15, wspace=0.28)
    for ax, (col, label, is_rate) in zip(axes, SWEEP_PANELS):
        for s in styles:
            if s.pending:
                continue
            d = df[df["tool"] == s.label].sort_values("min_p_bgc")
            flat = d[col].nunique() == 1
            ax.plot(d["min_p_bgc"], d[col], color=s.color, linewidth=2,
                    linestyle=(0, (4, 3)) if flat else "-",
                    marker=None if flat else "o", markersize=4.5,
                    markeredgecolor=SURFACE, markeredgewidth=1.2,
                    label=f"{s.label} (no score: flat)" if flat else s.label,
                    solid_capstyle="round", zorder=3)
        ax.set_title(label, pad=8)
        ax.set_xlabel("Minimum score kept (p_bgc)")
        ax.grid(axis="x", visible=False)
        top = pd.to_numeric(df[col], errors="coerce").max() or 1
        ax.set_ylim(0, top * 1.12)
        if is_rate:
            ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        else:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    handles, labels = axes[0].get_legend_handles_labels()
    for s in styles:
        if s.pending:
            handles.append(Patch(color=s.color)); labels.append(s.legend)
    fig.legend(handles, labels, loc="upper left", ncol=len(labels),
               bbox_to_anchor=(0.005, 0.86), handlelength=1.8, columnspacing=1.6)
    title(fig, "What a score cutoff buys",
          "Each point keeps only the calls scoring at least x · a tool with no "
          "score (antiSMASH) is the same at every cutoff, drawn dashed")
    return fig


# ═══════════════════════════════ raw ═══════════════════════════════════════

def _raw_scope(df: pd.DataFrame) -> str:
    row = df.iloc[0]
    return (f"{int(row['n_assemblies']):,} genomes · "
            f"{float(row['total_bp']) / 1e9:.1f} Gb · no ground truth")


def fig_raw_totals(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["raw_totals"]
    panels = [
        ("n_regions", "a   Regions called", compact),
        ("frac_bp_called", "b   Share of all sequence called", lambda v: pct(v, 1)),
        ("median_region_bp", "c   Median region length (kb)",
         lambda v: f"{v / 1e3:.1f} kb"),
    ]
    styles = styles_for(df, "label", [c for c, _, _ in panels], tool_col="tool")
    rows = df.set_index("label")
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.76, bottom=0.08, wspace=0.3)
    y = np.arange(len(styles))[::-1]
    for i, (ax, (col, label, fmt)) in enumerate(zip(axes, panels)):
        vals = [float(rows.loc[s.label, col]) for s in styles]
        for bar_y, v, s in zip(y, vals, styles):
            ax.barh(bar_y, v, height=0.62, color=s.color, hatch=s.hatch,
                    edgecolor=SURFACE, zorder=2)
        top = max(vals) or 1
        for bar_y, v, s in zip(y, vals, styles):
            ax.text(v + top * 0.02, bar_y, "pending" if s.pending else fmt(v),
                    va="center", fontsize=8.5, color=MUTED if s.pending else TEXT_2)
        ax.set_yticks(y, [s.legend for s in styles] if i == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(0, top * 1.3)
        ax.grid(axis="x"); ax.grid(axis="y", visible=False)
        ax.set_title(label, pad=8)
        ax.xaxis.set_major_formatter(FuncFormatter(
            lambda v, _, c=col: f"{100 * v:.0f}%" if c == "frac_bp_called"
            else (f"{v / 1e3:.0f}" if c == "median_region_bp" else compact(v))))
    title(fig, "How much each tool calls", _raw_scope(df))
    return fig


def fig_raw_by_class(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["raw_by_class"]
    styles = styles_for(df, "label", ["n_regions"], tool_col="tool")
    totals = df.groupby("class", sort=False)["n_regions"].sum()
    classes = [c for c in dict.fromkeys(df["class"]) if totals[c] > 0]
    values = {s.label: df[df["label"] == s.label].set_index("class")
              .reindex(classes)["n_regions"].astype(float).tolist()
              for s in styles}
    # One class can dwarf the rest (DeepBGC leaves most regions unclassified):
    # cap the axis at the largest *classified* bar and clip the outliers.
    classified = [v for s in styles for c, v in zip(classes, values[s.label])
                  if c != "Unclassified"]
    biggest = max(v for vals in values.values() for v in vals) if values else 0
    cap = max(classified) * 1.15 if classified and biggest > 1.6 * max(classified) else None

    fig, ax = plt.subplots(figsize=(11, 5.4))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.08)
    grouped_bars(ax, classes, styles, values, fmt=None, cap=cap)
    if cap is not None:
        ax.set_ylim(0, cap)
    ax.set_ylabel("Regions")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    legend_top(fig, ax, y=0.90)
    title(fig, "Regions called per BGC class",
          "Each tool's own class names mapped onto MiBIG's classes · "
          "Hybrid = more than one class · clipped bars show their true value")
    return fig


def fig_raw_by_genus(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["raw_by_genus"]
    col = "mean_regions_per_genome"
    styles = styles_for(df, "label", [col], tool_col="tool")
    genera = list(dict.fromkeys(df["genus"]))
    n = df.drop_duplicates("genus").set_index("genus")["n_assemblies"]
    values = {s.label: df[df["label"] == s.label].set_index("genus")
              .reindex(genera)[col].astype(float).tolist() for s in styles}
    fig, ax = plt.subplots(figsize=(12.5, 5.6))
    fig.subplots_adjust(left=0.06, right=0.99, top=0.84, bottom=0.14)
    grouped_bars(ax, genera, styles, values, fmt=None)
    ax.set_ylabel("Mean regions per genome")
    # Genus names are italic by convention; pooled/unknown groups are not genera.
    for tick, g in zip(ax.get_xticklabels(), genera):
        tick.set_fontstyle("normal" if g in ("Other genera", "Unknown") else "italic")
        tick.set_fontsize(8.5)
    under = blended_transform_factory(ax.transData, ax.transAxes)
    for x, g in enumerate(genera):
        ax.text(x, -0.06, f"n={int(n[g]):,}", transform=under, ha="center",
                va="top", fontsize=7.5, color=MUTED)
    legend_top(fig, ax, y=0.90)
    title(fig, "Regions called per genome, by genus",
          "The genera with the most genomes in the database; n = genomes")
    return fig


def fig_raw_agreement(t: dict[str, pd.DataFrame]) -> plt.Figure:
    df = t["raw_agreement"]
    order = list(dict.fromkeys(df["reference"]))
    grid = df.pivot(index="reference", columns="query", values="bp_agreement") \
        .reindex(index=order, columns=order).apply(pd.to_numeric, errors="coerce")
    regions = df.pivot(index="reference", columns="query", values="region_agreement") \
        .reindex(index=order, columns=order).apply(pd.to_numeric, errors="coerce")
    pending = {lbl for lbl in order
               if is_pending(df[df["reference"] == lbl], ["n_reference_regions"])}

    cmap = LinearSegmentedColormap.from_list("seq", SEQUENTIAL).with_extremes(bad=GRID)
    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    fig.subplots_adjust(left=0.2, right=0.93, top=0.80, bottom=0.17)
    values = grid.to_numpy(dtype=float)
    np.fill_diagonal(values, np.nan)          # a tool against itself says nothing
    im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=0, vmax=1,
                   aspect="auto")
    for i, ref in enumerate(order):
        for j, query in enumerate(order):
            v, r = grid.iloc[i, j], regions.iloc[i, j]
            if i == j:
                txt, color = "—", MUTED
            elif np.isnan(v):
                txt, color = "pending" if (ref in pending or query in pending) else "n/a", MUTED
            else:
                txt = f"{pct(v)}\n({pct(r)} of regions)"
                color = SURFACE if _luminance(cmap(v)) < 0.4 else TEXT
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.5, color=color)
    tick = [f"{lbl} (pending)" if lbl in pending else lbl for lbl in order]
    ax.set_xticks(range(len(order)), tick, rotation=25, ha="right",
                  rotation_mode="anchor")
    ax.set_yticks(range(len(order)), tick)
    ax.set_xlabel("…also called by this tool")
    ax.set_ylabel("Of the sequence this tool called…")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(len(order) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(order) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="both", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=8, length=0)
    title(fig, "Do the tools call the same regions?",
          "Share of the row tool's called sequence (bp) that the column tool "
          "also called;\nin brackets, the share of its regions the column tool "
          "covers by ≥50% with one call")
    return fig


# ═══════════════════════════════ driver ════════════════════════════════════

FIGURES: dict[str, tuple[str, Callable[[dict[str, pd.DataFrame]], plt.Figure]]] = {
    "mibig_metrics": ("mibig_metrics", fig_mibig_metrics),
    "mibig_by_class": ("mibig_by_class", fig_mibig_by_class),
    "mibig_threshold_sweep": ("mibig_threshold_sweep", fig_mibig_threshold_sweep),
    "raw_totals": ("raw_totals", fig_raw_totals),
    "raw_by_class": ("raw_by_class", fig_raw_by_class),
    "raw_by_genus": ("raw_by_genus", fig_raw_by_genus),
    "raw_agreement": ("raw_agreement", fig_raw_agreement),
}


def load_table(tables_dir: Path, name: str) -> pd.DataFrame | None:
    path = tables_dir / f"{name}.tsv"
    if not path.is_file():
        return None
    return pd.read_csv(path, sep="\t", keep_default_na=False,
                       na_values=[""], dtype={"status": str})


def run(tables_dir: Path, output_dir: Path, formats: list[str], dpi: int,
        only: list[str] | None = None) -> list[Path]:
    apply_theme()
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, (table, draw) in FIGURES.items():
        if only and name not in only:
            continue
        df = load_table(tables_dir, table)
        if df is None or df.empty:
            LOG.warning("%s: no %s.tsv in %s — figure skipped", name, table, tables_dir)
            continue
        fig = draw({table: df})
        for fmt in formats:
            path = output_dir / f"{name}.{fmt}"
            fig.savefig(path, dpi=dpi)
            written.append(path)
        plt.close(fig)
        LOG.info("%s → %s", name, ", ".join(f"{name}.{f}" for f in formats))
    return written


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--tables-dir", type=Path, required=True,
                   help="directory holding the build_report_tables.py TSVs")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="where the figures go")
    p.add_argument("--formats", nargs="+", default=["png", "svg"],
                   help="file formats to write (default: png svg)")
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    p.add_argument("--only", nargs="+", choices=sorted(FIGURES), default=None,
                   metavar="FIGURE", help=f"draw only these: {', '.join(FIGURES)}")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    written = run(args.tables_dir, args.output_dir, args.formats, args.dpi, args.only)
    if not written:
        raise SystemExit(f"no figures written — no tables found in {args.tables_dir}")


if __name__ == "__main__":
    main()
