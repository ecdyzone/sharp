#!/usr/bin/env python3
"""Benchmark artifacts → the tidy tables behind the report figures.

The figures in a report (wiki, paper, slides) are drawn by
`plot_report_figures.py` from the small TSVs this script writes. Splitting the
two is deliberate:

- **Numbers come from the tested metric code, not from plotting code.** Every
  value here is produced by `sharp.metrics.evaluate_predictions` (MiBIG) or by
  `summarize_predictions.summarize` (no ground truth) — the same functions
  behind `benchmark.json` and `summary_by_tool.tsv`, so a figure can never
  disagree with the result file it illustrates.
- **A tool with no results yet is still drawn.** Each `--placeholder` name (by
  default `SHARP`) gets all-zero rows with `status=placeholder`, so the
  figures already have its slot. Fill it either way:
    1. pass its predictions like any other tool (`SHARP=path.parquet`); the
       placeholder is dropped automatically, or
    2. type its numbers into the TSVs by hand and re-run only the plotter.

Two subcommands, one per comparison:

`mibig` — tools against MiBIG ground truth. Same inputs and same defaults as
`sharp.evaluate` (ground truth, `--contigs` scope, match criterion), so the
headline numbers reproduce that step's JSON exactly. Writes:

    mibig_metrics.tsv          one row per tool: the headline rates and counts
    mibig_by_class.tsv         clusters recovered per BGC class, per tool
    mibig_threshold_sweep.tsv  the rates at each p_bgc cutoff (scored tools move,
                               antiSMASH is flat — it has no score)

`raw` — tools against each other over a genome database with no ground truth.
Nothing here is recall or precision; see summarize_predictions.py for why raw
counts mislead and how these tables normalise for it. Writes:

    raw_totals.tsv     one row per series: regions, bp called, per-Mb rates
    raw_by_class.tsv   regions per BGC class, per series
    raw_by_genus.tsv   mean regions per genome for the most-sampled genera
    raw_agreement.tsv  every ordered pair of series: the fraction of the
                       reference's regions / bp the query also called

A *series* is one tool at one score cutoff. `--thresholds DeepBGC=0.5,0.8`
turns DeepBGC into two series ("DeepBGC ≥0.5", "DeepBGC ≥0.8"); a tool with no
score (every `p_bgc` identical, i.e. antiSMASH) is always a single series.

BGC classes are put on one shared vocabulary — PKS, NRPS, RiPP, Terpene,
Saccharide, Other, Hybrid (more than one of those), Unclassified (no class) —
from antiSMASH product types, DeepBGC `product_class` and MiBIG classes alike.
See CLASS_LOOKUP below.

Usage:
    # 1. MiBIG benchmark — the same ground truth and scope as sharp.evaluate
    python scripts/build_report_tables.py mibig \\
        --ground-truth data/interim/benchmark_set/benchmark_ground_truth.tsv \\
        --contigs data/interim/benchmark_set/analyzed_contigs.txt \\
        --predictions antiSMASH=data/interim/antismash_predictions_benchmark_set.parquet \\
                      DeepBGC=data/interim/deepbgc_predictions_benchmark_set.parquet \\
        --output-dir data/processed/report/tables

    # 2. Raw comparison over a genome database (no ground truth)
    python scripts/build_report_tables.py raw \\
        --manifest data/interim/actino_db/genomes.tsv \\
        --assemblies data/interim/actino_db/assemblies.tsv \\
        --predictions antiSMASH=data/interim/antismash_predictions_actino.parquet \\
                      DeepBGC=data/interim/deepbgc_predictions_actino.parquet \\
        --thresholds DeepBGC=0.5,0.8 \\
        --output-dir data/processed/report/tables

    # Once S(H)ARP has predictions, add it like any other tool — the default
    # SHARP placeholder is then dropped:
    #     --predictions ... SHARP=data/interim/sharptool_predictions_actino.parquet
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from build_genome_manifest import load_contig_index  # noqa: E402
from predictions_to_ground_truth import regions_to_clusters  # noqa: E402
from summarize_predictions import check_contigs, parse_spec, summarize  # noqa: E402
from sharp.io import (  # noqa: E402
    KnownCluster,
    PredictedRegion,
    load_contigs,
    load_ground_truth_tsv,
    load_predictions_parquet,
)
from sharp.metrics import (  # noqa: E402
    BenchmarkResult,
    MatchCriterion,
    evaluate_predictions,
    merge_intervals,
)

LOG = logging.getLogger("build_report_tables")

MEASURED = "measured"
PLACEHOLDER = "placeholder"

# The sweep grid. DeepBGC's own output cutoff is 0.5, so nothing moves below it.
DEFAULT_SWEEP = [round(0.5 + 0.05 * i, 2) for i in range(10)]   # 0.50 … 0.95


# ═════════════════════════ BGC class vocabulary ════════════════════════════
# Every tool-format assumption about class names lives in this block.
#
# Display categories, in plotting order. Hybrid = a region whose names resolve
# to more than one category; Unclassified = a region with no class at all
# (DeepBGC leaves ~65% of its candidates without one).
CLASS_ORDER = ["PKS", "NRPS", "RiPP", "Terpene", "Saccharide", "Other",
               "Hybrid", "Unclassified"]

# antiSMASH product types → the CATEGORY antiSMASH itself assigns each rule.
# Transcribed from antismash 8-0-stable,
# antismash/detection/hmm_detection/cluster_rules/{strict,relaxed,loose}.txt
# (every RULE's CATEGORY line). Using antiSMASH's own grouping, rather than a
# hand-made one, keeps the classification out of our hands.
_ANTISMASH_CATEGORIES = {
    "PKS": """arylpolyene benzoxazole hglE-KS HR-T2PKS ladderane PKS-like
        PpyS-KS prodigiosin PUFA T1PKS T2PKS T3PKS transAT-PKS
        transAT-PKS-like""",
    "NRPS": """CDPS fungal_CDPS isocyanide-nrp lysine mycosporine NAPAA
        NRP-metallophore NRPS NRPS-like RCDPS thioamide-NRP""",
    "RiPP": """archaeal-RiPP atropopeptide azole-containing-RiPP bottromycin
        crocagin cyanobactin cyclic-lactone-autoinducer darobactin epipeptide
        fungal-RiPP fungal-RiPP-like glycocin guanidinotides
        lanthipeptide-class-i lanthipeptide-class-ii lanthipeptide-class-iii
        lanthipeptide-class-iv lanthipeptide-class-v lassopeptide linaridin
        lipolanthine methanobactin microviridin proteusin ranthipeptide
        RaS-RiPP redox-cofactor RiPP-like RRE-containing sactipeptide
        spliceotide thioamitides triceptide""",
    "Terpene": "quinone_isoprenoid_chain terpene terpene-precursor",
    "Saccharide": "oligosaccharide saccharide",
    "Other": """2dos acyl_amino_acids amglyccycl aminocoumarin
        aminopolycarboxylic-acid azoxy-crosslink azoxy-dimer betalactone
        blactam butyrolactone cytokinin deazapurine ectoine fatty_acid furan
        halogenated hserlactone hydrogen-cyanide hydroxytropolone indole
        isocyanide leupeptin lincosamides melanin NAGGN NI-siderophore
        nitropropanoic_acid nucleoside opine-like-metallophore other PBDE
        phenazine phosphoglycolipid phosphonate phosphonate-like
        polyhalogenated-pyrrole polyyne pyrrolidine resorcinol
        tropodithietic-acid""",
}

# DeepBGC product_class parts (joined with "-", e.g. "NRP-Polyketide") and
# MiBIG classes (joined with "/", e.g. "NRPS/PKS"; "ribosomal" = RiPP).
_OTHER_VOCABULARIES = {
    "polyketide": "PKS", "nrp": "NRPS", "ripp": "RiPP", "alkaloid": "Other",
    "pks": "PKS", "ribosomal": "RiPP",
}

CLASS_LOOKUP: dict[str, str] = {
    name.lower(): category
    for category, names in _ANTISMASH_CATEGORIES.items()
    for name in names.split()
}
CLASS_LOOKUP.update(_OTHER_VOCABULARIES)
# Lower-cased, so DeepBGC's "Terpene"/"Saccharide"/"Other" and MiBIG's "NRPS"
# resolve through the antiSMASH entries of the same name.

_LIST_SEPARATORS = re.compile(r"[;/]")   # antiSMASH ";", MiBIG "/"
# ═══════════════════════════════════════════════════════════════════════════


def _name_categories(name: str, unknown: Counter | None) -> set[str]:
    """Categories one class name resolves to. Whole name first (antiSMASH
    names contain "-"), then DeepBGC-style "-" parts; anything unresolved is
    counted in `unknown` and filed under Other."""
    key = name.strip().lower()
    if not key:
        return set()
    if key in CLASS_LOOKUP:
        return {CLASS_LOOKUP[key]}
    parts = [p for p in key.split("-") if p]
    if parts and all(p in CLASS_LOOKUP for p in parts):
        return {CLASS_LOOKUP[p] for p in parts}
    if unknown is not None:
        unknown[name.strip()] += 1
    return {"Other"}


def normalize_class(raw: str | None, unknown: Counter | None = None) -> str:
    """Map any tool's class string onto CLASS_ORDER."""
    if raw is None:
        return "Unclassified"
    categories: set[str] = set()
    for name in _LIST_SEPARATORS.split(raw):
        categories |= _name_categories(name, unknown)
    if not categories:
        return "Unclassified"
    if len(categories) > 1:
        return "Hybrid"
    return categories.pop()


def _log_unknown(unknown: Counter, where: str) -> None:
    if unknown:
        LOG.warning("%s: %d class name(s) not in CLASS_LOOKUP, filed under "
                    "Other — add them if they are real: %s", where,
                    len(unknown), ", ".join(f"{k} ({v})" for k, v in unknown.most_common()))


# ═══════════════════════════════ shared ════════════════════════════════════

def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> int:
    """Write rows as TSV, columns in the given order. Floats keep 6 significant
    digits — plenty for a figure, short enough to edit by hand."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: _fmt(row.get(k)) for k in columns})
    LOG.info("wrote %d row(s) → %s", len(rows), path)
    return len(rows)


def _fmt(v: object) -> object:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


def load_tools(
    specs: list[str], placeholders: list[str]
) -> tuple[list[tuple[str, list[PredictedRegion]]], list[str]]:
    """Read every NAME=PATH spec; drop placeholders that now have data."""
    tools: list[tuple[str, list[PredictedRegion]]] = []
    for spec in specs:
        name, path = parse_spec(spec)
        if not path.is_file():
            raise SystemExit(f"no such predictions file: {path}")
        tools.append((name, load_predictions_parquet(path)))
        LOG.info("%s: %d region(s) from %s", name, len(tools[-1][1]), path)
    names = [n for n, _ in tools]
    if len(set(names)) != len(names):
        raise SystemExit(f"--predictions names must be unique, got {names}")
    pending = []
    for name in placeholders:
        if name in names:
            LOG.info("%s has predictions — its placeholder is dropped", name)
        elif name not in pending:
            pending.append(name)
    for name in pending:
        LOG.warning("%s: no predictions given — writing all-zero placeholder "
                    "rows (status=%s)", name, PLACEHOLDER)
    return tools, pending


# ═══════════════════════════════ mibig ═════════════════════════════════════

MIBIG_METRIC_COLUMNS = [
    "tool", "status",
    "detection_recall", "reciprocal_recall", "matched_prediction_frac",
    "nucleotide_recall", "nucleotide_precision", "median_prediction_coverage",
    "n_clusters", "n_recovered", "n_recovered_reciprocal",
    "n_predictions", "n_matched_predictions",
    "gt_bp", "predicted_bp", "intersect_bp",
    "n_clusters_recovered_by_union_only", "n_merged_predictions",
    "n_contigs", "min_p_bgc",
]
MIBIG_CLASS_COLUMNS = ["class", "n_clusters", "tool", "status",
                       "n_recovered", "recall"]
MIBIG_SWEEP_COLUMNS = [
    "tool", "status", "min_p_bgc", "n_predictions",
    "detection_recall", "reciprocal_recall", "matched_prediction_frac",
    "nucleotide_recall", "nucleotide_precision",
]


@dataclass(frozen=True)
class MibigSettings:
    """The knobs sharp.evaluate exposes, so the tables reproduce its JSON."""
    criterion: MatchCriterion = MatchCriterion()
    reciprocal_frac: float = 0.5
    min_p_bgc: float = 0.0


def _evaluate(
    preds: list[PredictedRegion], gt: list[KnownCluster], scope: set[str],
    s: MibigSettings, min_p_bgc: float | None = None,
) -> BenchmarkResult:
    return evaluate_predictions(
        preds, gt, scope=scope, scope_source="explicit",
        criterion=s.criterion, reciprocal_frac=s.reciprocal_frac,
        min_p_bgc=s.min_p_bgc if min_p_bgc is None else min_p_bgc,
        max_listed_ids=0,
    )


def mibig_metric_row(tool: str, r: BenchmarkResult) -> dict:
    """One tool's headline row, field for field from its BenchmarkResult."""
    return {
        "tool": tool, "status": MEASURED,
        "detection_recall": r.detection.recall,
        "reciprocal_recall": r.reciprocal.recall,
        "matched_prediction_frac": r.detection.matched_prediction_frac,
        "nucleotide_recall": r.nucleotide.recall,
        "nucleotide_precision": r.nucleotide.precision,
        "median_prediction_coverage": r.boundary.median_prediction_coverage,
        "n_clusters": r.detection.n_clusters,
        "n_recovered": r.detection.n_recovered,
        "n_recovered_reciprocal": r.reciprocal.n_recovered,
        "n_predictions": r.detection.n_predictions,
        "n_matched_predictions": r.detection.n_matched_predictions,
        "gt_bp": r.nucleotide.gt_bp,
        "predicted_bp": r.nucleotide.predicted_bp,
        "intersect_bp": r.nucleotide.intersect_bp,
        "n_clusters_recovered_by_union_only":
            r.boundary.n_clusters_recovered_by_union_only,
        "n_merged_predictions": r.boundary.n_merged_predictions,
        "n_contigs": r.scope.n_contigs,
        "min_p_bgc": r.scope.min_p_bgc,
    }


def mibig_placeholder_row(tool: str, gt_in_scope: list[KnownCluster],
                          n_contigs: int, min_p_bgc: float) -> dict:
    """All zeros, except what describes the ground truth rather than the tool."""
    row = {k: 0 for k in MIBIG_METRIC_COLUMNS}
    row.update({
        "tool": tool, "status": PLACEHOLDER,
        "n_clusters": len(gt_in_scope),
        "gt_bp": _union_bp((c.contig, c.start, c.end) for c in gt_in_scope),
        "n_contigs": n_contigs, "min_p_bgc": min_p_bgc,
    })
    return row


def _union_bp(intervals: Iterable[tuple[str, int, int]]) -> int:
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, start, end in intervals:
        by_contig.setdefault(contig, []).append((start, end))
    return sum(e - s for spans in by_contig.values()
               for s, e in merge_intervals(spans))


def gt_by_class(gt: list[KnownCluster]) -> dict[str, list[KnownCluster]]:
    """Ground truth bucketed by display class, in CLASS_ORDER."""
    unknown: Counter = Counter()
    out: dict[str, list[KnownCluster]] = {c: [] for c in CLASS_ORDER}
    for c in gt:
        out[normalize_class(c.cluster_class, unknown)].append(c)
    _log_unknown(unknown, "ground truth")
    return {k: v for k, v in out.items() if v}


def mibig_class_rows(
    tool: str, preds: list[PredictedRegion] | None,
    classes: dict[str, list[KnownCluster]], scope: set[str], s: MibigSettings,
) -> list[dict]:
    """Clusters recovered per class. Each class is scored as its own ground
    truth, so a cluster counts as recovered exactly when it does in the
    headline number — the per-class counts sum to `n_recovered`."""
    rows = []
    for cls, clusters in classes.items():
        n = len(clusters)
        if preds is None:
            rows.append({"class": cls, "n_clusters": n, "tool": tool,
                         "status": PLACEHOLDER, "n_recovered": 0, "recall": 0.0})
            continue
        det = _evaluate(preds, clusters, scope, s).detection
        rows.append({"class": cls, "n_clusters": n, "tool": tool,
                     "status": MEASURED, "n_recovered": det.n_recovered,
                     "recall": det.recall})
    return rows


def sweep_rows(
    tool: str, preds: list[PredictedRegion] | None, gt: list[KnownCluster],
    scope: set[str], s: MibigSettings, thresholds: list[float],
) -> list[dict]:
    """The headline rates at each p_bgc cutoff."""
    rows = []
    for t in thresholds:
        if preds is None:
            rows.append({"tool": tool, "status": PLACEHOLDER, "min_p_bgc": t,
                         **{k: 0 for k in MIBIG_SWEEP_COLUMNS[3:]}})
            continue
        r = _evaluate(preds, gt, scope, s, min_p_bgc=t)
        rows.append({
            "tool": tool, "status": MEASURED, "min_p_bgc": t,
            "n_predictions": r.detection.n_predictions,
            "detection_recall": r.detection.recall,
            "reciprocal_recall": r.reciprocal.recall,
            "matched_prediction_frac": r.detection.matched_prediction_frac,
            "nucleotide_recall": r.nucleotide.recall,
            "nucleotide_precision": r.nucleotide.precision,
        })
    return rows


def run_mibig(
    ground_truth: Path, contigs: Path, specs: list[str], placeholders: list[str],
    thresholds: list[float], output_dir: Path, s: MibigSettings,
) -> None:
    gt = load_ground_truth_tsv(ground_truth)
    scope = load_contigs(contigs)
    gt_in_scope = [c for c in gt if c.contig in scope]
    LOG.info("ground truth: %d cluster(s), %d on the %d contig(s) in scope",
             len(gt), len(gt_in_scope), len(scope))
    classes = gt_by_class(gt_in_scope)

    tools, pending = load_tools(specs, placeholders)
    metric_rows, class_rows, sweep = [], [], []
    for name, preds in tools:
        r = _evaluate(preds, gt, scope, s)
        LOG.info("%s: detection recall %.3f (%d/%d), %d prediction(s)", name,
                 r.detection.recall, r.detection.n_recovered,
                 r.detection.n_clusters, r.detection.n_predictions)
        metric_rows.append(mibig_metric_row(name, r))
        class_rows += mibig_class_rows(name, preds, classes, scope, s)
        sweep += sweep_rows(name, preds, gt, scope, s, thresholds)
    for name in pending:
        metric_rows.append(
            mibig_placeholder_row(name, gt_in_scope, len(scope), s.min_p_bgc))
        class_rows += mibig_class_rows(name, None, classes, scope, s)
        sweep += sweep_rows(name, None, gt, scope, s, thresholds)

    write_tsv(output_dir / "mibig_metrics.tsv", metric_rows, MIBIG_METRIC_COLUMNS)
    write_tsv(output_dir / "mibig_by_class.tsv", class_rows, MIBIG_CLASS_COLUMNS)
    write_tsv(output_dir / "mibig_threshold_sweep.tsv", sweep, MIBIG_SWEEP_COLUMNS)


# ════════════════════════════════ raw ══════════════════════════════════════

RAW_TOTAL_COLUMNS = [
    "label", "tool", "threshold", "status",
    "n_regions", "n_assemblies_called", "n_assemblies",
    "bp_called", "total_bp", "regions_per_mb", "frac_bp_called",
    "median_region_bp",
]
RAW_CLASS_COLUMNS = ["class", "label", "tool", "status", "n_regions",
                     "frac_regions"]
RAW_GENUS_COLUMNS = [
    "genus", "n_assemblies", "label", "tool", "status",
    "mean_regions_per_genome", "mean_regions_per_mb", "mean_frac_bp_called",
]
RAW_AGREEMENT_COLUMNS = [
    "reference", "query", "status", "n_reference_regions",
    "region_agreement", "bp_agreement",
]


@dataclass
class Series:
    """One tool at one score cutoff. `regions` is None for a placeholder."""
    label: str
    tool: str
    threshold: float
    regions: list[PredictedRegion] | None

    @property
    def status(self) -> str:
        return PLACEHOLDER if self.regions is None else MEASURED


def parse_thresholds(specs: list[str]) -> dict[str, list[float]]:
    """`["DeepBGC=0.5,0.8"]` -> {"DeepBGC": [0.5, 0.8]}."""
    out: dict[str, list[float]] = {}
    for spec in specs:
        name, sep, values = spec.partition("=")
        try:
            out[name] = sorted(float(v) for v in values.split(",") if v)
        except ValueError:
            out[name] = []
        if not sep or not name or not out[name]:
            raise SystemExit(f"--thresholds takes NAME=P[,P...], got {spec!r} "
                             "(e.g. DeepBGC=0.5,0.8)")
    return out


def series_label(tool: str, threshold: float, n_thresholds: int) -> str:
    """'DeepBGC' alone, or 'DeepBGC ≥0.8' once the cutoff is a choice."""
    if n_thresholds > 1 or threshold > 0:
        return f"{tool} ≥{threshold:g}"
    return tool


def make_series(tool: str, regions: list[PredictedRegion],
                thresholds: list[float]) -> list[Series]:
    """A tool's series. An unscored tool (one p_bgc value for every region,
    i.e. antiSMASH) gets one series: a cutoff cannot change it."""
    if len({r.p_bgc for r in regions}) <= 1:
        if len(thresholds) > 1:
            LOG.warning("%s has no score (every p_bgc identical); ignoring its "
                        "--thresholds", tool)
        return [Series(tool, tool, 0.0, regions)]
    return [
        Series(series_label(tool, t, len(thresholds)), tool, t,
               [r for r in regions if r.p_bgc >= t])
        for t in thresholds
    ]


def genus_of(organism: str) -> str:
    """Genus from an NCBI organism name.

    'Candidatus Nanopelagicus limnes' -> 'Candidatus Nanopelagicus'
    (Candidatus is a status, not a genus). Brackets are kept: NCBI writes
    '[Mycobacterium] stephanolepidis' precisely because the species is *not*
    in that genus any more, so folding it into Mycobacterium would be wrong.
    """
    words = organism.split()
    if not words:
        return "Unknown"
    if words[0] == "Candidatus" and len(words) > 1:
        return f"Candidatus {words[1]}"
    return words[0]


def load_genera(path: Path) -> dict[str, str]:
    """assemblies.tsv (build_genome_manifest.py) -> {assembly: genus}."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"assembly", "organism"} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} is not an assemblies.tsv "
                             f"(missing {sorted(missing)})")
        return {row["assembly"]: genus_of(row["organism"]) for row in reader}


def totals_rows(series: list[Series], index: dict[str, tuple[str, int]]):
    """raw_totals.tsv rows plus the per-assembly summaries behind them."""
    total_bp = sum(n for _, n in index.values())
    n_assemblies = len({a for a, _ in index.values()})
    rows, per_assembly = [], {}
    for s in series:
        if s.regions is None:
            rows.append({"label": s.label, "tool": s.tool, "threshold": s.threshold,
                         "status": PLACEHOLDER, "n_regions": 0,
                         "n_assemblies_called": 0, "n_assemblies": n_assemblies,
                         "bp_called": 0, "total_bp": total_bp,
                         "regions_per_mb": 0.0, "frac_bp_called": 0.0,
                         "median_region_bp": 0})
            continue
        # Regions are already cut at the threshold; summarize() must not cut again.
        overall, by_asm = summarize(s.label, s.regions, index, 0.0)
        per_assembly[s.label] = by_asm
        rows.append({
            "label": s.label, "tool": s.tool, "threshold": s.threshold,
            "status": MEASURED, "n_regions": overall.n_regions,
            "n_assemblies_called": overall.n_assemblies_called,
            "n_assemblies": n_assemblies,
            "bp_called": overall.bp_called, "total_bp": overall.total_bp,
            "regions_per_mb": overall.regions_per_mb,
            "frac_bp_called": overall.frac_bp_called,
            "median_region_bp": overall.median_region_bp,
        })
    return rows, per_assembly


def class_rows(series: list[Series]) -> list[dict]:
    rows = []
    for s in series:
        unknown: Counter = Counter()
        counts = Counter(normalize_class(r.predicted_class, unknown)
                         for r in (s.regions or []))
        _log_unknown(unknown, s.label)
        n = sum(counts.values())
        for cls in CLASS_ORDER:
            rows.append({"class": cls, "label": s.label, "tool": s.tool,
                         "status": s.status, "n_regions": counts.get(cls, 0),
                         "frac_regions": counts.get(cls, 0) / n if n else 0.0})
    return rows


def genus_rows(series: list[Series], per_assembly: dict, genera: dict[str, str],
               index: dict[str, tuple[str, int]], top_n: int) -> list[dict]:
    """Mean per-genome numbers for the `top_n` genera with the most assemblies;
    everything else pooled as 'Other genera'."""
    assemblies = sorted({a for a, _ in index.values()})
    genus = {a: genera.get(a, "Unknown") for a in assemblies}
    sizes = Counter(genus.values())
    top = [g for g, _ in sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))][:top_n]
    group = {a: (g if g in top else "Other genera") for a, g in genus.items()}
    order = top + (["Other genera"] if len(top) < len(sizes) else [])
    n_in = Counter(group.values())

    rows = []
    for s in series:
        acc: dict[str, list] = {g: [] for g in order}
        for a in per_assembly.get(s.label, []):
            acc[group[a.assembly]].append(a)
        for g in order:
            items = acc[g]
            rows.append({
                "genus": g, "n_assemblies": n_in[g], "label": s.label,
                "tool": s.tool, "status": s.status,
                "mean_regions_per_genome":
                    statistics.fmean(a.n_regions for a in items) if items else 0.0,
                "mean_regions_per_mb":
                    statistics.fmean(a.regions_per_mb for a in items) if items else 0.0,
                "mean_frac_bp_called":
                    statistics.fmean(a.frac_bp_called for a in items) if items else 0.0,
            })
    return rows


def agreement_rows(series: list[Series], scope: set[str]) -> list[dict]:
    """For every ordered pair: how much of the reference did the query call?

    The README's tool-vs-tool method, in memory: the reference's regions become
    the ground truth of an ordinary evaluate_predictions call, so
    `region_agreement` is detection.recall (reference regions the query covers
    by ≥ min_cluster_frac with a single region) and `bp_agreement` is
    nucleotide.recall (reference bp the query also called). Agreement with a
    tool, not with truth, and not symmetric.
    """
    rows = []
    for ref in series:
        ref_gt = (regions_to_clusters(ref.regions, prefix=ref.label)
                  if ref.regions is not None else None)
        for query in series:
            row = {"reference": ref.label, "query": query.label,
                   "n_reference_regions": len(ref.regions or [])}
            if ref_gt is None or query.regions is None:
                rows.append({**row, "status": PLACEHOLDER,
                             "region_agreement": None, "bp_agreement": None})
                continue
            if ref is query:
                rows.append({**row, "status": MEASURED,
                             "region_agreement": 1.0, "bp_agreement": 1.0})
                continue
            r = evaluate_predictions(query.regions, ref_gt, scope=scope,
                                     scope_source="explicit", max_listed_ids=1)
            LOG.info("%s covers %.1f%% of %s's regions and %.1f%% of its bp",
                     query.label, 100 * r.detection.recall, ref.label,
                     100 * r.nucleotide.recall)
            rows.append({**row, "status": MEASURED,
                         "region_agreement": r.detection.recall,
                         "bp_agreement": r.nucleotide.recall})
    return rows


def run_raw(
    manifest: Path, assemblies: Path | None, specs: list[str],
    placeholders: list[str], thresholds: dict[str, list[float]],
    top_genera: int, output_dir: Path,
) -> None:
    index = load_contig_index(manifest)
    LOG.info("manifest: %d contig(s) across %d assemblies, %.2f Gb",
             len(index), len({a for a, _ in index.values()}),
             sum(n for _, n in index.values()) / 1e9)

    tools, pending = load_tools(specs, placeholders)
    unused = set(thresholds) - {n for n, _ in tools}
    if unused:
        LOG.warning("--thresholds for tool(s) with no predictions: %s",
                    ", ".join(sorted(unused)))
    series: list[Series] = []
    for name, regions in tools:
        joined = check_contigs(name, regions, index)
        series += make_series(name, joined, thresholds.get(name, [0.0]))
    series += [Series(name, name, 0.0, None) for name in pending]

    totals, per_assembly = totals_rows(series, index)
    write_tsv(output_dir / "raw_totals.tsv", totals, RAW_TOTAL_COLUMNS)
    write_tsv(output_dir / "raw_by_class.tsv", class_rows(series), RAW_CLASS_COLUMNS)
    if assemblies is None:
        LOG.warning("no --assemblies given; raw_by_genus.tsv not written")
    else:
        write_tsv(output_dir / "raw_by_genus.tsv",
                  genus_rows(series, per_assembly, load_genera(assemblies),
                             index, top_genera),
                  RAW_GENUS_COLUMNS)
    write_tsv(output_dir / "raw_agreement.tsv",
              agreement_rows(series, set(index)), RAW_AGREEMENT_COLUMNS)


# ════════════════════════════════ CLI ══════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--predictions", nargs="+", default=[], metavar="NAME=PATH",
                        help="labelled predictions parquets; NAME is the label "
                             "drawn in the figures, e.g. antiSMASH=path.parquet")
        sp.add_argument("--placeholder", nargs="*", default=["SHARP"],
                        metavar="NAME",
                        help="tools with no results yet: written as all-zero rows "
                             "(default: SHARP; pass with no names to disable)")
        sp.add_argument("--output-dir", type=Path, required=True,
                        help="where the TSVs go")
        sp.add_argument("-v", "--verbose", action="store_true")

    m = sub.add_parser("mibig", help="tools against MiBIG ground truth")
    m.add_argument("--ground-truth", type=Path, required=True,
                   help="the scope's benchmark_ground_truth.tsv")
    m.add_argument("--contigs", type=Path, required=True,
                   help="the scope's analyzed_contigs.txt (the recall denominator)")
    m.add_argument("--thresholds", nargs="+", type=float, default=DEFAULT_SWEEP,
                   metavar="P", help="p_bgc cutoffs for the sweep table "
                                     "(default: 0.50 to 0.95 in steps of 0.05)")
    m.add_argument("--min-cluster-frac", type=float, default=0.5,
                   help="as sharp.evaluate (default 0.5)")
    m.add_argument("--min-prediction-frac", type=float, default=0.0,
                   help="as sharp.evaluate (default 0.0)")
    m.add_argument("--reciprocal-frac", type=float, default=0.5,
                   help="as sharp.evaluate (default 0.5)")
    m.add_argument("--min-p-bgc", type=float, default=0.0,
                   help="as sharp.evaluate, for the headline and per-class "
                        "tables (default 0.0)")
    common(m)

    r = sub.add_parser("raw", help="tools against each other, no ground truth")
    r.add_argument("--manifest", type=Path, required=True,
                   help="genomes.tsv from build_genome_manifest.py")
    r.add_argument("--assemblies", type=Path, default=None,
                   help="assemblies.tsv from build_genome_manifest.py "
                        "(for the per-genus table)")
    r.add_argument("--thresholds", nargs="+", default=[], metavar="NAME=P[,P...]",
                   help="score cutoffs per tool, e.g. DeepBGC=0.5,0.8 "
                        "(default: everything the tool reported)")
    r.add_argument("--top-genera", type=int, default=10,
                   help="genera shown individually in raw_by_genus.tsv (default 10)")
    common(r)
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    if not args.predictions and not args.placeholder:
        raise SystemExit("nothing to tabulate: give --predictions or --placeholder")
    if args.command == "mibig":
        run_mibig(
            ground_truth=args.ground_truth, contigs=args.contigs,
            specs=args.predictions, placeholders=args.placeholder,
            thresholds=sorted(set(args.thresholds)), output_dir=args.output_dir,
            s=MibigSettings(
                criterion=MatchCriterion(args.min_cluster_frac,
                                         args.min_prediction_frac),
                reciprocal_frac=args.reciprocal_frac,
                min_p_bgc=args.min_p_bgc,
            ),
        )
    else:
        run_raw(
            manifest=args.manifest, assemblies=args.assemblies,
            specs=args.predictions, placeholders=args.placeholder,
            thresholds=parse_thresholds(args.thresholds),
            top_genera=args.top_genera, output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
