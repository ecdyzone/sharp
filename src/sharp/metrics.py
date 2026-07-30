"""Metrics for the benchmark step.

Pure functions over PredictedRegion / KnownCluster collections.
No I/O, no logging — just math.

Half-open interval convention throughout: [start, end). end - start = length.

Three things this module insists on, each learned from a real run of
antiSMASH / DeepBGC / GECCO against MiBIG-derived ground truth:

1. **Scope.** Recall is measured only over ground-truth clusters that sit on a
   contig the tool was actually run on. Ground truth spans a database; a run
   spans one assembly. Dividing by the whole database makes recall meaningless
   (a real run scored 1/430 = 0.002 when the only evaluable cluster was found).

2. **Detection and boundary accuracy are different questions.** A symmetric
   reciprocal-overlap rule answers both at once and reports "not found" when
   only the bounds were loose. DeepBGC covered 100% of a true cluster inside a
   region 4.9x its length and scored recall 0.000. So the match rule is
   asymmetric by default: `min_cluster_frac` decides *found*,
   `min_prediction_frac` decides *tightly bounded*, and they are reported
   separately.

3. **Unmatched is not false.** Ground truth (MiBIG) is deliberately incomplete —
   ~53% of Streptomyces entries are dropped for missing coordinates, so absence
   from the ground truth carries no information. Nothing here is named
   "precision" or "false positive" at the region level; the honest quantity is
   `matched_prediction_frac`, a lower bound. Nucleotide-level `precision` keeps
   its name because there it is a plain bp ratio, not a claim of correctness.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median
from typing import Collection, Iterable, Sequence

from sharp.io import KnownCluster, PredictedRegion

Interval = tuple[int, int]


# ────────────────────────────── match rule ─────────────────────────────────

@dataclass(frozen=True)
class MatchCriterion:
    """How a (prediction, cluster) pair is judged to match.

    min_cluster_frac    — fraction of the CLUSTER that must be covered.
                          "Did the tool find this BGC?"
    min_prediction_frac — fraction of the PREDICTION that must be covered.
                          "Did it bound the BGC tightly, or call half a genome?"

    Setting both to the same value reproduces the symmetric reciprocal-overlap
    rule; the default leaves the prediction side off so detection is measured
    on its own.
    """
    min_cluster_frac: float = 0.5
    min_prediction_frac: float = 0.0


# ────────────────────────────── result types ───────────────────────────────

@dataclass(frozen=True)
class ScopeInfo:
    """What was actually evaluated, and what was excluded."""
    source: str                        # "explicit" | "inferred"
    n_contigs: int
    n_clusters_in_scope: int
    n_clusters_total: int
    n_predictions_in_scope: int
    n_predictions_total: int
    n_predictions_below_threshold: int
    min_p_bgc: float


@dataclass(frozen=True)
class MatchMetrics:
    """Region-level outcome under one match criterion.

    `matched_prediction_frac` is deliberately NOT called precision: an unmatched
    prediction may simply be absent from an incomplete ground truth. Treat it as
    a lower bound on precision.
    """
    recall: float
    n_recovered: int
    n_clusters: int
    matched_prediction_frac: float
    n_matched_predictions: int
    n_predictions: int
    f1: float


@dataclass(frozen=True)
class NucleotideMetrics:
    """Base-pair level agreement, computed over merged interval unions.

    Robust to boundary disagreement and to a cluster being split across several
    predicted regions, so it complements the region-level counts rather than
    restating them.
    """
    recall: float
    precision: float
    jaccard: float
    gt_bp: int
    predicted_bp: int
    intersect_bp: int


@dataclass(frozen=True)
class BoundaryMetrics:
    """How well-bounded the calls are — the axis detection recall ignores."""
    median_prediction_coverage: float       # ov / prediction length, matched preds
    median_cluster_coverage: float          # ov / cluster length, recovered clusters
    n_clusters_recovered_by_union_only: int  # split across ≥2 predictions
    n_merged_predictions: int                # one prediction matching ≥2 clusters


@dataclass(frozen=True)
class BenchmarkResult:
    """Output of evaluate_predictions(). Serializable to JSON via asdict()."""
    scope: ScopeInfo
    criterion: MatchCriterion
    reciprocal_frac: float
    detection: MatchMetrics
    reciprocal: MatchMetrics
    nucleotide: NucleotideMetrics
    boundary: BoundaryMetrics
    recovered_cluster_ids: list[str] = field(default_factory=list)
    missed_cluster_ids: list[str] = field(default_factory=list)
    matched_prediction_ids: list[str] = field(default_factory=list)
    unmatched_prediction_ids: list[str] = field(default_factory=list)
    ids_truncated: bool = False


# ────────────────────────────── overlap math ───────────────────────────────

def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Length of [a_start, a_end) ∩ [b_start, b_end). Zero if disjoint or
    if either interval is degenerate (start >= end)."""
    if a_start >= a_end or b_start >= b_end:
        return 0
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    """Sort and merge overlapping or abutting half-open intervals.

    [0, 100) and [100, 200) abut with no gap, so they merge into [0, 200).
    Degenerate intervals (start >= end) are dropped.
    """
    valid = sorted((s, e) for s, e in intervals if s < e)
    merged: list[Interval] = []
    for start, end in valid:
        if merged and start <= merged[-1][1]:
            last_start, last_end = merged[-1]
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def covered_bp(target: Interval, merged: Sequence[Interval]) -> int:
    """bp of `target` covered by `merged` — which must already be merged and
    sorted (see merge_intervals), so overlaps are never double-counted."""
    t_start, t_end = target
    return sum(overlap_bp(t_start, t_end, s, e) for s, e in merged)


def matches(p: PredictedRegion, c: KnownCluster, crit: MatchCriterion) -> bool:
    """True iff p and c are on the same contig and their overlap clears both
    sides of `crit`."""
    if p.contig != c.contig:
        return False
    ov = overlap_bp(p.start, p.end, c.start, c.end)
    if ov == 0:
        return False
    p_len = p.end - p.start
    c_len = c.end - c.start
    if p_len <= 0 or c_len <= 0:
        return False
    return (ov / c_len >= crit.min_cluster_frac
            and ov / p_len >= crit.min_prediction_frac)


def reciprocal_overlap(
    p: PredictedRegion, c: KnownCluster, min_frac: float
) -> bool:
    """True iff p and c overlap by ≥ min_frac of EACH interval's length.
    Symmetric in p and c — the strict rule, kept for labelling in train.py."""
    return matches(p, c, MatchCriterion(min_frac, min_frac))


# ────────────────────────────── aggregate ──────────────────────────────────

def _bucket_by_contig(items: Iterable) -> dict[str, list]:
    """Group predictions or clusters by contig.

    Matching only ever compares items sharing a contig, so bucketing first turns
    an O(P*C) sweep into O(P + C + pairs-per-contig). At BGC Atlas scale (204k
    ground-truth regions) the naive double loop is ~1e9 pair tests; MiBIG's
    busiest contig holds 15 clusters, so bucketed it is effectively linear.
    """
    out: dict[str, list] = {}
    for item in items:
        out.setdefault(item.contig, []).append(item)
    return out


def _match_metrics(
    predictions: list[PredictedRegion],
    ground_truth: list[KnownCluster],
    pred_by_contig: dict[str, list[PredictedRegion]],
    crit: MatchCriterion,
) -> tuple[MatchMetrics, set[str], set[str]]:
    """Region-level counts under one criterion.

    Returns the metrics plus the recovered-cluster and matched-prediction id
    sets, so callers can report the ids without matching twice.
    """
    recovered: set[str] = set()
    matched_preds: set[str] = set()

    for c in ground_truth:
        for p in pred_by_contig.get(c.contig, ()):
            if matches(p, c, crit):
                recovered.add(c.cluster_id)
                matched_preds.add(p.region_id)

    n_clusters = len(ground_truth)
    n_predictions = len(predictions)
    recall = len(recovered) / n_clusters if n_clusters else 0.0
    matched_frac = len(matched_preds) / n_predictions if n_predictions else 0.0
    f1 = (2 * matched_frac * recall / (matched_frac + recall)
          if (matched_frac + recall) else 0.0)

    return (
        MatchMetrics(
            recall=recall,
            n_recovered=len(recovered),
            n_clusters=n_clusters,
            matched_prediction_frac=matched_frac,
            n_matched_predictions=len(matched_preds),
            n_predictions=n_predictions,
            f1=f1,
        ),
        recovered,
        matched_preds,
    )


def evaluate_predictions(
    predictions: list[PredictedRegion],
    ground_truth: list[KnownCluster],
    *,
    scope: Collection[str] | None = None,
    scope_source: str = "inferred",
    criterion: MatchCriterion | None = None,
    reciprocal_frac: float = 0.5,
    min_p_bgc: float = 0.0,
    max_listed_ids: int = 1000,
) -> BenchmarkResult:
    """Score predictions against known clusters.

    scope           — contigs the tool was actually run on. None infers them from
                      the predictions themselves, which is optimistic: a contig
                      analyzed but not called on drops out of the denominator.
    criterion       — asymmetric match rule (see MatchCriterion).
    reciprocal_frac — threshold for the strict symmetric rule, reported alongside
                      so the two can be compared directly.
    min_p_bgc       — drop predictions scoring below this before anything else.
    max_listed_ids  — cap on each id list in the result; 0 means unlimited.

    A cluster is recovered if a SINGLE prediction clears the criterion — finding
    a BGC means finding it as a unit. Clusters covered only by several
    predictions together are counted in
    `boundary.n_clusters_recovered_by_union_only` rather than in recall.
    """
    crit = criterion if criterion is not None else MatchCriterion()

    n_predictions_total = len(predictions)
    n_clusters_total = len(ground_truth)

    scored = [p for p in predictions if p.p_bgc >= min_p_bgc]
    n_below_threshold = n_predictions_total - len(scored)

    resolved_scope = {p.contig for p in scored} if scope is None else set(scope)

    preds = [p for p in scored if p.contig in resolved_scope]
    clusters = [c for c in ground_truth if c.contig in resolved_scope]

    pred_by_contig = _bucket_by_contig(preds)

    detection, recovered, matched_preds = _match_metrics(
        preds, clusters, pred_by_contig, crit
    )
    reciprocal, _, _ = _match_metrics(
        preds, clusters, pred_by_contig,
        MatchCriterion(reciprocal_frac, reciprocal_frac),
    )

    nucleotide = _nucleotide_metrics(preds, clusters, pred_by_contig)
    boundary = _boundary_metrics(
        preds, clusters, pred_by_contig, crit, recovered, matched_preds
    )

    all_cluster_ids = {c.cluster_id for c in clusters}
    all_pred_ids = {p.region_id for p in preds}
    missed = sorted(all_cluster_ids - recovered)
    unmatched = sorted(all_pred_ids - matched_preds)

    lists = [sorted(recovered), missed, sorted(matched_preds), unmatched]
    truncated = bool(max_listed_ids) and any(len(x) > max_listed_ids for x in lists)
    if max_listed_ids:
        lists = [x[:max_listed_ids] for x in lists]

    return BenchmarkResult(
        scope=ScopeInfo(
            source=scope_source,
            n_contigs=len(resolved_scope),
            n_clusters_in_scope=len(clusters),
            n_clusters_total=n_clusters_total,
            n_predictions_in_scope=len(preds),
            n_predictions_total=n_predictions_total,
            n_predictions_below_threshold=n_below_threshold,
            min_p_bgc=min_p_bgc,
        ),
        criterion=crit,
        reciprocal_frac=reciprocal_frac,
        detection=detection,
        reciprocal=reciprocal,
        nucleotide=nucleotide,
        boundary=boundary,
        recovered_cluster_ids=lists[0],
        missed_cluster_ids=lists[1],
        matched_prediction_ids=lists[2],
        unmatched_prediction_ids=lists[3],
        ids_truncated=truncated,
    )


def _nucleotide_metrics(
    predictions: list[PredictedRegion],
    ground_truth: list[KnownCluster],
    pred_by_contig: dict[str, list[PredictedRegion]],
) -> NucleotideMetrics:
    """bp-level agreement between the union of predictions and the union of
    ground-truth clusters, restricted to contigs already in scope."""
    merged_preds = {
        contig: merge_intervals((p.start, p.end) for p in items)
        for contig, items in pred_by_contig.items()
    }
    gt_by_contig = _bucket_by_contig(ground_truth)
    merged_gt = {
        contig: merge_intervals((c.start, c.end) for c in items)
        for contig, items in gt_by_contig.items()
    }

    gt_bp = sum(e - s for spans in merged_gt.values() for s, e in spans)
    predicted_bp = sum(e - s for spans in merged_preds.values() for s, e in spans)
    intersect_bp = sum(
        covered_bp(span, merged_preds.get(contig, []))
        for contig, spans in merged_gt.items()
        for span in spans
    )
    union_bp = gt_bp + predicted_bp - intersect_bp

    return NucleotideMetrics(
        recall=intersect_bp / gt_bp if gt_bp else 0.0,
        precision=intersect_bp / predicted_bp if predicted_bp else 0.0,
        jaccard=intersect_bp / union_bp if union_bp else 0.0,
        gt_bp=gt_bp,
        predicted_bp=predicted_bp,
        intersect_bp=intersect_bp,
    )


def _boundary_metrics(
    predictions: list[PredictedRegion],
    ground_truth: list[KnownCluster],
    pred_by_contig: dict[str, list[PredictedRegion]],
    crit: MatchCriterion,
    recovered: set[str],
    matched_preds: set[str],
) -> BoundaryMetrics:
    """How tightly the matched calls bound their clusters, plus the split/merge
    diagnostics that explain a gap between detection and reciprocal recall."""
    gt_by_contig = _bucket_by_contig(ground_truth)

    # Tightness of each matched prediction: best ov / prediction length.
    pred_coverages: list[float] = []
    n_merged = 0
    for p in predictions:
        if p.region_id not in matched_preds:
            continue
        p_len = p.end - p.start
        overlaps = [overlap_bp(p.start, p.end, c.start, c.end)
                    for c in gt_by_contig.get(p.contig, ())]
        if p_len > 0 and overlaps:
            pred_coverages.append(max(overlaps) / p_len)
        n_matched_clusters = sum(
            1 for c in gt_by_contig.get(p.contig, ()) if matches(p, c, crit)
        )
        if n_matched_clusters > 1:
            n_merged += 1

    # How much of each recovered cluster was covered, and which clusters only
    # reach the threshold once several predictions are unioned together.
    cluster_coverages: list[float] = []
    n_union_only = 0
    for c in ground_truth:
        c_len = c.end - c.start
        if c_len <= 0:
            continue
        contig_preds = pred_by_contig.get(c.contig, [])
        if c.cluster_id in recovered:
            best = max(
                (overlap_bp(p.start, p.end, c.start, c.end) for p in contig_preds),
                default=0,
            )
            cluster_coverages.append(best / c_len)
        else:
            union = merge_intervals((p.start, p.end) for p in contig_preds)
            if covered_bp((c.start, c.end), union) / c_len >= crit.min_cluster_frac:
                n_union_only += 1

    return BoundaryMetrics(
        median_prediction_coverage=median(pred_coverages) if pred_coverages else 0.0,
        median_cluster_coverage=median(cluster_coverages) if cluster_coverages else 0.0,
        n_clusters_recovered_by_union_only=n_union_only,
        n_merged_predictions=n_merged,
    )
