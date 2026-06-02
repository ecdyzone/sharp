"""Metrics for the benchmark step.

Pure functions over PredictedRegion / KnownCluster collections.
No I/O, no logging — just math. Owned here:

  - overlap_bp:           bp of intersection between two half-open intervals
  - reciprocal_overlap:   does a prediction "match" a known cluster?
  - evaluate_predictions: aggregate precision / recall / F1 + per-item details

Half-open interval convention throughout: [start, end). end - start = length.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from sharp.io import KnownCluster, PredictedRegion


@dataclass(frozen=True)
class BenchmarkResult:
    """Output of evaluate_predictions(). Serializable to JSON."""
    n_predictions: int
    n_clusters: int
    n_recovered: int
    n_true_positives: int
    recall: float
    precision: float
    f1: float
    min_overlap_frac: float
    recovered_cluster_ids: list[str]      = field(default_factory=list)
    missed_cluster_ids: list[str]         = field(default_factory=list)
    matched_prediction_ids: list[str]     = field(default_factory=list)
    false_positive_prediction_ids: list[str] = field(default_factory=list)


# ────────────────────────────── overlap math ───────────────────────────────

def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Length of [a_start, a_end) ∩ [b_start, b_end). Zero if disjoint or
    if either interval is degenerate (start >= end)."""
    if a_start >= a_end or b_start >= b_end:
        return 0
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def reciprocal_overlap(
    p: PredictedRegion, c: KnownCluster, min_frac: float
) -> bool:
    """True iff p and c are on the same contig and their overlap is ≥ min_frac
    of EACH interval's length. Symmetric in p and c."""
    if p.contig != c.contig:
        return False
    ov = overlap_bp(p.start, p.end, c.start, c.end)
    if ov == 0:
        return False
    p_len = p.end - p.start
    c_len = c.end - c.start
    if p_len <= 0 or c_len <= 0:
        return False
    return (ov / p_len >= min_frac) and (ov / c_len >= min_frac)


# ────────────────────────────── aggregate ──────────────────────────────────

def evaluate_predictions(
    predictions: list[PredictedRegion],
    ground_truth: list[KnownCluster],
    min_overlap_frac: float = 0.5,
) -> BenchmarkResult:
    """Compute precision, recall, F1 for a set of predictions against a set
    of known clusters.

    A prediction is a true positive iff it matches at least one known cluster
    (reciprocal_overlap ≥ min_overlap_frac). A known cluster is recovered iff
    it is matched by at least one prediction. Each cluster contributes to
    recall at most once; each prediction contributes to precision at most
    once — set semantics throughout.
    """
    recovered: set[str] = set()
    matched_preds: set[str] = set()

    # O(P*C). For MVP scale (hundreds of each), this is fine. Speed up
    # later with an interval index per contig if it ever matters.
    for c in ground_truth:
        for p in predictions:
            if reciprocal_overlap(p, c, min_overlap_frac):
                recovered.add(c.cluster_id)
                matched_preds.add(p.region_id)

    n_pred = len(predictions)
    n_clus = len(ground_truth)
    n_rec = len(recovered)
    n_tp  = len(matched_preds)

    recall    = n_rec / n_clus if n_clus else 0.0
    precision = n_tp / n_pred if n_pred else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) else 0.0)

    all_cluster_ids = [c.cluster_id for c in ground_truth]
    all_pred_ids    = [p.region_id  for p in predictions]

    return BenchmarkResult(
        n_predictions=n_pred,
        n_clusters=n_clus,
        n_recovered=n_rec,
        n_true_positives=n_tp,
        recall=recall,
        precision=precision,
        f1=f1,
        min_overlap_frac=min_overlap_frac,
        recovered_cluster_ids=sorted(recovered),
        missed_cluster_ids=sorted(set(all_cluster_ids) - recovered),
        matched_prediction_ids=sorted(matched_preds),
        false_positive_prediction_ids=sorted(set(all_pred_ids) - matched_preds),
    )
