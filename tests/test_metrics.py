"""Tests for sharp.metrics — overlap math, matching, aggregation.

This is the file that has to be right. The metric definitions are the
methodological contract; if the math here is wrong, every benchmark
number downstream is meaningless.
"""
from __future__ import annotations

import pytest

from sharp.io import KnownCluster, PredictedRegion
from sharp.metrics import (
    BenchmarkResult,
    evaluate_predictions,
    overlap_bp,
    reciprocal_overlap,
)


# Helpers — short constructors keep the test bodies readable.
def pred(rid: str, contig: str, start: int, end: int, p: float = 0.9) -> PredictedRegion:
    return PredictedRegion(rid, contig, start, end, p)

def clu(cid: str, contig: str, start: int, end: int) -> KnownCluster:
    return KnownCluster(cid, contig, start, end)


# ────────────────────────────── overlap_bp ─────────────────────────────────

class TestOverlapBp:
    def test_disjoint_left(self) -> None:
        assert overlap_bp(0, 100, 200, 300) == 0

    def test_disjoint_right(self) -> None:
        assert overlap_bp(200, 300, 0, 100) == 0

    def test_touching_is_zero(self) -> None:
        # Half-open: [0, 100) and [100, 200) share no positions.
        assert overlap_bp(0, 100, 100, 200) == 0
        assert overlap_bp(100, 200, 0, 100) == 0

    def test_identical(self) -> None:
        assert overlap_bp(50, 150, 50, 150) == 100

    def test_full_containment(self) -> None:
        # B inside A
        assert overlap_bp(0, 200, 50, 150) == 100
        # A inside B (symmetric)
        assert overlap_bp(50, 150, 0, 200) == 100

    def test_partial_overlap(self) -> None:
        assert overlap_bp(0, 100, 50, 150) == 50
        assert overlap_bp(50, 150, 0, 100) == 50

    def test_one_bp_overlap(self) -> None:
        assert overlap_bp(0, 101, 100, 200) == 1

    def test_degenerate_zero_length(self) -> None:
        # A point interval has zero length → zero overlap.
        assert overlap_bp(50, 50, 0, 100) == 0
        assert overlap_bp(0, 100, 50, 50) == 0

    def test_inverted_interval_is_zero(self) -> None:
        # Defensive: a malformed [end, start) returns 0 instead of garbage.
        assert overlap_bp(100, 50, 0, 100) == 0
        assert overlap_bp(0, 100, 100, 50) == 0


# ────────────────────────────── reciprocal_overlap ─────────────────────────

class TestReciprocalOverlap:
    def test_different_contigs(self) -> None:
        p = pred("R1", "chr1", 0, 100)
        c = clu("C1", "chr2", 0, 100)
        assert not reciprocal_overlap(p, c, 0.5)

    def test_same_contig_disjoint(self) -> None:
        assert not reciprocal_overlap(pred("R1", "chr1", 0, 100),
                                      clu("C1", "chr1", 200, 300), 0.5)

    def test_identical_intervals(self) -> None:
        # 100% reciprocal — passes any threshold.
        assert reciprocal_overlap(pred("R1", "chr1", 0, 100),
                                  clu("C1", "chr1", 0, 100), 0.5)
        assert reciprocal_overlap(pred("R1", "chr1", 0, 100),
                                  clu("C1", "chr1", 0, 100), 1.0)

    def test_below_threshold_symmetric(self) -> None:
        # 30% overlap on each side — fails 0.5 threshold.
        p = pred("R1", "chr1", 0, 100)
        c = clu("C1", "chr1", 70, 170)         # overlap = 30
        assert not reciprocal_overlap(p, c, 0.5)
        # Lowering the threshold lets it through.
        assert reciprocal_overlap(p, c, 0.3)
        assert reciprocal_overlap(p, c, 0.29)

    def test_above_threshold_one_side_only(self) -> None:
        # Prediction is 10x larger than cluster. Cluster is 100% covered,
        # but prediction is only 10% covered → reciprocal fails 0.5.
        p = pred("R1", "chr1", 0, 1000)
        c = clu("C1", "chr1", 0, 100)
        assert not reciprocal_overlap(p, c, 0.5)
        assert reciprocal_overlap(p, c, 0.1)   # generous threshold lets it pass

    def test_exact_threshold_boundary(self) -> None:
        # Both sides exactly 50% — should pass at min_frac=0.5 (≥, not >).
        p = pred("R1", "chr1", 0, 200)
        c = clu("C1", "chr1", 100, 300)        # overlap = 100
        # p_len=200, c_len=200, ov=100 → 0.5 on both sides
        assert reciprocal_overlap(p, c, 0.5)
        assert not reciprocal_overlap(p, c, 0.5001)

    def test_zero_length_intervals(self) -> None:
        assert not reciprocal_overlap(pred("R1", "chr1", 50, 50),
                                      clu("C1", "chr1", 0, 100), 0.5)
        assert not reciprocal_overlap(pred("R1", "chr1", 0, 100),
                                      clu("C1", "chr1", 50, 50), 0.5)


# ────────────────────────────── evaluate_predictions ───────────────────────

class TestEvaluatePredictions:
    def test_empty_both(self) -> None:
        result = evaluate_predictions([], [], min_overlap_frac=0.5)
        assert result.precision == 0.0
        assert result.recall == 0.0
        assert result.f1 == 0.0
        assert result.n_predictions == 0
        assert result.n_clusters == 0

    def test_empty_predictions_with_clusters(self) -> None:
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr1", 200, 300)]
        result = evaluate_predictions([], clusters, min_overlap_frac=0.5)
        assert result.recall == 0.0
        assert result.precision == 0.0          # 0/0 → 0 by convention here
        assert result.missed_cluster_ids == ["C1", "C2"]

    def test_empty_clusters_with_predictions(self) -> None:
        preds = [pred("R1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, [], min_overlap_frac=0.5)
        assert result.recall == 0.0
        assert result.precision == 0.0          # no clusters → nothing can be a TP
        assert result.false_positive_prediction_ids == ["R1"]

    def test_perfect_one_to_one_match(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100),  pred("R2", "chr1", 200, 300)]
        clusters = [clu("C1", "chr1", 0, 100),    clu("C2", "chr1", 200, 300)]
        result = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        assert result.precision == 1.0
        assert result.recall == 1.0
        assert result.f1 == 1.0
        assert sorted(result.recovered_cluster_ids) == ["C1", "C2"]
        assert result.missed_cluster_ids == []
        assert result.false_positive_prediction_ids == []

    def test_one_cluster_two_overlapping_predictions(self) -> None:
        # Two predictions both well-overlap the same cluster.
        # Recovered should be {C1} (1 unique cluster), TPs should be {R1, R2}.
        # → recall = 1/1, precision = 2/2 = 1.0
        preds    = [pred("R1", "chr1", 0, 100), pred("R2", "chr1", 10, 110)]
        clusters = [clu("C1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        assert result.recall == 1.0
        assert result.precision == 1.0
        assert result.n_recovered == 1
        assert result.n_true_positives == 2

    def test_half_recall_no_false_positives(self) -> None:
        # 1 of 2 clusters recovered, 0 false positives.
        # → precision=1.0, recall=0.5, F1 = 2*1*0.5/1.5 = 2/3
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr1", 500, 600)]
        result = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        assert result.precision == 1.0
        assert result.recall == 0.5
        assert result.f1 == pytest.approx(2 / 3, abs=1e-9)
        assert result.recovered_cluster_ids == ["C1"]
        assert result.missed_cluster_ids == ["C2"]

    def test_one_true_positive_one_false_positive(self) -> None:
        # R1 matches C1; R2 matches nothing.
        # → recall = 1/1, precision = 1/2, F1 = 2/3
        preds = [
            pred("R1", "chr1", 0, 100),
            pred("R2", "chr1", 500, 600),
        ]
        clusters = [clu("C1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        assert result.precision == 0.5
        assert result.recall == 1.0
        assert result.f1 == pytest.approx(2 / 3, abs=1e-9)
        assert result.false_positive_prediction_ids == ["R2"]

    def test_overlapping_different_contigs_no_match(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr2", 0, 100)]
        result = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        assert result.recall == 0.0
        assert result.precision == 0.0

    def test_threshold_changes_result(self) -> None:
        # 30% reciprocal overlap.
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 70, 170)]
        strict = evaluate_predictions(preds, clusters, min_overlap_frac=0.5)
        lenient = evaluate_predictions(preds, clusters, min_overlap_frac=0.3)
        assert strict.recall == 0.0
        assert lenient.recall == 1.0

    def test_threshold_recorded_in_output(self) -> None:
        result = evaluate_predictions([], [], min_overlap_frac=0.42)
        assert result.min_overlap_frac == 0.42

    def test_result_is_dataclass(self) -> None:
        result = evaluate_predictions([], [], min_overlap_frac=0.5)
        assert isinstance(result, BenchmarkResult)
