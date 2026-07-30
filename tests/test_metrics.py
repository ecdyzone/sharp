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
    MatchCriterion,
    covered_bp,
    evaluate_predictions,
    matches,
    merge_intervals,
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


# ────────────────────────────── interval union ─────────────────────────────

class TestMergeIntervals:
    def test_empty(self) -> None:
        assert merge_intervals([]) == []

    def test_disjoint_preserved_and_sorted(self) -> None:
        assert merge_intervals([(200, 300), (0, 100)]) == [(0, 100), (200, 300)]

    def test_overlapping_merged(self) -> None:
        assert merge_intervals([(0, 100), (50, 150)]) == [(0, 150)]

    def test_abutting_merged(self) -> None:
        # Half-open [0,100) and [100,200) cover [0,200) with no gap.
        assert merge_intervals([(0, 100), (100, 200)]) == [(0, 200)]

    def test_nested_absorbed(self) -> None:
        assert merge_intervals([(0, 500), (100, 200)]) == [(0, 500)]

    def test_degenerate_dropped(self) -> None:
        assert merge_intervals([(50, 50), (0, 100)]) == [(0, 100)]

    def test_chain(self) -> None:
        assert merge_intervals([(0, 10), (5, 20), (18, 30), (100, 110)]) == [
            (0, 30), (100, 110),
        ]


class TestCoveredBp:
    def test_no_coverage(self) -> None:
        assert covered_bp((0, 100), [(200, 300)]) == 0

    def test_full_coverage(self) -> None:
        assert covered_bp((0, 100), [(0, 100)]) == 100

    def test_partial_from_two_spans(self) -> None:
        # Two disjoint predictions each covering 25 bp of the target.
        assert covered_bp((0, 100), [(0, 25), (75, 100)]) == 50

    def test_no_double_counting_when_merged(self) -> None:
        # merge_intervals is the contract; merged input can't double-count.
        merged = merge_intervals([(0, 60), (40, 100)])
        assert covered_bp((0, 100), merged) == 100


# ────────────────────────────── match rules ────────────────────────────────

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

    def test_exact_threshold_boundary(self) -> None:
        # Both sides exactly 50% — should pass at min_frac=0.5 (≥, not >).
        p = pred("R1", "chr1", 0, 200)
        c = clu("C1", "chr1", 100, 300)        # overlap = 100
        assert reciprocal_overlap(p, c, 0.5)
        assert not reciprocal_overlap(p, c, 0.5001)

    def test_zero_length_intervals(self) -> None:
        assert not reciprocal_overlap(pred("R1", "chr1", 50, 50),
                                      clu("C1", "chr1", 0, 100), 0.5)
        assert not reciprocal_overlap(pred("R1", "chr1", 0, 100),
                                      clu("C1", "chr1", 50, 50), 0.5)


class TestAsymmetricMatching:
    """The DeepBGC case: a prediction that fully contains a cluster but is far
    wider than it. Detection and boundary accuracy are separate questions."""

    # Prediction 10x the cluster: cluster 100% covered, prediction 10% covered.
    WIDE_PRED = pred("R1", "chr1", 0, 1000)
    SMALL_CLU = clu("C1", "chr1", 0, 100)

    def test_symmetric_rule_rejects_a_cluster_it_fully_contains(self) -> None:
        # The old behaviour, kept as the strict criterion.
        assert not reciprocal_overlap(self.WIDE_PRED, self.SMALL_CLU, 0.5)
        assert reciprocal_overlap(self.WIDE_PRED, self.SMALL_CLU, 0.1)

    def test_default_criterion_counts_it_as_found(self) -> None:
        # Regression test for the real DeepBGC result: a 94 kb region covering
        # 100% of a 19 kb cluster used to score recall 0.000.
        assert matches(self.WIDE_PRED, self.SMALL_CLU, MatchCriterion())

    def test_prediction_side_can_be_re_enabled(self) -> None:
        assert not matches(self.WIDE_PRED, self.SMALL_CLU,
                           MatchCriterion(min_cluster_frac=0.5,
                                          min_prediction_frac=0.5))

    def test_detection_recall_unaffected_by_prediction_width(self) -> None:
        result = evaluate_predictions([self.WIDE_PRED], [self.SMALL_CLU])
        assert result.detection.recall == 1.0
        assert result.reciprocal.recall == 0.0
        # Boundary metrics carry the width penalty instead.
        assert result.boundary.median_prediction_coverage == pytest.approx(0.1)
        assert result.boundary.median_cluster_coverage == pytest.approx(1.0)

    def test_different_contigs_never_match(self) -> None:
        assert not matches(pred("R1", "chr1", 0, 100),
                           clu("C1", "chr2", 0, 100), MatchCriterion())


# ────────────────────────────── scope ──────────────────────────────────────

class TestScope:
    """Recall's denominator must exclude contigs that were never analyzed."""

    def test_inferred_scope_excludes_unanalyzed_contigs(self) -> None:
        # The real shape: ground truth spans a database, the run spans one contig.
        preds = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100)] + [
            clu(f"C{i}", f"other{i}", 0, 100) for i in range(2, 50)
        ]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 1.0          # not 1/49
        assert result.scope.n_clusters_in_scope == 1
        assert result.scope.n_clusters_total == 49
        assert result.scope.source == "inferred"
        assert result.missed_cluster_ids == []

    def test_explicit_scope_keeps_uncalled_contigs_in_denominator(self) -> None:
        # chr2 was analyzed but produced nothing — that is a real miss.
        preds = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr2", 0, 100)]
        result = evaluate_predictions(
            preds, clusters, scope={"chr1", "chr2"}, scope_source="explicit"
        )
        assert result.detection.recall == 0.5
        assert result.scope.n_contigs == 2
        assert result.scope.source == "explicit"
        assert result.missed_cluster_ids == ["C2"]

    def test_inferred_scope_is_optimistic_relative_to_explicit(self) -> None:
        # Documents the warning evaluate.py emits: same inputs, higher recall.
        preds = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr2", 0, 100)]
        inferred = evaluate_predictions(preds, clusters)
        explicit = evaluate_predictions(preds, clusters, scope={"chr1", "chr2"})
        assert inferred.detection.recall > explicit.detection.recall

    def test_predictions_outside_scope_are_excluded(self) -> None:
        preds = [pred("R1", "chr1", 0, 100), pred("R2", "chr9", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, clusters, scope={"chr1"})
        assert result.scope.n_predictions_in_scope == 1
        assert result.scope.n_predictions_total == 2
        assert result.unmatched_prediction_ids == []

    def test_empty_scope_yields_zero_metrics(self) -> None:
        result = evaluate_predictions(
            [pred("R1", "chr1", 0, 100)], [clu("C1", "chr1", 0, 100)], scope=set()
        )
        assert result.detection.recall == 0.0
        assert result.scope.n_clusters_in_scope == 0


# ────────────────────────────── p_bgc filtering ────────────────────────────

class TestMinPBgc:
    def test_default_keeps_everything(self) -> None:
        preds = [pred("R1", "chr1", 0, 100, 0.1)]
        result = evaluate_predictions(preds, [clu("C1", "chr1", 0, 100)])
        assert result.scope.n_predictions_below_threshold == 0
        assert result.detection.recall == 1.0

    def test_low_scoring_predictions_dropped(self) -> None:
        preds = [pred("R1", "chr1", 0, 100, 0.9), pred("R2", "chr1", 500, 600, 0.2)]
        result = evaluate_predictions(
            preds, [clu("C1", "chr1", 0, 100)], scope={"chr1"}, min_p_bgc=0.5
        )
        assert result.scope.n_predictions_below_threshold == 1
        assert result.scope.n_predictions_in_scope == 1
        assert result.unmatched_prediction_ids == []

    def test_threshold_can_remove_the_only_match(self) -> None:
        preds = [pred("R1", "chr1", 0, 100, 0.4)]
        result = evaluate_predictions(
            preds, [clu("C1", "chr1", 0, 100)], scope={"chr1"}, min_p_bgc=0.5
        )
        assert result.detection.recall == 0.0

    def test_threshold_recorded(self) -> None:
        result = evaluate_predictions(
            [], [clu("C1", "chr1", 0, 100)], scope={"chr1"}, min_p_bgc=0.42
        )
        assert result.scope.min_p_bgc == 0.42


# ────────────────────────────── evaluate_predictions ───────────────────────

class TestEvaluatePredictions:
    def test_empty_both(self) -> None:
        result = evaluate_predictions([], [])
        assert result.detection.recall == 0.0
        assert result.detection.matched_prediction_frac == 0.0
        assert result.detection.f1 == 0.0
        assert result.scope.n_predictions_total == 0
        assert result.scope.n_clusters_total == 0

    def test_empty_predictions_with_clusters(self) -> None:
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr1", 200, 300)]
        result = evaluate_predictions([], clusters, scope={"chr1"})
        assert result.detection.recall == 0.0
        assert result.missed_cluster_ids == ["C1", "C2"]

    def test_empty_clusters_with_predictions(self) -> None:
        preds = [pred("R1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, [])
        assert result.detection.recall == 0.0
        assert result.unmatched_prediction_ids == ["R1"]

    def test_perfect_one_to_one_match(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100),  pred("R2", "chr1", 200, 300)]
        clusters = [clu("C1", "chr1", 0, 100),    clu("C2", "chr1", 200, 300)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 1.0
        assert result.detection.matched_prediction_frac == 1.0
        assert result.detection.f1 == 1.0
        assert sorted(result.recovered_cluster_ids) == ["C1", "C2"]
        assert result.missed_cluster_ids == []
        assert result.unmatched_prediction_ids == []

    def test_one_cluster_two_overlapping_predictions(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100), pred("R2", "chr1", 10, 110)]
        clusters = [clu("C1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 1.0
        assert result.detection.n_recovered == 1
        assert result.detection.n_matched_predictions == 2

    def test_half_recall_no_unmatched(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr1", 500, 600)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.matched_prediction_frac == 1.0
        assert result.detection.recall == 0.5
        assert result.detection.f1 == pytest.approx(2 / 3, abs=1e-9)
        assert result.recovered_cluster_ids == ["C1"]
        assert result.missed_cluster_ids == ["C2"]

    def test_unmatched_prediction_is_not_called_false(self) -> None:
        # Naming matters: the ground truth is incomplete, so an unmatched
        # prediction is unvalidated rather than wrong.
        preds = [pred("R1", "chr1", 0, 100), pred("R2", "chr1", 500, 600)]
        result = evaluate_predictions(preds, [clu("C1", "chr1", 0, 100)])
        assert result.unmatched_prediction_ids == ["R2"]
        assert not hasattr(result, "false_positive_prediction_ids")
        assert not hasattr(result.detection, "precision")

    def test_overlapping_different_contigs_no_match(self) -> None:
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr2", 0, 100)]
        result = evaluate_predictions(preds, clusters, scope={"chr1", "chr2"})
        assert result.detection.recall == 0.0

    def test_threshold_changes_result(self) -> None:
        # Cluster 30% covered by the prediction.
        preds    = [pred("R1", "chr1", 0, 100)]
        clusters = [clu("C1", "chr1", 70, 170)]
        strict = evaluate_predictions(preds, clusters,
                                      criterion=MatchCriterion(0.5))
        lenient = evaluate_predictions(preds, clusters,
                                       criterion=MatchCriterion(0.3))
        assert strict.detection.recall == 0.0
        assert lenient.detection.recall == 1.0

    def test_criterion_recorded_in_output(self) -> None:
        result = evaluate_predictions([], [], criterion=MatchCriterion(0.42, 0.1),
                                      reciprocal_frac=0.7)
        assert result.criterion.min_cluster_frac == 0.42
        assert result.criterion.min_prediction_frac == 0.1
        assert result.reciprocal_frac == 0.7

    def test_result_is_dataclass(self) -> None:
        assert isinstance(evaluate_predictions([], []), BenchmarkResult)

    def test_contig_bucketing_does_not_change_results(self) -> None:
        # Same intervals repeated on many contigs — bucketing must keep them apart.
        preds = [pred(f"R{i}", f"chr{i}", 0, 100) for i in range(20)]
        clusters = [clu(f"C{i}", f"chr{i}", 0, 100) for i in range(20)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 1.0
        assert result.detection.n_recovered == 20
        assert result.boundary.n_merged_predictions == 0


# ────────────────────────────── nucleotide metrics ─────────────────────────

class TestNucleotideMetrics:
    def test_exact_match(self) -> None:
        result = evaluate_predictions([pred("R1", "chr1", 0, 100)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.nucleotide.recall == 1.0
        assert result.nucleotide.precision == 1.0
        assert result.nucleotide.jaccard == 1.0

    def test_over_prediction_lowers_precision_not_recall(self) -> None:
        # Prediction 10x too wide: all of the cluster found, most bp spurious.
        result = evaluate_predictions([pred("R1", "chr1", 0, 1000)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.nucleotide.recall == 1.0
        assert result.nucleotide.precision == pytest.approx(0.1)
        assert result.nucleotide.gt_bp == 100
        assert result.nucleotide.predicted_bp == 1000
        assert result.nucleotide.intersect_bp == 100

    def test_partial_coverage(self) -> None:
        result = evaluate_predictions([pred("R1", "chr1", 0, 50)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.nucleotide.recall == pytest.approx(0.5)
        assert result.nucleotide.precision == 1.0
        assert result.nucleotide.jaccard == pytest.approx(0.5)

    def test_overlapping_predictions_not_double_counted(self) -> None:
        # Two predictions overlapping each other cover 150 bp between them.
        preds = [pred("R1", "chr1", 0, 100), pred("R2", "chr1", 50, 150)]
        result = evaluate_predictions(preds, [clu("C1", "chr1", 0, 150)])
        assert result.nucleotide.predicted_bp == 150
        assert result.nucleotide.recall == 1.0

    def test_no_overlap(self) -> None:
        result = evaluate_predictions([pred("R1", "chr1", 500, 600)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.nucleotide.recall == 0.0
        assert result.nucleotide.jaccard == 0.0


# ────────────────────────────── boundary diagnostics ───────────────────────

class TestBoundaryMetrics:
    def test_split_cluster_not_counted_in_recall(self) -> None:
        # Two adjacent predictions each cover 40% of the cluster — 80% together,
        # but neither alone clears 0.5. Finding a BGC means finding it as a unit,
        # so recall stays 0 and the split is reported separately.
        preds = [pred("R1", "chr1", 0, 40), pred("R2", "chr1", 40, 80)]
        clusters = [clu("C1", "chr1", 0, 100)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 0.0
        assert result.boundary.n_clusters_recovered_by_union_only == 1
        assert result.nucleotide.recall == pytest.approx(0.8)

    def test_merged_prediction_spanning_two_clusters(self) -> None:
        preds = [pred("R1", "chr1", 0, 1000)]
        clusters = [clu("C1", "chr1", 0, 100), clu("C2", "chr1", 500, 600)]
        result = evaluate_predictions(preds, clusters)
        assert result.detection.recall == 1.0
        assert result.boundary.n_merged_predictions == 1

    def test_tight_calls_score_high_prediction_coverage(self) -> None:
        result = evaluate_predictions([pred("R1", "chr1", 0, 100)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.boundary.median_prediction_coverage == pytest.approx(1.0)

    def test_no_matches_leaves_medians_at_zero(self) -> None:
        result = evaluate_predictions([pred("R1", "chr1", 500, 600)],
                                      [clu("C1", "chr1", 0, 100)])
        assert result.boundary.median_prediction_coverage == 0.0
        assert result.boundary.median_cluster_coverage == 0.0


# ────────────────────────────── id list capping ────────────────────────────

class TestIdListCapping:
    def test_lists_capped(self) -> None:
        clusters = [clu(f"C{i:04d}", "chr1", i * 1000, i * 1000 + 100)
                    for i in range(50)]
        result = evaluate_predictions([], clusters, scope={"chr1"},
                                      max_listed_ids=10)
        assert len(result.missed_cluster_ids) == 10
        assert result.ids_truncated is True

    def test_zero_means_unlimited(self) -> None:
        clusters = [clu(f"C{i:04d}", "chr1", i * 1000, i * 1000 + 100)
                    for i in range(50)]
        result = evaluate_predictions([], clusters, scope={"chr1"},
                                      max_listed_ids=0)
        assert len(result.missed_cluster_ids) == 50
        assert result.ids_truncated is False

    def test_not_truncated_when_under_cap(self) -> None:
        result = evaluate_predictions([], [clu("C1", "chr1", 0, 100)],
                                      scope={"chr1"}, max_listed_ids=10)
        assert result.ids_truncated is False

    def test_counts_are_unaffected_by_capping(self) -> None:
        clusters = [clu(f"C{i:04d}", "chr1", i * 1000, i * 1000 + 100)
                    for i in range(50)]
        result = evaluate_predictions([], clusters, scope={"chr1"},
                                      max_listed_ids=5)
        assert result.detection.n_clusters == 50
