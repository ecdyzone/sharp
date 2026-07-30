"""Benchmark step: compare predictions against a known ground truth.

Reads:  predictions.parquet (from predict.py, or from a baseline converter)
        ground_truth.tsv    (curated from MiBIG, etc.)
        contigs.txt         (optional — the contigs the tool was run on)
Writes: benchmark.json      (scope, detection, reciprocal, nucleotide, boundary)

Delegates:
  - I/O                  → sharp.io
  - Metric math          → sharp.metrics
  - Configuration        → sharp.config

Two things to know before reading the numbers:

  * **Scope.** Recall counts only ground-truth clusters on contigs that were
    actually analyzed. Pass `--contigs` with the assembly you ran the tool on.
    Without it the scope is inferred from the predictions, which is optimistic —
    a contig analyzed but not called on silently leaves the denominator — so
    every tool in a comparison should be given the same `--contigs` file.

  * **Unmatched is not false.** Ground truth is incomplete by construction, so
    predictions absent from it are unvalidated rather than wrong. The output
    reports `matched_prediction_frac` (a lower bound on precision), never a
    region-level "precision" or "false positive" count.

Usage (after `pip install -e .`):
    python -m sharp.evaluate \\
        --predictions data/interim/predictions.parquet \\
        --ground-truth data/raw/mibig_ground_truth.tsv \\
        --output data/processed/benchmark.json \\
        --contigs data/raw/analyzed_contigs.txt \\
        --min-cluster-frac 0.5
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from sharp.config import EvaluateConfig
from sharp.io import (
    load_contigs,
    load_ground_truth_tsv,
    load_predictions_parquet,
    write_benchmark_json,
)
from sharp.metrics import BenchmarkResult, MatchCriterion, evaluate_predictions

LOG = logging.getLogger("evaluate")

_SAMPLE = 5   # contig names shown per side when scope and ground truth disjoin


# ────────────────────────────── orchestration ──────────────────────────────

def _log_summary(result: BenchmarkResult) -> None:
    """Report the honest headline, including why detection and reciprocal
    disagree when they do."""
    s, d, r = result.scope, result.detection, result.reciprocal

    origin = "explicit" if s.source == "explicit" else "INFERRED from predictions"
    LOG.info(
        "scope: %d contig(s) (%s) — %d of %d ground-truth clusters evaluable",
        s.n_contigs, origin, s.n_clusters_in_scope, s.n_clusters_total,
    )
    if s.n_predictions_below_threshold:
        LOG.info(
            "       %d prediction(s) dropped below p_bgc %.2f",
            s.n_predictions_below_threshold, s.min_p_bgc,
        )

    LOG.info(
        "detection : recall=%.3f (%d/%d)   matched %d/%d predictions",
        d.recall, d.n_recovered, d.n_clusters,
        d.n_matched_predictions, d.n_predictions,
    )
    LOG.info("reciprocal: recall=%.3f @%.2f", r.recall, result.reciprocal_frac)
    if d.n_recovered > r.n_recovered:
        LOG.info(
            "            ← %d cluster(s) found, but bounds too loose to pass "
            "the symmetric match", d.n_recovered - r.n_recovered,
        )

    n = result.nucleotide
    LOG.info(
        "nucleotide: recall=%.3f  precision=%.3f  jaccard=%.3f",
        n.recall, n.precision, n.jaccard,
    )
    b = result.boundary
    LOG.info(
        "boundary  : median prediction coverage %.3f, median cluster coverage %.3f",
        b.median_prediction_coverage, b.median_cluster_coverage,
    )
    if b.n_clusters_recovered_by_union_only:
        LOG.info(
            "            %d cluster(s) covered only by several predictions "
            "together (not counted in recall)",
            b.n_clusters_recovered_by_union_only,
        )
    if b.n_merged_predictions:
        LOG.info(
            "            %d prediction(s) span more than one cluster",
            b.n_merged_predictions,
        )
    LOG.info(
        "note: %d unmatched prediction(s) are unvalidated, not false — the "
        "ground truth is incomplete", len(result.unmatched_prediction_ids),
    )


def _report_contig_mismatch(predictions, ground_truth, scope) -> None:
    """No cluster fell in scope. Usually the contig names simply don't match
    between the two files (assembly-qualified vs. bare accession), so show both
    sides rather than reporting recall over nothing."""
    LOG.error(
        "no ground-truth clusters on any contig in scope — nothing to evaluate"
    )
    gt_contigs = sorted({c.contig for c in ground_truth})
    pred_contigs = sorted({p.contig for p in predictions})
    LOG.error(
        "  scope (%d contig(s)):      %s",
        len(scope), ", ".join(sorted(scope)[:_SAMPLE]) or "(empty)",
    )
    LOG.error(
        "  prediction contigs (%d):   %s",
        len(pred_contigs), ", ".join(pred_contigs[:_SAMPLE]) or "(none)",
    )
    LOG.error(
        "  ground-truth contigs (%d): %s",
        len(gt_contigs), ", ".join(gt_contigs[:_SAMPLE]) or "(none)",
    )
    LOG.error(
        "  if these look like the same contigs under different names, the "
        "ground truth and the predictions disagree on naming — fix that before "
        "reading any benchmark number"
    )


def run(cfg: EvaluateConfig) -> None:
    LOG.info("reading %s", cfg.predictions_path)
    predictions = load_predictions_parquet(cfg.predictions_path)

    LOG.info("reading %s", cfg.ground_truth_path)
    ground_truth = load_ground_truth_tsv(cfg.ground_truth_path)

    if not ground_truth:
        LOG.error("ground truth is empty — nothing to evaluate against")
        sys.exit(1)

    # Scope: what the tool was actually run on. Explicit beats inferred, and an
    # inferred scope is worth a warning because it flatters a silent tool.
    if cfg.contigs_path is not None:
        LOG.info("reading %s", cfg.contigs_path)
        scope = load_contigs(cfg.contigs_path)
        scope_source = "explicit"
        if not scope:
            LOG.error("contigs file %s is empty", cfg.contigs_path)
            sys.exit(1)
    else:
        scope = {p.contig for p in predictions}
        scope_source = "inferred"
        LOG.warning(
            "no --contigs given: scope inferred from the %d contig(s) present in "
            "the predictions. Recall may be optimistic — a contig that was "
            "analyzed but produced no prediction drops out of the denominator. "
            "Pass --contigs (the same file for every tool) for a fair comparison.",
            len(scope),
        )

    if not any(c.contig in scope for c in ground_truth):
        _report_contig_mismatch(predictions, ground_truth, scope)
        sys.exit(1)

    result = evaluate_predictions(
        predictions,
        ground_truth,
        scope=scope,
        scope_source=scope_source,
        criterion=MatchCriterion(cfg.min_cluster_frac, cfg.min_prediction_frac),
        reciprocal_frac=cfg.reciprocal_frac,
        min_p_bgc=cfg.min_p_bgc,
        max_listed_ids=cfg.max_listed_ids,
    )

    _log_summary(result)

    write_benchmark_json(cfg.output_path, result)
    LOG.info("wrote %s", cfg.output_path)


# ────────────────────────────── cli ────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--predictions", type=Path, required=True,
                   help="predictions.parquet from predict.py or a converter")
    p.add_argument("--ground-truth", type=Path, required=True,
                   help="TSV with columns: cluster_id, contig, start, end[, class]")
    p.add_argument("--output", type=Path, required=True,
                   help="output benchmark.json")
    p.add_argument("--contigs", type=Path, default=None,
                   help="contigs the tool was run on (one per line, or a .fai). "
                        "Omitted: inferred from the predictions, with a warning")
    p.add_argument("--min-cluster-frac", type=float, default=0.5,
                   help="fraction of a cluster that must be covered to count it "
                        "as found (default: %(default)s)")
    p.add_argument("--min-prediction-frac", type=float, default=0.0,
                   help="fraction of a prediction that must be covered for it to "
                        "count as tightly bounded; 0 keeps detection and boundary "
                        "accuracy separate (default: %(default)s)")
    p.add_argument("--reciprocal-frac", type=float, default=0.5,
                   help="threshold for the strict symmetric rule reported "
                        "alongside (default: %(default)s)")
    p.add_argument("--min-p-bgc", type=float, default=0.0,
                   help="drop predictions scoring below this before scoring "
                        "(default: %(default)s)")
    p.add_argument("--max-listed-ids", type=int, default=1000,
                   help="cap on each id list in the JSON; 0 = unlimited "
                        "(default: %(default)s)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    if not 0.0 < args.min_cluster_frac <= 1.0:
        LOG.error("--min-cluster-frac must be in (0, 1], got %s", args.min_cluster_frac)
        sys.exit(2)
    if not 0.0 <= args.min_prediction_frac <= 1.0:
        LOG.error("--min-prediction-frac must be in [0, 1], got %s",
                  args.min_prediction_frac)
        sys.exit(2)
    if not 0.0 < args.reciprocal_frac <= 1.0:
        LOG.error("--reciprocal-frac must be in (0, 1], got %s", args.reciprocal_frac)
        sys.exit(2)
    if args.max_listed_ids < 0:
        LOG.error("--max-listed-ids must be >= 0, got %s", args.max_listed_ids)
        sys.exit(2)

    run(EvaluateConfig(
        predictions_path=args.predictions,
        ground_truth_path=args.ground_truth,
        output_path=args.output,
        contigs_path=args.contigs,
        min_cluster_frac=args.min_cluster_frac,
        min_prediction_frac=args.min_prediction_frac,
        reciprocal_frac=args.reciprocal_frac,
        min_p_bgc=args.min_p_bgc,
        max_listed_ids=args.max_listed_ids,
    ))


if __name__ == "__main__":
    main()
