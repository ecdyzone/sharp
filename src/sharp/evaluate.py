"""Benchmark step: compare predictions against a known ground truth.

Reads:  predictions.parquet (from predict.py)
        ground_truth.tsv    (curated from MiBIG, etc.)
Writes: benchmark.json      (precision, recall, F1 + per-item details)

Delegates:
  - I/O                  → sharp.io
  - Metric math          → sharp.metrics
  - Configuration        → sharp.config

Usage (after `pip install -e .`):
    python -m sharp.evaluate \\
        --predictions data/interim/predictions.parquet \\
        --ground-truth data/raw/mibig_ground_truth.tsv \\
        --output data/processed/benchmark.json \\
        --min-overlap-frac 0.5
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from sharp.config import EvaluateConfig
from sharp.io import (
    load_ground_truth_tsv,
    load_predictions_parquet,
    write_benchmark_json,
)
from sharp.metrics import evaluate_predictions

LOG = logging.getLogger("evaluate")


# ────────────────────────────── orchestration ──────────────────────────────

def run(cfg: EvaluateConfig) -> None:
    LOG.info("reading %s", cfg.predictions_path)
    predictions = load_predictions_parquet(cfg.predictions_path)

    LOG.info("reading %s", cfg.ground_truth_path)
    ground_truth = load_ground_truth_tsv(cfg.ground_truth_path)

    if not ground_truth:
        LOG.error("ground truth is empty — nothing to evaluate against")
        sys.exit(1)

    LOG.info(
        "evaluating %d predictions against %d known clusters (min_overlap_frac=%.2f)",
        len(predictions), len(ground_truth), cfg.min_overlap_frac,
    )
    result = evaluate_predictions(predictions, ground_truth, cfg.min_overlap_frac)

    LOG.info(
        "→ precision=%.3f  recall=%.3f  F1=%.3f   "
        "(TP=%d, recovered=%d/%d, FP=%d)",
        result.precision, result.recall, result.f1,
        result.n_true_positives, result.n_recovered, result.n_clusters,
        len(result.false_positive_prediction_ids),
    )

    write_benchmark_json(cfg.output_path, result)
    LOG.info("wrote %s", cfg.output_path)


# ────────────────────────────── cli ────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--predictions", type=Path, required=True,
                   help="predictions.parquet from predict.py")
    p.add_argument("--ground-truth", type=Path, required=True,
                   help="TSV with columns: cluster_id, contig, start, end[, class]")
    p.add_argument("--output", type=Path, required=True,
                   help="output benchmark.json")
    p.add_argument("--min-overlap-frac", type=float, default=0.5,
                   help="reciprocal overlap threshold for a match (default: %(default)s)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    if not 0.0 < args.min_overlap_frac <= 1.0:
        LOG.error("--min-overlap-frac must be in (0, 1], got %s", args.min_overlap_frac)
        sys.exit(2)
    run(EvaluateConfig(
        predictions_path=args.predictions,
        ground_truth_path=args.ground_truth,
        output_path=args.output,
        min_overlap_frac=args.min_overlap_frac,
    ))


if __name__ == "__main__":
    main()
