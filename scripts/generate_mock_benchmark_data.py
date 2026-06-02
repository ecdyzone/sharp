#!/usr/bin/env python3
"""Generate synthetic predictions + ground truth for smoke-testing the
benchmark step end-to-end.

The two outputs are correlated by construction:
  - N clusters are placed at random positions across a few mock contigs
  - Of those, `recall_rate` fraction get a high-confidence prediction
    that overlaps them well (true positives)
  - The rest are "missed" (no overlapping prediction)
  - An additional `n_false_positives` random predictions are added that
    do NOT overlap any cluster (false positives)

This means you can predict the resulting metrics before running evaluate.py
— useful for verifying the pipeline is wired up correctly.

Usage:
    python scripts/generate_mock_benchmark_data.py
    python scripts/generate_mock_benchmark_data.py --n-clusters 30 --recall-rate 0.7
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

from sharp.config import MOCK_DIR
from sharp.io import (
    KnownCluster,
    PredictedRegion,
    write_ground_truth_tsv,
    write_predictions_parquet,
)

LOG = logging.getLogger("generate_mock_benchmark_data")

CONTIG_LENGTH = 5_000_000   # 5 Mb per mock contig
CLUSTER_LEN_MEAN = 30_000   # ~30 kb is a typical BGC size
CLUSTER_LEN_STD  = 10_000


def _place_intervals(rng: np.random.Generator, n: int, contigs: list[str]) -> list[tuple[str, int, int]]:
    """Place n non-overlapping intervals across the given contigs."""
    out: list[tuple[str, int, int]] = []
    per_contig: dict[str, list[tuple[int, int]]] = {c: [] for c in contigs}
    attempts = 0
    while len(out) < n and attempts < n * 50:
        attempts += 1
        contig = contigs[rng.integers(0, len(contigs))]
        length = int(np.clip(rng.normal(CLUSTER_LEN_MEAN, CLUSTER_LEN_STD), 5_000, 80_000))
        start = int(rng.integers(0, CONTIG_LENGTH - length))
        end = start + length
        # reject if overlaps an existing one on this contig
        if any(not (end <= s or start >= e) for s, e in per_contig[contig]):
            continue
        per_contig[contig].append((start, end))
        out.append((contig, start, end))
    if len(out) < n:
        LOG.warning("only placed %d / %d non-overlapping intervals", len(out), n)
    return out


def generate_mock_benchmark_data(
    n_clusters: int = 20,
    recall_rate: float = 0.7,
    n_false_positives: int = 5,
    n_contigs: int = 2,
    seed: int = 42,
) -> tuple[list[KnownCluster], list[PredictedRegion]]:
    """Build correlated ground truth and predictions.

    Returns (clusters, predictions) such that:
      - len(clusters) == n_clusters
      - ~recall_rate * n_clusters predictions overlap a cluster well (TPs)
      - n_false_positives predictions overlap no cluster (FPs)

    Note: the exact recall depends on how many cluster placements succeed
    and may be slightly below `recall_rate` for crowded inputs.
    """
    rng = np.random.default_rng(seed)
    contigs = [f"chr{i+1}" for i in range(n_contigs)]

    # 1. Place ground-truth clusters
    cluster_intervals = _place_intervals(rng, n_clusters, contigs)
    clusters = [
        KnownCluster(
            cluster_id=f"BGC{i+1:04d}",
            contig=ct,
            start=s,
            end=e,
            cluster_class=rng.choice(["T1PKS", "T2PKS", "NRPS", "terpene", "RiPP"]),
        )
        for i, (ct, s, e) in enumerate(cluster_intervals)
    ]

    # 2. True-positive predictions: shift each chosen cluster's coords slightly
    #    so the prediction isn't an exact copy, but still passes ≥50% reciprocal.
    n_tp = int(round(recall_rate * len(clusters)))
    tp_indices = rng.choice(len(clusters), size=n_tp, replace=False)
    predictions: list[PredictedRegion] = []
    for j, idx in enumerate(tp_indices):
        c = clusters[idx]
        length = c.end - c.start
        shift = int(rng.integers(-length // 5, length // 5 + 1))  # up to ±20%
        predictions.append(PredictedRegion(
            region_id=f"R{j+1:04d}",
            contig=c.contig,
            start=max(0, c.start + shift),
            end=c.end + shift,
            p_bgc=float(rng.uniform(0.7, 0.99)),
            predicted_class=c.cluster_class,
        ))

    # 3. False-positive predictions: place in regions not overlapping any cluster
    fp_intervals = _place_intervals(rng, n_false_positives, contigs)
    # Filter any that happen to overlap a real cluster
    cluster_by_contig: dict[str, list[tuple[int, int]]] = {ct: [] for ct in contigs}
    for c in clusters:
        cluster_by_contig.setdefault(c.contig, []).append((c.start, c.end))
    fp_predictions: list[PredictedRegion] = []
    for ct, s, e in fp_intervals:
        if any(not (e <= cs or s >= ce) for cs, ce in cluster_by_contig.get(ct, [])):
            continue
        fp_predictions.append(PredictedRegion(
            region_id=f"R{len(predictions) + len(fp_predictions) + 1:04d}",
            contig=ct,
            start=s,
            end=e,
            p_bgc=float(rng.uniform(0.5, 0.7)),
            predicted_class=None,
        ))
    predictions.extend(fp_predictions)

    return clusters, predictions


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--n-clusters", type=int, default=20)
    p.add_argument("--recall-rate", type=float, default=0.7,
                   help="fraction of clusters that get a matching prediction")
    p.add_argument("--n-false-positives", type=int, default=5)
    p.add_argument("--n-contigs", type=int, default=2)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--ground-truth-out", type=Path,
                   default=MOCK_DIR / "ground_truth.tsv")
    p.add_argument("--predictions-out", type=Path,
                   default=MOCK_DIR / "predictions.parquet")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )

    clusters, predictions = generate_mock_benchmark_data(
        n_clusters=args.n_clusters,
        recall_rate=args.recall_rate,
        n_false_positives=args.n_false_positives,
        n_contigs=args.n_contigs,
        seed=args.seed,
    )

    n_c = write_ground_truth_tsv(args.ground_truth_out, clusters)
    n_p = write_predictions_parquet(args.predictions_out, predictions)
    LOG.info("wrote %d clusters → %s", n_c, args.ground_truth_out)
    LOG.info("wrote %d predictions → %s", n_p, args.predictions_out)
    LOG.info(
        "expected metrics: recall ≈ %.2f, precision ≈ %.2f",
        args.recall_rate,
        (n_p - args.n_false_positives) / n_p if n_p else 0.0,
    )


if __name__ == "__main__":
    main()
