#!/usr/bin/env python3
"""Turn one tool's predictions into a ground-truth table, for tool-vs-tool agreement.

`sharp.evaluate` scores predictions against a ground truth. With no MiBiG
labels there is none — but the machinery does not care what the reference
*means*, only that it is a set of intervals. Converting tool B's predictions
into ground-truth shape therefore lets `evaluate.py` answer "how much of what B
called did A also call?" with no new metric code.

**Read the resulting JSON carefully.** The field names still say recall and
ground truth, and they now mean something weaker:

    detection.recall      -> fraction of B's regions that A also called
                             (AGREEMENT with a tool, not with truth)
    scope.n_clusters      -> number of B regions in scope
    matched_prediction_frac -> fraction of A's regions that land on a B region

Neither direction is correctness. Run it both ways: A-vs-B and B-vs-A are
different numbers, and the interesting one for a new tool is usually the
asymmetry — what it calls that the established tool does not.

`--min-p-bgc` filters before conversion, because the reference's size is a
threshold choice exactly like the predictions' (see summarize_predictions.py).

Usage:
    # antiSMASH regions as the reference
    python scripts/predictions_to_ground_truth.py \\
        --input data/interim/antismash_predictions_actino.parquet \\
        --output data/interim/agreement/antismash_as_gt.tsv

    # then: how much of antiSMASH does DeepBGC agree with?
    python -m sharp.evaluate \\
        --predictions data/interim/deepbgc_predictions_actino.parquet \\
        --ground-truth data/interim/agreement/antismash_as_gt.tsv \\
        --contigs data/interim/actino_db/analyzed_contigs.txt \\
        --output data/processed/actino_comparison/agreement_deepbgc_vs_antismash.json
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from sharp.io import (
    KnownCluster,
    load_predictions_parquet,
    write_ground_truth_tsv,
)

LOG = logging.getLogger("predictions_to_ground_truth")


def to_clusters(
    path: Path, min_p_bgc: float = 0.0, prefix: str | None = None
) -> list[KnownCluster]:
    """Load predictions and re-label them as KnownClusters.

    Coordinates pass through untouched — both types are 0-based half-open, so
    there is nothing to convert and nothing to get wrong.
    """
    regions = [r for r in load_predictions_parquet(path) if r.p_bgc >= min_p_bgc]
    return [
        KnownCluster(
            # Prefixing keeps the origin visible in a merged table and stops
            # two tools' identically-named regions colliding.
            cluster_id=f"{prefix}:{r.region_id}" if prefix else r.region_id,
            contig=r.contig,
            start=r.start,
            end=r.end,
            cluster_class=r.predicted_class,
        )
        for r in regions
    ]


def run(
    input_path: Path, output_path: Path, min_p_bgc: float, prefix: str | None
) -> None:
    if not input_path.is_file():
        raise SystemExit(f"no such predictions file: {input_path}")
    clusters = to_clusters(input_path, min_p_bgc, prefix)
    if not clusters:
        raise SystemExit(
            f"{input_path} yielded no regions at --min-p-bgc {min_p_bgc} — "
            "an empty reference would make every agreement number meaningless"
        )
    n = write_ground_truth_tsv(output_path, clusters)
    LOG.info("wrote %d reference interval(s) → %s", n, output_path)
    LOG.warning("this is a TOOL, not truth: evaluate.py's 'recall' against it "
                "means agreement with that tool")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", type=Path, required=True,
                   help="predictions parquet to use as the reference")
    p.add_argument("--output", type=Path, required=True,
                   help="ground-truth-shaped TSV to write")
    p.add_argument("--min-p-bgc", type=float, default=0.0,
                   help="drop reference regions scoring below this (default 0)")
    p.add_argument("--prefix", default=None,
                   help="prepend '<prefix>:' to every cluster_id, e.g. the tool name")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    run(args.input, args.output, args.min_p_bgc, args.prefix)


if __name__ == "__main__":
    main()
