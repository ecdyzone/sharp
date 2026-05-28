#!/usr/bin/env python3
"""Generate synthetic protein records for smoke-testing the pipeline.

Writes a FASTA file conforming to the same schema as the real
neighborhood_proteins.faa (region_id in each header), so downstream
steps consume it without knowing it's synthetic.

Usage:
    python scripts/generate_mock_data.py --n 100
    python scripts/generate_mock_data.py --n 50 --output data/mock/small.faa
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

from sharp.config import MOCK_DIR
from sharp.io import ProteinRecord, write_fasta

LOG = logging.getLogger("generate_mock_data")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def generate_mock_proteins(n: int, seed: int = 42) -> list[ProteinRecord]:
    """Synthetic proteins. Lengths ~N(300, 100) clipped to [50, 800].
    Proteins are distributed across ~n/10 regions."""
    rng = np.random.default_rng(seed)
    lengths = np.clip(rng.normal(300, 100, n), 50, 800).astype(int)
    n_regions = max(1, n // 10)
    region_ids = [f"R{i:03d}" for i in rng.integers(0, n_regions, n)]
    aa = np.array(list(AMINO_ACIDS))
    return [
        ProteinRecord(
            protein_id=f"mock_{i:05d}",
            region_id=rid,
            sequence="".join(rng.choice(aa, L).tolist()),
        )
        for i, (L, rid) in enumerate(zip(lengths, region_ids))
    ]


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--n", type=int, default=50,
                   help="number of proteins (default: %(default)s)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", type=Path,
                   default=MOCK_DIR / "neighborhood_proteins.faa",
                   help="output FASTA (default: %(default)s)")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    records = generate_mock_proteins(args.n, args.seed)
    n = write_fasta(args.output, records)
    LOG.info("wrote %d mock proteins to %s", n, args.output)


if __name__ == "__main__":
    main()
