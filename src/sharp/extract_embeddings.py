"""ESM-2 embeddings per protein for S(H)ARP.

Reads:  neighborhood_proteins.faa  (FASTA; headers carry region_id=<id>)
Writes: embeddings.parquet         (protein_id, region_id, embedding[D])

This module orchestrates the embedding extraction step only. It delegates:
  - I/O                  → sharp.io
  - Model lifecycle      → sharp.model_management
  - Configuration        → sharp.config

Usage (after `pip install -e .`):
    python -m sharp.extract_embeddings \\
        --input data/interim/neighborhood_proteins.faa \\
        --output data/interim/embeddings.parquet
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Iterator

import numpy as np
import pyarrow.parquet as pq

from sharp.config import EmbeddingConfig
from sharp.io import ProteinRecord, parse_fasta, write_embeddings_parquet
from sharp.model_management import (
    MODEL_REGISTRY,
    Embedder,
    ensure_model_available,
    select_device,
)

LOG = logging.getLogger("extract_embeddings")


# ────────────────────────────── batching ───────────────────────────────────

def length_bucketed_batches(
    records: list[ProteinRecord], batch_size: int
) -> Iterator[list[ProteinRecord]]:
    """Sort proteins by length, chunk contiguously. Padding waste within
    each batch is bounded by the spread of `batch_size` adjacent proteins."""
    by_length = sorted(records, key=lambda r: r.length)
    for i in range(0, len(by_length), batch_size):
        yield by_length[i : i + batch_size]


# ────────────────────────────── orchestration ──────────────────────────────

def _report_input_stats(records: list[ProteinRecord], max_length: int) -> None:
    lengths = np.fromiter(
        (r.length for r in records), dtype=np.int32, count=len(records)
    )
    LOG.info(
        "loaded %d proteins  min=%d  median=%d  max=%d",
        len(records), lengths.min(), int(np.median(lengths)), lengths.max(),
    )
    n_trunc = int((lengths > max_length).sum())
    if n_trunc:
        LOG.warning(
            "%d / %d proteins exceed max_length=%d and will be truncated",
            n_trunc, len(records), max_length,
        )


def run(cfg: EmbeddingConfig) -> None:
    LOG.info("reading %s", cfg.input_path)
    records = list(parse_fasta(cfg.input_path))
    if not records:
        LOG.error("no records to process")
        sys.exit(1)
    _report_input_stats(records, cfg.max_length)

    # Ensure model is local before doing anything compute-intensive.
    ensure_model_available(cfg.model_name)

    device = select_device(cfg.device)
    embedder = Embedder(cfg.model_name, device, cfg.max_length)

    def iter_embeddings() -> Iterator[tuple[list[ProteinRecord], np.ndarray]]:
        t0 = time.time()
        done = 0
        for batch in length_bucketed_batches(records, cfg.batch_size):
            vecs = embedder.embed_batch([r.sequence for r in batch])
            done += len(batch)
            if (done // cfg.log_every) > ((done - len(batch)) // cfg.log_every):
                rate = done / max(time.time() - t0, 1e-6)
                LOG.info("  %d / %d  (%.1f prot/s)", done, len(records), rate)
            yield batch, vecs

    n = write_embeddings_parquet(cfg.output_path, iter_embeddings(), embedder.dim)

    meta = pq.read_metadata(cfg.output_path)
    size_kb = cfg.output_path.stat().st_size / 1024
    LOG.info(
        "wrote %d rows × %d cols → %s  (%.1f KB on disk)",
        meta.num_rows, meta.num_columns, cfg.output_path, size_kb,
    )
    assert meta.num_rows == n


# ────────────────────────────── cli ────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", type=Path, required=True,
                   help="neighborhood_proteins.faa (FASTA with region_id in header)")
    p.add_argument("--output", type=Path, required=True,
                   help="output embeddings.parquet")
    p.add_argument("--model", default="esm2_t6_8M_UR50D",
                   choices=list(MODEL_REGISTRY),
                   help="ESM-2 variant (default: %(default)s)")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--max-length", type=int, default=1024,
                   help="truncate proteins longer than this (residues)")
    p.add_argument("--device", default="auto", help="auto | cpu | cuda | mps")
    p.add_argument("--log-every", type=int, default=50)
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    run(EmbeddingConfig(
        input_path=args.input,
        output_path=args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        max_length=args.max_length,
        device=args.device,
        log_every=args.log_every,
    ))


if __name__ == "__main__":
    main()
