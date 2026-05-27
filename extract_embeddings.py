#!/usr/bin/env python3
"""
extract_embeddings.py — ESM-2 embeddings per protein for S(H)ARP.

Reads:  neighborhood_proteins.faa  (FASTA; headers carry `region_id=<id>`)
Writes: embeddings.parquet         (columns: protein_id, region_id, embedding)

MVP scope:
  - single small ESM-2 model (8M params by default; registry supports up to 650M)
  - length-bucketed batching to minimize padding waste on CPU
  - mean-pool over residue tokens, masking pad / CLS / EOS
  - streaming parquet write (constant memory regardless of input size)
  - CPU / MPS / CUDA auto-detection
  - synthetic mock data for end-to-end smoke testing without real input

Usage:
    # Real data
    python extract_embeddings.py \\
        --input neighborhood_proteins.faa \\
        --output embeddings.parquet \\
        --model esm2_t6_8M_UR50D

    # Smoke test with synthetic proteins (no input file needed)
    python extract_embeddings.py --mock 50 --output mock_embeddings.parquet

FASTA header format expected:
    >PROTEIN_ID region_id=R001 [...other tokens ignored]
    MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ...

Requires Python 3.10+, torch, transformers, pyarrow, numpy.
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import torch
from transformers import AutoModel, AutoTokenizer


# ────────────────────────────── constants ───────────────────────────────────

# Short name → (HuggingFace hub id, hidden dimension).
# Bigger models give richer embeddings but cost more memory / time.
MODEL_REGISTRY: dict[str, tuple[str, int]] = {
    "esm2_t6_8M_UR50D":    ("facebook/esm2_t6_8M_UR50D",    320),
    "esm2_t12_35M_UR50D":  ("facebook/esm2_t12_35M_UR50D",  480),
    "esm2_t30_150M_UR50D": ("facebook/esm2_t30_150M_UR50D", 640),
    "esm2_t33_650M_UR50D": ("facebook/esm2_t33_650M_UR50D", 1280),
}

DEFAULT_MODEL = "esm2_t6_8M_UR50D"
REGION_ID_RE = re.compile(r"region_id=([^\s]+)")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
LOG = logging.getLogger("extract_embeddings")


# ────────────────────────────── config ──────────────────────────────────────

@dataclass(frozen=True)
class Config:
    input_path: Path | None
    output_path: Path
    model_name: str = DEFAULT_MODEL
    batch_size: int = 8
    max_length: int = 1024
    device: str = "auto"
    mock_n: int | None = None
    mock_seed: int = 42
    log_every: int = 50


# ────────────────────────────── data model ──────────────────────────────────

@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    region_id: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def parse_fasta(path: Path) -> Iterator[ProteinRecord]:
    """Stream-parse FASTA. Records missing region_id or with empty sequence
    are skipped with a warning."""
    pid: str | None = None
    region: str | None = None
    seq_chunks: list[str] = []

    def maybe_emit() -> Iterator[ProteinRecord]:
        if pid is None:
            return
        seq = "".join(seq_chunks)
        if not seq:
            LOG.warning("skip %s — empty sequence", pid)
        elif region is None:
            LOG.warning("skip %s — no region_id in header", pid)
        else:
            yield ProteinRecord(pid, region, seq)

    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                yield from maybe_emit()
                header = line[1:]
                pid = header.split(None, 1)[0]
                m = REGION_ID_RE.search(header)
                region = m.group(1) if m else None
                seq_chunks = []
            else:
                seq_chunks.append(line)
    yield from maybe_emit()


def generate_mock_records(n: int, seed: int) -> list[ProteinRecord]:
    """Synthetic proteins for smoke testing. Lengths ~N(300, 100) clipped
    to [50, 800]. Proteins assigned to ~n/10 mock regions."""
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


# ────────────────────────────── embedder ────────────────────────────────────

class Embedder:
    """Wraps an ESM-2 model. Produces one fixed-dim vector per protein by
    mean-pooling per-residue embeddings, with CLS / EOS / pad tokens masked
    out of the mean."""

    def __init__(self, model_name: str, device: torch.device, max_length: int):
        hub_id, self.dim = MODEL_REGISTRY[model_name]
        LOG.info("loading %s (%d-dim) on %s", hub_id, self.dim, device)
        self.tokenizer = AutoTokenizer.from_pretrained(hub_id)
        self.model = AutoModel.from_pretrained(hub_id).to(device).eval()
        self.device = device
        self.max_length = max_length

    @torch.inference_mode()
    def embed_batch(self, sequences: list[str]) -> np.ndarray:
        tok = self.tokenizer(
            sequences,
            padding=True,
            truncation=True,
            max_length=self.max_length,
            return_tensors="pt",
        ).to(self.device)

        # (B, L, D) — last hidden state for every token position
        hidden = self.model(**tok).last_hidden_state

        # Build a residue-only mask. ESM-2's attention_mask covers CLS, residues
        # and EOS; we want only the residues for the protein-level mean.
        mask = tok.attention_mask.clone()                              # (B, L)
        mask[:, 0] = 0                                                  # zero CLS
        last_idx = tok.attention_mask.sum(dim=1) - 1                   # EOS index per row
        mask[torch.arange(mask.size(0), device=self.device), last_idx] = 0

        m = mask.unsqueeze(-1).to(hidden.dtype)                         # (B, L, 1)
        pooled = (hidden * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)    # (B, D)
        return pooled.float().cpu().numpy()


# ────────────────────────────── batching ────────────────────────────────────

def length_bucketed_batches(
    records: list[ProteinRecord], batch_size: int
) -> Iterator[list[ProteinRecord]]:
    """Sort by length, then chunk contiguously. Padding waste within a batch
    is bounded by the length spread of `batch_size` adjacent proteins."""
    by_length = sorted(records, key=lambda r: r.length)
    for i in range(0, len(by_length), batch_size):
        yield by_length[i : i + batch_size]


# ────────────────────────────── writer ──────────────────────────────────────

def write_parquet_streaming(
    output_path: Path,
    batches: Iterable[tuple[list[ProteinRecord], np.ndarray]],
    embedding_dim: int,
) -> int:
    """Write embeddings batch by batch. Uses fixed-size lists for the embedding
    column so downstream tools can lay it out as a 2-D matrix efficiently."""
    schema = pa.schema([
        ("protein_id", pa.string()),
        ("region_id",  pa.string()),
        ("embedding",  pa.list_(pa.float32(), embedding_dim)),
    ])

    n_written = 0
    with pq.ParquetWriter(output_path, schema, compression="zstd") as writer:
        for records, vecs in batches:
            assert vecs.shape == (len(records), embedding_dim), \
                f"unexpected vec shape {vecs.shape}, expected ({len(records)}, {embedding_dim})"
            # Build the fixed-size-list array from a flat float32 buffer —
            # avoids materializing Python-side nested lists.
            flat = pa.array(vecs.reshape(-1).astype(np.float32, copy=False),
                            type=pa.float32())
            emb_arr = pa.FixedSizeListArray.from_arrays(flat, embedding_dim)
            table = pa.Table.from_arrays(
                [
                    pa.array([r.protein_id for r in records], type=pa.string()),
                    pa.array([r.region_id  for r in records], type=pa.string()),
                    emb_arr,
                ],
                schema=schema,
            )
            writer.write_table(table)
            n_written += len(records)
    return n_written


# ────────────────────────────── device ──────────────────────────────────────

def select_device(spec: str) -> torch.device:
    if spec != "auto":
        return torch.device(spec)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# ────────────────────────────── pipeline ────────────────────────────────────

def load_records(cfg: Config) -> list[ProteinRecord]:
    if cfg.mock_n is not None:
        LOG.info("generating %d mock proteins (seed=%d)", cfg.mock_n, cfg.mock_seed)
        return generate_mock_records(cfg.mock_n, cfg.mock_seed)
    assert cfg.input_path is not None, "either input_path or mock_n must be set"
    LOG.info("reading %s", cfg.input_path)
    return list(parse_fasta(cfg.input_path))


def report_input_stats(records: list[ProteinRecord], max_length: int) -> None:
    lengths = np.fromiter((r.length for r in records), dtype=np.int32, count=len(records))
    LOG.info("loaded %d proteins  min=%d  median=%d  max=%d",
             len(records), lengths.min(), int(np.median(lengths)), lengths.max())
    n_trunc = int((lengths > max_length).sum())
    if n_trunc:
        LOG.warning("%d / %d proteins exceed max_length=%d and will be truncated",
                    n_trunc, len(records), max_length)


def run(cfg: Config) -> None:
    records = load_records(cfg)
    if not records:
        LOG.error("no records to process")
        sys.exit(1)
    report_input_stats(records, cfg.max_length)

    device = select_device(cfg.device)
    embedder = Embedder(cfg.model_name, device, cfg.max_length)

    def iter_embeddings() -> Iterator[tuple[list[ProteinRecord], np.ndarray]]:
        t0 = time.time()
        done = 0
        for batch in length_bucketed_batches(records, cfg.batch_size):
            vecs = embedder.embed_batch([r.sequence for r in batch])
            done += len(batch)
            # Log roughly every `log_every` proteins.
            if (done // cfg.log_every) > ((done - len(batch)) // cfg.log_every):
                rate = done / max(time.time() - t0, 1e-6)
                LOG.info("  %d / %d  (%.1f prot/s)", done, len(records), rate)
            yield batch, vecs

    cfg.output_path.parent.mkdir(parents=True, exist_ok=True)
    n = write_parquet_streaming(cfg.output_path, iter_embeddings(), embedder.dim)

    # Verify by reading back metadata (cheap; doesn't load the embeddings).
    meta = pq.read_metadata(cfg.output_path)
    size_kb = cfg.output_path.stat().st_size / 1024
    LOG.info("wrote %d rows × %d cols → %s  (%.1f KB on disk)",
             meta.num_rows, meta.num_columns, cfg.output_path, size_kb)
    assert meta.num_rows == n, f"row mismatch: wrote {n}, parquet has {meta.num_rows}"


# ────────────────────────────── cli ─────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", type=Path, help="neighborhood_proteins.faa (FASTA)")
    p.add_argument("--output", type=Path, required=True, help="embeddings.parquet")
    p.add_argument("--model", default=DEFAULT_MODEL, choices=list(MODEL_REGISTRY),
                   help=f"ESM-2 variant (default: {DEFAULT_MODEL})")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--max-length", type=int, default=1024,
                   help="truncate proteins longer than this (residues)")
    p.add_argument("--device", default="auto", help="auto | cpu | cuda | mps")
    p.add_argument("--mock", type=int, default=None, metavar="N",
                   help="generate N synthetic proteins instead of reading --input")
    p.add_argument("--mock-seed", type=int, default=42)
    p.add_argument("--log-every", type=int, default=50)
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    if args.mock is None and args.input is None:
        LOG.error("either --input or --mock N must be provided")
        sys.exit(2)
    run(Config(
        input_path=args.input,
        output_path=args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        max_length=args.max_length,
        device=args.device,
        mock_n=args.mock,
        mock_seed=args.mock_seed,
        log_every=args.log_every,
    ))


if __name__ == "__main__":
    main()
