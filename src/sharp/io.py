"""I/O for the pipeline.

Owns the `ProteinRecord` data type that flows between steps, and the
FASTA / parquet read+write functions. Everything that touches disk for
the embedding step routes through here.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

LOG = logging.getLogger(__name__)
REGION_ID_RE = re.compile(r"region_id=([^\s]+)")


@dataclass(frozen=True)
class ProteinRecord:
    """A protein with its identifier and the BGC-candidate region it belongs to."""
    protein_id: str
    region_id: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


# ────────────────────────────── FASTA ──────────────────────────────────────

def parse_fasta(path: Path) -> Iterator[ProteinRecord]:
    """Stream-parse a FASTA file into ProteinRecord objects.

    Records missing `region_id=` in their header or with empty sequence
    are skipped with a warning.

    Header format expected:
        >PROTEIN_ID region_id=R001 [other tokens are ignored]
    """
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


def write_fasta(
    path: Path, records: Iterable[ProteinRecord], line_width: int = 60
) -> int:
    """Write protein records to FASTA. Returns the number of records written."""
    n = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for r in records:
            fh.write(f">{r.protein_id} region_id={r.region_id}\n")
            for i in range(0, len(r.sequence), line_width):
                fh.write(r.sequence[i : i + line_width] + "\n")
            n += 1
    return n


# ────────────────────────────── parquet ────────────────────────────────────

def write_embeddings_parquet(
    output_path: Path,
    batches: Iterable[tuple[list[ProteinRecord], np.ndarray]],
    embedding_dim: int,
) -> int:
    """Stream embeddings into a parquet file batch by batch.

    Uses a fixed-size list type for the embedding column so downstream tools
    can lay it out as a 2-D matrix efficiently. Holds one batch in memory
    at a time — constant memory regardless of input size.

    Returns the total number of rows written.
    """
    schema = pa.schema([
        ("protein_id", pa.string()),
        ("region_id",  pa.string()),
        ("embedding",  pa.list_(pa.float32(), embedding_dim)),
    ])

    n_written = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pq.ParquetWriter(output_path, schema, compression="zstd") as writer:
        for records, vecs in batches:
            assert vecs.shape == (len(records), embedding_dim), (
                f"unexpected vec shape {vecs.shape}, "
                f"expected ({len(records)}, {embedding_dim})"
            )
            # Build the fixed-size-list array from a flat float32 buffer —
            # avoids materializing Python-side nested lists.
            flat = pa.array(
                vecs.reshape(-1).astype(np.float32, copy=False),
                type=pa.float32(),
            )
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
