"""I/O for the pipeline.

Owns the data types that flow between steps (`ProteinRecord`,
`PredictedRegion`, `KnownCluster`) and the read/write functions for the
file formats produced and consumed by the pipeline.
"""
from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

if TYPE_CHECKING:
    from sharp.metrics import BenchmarkResult

LOG = logging.getLogger(__name__)
REGION_ID_RE = re.compile(r"region_id=([^\s]+)")


# ────────────────────────────── data types ─────────────────────────────────

@dataclass(frozen=True)
class ProteinRecord:
    """A protein with its identifier and the BGC-candidate region it belongs to."""
    protein_id: str
    region_id: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class PredictedRegion:
    """A region the pipeline predicted as a candidate BGC, with the model's
    probability score. Coordinates are half-open: [start, end)."""
    region_id: str
    contig: str
    start: int
    end: int
    p_bgc: float
    predicted_class: str | None = None


@dataclass(frozen=True)
class KnownCluster:
    """A known BGC from a ground-truth source (e.g. MiBIG). Coordinates are
    half-open: [start, end)."""
    cluster_id: str
    contig: str
    start: int
    end: int
    cluster_class: str | None = None


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


# ────────────────────────────── parquet (embeddings) ───────────────────────

def write_embeddings_parquet(
    output_path: Path,
    batches: Iterable[tuple[list[ProteinRecord], np.ndarray]],
    embedding_dim: int,
) -> int:
    """Stream embeddings into a parquet file batch by batch."""
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


# ────────────────────────────── parquet (predictions) ──────────────────────

PREDICTIONS_SCHEMA = pa.schema([
    ("region_id",       pa.string()),
    ("contig",          pa.string()),
    ("start",           pa.int64()),
    ("end",             pa.int64()),
    ("p_bgc",           pa.float32()),
    ("predicted_class", pa.string()),   # nullable
])


def load_predictions_parquet(path: Path) -> list[PredictedRegion]:
    """Read a predictions.parquet file into PredictedRegion objects."""
    table = pq.read_table(path)
    cols = {name: table.column(name).to_pylist() for name in table.column_names}
    n = table.num_rows
    classes = cols.get("predicted_class", [None] * n)
    return [
        PredictedRegion(
            region_id=cols["region_id"][i],
            contig=cols["contig"][i],
            start=int(cols["start"][i]),
            end=int(cols["end"][i]),
            p_bgc=float(cols["p_bgc"][i]),
            predicted_class=classes[i],
        )
        for i in range(n)
    ]


def write_predictions_parquet(
    path: Path, predictions: list[PredictedRegion]
) -> int:
    """Write PredictedRegion objects to parquet. Used by mock data generation
    and (eventually) by predict.py."""
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pydict(
        {
            "region_id":       [p.region_id for p in predictions],
            "contig":          [p.contig for p in predictions],
            "start":           [p.start for p in predictions],
            "end":             [p.end for p in predictions],
            "p_bgc":           [p.p_bgc for p in predictions],
            "predicted_class": [p.predicted_class for p in predictions],
        },
        schema=PREDICTIONS_SCHEMA,
    )
    pq.write_table(table, path, compression="zstd")
    return len(predictions)


# ────────────────────────────── TSV (ground truth) ─────────────────────────

GROUND_TRUTH_COLUMNS = ["cluster_id", "contig", "start", "end", "class"]


def load_ground_truth_tsv(path: Path) -> list[KnownCluster]:
    """Read a ground truth TSV. Required columns: cluster_id, contig, start,
    end. Optional: class. Extra columns are ignored.

    Rows with a degenerate interval (`end <= start`) are skipped with a warning.
    A zero- or negative-length cluster can never be matched and silently drags
    every recall number down, so it is dropped at the door rather than counted.
    """
    out: list[KnownCluster] = []
    n_degenerate = 0
    with Path(path).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"cluster_id", "contig", "start", "end"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"ground truth TSV missing columns: {sorted(missing)}")
        for row in reader:
            start, end = int(row["start"]), int(row["end"])
            if end <= start:
                n_degenerate += 1
                LOG.warning(
                    "skip %s — degenerate interval [%d, %d)",
                    row["cluster_id"], start, end,
                )
                continue
            out.append(KnownCluster(
                cluster_id=row["cluster_id"],
                contig=row["contig"],
                start=start,
                end=end,
                cluster_class=row.get("class") or None,
            ))
    if n_degenerate:
        LOG.warning("%d ground-truth rows dropped for degenerate coordinates", n_degenerate)
    return out


def write_ground_truth_tsv(path: Path, clusters: list[KnownCluster]) -> int:
    """Write KnownCluster objects to TSV. Used by mock data generation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=GROUND_TRUTH_COLUMNS, delimiter="\t")
        writer.writeheader()
        for c in clusters:
            writer.writerow({
                "cluster_id": c.cluster_id,
                "contig":     c.contig,
                "start":      c.start,
                "end":        c.end,
                "class":      c.cluster_class or "",
            })
    return len(clusters)


# ────────────────────────────── contig scope list ──────────────────────────

def load_contigs(path: Path) -> set[str]:
    """Read the set of contigs a tool was actually run on.

    This is the denominator of recall: ground-truth clusters on contigs that
    were never analyzed cannot be found and must not be counted as missed.

    Accepts either one contig name per line, or any tab-delimited file whose
    first field is the contig name — so a samtools `.fai` index, or the first
    two columns of a genome file, can be passed unmodified. Blank lines and
    `#` comments are ignored.
    """
    contigs: set[str] = set()
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            contigs.add(line.split("\t")[0].strip())
    return contigs


# ────────────────────────────── JSON (benchmark) ───────────────────────────

def write_benchmark_json(path: Path, result: "BenchmarkResult") -> None:
    """Serialize a BenchmarkResult to JSON. Stable, human-readable layout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(result), indent=2, sort_keys=False) + "\n")
