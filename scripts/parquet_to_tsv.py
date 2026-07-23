#!/usr/bin/env python3
"""Convert any S(H)ARP parquet file into a TSV.

Generic dump: reads a parquet file (predictions.parquet, embeddings.parquet,
kg_features.parquet, ...) and writes one TSV with the same columns, in
row-group batches so large files are never pulled fully into memory.

List-typed columns (e.g. embeddings.parquet's `embedding` column, a
fixed-size list of floats per protein) have no TSV representation. This
script joins them into a single cell with `,` so no value is silently
dropped, but the conversion is one-way — a TSV produced from a list-typed
column cannot be read back into the original parquet schema.

Usage:
    # Inspect a file's schema first (recommended, especially for anything
    # other than predictions.parquet)
    python scripts/parquet_to_tsv.py --inspect data/interim/embeddings.parquet

    # Convert
    python scripts/parquet_to_tsv.py \\
        --input data/interim/predictions.parquet \\
        --output data/interim/predictions.tsv
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import pyarrow.parquet as pq

LOG = logging.getLogger("parquet_to_tsv")


def cell(value: object) -> object:
    """Render one parquet value as a TSV-safe cell. Lists (e.g. embedding
    vectors) are joined with ',' since TSV has no native list type."""
    if isinstance(value, list):
        return ",".join(str(v) for v in value)
    return value


def convert(input_path: Path, output_path: Path) -> int:
    reader = pq.ParquetFile(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with output_path.open("w", newline="") as fh:
        writer = None
        for batch in reader.iter_batches():
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=batch.schema.names, delimiter="\t")
                writer.writeheader()
            for row in batch.to_pylist():
                writer.writerow({k: cell(v) for k, v in row.items()})
                n_written += 1

    LOG.info("wrote %d row(s) → %s", n_written, output_path)
    return n_written


def inspect(input_path: Path) -> None:
    """Print the schema of a parquet file so list-typed columns (which get
    joined into a single TSV cell) can be spotted before converting."""
    pf = pq.ParquetFile(input_path)
    print(f"\n{'='*70}\nFILE: {input_path}")
    print(f"{'='*70}")
    print(f"n_rows: {pf.metadata.num_rows}")
    print("schema:")
    for field in pf.schema_arrow:
        print(f"   {field.name}: {field.type}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="PATH",
                   help="print the schema of a parquet file and exit")
    p.add_argument("--input", type=Path, help="input parquet file")
    p.add_argument("--output", type=Path, help="output .tsv path")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> None:
    p = build_parser()
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.inspect is not None:
        inspect(args.inspect)
        return

    if args.input is None or args.output is None:
        p.error("either --inspect PATH, or both --input and --output, are required")
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
