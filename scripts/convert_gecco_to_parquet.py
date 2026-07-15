#!/usr/bin/env python3
"""Convert a GECCO run's output into predictions.parquet.

Parses GECCO's `<genome>.clusters.tsv` (one row per predicted cluster) into
the `PredictedRegion` schema used throughout S(H)ARP. S(H)ARP never invokes
GECCO itself — the user runs it in its own isolated env (see
`scripts/setup_gecco.sh`); this script only parses the output file it leaves
behind.

IMPORTANT — verify the schema first:
    Run --inspect on a real GECCO output directory (or its `.clusters.tsv`
    file directly) before trusting the parser:

        python scripts/convert_gecco_to_parquet.py --inspect <output dir>

    This prints the column list and what each accessor resolves to. If
    GECCO's TSV layout differs from what's encoded here, the FIELD PATHS
    section below is the only thing that needs editing.

Coordinate convention (verified 2026-07-15 against a real GECCO 0.10.3 run —
see CLAUDE.md "Baseline integration" for full evidence): `start`/`end` are
**1-based inclusive** — the one baseline tool where this holds (antiSMASH and
DeepBGC are both already 0-based half-open). Confirmed across all 5 rows of
the verification run: matching `.gbk` LOCUS bp length = `end - start + 1`
every time. Converted here to S(H)ARP's 0-based half-open with `start - 1`,
`end` unchanged (same pattern as MiBIG ingest).

Usage:
    # Inspect the real schema (do this first)
    python scripts/convert_gecco_to_parquet.py --inspect <gecco output dir>

    # Convert
    python scripts/convert_gecco_to_parquet.py \\
        --input <gecco output dir or <genome>.clusters.tsv> \\
        --output data/interim/gecco_predictions.parquet
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

from sharp.io import PredictedRegion, write_predictions_parquet

LOG = logging.getLogger("convert_gecco_to_parquet")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# These helpers isolate every assumption about GECCO's `.clusters.tsv`
# layout in ONE place. If --inspect shows a different layout, edit only this
# section.
#
# As of GECCO 0.10.3's `<genome>.clusters.tsv`, the columns relevant here are
# (15 columns total; unlisted ones — per-class probabilities, protein/domain
# id lists — are not used):
#   sequence_id  — contig
#   cluster_id   — unique per-cluster id, e.g. "AL589148.1_cluster_1"
#   start        — 1-based inclusive (verified; convert with start-1)
#   end          — 1-based inclusive (verified; unchanged after conversion)
#   average_p    — the model's cluster probability
#   type         — predicted BGC class; frequently "Unknown" in practice
# ════════════════════════════════════════════════════════════════════════════

def get_contig(row: dict[str, str]) -> str | None:
    return row.get("sequence_id") or None


def get_region_id(row: dict[str, str]) -> str | None:
    return row.get("cluster_id") or None


def get_coords(row: dict[str, str]) -> tuple[int, int] | None:
    """Parse (start, end) and convert 1-based inclusive -> 0-based
    half-open: start-1, end unchanged."""
    start, end = row.get("start"), row.get("end")
    if start is None or end is None:
        return None
    try:
        start_i, end_i = int(start), int(end)
    except ValueError:
        return None
    return max(0, start_i - 1), end_i


def get_p_bgc(row: dict[str, str]) -> float | None:
    p = row.get("average_p")
    if p is None:
        return None
    try:
        return float(p)
    except ValueError:
        return None


def get_predicted_class(row: dict[str, str]) -> str | None:
    cls = row.get("type")
    return cls if cls else None


# ══════════════════════════════ parsing ════════════════════════════════════

def row_to_region(row: dict[str, str]) -> PredictedRegion | None:
    """Convert one `.clusters.tsv` row into a PredictedRegion. Returns None
    if the row is missing a required field or has invalid coordinates."""
    contig = get_contig(row)
    region_id = get_region_id(row)
    coords = get_coords(row)
    p_bgc = get_p_bgc(row)
    if contig is None or region_id is None or coords is None or p_bgc is None:
        return None
    start, end = coords
    if end <= start:
        return None

    return PredictedRegion(
        region_id=region_id,
        contig=contig,
        start=start,
        end=end,
        p_bgc=p_bgc,
        predicted_class=get_predicted_class(row),
    )


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rows_to_regions(rows: list[dict[str, str]]) -> list[PredictedRegion]:
    out: list[PredictedRegion] = []
    for row in rows:
        region = row_to_region(row)
        if region is None:
            LOG.warning("skip unparseable row: cluster_id=%r", row.get("cluster_id"))
            continue
        out.append(region)
    return out


def resolve_clusters_tsv_path(input_path: Path) -> Path:
    """Accept either a direct path to <genome>.clusters.tsv or the GECCO
    output directory containing it (GECCO names the file after the input
    genome, e.g. sequence.clusters.tsv for sequence.fasta — not a fixed
    name)."""
    if input_path.is_file():
        return input_path
    if not input_path.is_dir():
        raise FileNotFoundError(f"no such file or directory: {input_path}")
    candidates = sorted(input_path.glob("*.clusters.tsv"))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"no *.clusters.tsv found in {input_path}")
    raise ValueError(
        f"multiple *.clusters.tsv files in {input_path}, pass one directly via "
        f"--input: {[c.name for c in candidates]}"
    )


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(input_path: Path) -> None:
    """Print the structure of the .clusters.tsv so the user can verify the
    field paths this script reads against the actual schema on disk."""
    tsv_path = resolve_clusters_tsv_path(input_path)
    rows = load_rows(tsv_path)

    print(f"\n{'='*70}\nFILE: {tsv_path}")
    print(f"{'='*70}")
    print(f"n_rows: {len(rows)}")
    if rows:
        print("columns:", list(rows[0].keys()))

    regions = rows_to_regions(rows)
    print(f"\n→ would produce {len(regions)} region row(s):")
    for r in regions[:10]:
        print(f"   {r}")


# ══════════════════════════════ orchestration ══════════════════════════════

def convert(input_path: Path, output_path: Path) -> int:
    tsv_path = resolve_clusters_tsv_path(input_path)
    rows = load_rows(tsv_path)
    regions = rows_to_regions(rows)

    if not regions:
        LOG.warning("no usable rows found in %s", tsv_path)

    n = write_predictions_parquet(output_path, regions)
    LOG.info("wrote %d region rows → %s", n, output_path)
    return n


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="PATH",
                   help="print the structure of a real GECCO output "
                        "(dir or .clusters.tsv) and exit — use this to "
                        "verify the schema before converting")
    p.add_argument("--input", type=Path,
                   help="GECCO output directory, or its <genome>.clusters.tsv directly")
    p.add_argument("--output", type=Path,
                   help="output predictions.parquet path")
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
