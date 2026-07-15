#!/usr/bin/env python3
"""Convert a DeepBGC run's output into predictions.parquet.

Parses DeepBGC's `<prefix>.bgc.tsv` (one row per predicted BGC candidate)
into the `PredictedRegion` schema used throughout S(H)ARP. S(H)ARP never
invokes DeepBGC itself — the user runs it in its own isolated env (see
`scripts/setup_deepbgc.sh`); this script only parses the output file it
leaves behind.

IMPORTANT — verify the schema first:
    Run --inspect on a real DeepBGC output directory (or its `.bgc.tsv`
    file directly) before trusting the parser:

        python scripts/convert_deepbgc_to_parquet.py --inspect <output dir>

    This prints the column list and what each accessor resolves to. If
    DeepBGC's TSV layout differs from what's encoded here, the FIELD PATHS
    section below is the only thing that needs editing.

Coordinate convention (verified 2026-07-15 against a real DeepBGC run — see
CLAUDE.md "Baseline integration" for full evidence): `nucl_start`/`nucl_end`
are already 0-based half-open, matching S(H)ARP's internal convention. No
conversion is applied. Note these are NOT named `start`/`end` (unlike GECCO,
which uses those names for a different, 1-based, convention).

Usage:
    # Inspect the real schema (do this first)
    python scripts/convert_deepbgc_to_parquet.py --inspect <deepbgc output dir>

    # Convert
    python scripts/convert_deepbgc_to_parquet.py \\
        --input <deepbgc output dir or <prefix>.bgc.tsv> \\
        --output data/interim/deepbgc_predictions.parquet
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any

from sharp.io import PredictedRegion, write_predictions_parquet

LOG = logging.getLogger("convert_deepbgc_to_parquet")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# These helpers isolate every assumption about DeepBGC's `.bgc.tsv` layout in
# ONE place. If --inspect shows a different layout, edit only this section.
#
# As of DeepBGC 0.1.0's `<prefix>.bgc.tsv`, the columns relevant here are
# (28 columns total; unlisted ones — detector metadata, per-activity /
# per-class probabilities, protein/pfam id lists — are not used):
#   sequence_id       — contig
#   bgc_candidate_id  — unique per-candidate id, e.g. "AL589148.1_31460-41750.1"
#   nucl_start        — 0-based half-open (verified; NOT named "start")
#   nucl_end          — 0-based half-open (verified; NOT named "end")
#   deepbgc_score     — the model's BGC probability
#   product_class     — often EMPTY STRING when DeepBGC has no confident
#                        product-class call — treat blank as no class, not
#                        as a missing/invalid row.
# ════════════════════════════════════════════════════════════════════════════

def get_contig(row: dict[str, str]) -> str | None:
    return row.get("sequence_id") or None


def get_region_id(row: dict[str, str]) -> str | None:
    return row.get("bgc_candidate_id") or None


def get_coords(row: dict[str, str]) -> tuple[int, int] | None:
    start, end = row.get("nucl_start"), row.get("nucl_end")
    if start is None or end is None:
        return None
    try:
        return int(start), int(end)
    except ValueError:
        return None


def get_p_bgc(row: dict[str, str]) -> float | None:
    score = row.get("deepbgc_score")
    if score is None:
        return None
    try:
        return float(score)
    except ValueError:
        return None


def get_predicted_class(row: dict[str, str]) -> str | None:
    cls = row.get("product_class")
    return cls if cls else None


# ══════════════════════════════ parsing ════════════════════════════════════

def row_to_region(row: dict[str, str]) -> PredictedRegion | None:
    """Convert one `.bgc.tsv` row into a PredictedRegion. Returns None if the
    row is missing a required field or has invalid coordinates."""
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
            LOG.warning("skip unparseable row: bgc_candidate_id=%r",
                        row.get("bgc_candidate_id"))
            continue
        out.append(region)
    return out


def resolve_bgc_tsv_path(input_path: Path) -> Path:
    """Accept either a direct path to <prefix>.bgc.tsv or the DeepBGC output
    directory containing it (DeepBGC names files after the --output prefix
    passed to `deepbgc pipeline`, e.g. out.bgc.tsv for --output out — not a
    fixed name)."""
    if input_path.is_file():
        return input_path
    if not input_path.is_dir():
        raise FileNotFoundError(f"no such file or directory: {input_path}")
    candidates = sorted(input_path.glob("*.bgc.tsv"))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"no *.bgc.tsv found in {input_path}")
    raise ValueError(
        f"multiple *.bgc.tsv files in {input_path}, pass one directly via --input: "
        f"{[c.name for c in candidates]}"
    )


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(input_path: Path) -> None:
    """Print the structure of the .bgc.tsv so the user can verify the field
    paths this script reads against the actual schema on disk."""
    tsv_path = resolve_bgc_tsv_path(input_path)
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
    tsv_path = resolve_bgc_tsv_path(input_path)
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
                   help="print the structure of a real DeepBGC output "
                        "(dir or .bgc.tsv) and exit — use this to verify "
                        "the schema before converting")
    p.add_argument("--input", type=Path,
                   help="DeepBGC output directory, or its <prefix>.bgc.tsv directly")
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
