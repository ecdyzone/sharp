#!/usr/bin/env python3
"""Convert an antiSMASH run's output into predictions.parquet.

Parses the antiSMASH JSON summary (named after the input FASTA, e.g.
`sequence.json` for `sequence.fasta`) into the `PredictedRegion` schema used
throughout S(H)ARP. S(H)ARP never invokes antiSMASH itself — the user runs it
in its own isolated env (see `scripts/setup_antismash.sh`); this script only
parses the output files it leaves behind.

IMPORTANT — verify the schema first:
    Run --inspect on a real antiSMASH output directory (or its JSON file)
    before trusting the parser:

        python scripts/convert_antismash_to_parquet.py --inspect <output dir>

    This prints the record/feature structure and what each accessor resolves
    to. If antiSMASH's JSON layout differs from what's encoded here, the
    FIELD PATHS section below is the only thing that needs editing.

Coordinate convention (verified 2026-07-15 against a real antismash 8.0.4
run — see CLAUDE.md "Baseline integration" for full evidence): the `location`
string is already 0-based half-open, matching S(H)ARP's internal convention.
No conversion is applied.

Usage:
    # Inspect the real schema (do this first)
    python scripts/convert_antismash_to_parquet.py --inspect <antismash output dir>

    # Convert
    python scripts/convert_antismash_to_parquet.py \\
        --input <antismash output dir or sequence.json> \\
        --output data/interim/antismash_predictions.parquet
"""
from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Any

from sharp.io import PredictedRegion, write_predictions_parquet

LOG = logging.getLogger("convert_antismash_to_parquet")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# These helpers isolate every assumption about antiSMASH's JSON summary
# layout in ONE place. If --inspect shows a different layout, edit only this
# section.
#
# As of antiSMASH 8.0.4, the summary JSON looks like (abbreviated):
#   {
#     "records": [
#       {
#         "id": "AL589148.1",
#         "features": [
#           {
#             "type": "region",
#             "location": "[201195:222794](+)",   # string, NOT plain ints
#             "qualifiers": {
#               "region_number": ["1"],
#               "product": ["terpene"],            # can list >1 for hybrids
#               ...
#             }
#           },
#           ... (other feature types: CDS, gene, source, ...)
#         ]
#       }
#     ]
#   }
# ════════════════════════════════════════════════════════════════════════════

LOCATION_RE = re.compile(r"\[(\d+):(\d+)\]")


def get_location(feature: dict[str, Any]) -> tuple[int, int] | None:
    """Parse a feature's `location` string into (start, end), 0-based
    half-open. Returns None if the field is missing or unparseable."""
    loc = feature.get("location")
    if not isinstance(loc, str):
        return None
    m = LOCATION_RE.search(loc)
    if m is None:
        return None
    return int(m.group(1)), int(m.group(2))


def get_region_number(feature: dict[str, Any]) -> int | None:
    numbers = feature.get("qualifiers", {}).get("region_number")
    if not numbers:
        return None
    try:
        return int(numbers[0])
    except (TypeError, ValueError):
        return None


def get_products(feature: dict[str, Any]) -> list[str]:
    products = feature.get("qualifiers", {}).get("product")
    return [str(p) for p in products] if isinstance(products, list) else []


# ══════════════════════════════ parsing ════════════════════════════════════

def feature_to_region(feature: dict[str, Any], contig: str) -> PredictedRegion | None:
    """Convert one `type == "region"` feature into a PredictedRegion. Returns
    None if the feature is missing coordinates or a region number."""
    coords = get_location(feature)
    region_number = get_region_number(feature)
    if coords is None or region_number is None:
        return None
    start, end = coords
    if end <= start:
        return None

    products = get_products(feature)
    predicted_class = ";".join(products) if products else None

    return PredictedRegion(
        region_id=f"{contig}.region{region_number:03d}",
        contig=contig,
        start=start,
        end=end,
        p_bgc=1.0,
        predicted_class=predicted_class,
    )


def record_to_regions(record: dict[str, Any]) -> list[PredictedRegion]:
    """Convert one record's region features into PredictedRegion rows."""
    contig = record.get("id")
    if not contig:
        return []
    out: list[PredictedRegion] = []
    for feature in record.get("features", []):
        if feature.get("type") != "region":
            continue
        region = feature_to_region(feature, contig)
        if region is None:
            LOG.warning("skip unparseable region feature in %s", contig)
            continue
        out.append(region)
    return out


def resolve_summary_path(input_path: Path) -> Path:
    """Accept either a direct path to the summary JSON or the antiSMASH
    output directory containing it (antiSMASH names the file after the input
    FASTA, e.g. sequence.json for sequence.fasta — not a fixed name)."""
    if input_path.is_file():
        return input_path
    if not input_path.is_dir():
        raise FileNotFoundError(f"no such file or directory: {input_path}")
    candidates = sorted(input_path.glob("*.json"))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"no .json summary found in {input_path}")
    raise ValueError(
        f"multiple .json files in {input_path}, pass one directly via --input: "
        f"{[c.name for c in candidates]}"
    )


def load_summary(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(input_path: Path) -> None:
    """Print the structure of the summary JSON so the user can verify the
    field paths this script reads against the actual schema on disk."""
    summary_path = resolve_summary_path(input_path)
    data = load_summary(summary_path)

    print(f"\n{'='*70}\nFILE: {summary_path}")
    print(f"{'='*70}")
    print("top-level keys:", sorted(data.keys()))
    records = data.get("records", [])
    print(f"n_records: {len(records)}")

    for rec in records:
        contig = rec.get("id")
        feats = rec.get("features", [])
        types: dict[str, int] = {}
        for f in feats:
            t = f.get("type", "?")
            types[t] = types.get(t, 0) + 1
        print(f"\nrecord id={contig!r}  n_features={len(feats)}  types={types}")

        regions = record_to_regions(rec)
        print(f"→ would produce {len(regions)} region row(s):")
        for r in regions:
            print(f"   {r}")


# ══════════════════════════════ orchestration ══════════════════════════════

def convert(input_path: Path, output_path: Path) -> int:
    summary_path = resolve_summary_path(input_path)
    data = load_summary(summary_path)

    regions: list[PredictedRegion] = []
    for record in data.get("records", []):
        regions.extend(record_to_regions(record))

    if not regions:
        LOG.warning("no region features found in %s", summary_path)

    n = write_predictions_parquet(output_path, regions)
    LOG.info("wrote %d region rows → %s", n, output_path)
    return n


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="PATH",
                   help="print the structure of a real antiSMASH output "
                        "(dir or summary .json) and exit — use this to "
                        "verify the schema before converting")
    p.add_argument("--input", type=Path,
                   help="antiSMASH output directory, or its summary .json directly")
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
