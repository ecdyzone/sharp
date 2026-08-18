#!/usr/bin/env python3
"""Merge per-genome baseline output into one predictions.parquet.

The array jobs (`run_<tool>_array.sbatch`) leave one output directory per
genome. `evaluate.py` wants a single predictions file covering the whole
benchmark scope, so this walks those directories, runs the tool's existing
converter over each, and concatenates the result.

This deliberately reuses `convert_<tool>_to_parquet.py` rather than reparsing
anything: every assumption about column names and coordinate bases stays in the
one converter per tool that is tested against a real output fixture. Adding a
tool here is a dispatch-table entry, not a parser.

A genome that produced no output is reported rather than skipped silently — a
missing genome shrinks the predictions but *not* the recall denominator (that
comes from `--contigs`), so it would otherwise look like the tool simply found
nothing there.

Usage:
    # Inspect what would be merged, without writing
    python scripts/merge_predictions.py --tool antismash \\
        --input-dir ~/projects/antismash/out_benchmark --inspect

    python scripts/merge_predictions.py --tool antismash \\
        --input-dir ~/projects/antismash/out_benchmark \\
        --output data/interim/antismash_predictions.parquet

    python scripts/merge_predictions.py --tool deepbgc \\
        --input-dir ~/projects/deepbgc/out_benchmark \\
        --output data/interim/deepbgc_predictions.parquet

    # Warn about genomes in the scope file that produced no output at all
    python scripts/merge_predictions.py --tool deepbgc \\
        --input-dir ~/projects/deepbgc/out_benchmark \\
        --contigs data/interim/benchmark_set/analyzed_contigs.txt \\
        --output data/interim/deepbgc_predictions.parquet
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Callable

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import convert_antismash_to_parquet as _antismash  # noqa: E402
import convert_deepbgc_to_parquet as _deepbgc  # noqa: E402
import convert_gecco_to_parquet as _gecco  # noqa: E402
from sharp.config import INTERIM_DIR  # noqa: E402
from sharp.io import PredictedRegion, write_predictions_parquet  # noqa: E402

LOG = logging.getLogger("merge_predictions")


# ═══════════════════════════ per-tool dispatch ═════════════════════════════
# Each entry turns one genome's output directory into regions, by calling the
# same functions that tool's own converter uses. Nothing here parses a tool
# format directly — that stays in the converters, which have the fixtures.
# ═══════════════════════════════════════════════════════════════════════════

def _antismash_regions(path: Path) -> list[PredictedRegion]:
    data = _antismash.load_summary(_antismash.resolve_summary_path(path))
    regions: list[PredictedRegion] = []
    for record in data.get("records", []):
        regions.extend(_antismash.record_to_regions(record))
    return regions


def _deepbgc_regions(path: Path) -> list[PredictedRegion]:
    rows = _deepbgc.load_rows(_deepbgc.resolve_bgc_tsv_path(path))
    return _deepbgc.rows_to_regions(rows)


def _gecco_regions(path: Path) -> list[PredictedRegion]:
    rows = _gecco.load_rows(_gecco.resolve_clusters_tsv_path(path))
    return _gecco.rows_to_regions(rows)


TOOLS: dict[str, Callable[[Path], list[PredictedRegion]]] = {
    "antismash": _antismash_regions,
    "deepbgc": _deepbgc_regions,
    "gecco": _gecco_regions,
}


# ══════════════════════════════ merging ════════════════════════════════════

def genome_dirs(input_dir: Path) -> list[Path]:
    """The per-genome output directories, one per array task, sorted by name."""
    return sorted(d for d in input_dir.iterdir() if d.is_dir())


def load_contigs(path: Path) -> list[str]:
    """Read the `--contigs` scope file: one accession per line."""
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def merge(
    input_dir: Path, tool: str
) -> tuple[list[PredictedRegion], dict[str, int], list[str]]:
    """Convert every per-genome output under input_dir.

    Returns (regions, per_genome_counts, failures) where `failures` names the
    directories whose output could not be read at all.
    """
    parse = TOOLS[tool]
    regions: list[PredictedRegion] = []
    counts: dict[str, int] = {}
    failures: list[str] = []

    for d in genome_dirs(input_dir):
        try:
            got = parse(d)
        except (FileNotFoundError, ValueError, KeyError) as e:
            LOG.warning("%s — could not read output: %s", d.name, e)
            failures.append(d.name)
            continue
        counts[d.name] = len(got)
        regions.extend(got)
        LOG.debug("%s — %d regions", d.name, len(got))

    return regions, counts, failures


def report_missing(counts: dict[str, int], contigs: list[str]) -> list[str]:
    """Accessions in the scope file with no output directory at all.

    Worth surfacing loudly: the recall denominator comes from `--contigs`, not
    from the predictions, so a genome whose run never completed is scored as
    though the tool looked and found nothing.
    """
    return [c for c in contigs if c not in counts]


# ══════════════════════════════ orchestration ══════════════════════════════

def run(
    input_dir: Path,
    tool: str,
    output_path: Path | None,
    contigs_path: Path | None,
    do_inspect: bool,
) -> None:
    if not input_dir.is_dir():
        raise SystemExit(f"not a directory: {input_dir}")

    regions, counts, failures = merge(input_dir, tool)

    contigs = load_contigs(contigs_path) if contigs_path else []
    missing = report_missing(counts, contigs) if contigs else []

    LOG.info("%s: %d genome dir(s), %d region(s) total",
             tool, len(counts), len(regions))
    if failures:
        LOG.warning("%d genome(s) had unreadable output: %s",
                    len(failures), ", ".join(failures[:10])
                    + (" ..." if len(failures) > 10 else ""))
    if missing:
        LOG.warning("%d genome(s) in the scope file produced NO output — they "
                    "will be scored as though the tool found nothing: %s",
                    len(missing), ", ".join(missing[:10])
                    + (" ..." if len(missing) > 10 else ""))
    empty = [g for g, n in counts.items() if n == 0]
    if empty:
        LOG.info("%d genome(s) ran but predicted nothing: %s",
                 len(empty), ", ".join(empty[:10])
                 + (" ..." if len(empty) > 10 else ""))

    if do_inspect:
        print(f"\n{'='*68}\nMERGE PREVIEW — {tool}\n{'='*68}")
        print(f"{'genome':<24} {'regions':>8}")
        for genome in sorted(counts):
            print(f"{genome:<24} {counts[genome]:>8}")
        print(f"{'-'*34}\n{'TOTAL':<24} {len(regions):>8}")
        if missing:
            print(f"\nno output at all ({len(missing)}): {', '.join(missing)}")
        if failures:
            print(f"\nunreadable ({len(failures)}): {', '.join(failures)}")
        return

    if output_path is None:
        raise SystemExit("--output is required unless --inspect is given")
    n = write_predictions_parquet(output_path, regions)
    LOG.info("wrote %d region rows → %s", n, output_path)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--tool", required=True, choices=sorted(TOOLS),
                   help="which baseline produced these outputs")
    p.add_argument("--input-dir", type=Path, required=True,
                   help="directory holding one output subdirectory per genome")
    p.add_argument("--output", type=Path, default=None,
                   help="predictions parquet to write "
                        f"(e.g. {INTERIM_DIR / '<tool>_predictions.parquet'})")
    p.add_argument("--contigs", type=Path, default=None,
                   help="scope file; genomes listed here with no output are "
                        "reported, since they would otherwise look like the "
                        "tool found nothing")
    p.add_argument("--inspect", action="store_true",
                   help="print the per-genome region counts and write nothing")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    run(
        input_dir=args.input_dir,
        tool=args.tool,
        output_path=args.output,
        contigs_path=args.contigs,
        do_inspect=args.inspect,
    )


if __name__ == "__main__":
    main()
