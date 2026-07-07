#!/usr/bin/env python3
"""Build a ground-truth TSV from a BGC Atlas ``complete-bgcs`` dump.

BGC Atlas (https://bgc-atlas.cs.uni-tuebingen.de) distributes complete BGCs as
one antiSMASH-produced GenBank file per region, e.g.::

    MGYA00004361_contig00011.region001.gbk

Each file is a single antiSMASH ``region`` extracted from a larger metagenome
assembly. This script walks such a directory and emits a TSV in the format
``sharp.io.load_ground_truth_tsv`` expects — the SAME schema as
``mibig_ground_truth.tsv`` — so BGC Atlas flows through ``evaluate.py`` unchanged::

    cluster_id  contig  start  end  class

⚠️  BGC Atlas is a SECONDARY, noisy ground truth: its positive labels are
    themselves antiSMASH predictions with no manual curation. Benchmark numbers
    against it are systematically optimistic. Always report alongside MiBIG.
    See CLAUDE.md → "Ground truth sources".

Schema facts (verified against the real 204k-file dump, 2026-07-07):

  * Genomic coordinates are NOT the LOCUS coordinates (those are region-local,
    ``1..N``). The region's position on its original contig lives in the
    antiSMASH structured COMMENT as ``Orig. start`` / ``Orig. end``.
  * ``Orig. start`` / ``Orig. end`` are ALREADY 0-based half-open — verified by
    ``Orig. end - Orig. start == len(seq)`` across thousands of files. Unlike
    MiBIG (1-based inclusive), NO coordinate conversion is applied here.
  * ``rec.id`` (the contig) is NOT globally unique — it repeats across
    assemblies, and a single contig can carry multiple regions (region001,
    region002, ...). So:
      - ``cluster_id`` = the filename stem (unique; includes assembly + contig
        + ``.regionNNN``).
      - ``contig``     = ``<MGYA assembly>_<rec.id>`` (assembly-qualified, so it
        is globally unique and safe for reciprocal-overlap matching).
  * Exactly one ``region`` feature per file. Its ``/product`` qualifier may hold
    multiple values (e.g. ``["thioamitides", "thiopeptide"]``), joined with "/".

Usage:
    # Inspect a few real files first (verify the schema on disk)
    python scripts/prepare_bgcatlas_ground_truth.py --inspect data/raw/complete-bgcs

    # Build the TSV (streams over all files; ~204k rows)
    python scripts/prepare_bgcatlas_ground_truth.py \\
        --input-dir data/raw/complete-bgcs \\
        --output data/raw/bgcatlas_ground_truth.tsv

    # Develop / test against a small subset without walking 10 GB
    python scripts/prepare_bgcatlas_ground_truth.py \\
        --input-dir data/raw/complete-bgcs \\
        --output /tmp/bgcatlas_sample.tsv --limit 100
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sharp.config import RAW_DIR
from sharp.io import KnownCluster, write_ground_truth_tsv

LOG = logging.getLogger("prepare_bgcatlas_ground_truth")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# Every assumption about the BGC Atlas .gbk layout is isolated here. If
# `--inspect` shows a different structure on disk, this is the only section
# that needs editing.
#
# A typical file (abbreviated):
#   LOCUS       contig00011   25524 bp ...
#   ACCESSION   contig00011                          <- rec.id (contig)
#   COMMENT     ##antiSMASH-Data-START##
#               Orig. start  :: 36256                <- genomic start (0-based)
#               Orig. end    :: 61780                <- genomic end (0-based, half-open)
#   FEATURES
#        region  1..25524
#                /product="thioamitides"             <- class(es)
#                /product="thiopeptide"
#                /region_number="1"
#
# Filename: <MGYA assembly>_<contig>.region<NNN>.gbk
# ════════════════════════════════════════════════════════════════════════════

# Assembly accession = everything up to the first underscore in the filename.
_ASSEMBLY_RE = re.compile(r"^([^_]+)_")

# antiSMASH stores the original coordinates in the structured comment under
# this section key.
_ANTISMASH_COMMENT_KEY = "antiSMASH-Data"


def get_assembly_id(filename: str) -> str | None:
    """Parse the MGnify assembly accession (the ``MGYA...`` prefix) from a
    BGC Atlas filename. Returns None if the filename has no ``_`` separator."""
    m = _ASSEMBLY_RE.match(filename)
    return m.group(1) if m else None


def get_cluster_id(filename: str) -> str:
    """The globally-unique cluster id is the filename stem — it already carries
    assembly + contig + ``.regionNNN``, so it distinguishes region001 from
    region002 on the same contig."""
    return Path(filename).stem


def get_contig(record: SeqRecord, filename: str) -> str | None:
    """Assembly-qualified contig: ``<assembly>_<rec.id>``.

    ``rec.id`` alone is not globally unique (it repeats across assemblies), so
    it is prefixed with the assembly accession from the filename. Falls back to
    the raw ``rec.id`` if the filename has no parseable assembly prefix."""
    contig = record.id
    if not contig or contig == "<unknown id>":
        return None
    assembly = get_assembly_id(filename)
    return f"{assembly}_{contig}" if assembly else contig


def get_orig_coords(record: SeqRecord) -> tuple[int | None, int | None]:
    """Extract (start, end) genomic coordinates from the antiSMASH structured
    comment. These are ALREADY 0-based half-open — returned as-is, no conversion.

    Returns (None, None) if the comment or its coordinate fields are absent."""
    sc = record.annotations.get("structured_comment")
    if not isinstance(sc, dict):
        return None, None
    data = sc.get(_ANTISMASH_COMMENT_KEY)
    if not isinstance(data, dict):
        return None, None
    start = data.get("Orig. start")
    end = data.get("Orig. end")
    try:
        start = int(start) if start is not None else None
        end = int(end) if end is not None else None
    except (TypeError, ValueError):
        return None, None
    return start, end


def get_region_class(record: SeqRecord) -> str | None:
    """Join the ``/product`` qualifier(s) of the (single) ``region`` feature
    with '/'. Returns None if there is no region feature or no product."""
    for feat in record.features:
        if feat.type == "region":
            products = feat.qualifiers.get("product")
            if products:
                return "/".join(str(p) for p in products)
    return None


# ══════════════════════════════ parsing ════════════════════════════════════

def record_to_cluster(record: SeqRecord, filename: str) -> KnownCluster | None:
    """Convert one parsed BGC Atlas GenBank record into a KnownCluster, or None
    if it lacks usable coordinates / contig.

    Coordinates come from the antiSMASH ``Orig. start``/``Orig. end`` comment
    and are already 0-based half-open — no conversion, unlike MiBIG."""
    contig = get_contig(record, filename)
    start, end = get_orig_coords(record)
    if contig is None or start is None or end is None:
        return None
    if end <= start:   # degenerate / inverted — skip defensively
        return None
    return KnownCluster(
        cluster_id=get_cluster_id(filename),
        contig=contig,
        start=start,
        end=end,
        cluster_class=get_region_class(record),
    )


def iter_gbk_files(input_dir: Path, limit: int | None = None) -> Iterator[Path]:
    """Yield ``.gbk`` paths under ``input_dir`` (recursively), in sorted order,
    stopping after ``limit`` files if given.

    Sorted so ``--limit`` selects a deterministic subset (and so region001 and
    region002 of a contig stay adjacent). The directory holds ~204k files, so
    this is a lazy generator — it never materializes the full list unless the
    caller consumes it all."""
    n = 0
    for path in sorted(input_dir.rglob("*.gbk")):
        yield path
        n += 1
        if limit is not None and n >= limit:
            return


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(input_dir: Path, n: int = 3) -> None:
    """Print the structure of the first n files so the schema accessors can be
    verified against the real data on disk (mirrors the MiBIG script)."""
    files = list(iter_gbk_files(input_dir, limit=n))
    if not files:
        LOG.error("no .gbk files found under %s", input_dir)
        return

    LOG.info("inspecting first %d .gbk file(s) under %s", len(files), input_dir)
    for f in files:
        print(f"\n{'='*70}\nFILE: {f.name}")
        print(f"{'='*70}")
        try:
            record = SeqIO.read(f, "genbank")
        except (ValueError, OSError) as e:
            print(f"  !! failed to parse: {e}")
            continue

        sc = record.annotations.get("structured_comment", {})
        print("structured_comment sections:", list(sc.keys()) if sc else "(none)")
        print("region features:",
              sum(1 for feat in record.features if feat.type == "region"))

        print("\n→ accessors resolve to:")
        print(f"   assembly_id   = {get_assembly_id(f.name)!r}")
        print(f"   cluster_id    = {get_cluster_id(f.name)!r}")
        print(f"   rec.id        = {record.id!r}")
        print(f"   contig        = {get_contig(record, f.name)!r}")
        start, end = get_orig_coords(record)
        print(f"   orig coords   = ({start}, {end})   len(seq)={len(record.seq)}")
        if start is not None and end is not None:
            match = (end - start == len(record.seq))
            print(f"   0-based check = end-start == len(seq)? {match}")
        print(f"   class         = {get_region_class(record)!r}")

        cluster = record_to_cluster(record, f.name)
        print(f"\n→ would produce: {cluster}")


# ══════════════════════════════ orchestration ══════════════════════════════

def build_ground_truth(
    input_dir: Path,
    output_path: Path,
    limit: int | None = None,
    log_every: int = 5000,
) -> int:
    """Parse every ``.gbk`` under ``input_dir`` into ground-truth rows and write
    the TSV. Returns the number of rows written.

    One file → at most one row (exactly one region per file). Files that fail to
    parse or lack usable coordinates are skipped with a count logged at the end.
    """
    clusters: list[KnownCluster] = []
    n_seen = 0
    n_parse_fail = 0
    n_no_coords = 0

    for path in iter_gbk_files(input_dir, limit=limit):
        n_seen += 1
        try:
            record = SeqIO.read(path, "genbank")
        except (ValueError, OSError) as e:
            n_parse_fail += 1
            LOG.warning("skip %s — parse error: %s", path.name, e)
            continue
        cluster = record_to_cluster(record, path.name)
        if cluster is None:
            n_no_coords += 1
            continue
        clusters.append(cluster)
        if log_every and n_seen % log_every == 0:
            LOG.info("processed %d files (%d rows so far)", n_seen, len(clusters))

    if n_seen == 0:
        LOG.error("no .gbk files found under %s", input_dir)
        return 0
    if n_parse_fail:
        LOG.warning("%d file(s) failed to parse", n_parse_fail)
    if n_no_coords:
        LOG.warning("%d file(s) had no usable coordinates/contig", n_no_coords)

    n = write_ground_truth_tsv(output_path, clusters)
    LOG.info("wrote %d cluster rows (from %d files) → %s",
             n, n_seen, output_path)
    return n


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="DIR",
                   help="print the structure of a few .gbk files and exit "
                        "(use this to verify the schema before building)")
    p.add_argument("--input-dir", type=Path,
                   help="directory of BGC Atlas .gbk files (complete-bgcs)")
    p.add_argument("--output", type=Path,
                   default=RAW_DIR / "bgcatlas_ground_truth.tsv",
                   help="output TSV (default: %(default)s)")
    p.add_argument("--limit", type=int, default=None,
                   help="process only the first N files (for dev/tests; "
                        "default: all)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.inspect is not None:
        inspect(args.inspect)
        return

    if args.input_dir is None:
        p.error("either --inspect DIR or --input-dir DIR is required")
    build_ground_truth(args.input_dir, args.output, limit=args.limit)


if __name__ == "__main__":
    main()
