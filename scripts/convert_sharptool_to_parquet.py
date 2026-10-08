#!/usr/bin/env python3
"""Convert the team's SHARP tool output into predictions.parquet.

"sharptool" is the team's SHARP BGC finder, named here to keep it apart from
this repo (also called sharp). Like the other baselines it runs outside this
env; this script only parses the `neighborhoods.tsv` it leaves behind.

That file has one row per protein, not per region. SHARP places anchors
(SARP proteins, afsR-box heptarepeats, core biosynthetic domains), draws a
window around each, and lists the genes inside it under one `block_id`. A
**block is one predicted region**: rows are grouped by block and collapsed.

Region extent (`--extent`):
    genes  (default) — first gene start to last gene end of the block, clipped
           to the block window. This is what SHARP reports as the cluster's
           contents, and how its author measures blocks against MiBIG.
    block  — the `block_id` window itself, which pads ~6 kb (median) of
           gene-less flank on each side.

Both run outputs are handled the same way:
    * per-locus MiBIG runs annotate with Bakta, which renames every record
      `contig_1`, so `nucleotide` no longer names the sequence. Rows whose
      contig is Bakta's `contig_N` are mapped back to `sample` (the locus
      accession the run was given); a sample with two Bakta contigs is refused,
      since only one of them can be the accession.
    * genome-database runs keep the real contig accession in `nucleotide`.

SHARP has no score and no BGC class: `p_bgc` is 1.0 (as for antiSMASH) and
`predicted_class` is left empty. `query_source` names the block's anchors, not
a product class, so it is deliberately not used as one.

IMPORTANT — verify the schema first:
        python scripts/convert_sharptool_to_parquet.py --inspect <neighborhoods.tsv>

    This prints the columns, the contig naming, the coordinate check below and
    the blocks it would produce. If the layout changed, the FIELD PATHS section
    is the only thing that needs editing.

Coordinate convention (verified 2026-10-08 on full real runs): gene `start`/
`end` and the `block_id` window are **1-based inclusive**. Evidence: CDS span
`end - start == 3 * (plen + 1) - 1` on 302,580 of 302,586 CDS rows of the
MiBIG-loci run and 2,080,203 of 2,091,002 of the genome-database run (never
`3 * (plen + 1)`); no block starts at 0, many start at exactly 1, and block
ends reach the contig length (`nlen`) but never pass it. Converted with
`start - 1`, `end` unchanged, as for GECCO.

Windows on circular records can wrap the origin (`block_id` end < start). Such
a block becomes two regions, one each side of the origin, suffixed
`.pre_origin` / `.post_origin`.

Genes are clipped to their block window before they enter a span. On complete
circular chromosomes a gene crossing the origin is written as `start=1,
end=<contig length>`, so unclipped, a block starting at position 1 inherits a
whole-chromosome span (7 regions of 2.5-10 Mb in the 2026-10-08 genome-database
run, against windows of at most 440 kb). Genes that only overhang the window
edge by a few kb are trimmed the same way.

Overlapping blocks (~15% of them) are kept as separate regions, since that is
what SHARP outputs. `--merge-overlapping` merges them per contig instead.

Usage:
    # Inspect the real schema (do this first)
    python scripts/convert_sharptool_to_parquet.py --inspect <neighborhoods.tsv>

    # MiBIG loci run → score with sharp.evaluate on the pool scope
    python scripts/convert_sharptool_to_parquet.py \\
        --input <sharp_loci_fasta_run/all_neighborhood.tsv> \\
        --output data/interim/sharptool_predictions_pool_bact_plus130.parquet

    # Genome-database run → raw tool-vs-tool tables
    python scripts/convert_sharptool_to_parquet.py \\
        --input <sharp_batch/neighborhoods.tsv> \\
        --output data/interim/sharptool_predictions_actino.parquet

    # Variants: the padded window instead of the gene span; merge overlaps
    python scripts/convert_sharptool_to_parquet.py --input ... --output ... \\
        --extent block --merge-overlapping
"""
from __future__ import annotations

import argparse
import collections
import csv
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from sharp.io import PredictedRegion, write_predictions_parquet
from sharp.metrics import merge_intervals

LOG = logging.getLogger("convert_sharptool_to_parquet")

EXTENTS = ("genes", "block")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# These helpers isolate every assumption about sharptool's neighborhoods.tsv
# in ONE place. If --inspect shows a different layout, edit only this section.
#
# 30 columns; the ones read here:
#   sample       — run input name: the locus accession (MiBIG loci run) or
#                  the assembly directory, e.g. GCF_002150765.1_ASM215076v1_genomic
#   nucleotide   — contig; Bakta's "contig_N" in the loci run
#   start, end   — the row's gene, 1-based inclusive
#   nlen         — contig length
#   block_id     — "<nucleotide>:<start>-<end>", the block window, 1-based
#                  inclusive; end < start when it wraps the origin
#   type         — feature type: CDS, tRNA, ncRNA, PSE, gap, assembly_gap, ...
#   plen         — protein length in aa (used by --inspect's coordinate check)
#   query_source — the block's anchor set, e.g. "Heptarepeat+amp_binding_nrps"
# ════════════════════════════════════════════════════════════════════════════

BAKTA_CONTIG = re.compile(r"^contig_\d+$")
BLOCK_ID = re.compile(r"^(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+)$")

# Assembly gaps are listed as rows but are not genes; they must not stretch a
# gene span.
NON_GENE_TYPES = frozenset({"gap", "assembly_gap"})


def get_sample(row: dict[str, str]) -> str | None:
    return row.get("sample") or None


def get_contig(row: dict[str, str]) -> str | None:
    """The real contig name. Bakta's `contig_N` is mapped back to `sample`."""
    nucleotide = row.get("nucleotide") or None
    if nucleotide is not None and BAKTA_CONTIG.match(nucleotide):
        return get_sample(row)
    return nucleotide


def is_bakta_renamed(row: dict[str, str]) -> bool:
    return bool(BAKTA_CONTIG.match(row.get("nucleotide") or ""))


def parse_block_id(block_id: str | None) -> tuple[str, int, int] | None:
    """`"contig:start-end"` → (contig, start, end), still 1-based inclusive."""
    m = BLOCK_ID.match(block_id or "")
    if m is None:
        return None
    return m["contig"], int(m["start"]), int(m["end"])


def get_gene_coords(row: dict[str, str]) -> tuple[int, int] | None:
    """The row's gene as (start, end), still 1-based inclusive."""
    try:
        return int(row["start"]), int(row["end"])
    except (KeyError, TypeError, ValueError):
        return None


def get_contig_length(row: dict[str, str]) -> int | None:
    try:
        return int(row["nlen"])
    except (KeyError, TypeError, ValueError):
        return None


def is_gene(row: dict[str, str]) -> bool:
    return row.get("type", "") not in NON_GENE_TYPES


def get_anchor_set(row: dict[str, str]) -> str:
    return row.get("query_source") or ""


# ══════════════════════════════ blocks ═════════════════════════════════════

PRE, POST = "pre_origin", "post_origin"


@dataclass
class Block:
    """One SHARP block, accumulated over its rows. Coordinates are 1-based
    inclusive, as in the file; conversion happens in `block_to_regions`."""
    contig: str
    start: int
    end: int
    contig_length: int | None = None
    anchors: str = ""
    # gene span per side of the origin; a block that does not wrap only uses PRE
    spans: dict[str, list[int]] = field(default_factory=dict)

    @property
    def wraps(self) -> bool:
        return self.end < self.start

    @property
    def block_id(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}"

    def add_gene(self, start: int, end: int) -> bool:
        """Extend the span with one gene, clipped to the window. Returns True
        if the gene reached past the window (clipped, or dropped if nothing of
        it is left inside)."""
        # In a wrapping block, a gene at or past the window start sits before
        # the origin; anything else sits after it.
        if not self.wraps:
            side, lo, hi = PRE, self.start, self.end
        elif start >= self.start:
            side, lo, hi = PRE, self.start, self.contig_length or end
        else:
            side, lo, hi = POST, 1, self.end
        clipped = start < lo or end > hi
        start, end = max(start, lo), min(end, hi)
        if end < start:
            return clipped
        span = self.spans.get(side)
        if span is None:
            self.spans[side] = [start, end]
        else:
            span[0] = min(span[0], start)
            span[1] = max(span[1], end)
        return clipped


def block_to_regions(block: Block, extent: str = "genes") -> list[PredictedRegion]:
    """A block's region(s), converted to 0-based half-open (start - 1).

    A wrapping block yields one region per side of the origin. With
    extent="genes", a side with no genes yields nothing.
    """
    if extent not in EXTENTS:
        raise ValueError(f"extent must be one of {EXTENTS}, got {extent!r}")

    if extent == "block":
        if not block.wraps:
            intervals = {None: (block.start, block.end)}
        elif block.contig_length is None:
            LOG.warning("skip wrapping block %s: no contig length", block.block_id)
            return []
        else:
            intervals = {PRE: (block.start, block.contig_length), POST: (1, block.end)}
    else:
        intervals = {
            (side if block.wraps else None): tuple(span)
            for side, span in block.spans.items()
        }

    regions = []
    for side, (start, end) in sorted(intervals.items(), key=lambda kv: kv[1]):
        if end < start:
            continue
        region_id = block.block_id if side is None else f"{block.block_id}.{side}"
        regions.append(PredictedRegion(
            region_id=region_id,
            contig=block.contig,
            start=start - 1,
            end=end,
            p_bgc=1.0,
            predicted_class=None,
        ))
    return regions


@dataclass
class ParseStats:
    rows: int = 0
    skipped_rows: int = 0
    origin_genes: int = 0       # a gene that itself spans the origin (end < start)
    clipped_genes: int = 0      # a gene reaching past its block window
    bakta_rows: int = 0
    samples: set[str] = field(default_factory=set)
    bakta_contigs: dict[str, set[str]] = field(default_factory=lambda: collections.defaultdict(set))


def accumulate_blocks(
    rows: Iterable[dict[str, str]], stats: ParseStats | None = None,
) -> dict[tuple[str, int, int], Block]:
    """Group rows into blocks keyed by (contig, start, end).

    Raises ValueError if a sample carries two Bakta contigs: `contig_N` can then
    no longer be mapped back to the one accession the sample names.
    """
    stats = stats if stats is not None else ParseStats()
    blocks: dict[tuple[str, int, int], Block] = {}

    for row in rows:
        stats.rows += 1
        sample = get_sample(row)
        contig = get_contig(row)
        parsed = parse_block_id(row.get("block_id"))
        # block_id names the row's own contig; anything else means the layout
        # changed under us.
        if (sample is None or contig is None or parsed is None
                or parsed[0] != row.get("nucleotide")):
            stats.skipped_rows += 1
            continue
        stats.samples.add(sample)
        if is_bakta_renamed(row):
            stats.bakta_rows += 1
            stats.bakta_contigs[sample].add(row["nucleotide"])

        _, start, end = parsed
        key = (contig, start, end)
        block = blocks.get(key)
        if block is None:
            block = blocks[key] = Block(
                contig=contig, start=start, end=end,
                contig_length=get_contig_length(row),
                anchors=get_anchor_set(row),
            )

        coords = get_gene_coords(row)
        if coords is None or not is_gene(row):
            continue
        if coords[1] < coords[0]:
            stats.origin_genes += 1
            continue
        stats.clipped_genes += block.add_gene(*coords)

    ambiguous = {s: c for s, c in stats.bakta_contigs.items() if len(c) > 1}
    if ambiguous:
        example = next(iter(ambiguous.items()))
        raise ValueError(
            f"{len(ambiguous)} sample(s) have more than one Bakta contig, so "
            f"contig_N cannot be mapped back to the sample accession "
            f"(e.g. {example[0]}: {sorted(example[1])})"
        )
    return blocks


def merge_overlapping(regions: list[PredictedRegion]) -> list[PredictedRegion]:
    """Merge overlapping or abutting regions per contig into one each."""
    by_contig: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for r in regions:
        by_contig[r.contig].append((r.start, r.end))
    merged = []
    for contig in sorted(by_contig):
        for start, end in merge_intervals(by_contig[contig]):
            merged.append(PredictedRegion(
                region_id=f"{contig}:merged:{start + 1}-{end}",
                contig=contig, start=start, end=end,
                p_bgc=1.0, predicted_class=None,
            ))
    return merged


def blocks_to_regions(
    blocks: dict[tuple[str, int, int], Block],
    extent: str = "genes",
    merge: bool = False,
) -> list[PredictedRegion]:
    regions: list[PredictedRegion] = []
    n_empty = 0
    for key in sorted(blocks):
        found = block_to_regions(blocks[key], extent)
        n_empty += not found
        regions.extend(found)
    if n_empty:
        LOG.warning("%d block(s) produced no region (no gene rows)", n_empty)
    return merge_overlapping(regions) if merge else regions


# ══════════════════════════════ I/O ════════════════════════════════════════

def iter_rows(path: Path) -> Iterator[dict[str, str]]:
    """Stream rows; the genome-database run is ~1.6 GB, so never load it."""
    # The `sequence` column holds whole proteins; lift csv's 128 KiB field cap.
    csv.field_size_limit(2**31 - 1)
    with path.open(newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def read_header(path: Path) -> list[str]:
    with path.open(newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"), [])


# ══════════════════════════════ inspect mode ═══════════════════════════════

def _cds_span_convention(row: dict[str, str]) -> str | None:
    """Classify one CDS row's span against its protein length (+1 stop codon):
    3*(plen+1) - 1 is 1-based inclusive, 3*(plen+1) is 0-based half-open."""
    if row.get("type") != "CDS":
        return None
    coords = get_gene_coords(row)
    try:
        nt = 3 * (int(row["plen"]) + 1)
    except (KeyError, TypeError, ValueError):
        return None
    if coords is None:
        return None
    span = coords[1] - coords[0]
    return {nt - 1: "1-based inclusive", nt: "0-based half-open"}.get(span, "other")


def inspect(path: Path, extent: str = "genes") -> None:
    """Print the file's structure so the field paths and the coordinate
    convention encoded above can be checked against real output."""
    print(f"\n{'='*70}\nFILE: {path}\n{'='*70}")
    print("columns:", read_header(path))

    conventions: collections.Counter[str] = collections.Counter()
    types: collections.Counter[str] = collections.Counter()

    def tap(rows: Iterable[dict[str, str]]) -> Iterator[dict[str, str]]:
        for row in rows:
            types[row.get("type", "")] += 1
            c = _cds_span_convention(row)
            if c is not None:
                conventions[c] += 1
            yield row

    stats = ParseStats()
    blocks = accumulate_blocks(tap(iter_rows(path)), stats)

    print(f"n_rows: {stats.rows}  (skipped, unparseable: {stats.skipped_rows})")
    print(f"n_samples: {len(stats.samples)}")
    print(f"rows with a Bakta contig_N (mapped to sample): {stats.bakta_rows}")
    print(f"feature types: {dict(types.most_common())}")
    print(f"CDS span vs protein length: {dict(conventions.most_common())}")
    print(f"genes spanning the origin (left out of spans): {stats.origin_genes}")
    print(f"genes reaching past their block window (clipped): {stats.clipped_genes}")

    n_wrap = sum(b.wraps for b in blocks.values())
    print(f"\nn_blocks: {len(blocks)}  (wrapping the origin: {n_wrap})")
    starts_at_0 = sum(b.start == 0 for b in blocks.values())
    starts_at_1 = sum(b.start == 1 for b in blocks.values())
    past_end = sum(
        b.contig_length is not None and max(b.start, b.end) > b.contig_length
        for b in blocks.values()
    )
    print(f"block starts at 0: {starts_at_0}, at 1: {starts_at_1}; "
          f"ends past the contig: {past_end}  (1-based inclusive: 0 / >0 / 0)")
    anchors = collections.Counter(
        a for b in blocks.values() for a in b.anchors.split("+") if a
    )
    print(f"anchor types (blocks carrying each): {dict(anchors.most_common())}")

    regions = blocks_to_regions(blocks, extent)
    print(f"\n→ would produce {len(regions)} region row(s) with --extent {extent}:")
    for r in regions[:10]:
        print(f"   {r}")


# ══════════════════════════════ orchestration ══════════════════════════════

def convert(
    input_path: Path, output_path: Path,
    extent: str = "genes", merge: bool = False,
) -> int:
    if not input_path.is_file():
        raise FileNotFoundError(f"no such file: {input_path}")

    stats = ParseStats()
    blocks = accumulate_blocks(iter_rows(input_path), stats)
    if stats.skipped_rows:
        LOG.warning("skipped %d unparseable row(s)", stats.skipped_rows)
    if stats.origin_genes:
        LOG.info("%d gene(s) span the origin; left out of gene spans", stats.origin_genes)
    if stats.clipped_genes:
        LOG.info("%d gene(s) reach past their block window; clipped to it", stats.clipped_genes)
    LOG.info("%d rows → %d blocks over %d samples (%d rows Bakta-renamed)",
             stats.rows, len(blocks), len(stats.samples), stats.bakta_rows)

    regions = blocks_to_regions(blocks, extent, merge)
    if not regions:
        LOG.warning("no usable blocks found in %s", input_path)

    n = write_predictions_parquet(output_path, regions)
    LOG.info("wrote %d region rows (extent=%s%s) → %s",
             n, extent, ", merged" if merge else "", output_path)
    return n


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="PATH",
                   help="print the structure of a real neighborhoods.tsv and "
                        "exit — use this to verify the schema before converting")
    p.add_argument("--input", type=Path,
                   help="sharptool neighborhoods.tsv (one row per protein)")
    p.add_argument("--output", type=Path,
                   help="output predictions.parquet path")
    p.add_argument("--extent", choices=EXTENTS, default="genes",
                   help="region extent: the block's gene span (default) or "
                        "its padded block_id window")
    p.add_argument("--merge-overlapping", action="store_true",
                   help="merge overlapping blocks per contig into one region")
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
        inspect(args.inspect, args.extent)
        return

    if args.input is None or args.output is None:
        p.error("either --inspect PATH, or both --input and --output, are required")
    convert(args.input, args.output, args.extent, args.merge_overlapping)


if __name__ == "__main__":
    main()
