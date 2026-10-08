"""Tests for scripts.convert_sharptool_to_parquet.

Two real, unmodified fixtures, each every row of the first 3 blocks of a full
2026-10-08 run (one row per protein, 30 columns):
    sharptool_loci_neighborhoods.tsv — per-locus MiBIG run (JOJE01000020.1),
        Bakta-annotated, so `nucleotide` is "contig_1" for every row
    sharptool_batch_neighborhoods.tsv — genome-database run
        (GCF_002150765.1), real contig accessions, two contigs

What is under test is what the diagnosis found: blocks are the regions, gene
and block coordinates are 1-based inclusive (start - 1 on conversion), Bakta's
contig_N maps back to the sample, and origin-wrapping blocks split in two.
Wrapping, gaps and ambiguous Bakta contigs do not occur in the fixtures' first
3 blocks, so they are covered with synthetic rows below.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from convert_sharptool_to_parquet import (  # noqa: E402
    Block,
    ParseStats,
    accumulate_blocks,
    block_to_regions,
    blocks_to_regions,
    convert,
    get_contig,
    get_gene_coords,
    inspect,
    merge_overlapping,
    parse_block_id,
)
from sharp.io import PredictedRegion, load_predictions_parquet  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"
LOCI = FIXTURES / "sharptool_loci_neighborhoods.tsv"
BATCH = FIXTURES / "sharptool_batch_neighborhoods.tsv"


def gene_row(
    sample="ACC.1", nucleotide="contig_1", start="101", end="400",
    nlen="10000", block_id="contig_1:1-1000", feature_type="CDS",
    plen="99", query_source="btad_sarp",
):
    return {
        "sample": sample, "nucleotide": nucleotide, "start": start, "end": end,
        "nlen": nlen, "block_id": block_id, "type": feature_type, "plen": plen,
        "query_source": query_source,
    }


def region(contig: str, start: int, end: int) -> PredictedRegion:
    return PredictedRegion(f"{contig}:{start}", contig, start, end, 1.0)


# ────────────────────────────── accessors ───────────────────────────────────

class TestAccessors:
    def test_bakta_contig_maps_to_sample(self) -> None:
        assert get_contig(gene_row(sample="JOJE01000020.1", nucleotide="contig_1")) == "JOJE01000020.1"

    def test_real_contig_is_kept(self) -> None:
        row = gene_row(sample="GCF_002150765.1_ASM215076v1_genomic",
                       nucleotide="NZ_MUYM01000001.1")
        assert get_contig(row) == "NZ_MUYM01000001.1"

    def test_contig_named_like_bakta_only_when_exact(self) -> None:
        assert get_contig(gene_row(nucleotide="contig_1_extra")) == "contig_1_extra"

    def test_parse_block_id(self) -> None:
        assert parse_block_id("NZ_MUYM01000001.1:37736-55163") == ("NZ_MUYM01000001.1", 37736, 55163)

    def test_parse_block_id_wrapping_keeps_order(self) -> None:
        assert parse_block_id("CP003238.1:640689-6419") == ("CP003238.1", 640689, 6419)

    def test_parse_block_id_contig_with_colon(self) -> None:
        assert parse_block_id("gnl|X:ab:5-9") == ("gnl|X:ab", 5, 9)

    @pytest.mark.parametrize("bad", ["", None, "contig_1", "contig_1:5", "c:a-b"])
    def test_parse_block_id_unparseable(self, bad) -> None:
        assert parse_block_id(bad) is None

    def test_get_gene_coords(self) -> None:
        assert get_gene_coords(gene_row(start="55419", end="55898")) == (55419, 55898)

    def test_get_gene_coords_non_numeric(self) -> None:
        assert get_gene_coords(gene_row(start="")) is None


# ────────────────────────────── blocks → regions ────────────────────────────

class TestBlockToRegions:
    def test_gene_span_converts_1based_inclusive(self) -> None:
        b = Block("C.1", 100, 900)
        b.add_gene(201, 300)
        b.add_gene(501, 700)
        [r] = block_to_regions(b, "genes")
        assert (r.region_id, r.contig, r.start, r.end) == ("C.1:100-900", "C.1", 200, 700)
        assert r.p_bgc == 1.0 and r.predicted_class is None

    def test_block_extent_converts_1based_inclusive(self) -> None:
        b = Block("C.1", 100, 900)
        b.add_gene(201, 300)
        [r] = block_to_regions(b, "block")
        assert (r.start, r.end) == (99, 900)

    def test_block_starting_at_1_becomes_0(self) -> None:
        [r] = block_to_regions(Block("C.1", 1, 50), "block")
        assert (r.start, r.end) == (0, 50)

    def test_wrapping_block_extent_splits_at_origin(self) -> None:
        # Real case: CP003238.1:640689-6419 on a 690,188 bp record.
        b = Block("CP003238.1", 640689, 6419, contig_length=690188)
        pre, post = sorted(block_to_regions(b, "block"), key=lambda r: -r.start)
        assert (pre.region_id, pre.start, pre.end) == ("CP003238.1:640689-6419.pre_origin", 640688, 690188)
        assert (post.region_id, post.start, post.end) == ("CP003238.1:640689-6419.post_origin", 0, 6419)

    def test_wrapping_block_gene_spans_one_per_side(self) -> None:
        b = Block("C.1", 9000, 500, contig_length=10000)
        b.add_gene(9101, 9900)   # before the origin
        b.add_gene(9500, 9800)
        b.add_gene(11, 400)      # after it
        regions = {r.region_id: (r.start, r.end) for r in block_to_regions(b, "genes")}
        assert regions == {
            "C.1:9000-500.pre_origin": (9100, 9900),
            "C.1:9000-500.post_origin": (10, 400),
        }

    def test_wrapping_block_side_without_genes_yields_nothing(self) -> None:
        # Real case: LC361337.1:22928-20275 lists genes only after the origin.
        b = Block("LC361337.1", 22928, 20275, contig_length=28777)
        b.add_gene(2235, 14722)
        [r] = block_to_regions(b, "genes")
        assert (r.region_id, r.start, r.end) == ("LC361337.1:22928-20275.post_origin", 2234, 14722)

    def test_origin_gene_written_full_length_is_clipped_to_window(self) -> None:
        # Real case: NZ_CP016793.1:1-6575 on a ~10 Mb chromosome. A gene
        # crossing the origin is written start=1, end=<contig length>;
        # unclipped, the region was 9,997,872 bp long.
        b = Block("NZ_CP016793.1", 1, 6575, contig_length=9997872)
        assert b.add_gene(1, 9997872) is True
        assert b.add_gene(2001, 3000) is False
        [r] = block_to_regions(b, "genes")
        assert (r.start, r.end) == (0, 6575)

    def test_gene_overhanging_window_end_is_trimmed(self) -> None:
        # Real case: LC208006.1:165-139022, last gene ending at 144335.
        b = Block("LC208006.1", 165, 139022)
        b.add_gene(10001, 20000)
        assert b.add_gene(130001, 144335) is True
        [r] = block_to_regions(b, "genes")
        assert (r.start, r.end) == (10000, 139022)

    def test_gene_wholly_outside_window_is_dropped(self) -> None:
        b = Block("C.1", 100, 900)
        assert b.add_gene(1001, 1200) is True
        assert b.spans == {}

    def test_wrapping_pre_origin_side_is_clipped_to_contig_length(self) -> None:
        b = Block("C.1", 9000, 500, contig_length=10000)
        assert b.add_gene(9501, 10300) is True
        assert b.spans == {"pre_origin": [9501, 10000]}

    def test_wrapping_block_without_length_is_skipped(self) -> None:
        assert block_to_regions(Block("C.1", 9000, 500), "block") == []

    def test_block_without_genes_yields_nothing_for_gene_extent(self) -> None:
        assert block_to_regions(Block("C.1", 100, 900), "genes") == []

    def test_unknown_extent_raises(self) -> None:
        with pytest.raises(ValueError):
            block_to_regions(Block("C.1", 100, 900), "window")


# ────────────────────────────── accumulate_blocks ───────────────────────────

class TestAccumulateBlocks:
    def test_rows_group_into_one_block_per_block_id(self) -> None:
        rows = [
            gene_row(start="101", end="400"),
            gene_row(start="601", end="900"),
            gene_row(start="2101", end="2400", block_id="contig_1:2000-3000"),
        ]
        blocks = accumulate_blocks(rows)
        assert set(blocks) == {("ACC.1", 1, 1000), ("ACC.1", 2000, 3000)}
        assert blocks[("ACC.1", 1, 1000)].spans == {"pre_origin": [101, 900]}

    def test_gap_rows_do_not_stretch_the_span(self) -> None:
        rows = [
            gene_row(start="101", end="400"),
            gene_row(start="401", end="990", feature_type="assembly_gap"),
            gene_row(start="5", end="50", feature_type="gap"),
        ]
        [block] = accumulate_blocks(rows).values()
        assert block.spans == {"pre_origin": [101, 400]}

    def test_gene_spanning_origin_is_left_out_and_counted(self) -> None:
        stats = ParseStats()
        rows = [gene_row(start="101", end="400"), gene_row(start="9900", end="20")]
        [block] = accumulate_blocks(rows, stats).values()
        assert block.spans == {"pre_origin": [101, 400]}
        assert stats.origin_genes == 1

    def test_genes_past_the_window_are_counted(self) -> None:
        stats = ParseStats()
        rows = [gene_row(start="101", end="400"), gene_row(start="901", end="1100")]
        [block] = accumulate_blocks(rows, stats).values()
        assert block.spans == {"pre_origin": [101, 1000]}
        assert stats.clipped_genes == 1

    def test_sample_with_two_bakta_contigs_is_refused(self) -> None:
        rows = [gene_row(), gene_row(nucleotide="contig_2", block_id="contig_2:1-1000")]
        with pytest.raises(ValueError, match="more than one Bakta contig"):
            accumulate_blocks(rows)

    def test_unparseable_and_mismatched_rows_are_skipped(self) -> None:
        stats = ParseStats()
        rows = [
            gene_row(),
            gene_row(block_id="garbage"),
            gene_row(block_id="other_contig:1-1000"),  # names a different contig
            gene_row(sample=""),
        ]
        blocks = accumulate_blocks(rows, stats)
        assert len(blocks) == 1
        assert stats.skipped_rows == 3

    def test_records_anchor_set_and_contig_length(self) -> None:
        [block] = accumulate_blocks([gene_row(query_source="Heptarepeat+btad_sarp")]).values()
        assert block.anchors == "Heptarepeat+btad_sarp"
        assert block.contig_length == 10000


# ────────────────────────────── merge_overlapping ───────────────────────────

class TestMergeOverlapping:
    def test_overlapping_regions_merge_per_contig(self) -> None:
        merged = merge_overlapping([
            region("A", 0, 100), region("A", 50, 200), region("A", 500, 600),
            region("B", 60, 90),
        ])
        assert [(r.contig, r.start, r.end) for r in merged] == [
            ("A", 0, 200), ("A", 500, 600), ("B", 60, 90),
        ]
        assert merged[0].region_id == "A:merged:1-200"

    def test_blocks_to_regions_merge_flag(self) -> None:
        blocks = accumulate_blocks([
            gene_row(start="101", end="600"),
            gene_row(start="501", end="900", block_id="contig_1:400-1200"),
        ])
        assert len(blocks_to_regions(blocks, "genes")) == 2
        [merged] = blocks_to_regions(blocks, "genes", merge=True)
        assert (merged.start, merged.end) == (100, 900)


# ────────────────────────────── convert (I/O, real fixtures) ────────────────

class TestConvertEndToEnd:
    def test_loci_run_maps_bakta_contig_and_uses_gene_span(self, tmp_path: Path) -> None:
        out = tmp_path / "sharptool_predictions.parquet"
        assert convert(LOCI, out) == 3

        by_id = {r.region_id: r for r in load_predictions_parquet(out)}
        assert set(by_id) == {
            "JOJE01000020.1:52340-79117",     # btad_sarp anchor
            "JOJE01000020.1:196998-219036",   # amp_binding_nrps anchor
            "JOJE01000020.1:233244-254227",   # Heptarepeat anchor
        }
        assert {r.contig for r in by_id.values()} == {"JOJE01000020.1"}
        # First gene 55419, last gene end 76918 (1-based inclusive) → [55418, 76918)
        r1 = by_id["JOJE01000020.1:52340-79117"]
        assert (r1.start, r1.end) == (55418, 76918)
        assert r1.p_bgc == 1.0 and r1.predicted_class is None

    def test_loci_run_block_extent(self, tmp_path: Path) -> None:
        out = tmp_path / "p.parquet"
        convert(LOCI, out, extent="block")
        by_id = {r.region_id: r for r in load_predictions_parquet(out)}
        r1 = by_id["JOJE01000020.1:52340-79117"]
        assert (r1.start, r1.end) == (52339, 79117)

    def test_batch_run_keeps_real_contigs(self, tmp_path: Path) -> None:
        out = tmp_path / "p.parquet"
        assert convert(BATCH, out) == 3
        spans = {(r.contig, r.start, r.end) for r in load_predictions_parquet(out)}
        assert spans == {
            ("NZ_MUYM01000001.1", 41320, 53045),
            ("NZ_MUYM01000001.1", 117197, 127909),
            ("NZ_MUYM01000002.1", 34040, 46291),
        }

    def test_missing_input_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            convert(tmp_path / "nope.tsv", tmp_path / "p.parquet")


class TestInspect:
    def test_reports_coordinate_evidence(self, capsys) -> None:
        inspect(LOCI)
        out = capsys.readouterr().out
        assert "CDS span vs protein length: {'1-based inclusive': 38}" in out
        assert "rows with a Bakta contig_N (mapped to sample): 38" in out
        assert "n_blocks: 3" in out
