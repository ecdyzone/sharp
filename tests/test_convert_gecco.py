"""Tests for scripts.convert_gecco_to_parquet.

`tests/fixtures/gecco_sequence.clusters.tsv` is the real, unmodified
`sequence.clusters.tsv` from a GECCO 0.10.3 run on AL589148.1 (5 cluster
rows, small enough to check in verbatim). All 5 rows have `type == "Unknown"`
— this is the real distribution, not an edge case we invented.

Coordinate conversion is the key thing under test here: GECCO's start/end
are 1-based inclusive (verified against matching .gbk LOCUS bp lengths — see
CLAUDE.md "Baseline integration"), the one baseline tool where a conversion
is actually needed (antiSMASH and DeepBGC are both already 0-based
half-open). If GECCO's TSV layout or coordinate convention changes,
--inspect on a fresh real output is what catches it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from convert_gecco_to_parquet import (  # noqa: E402
    convert,
    get_coords,
    get_contig,
    get_p_bgc,
    get_predicted_class,
    get_region_id,
    resolve_clusters_tsv_path,
    row_to_region,
    rows_to_regions,
)
from sharp.io import load_predictions_parquet  # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "gecco_sequence.clusters.tsv"


# ────────────────────────────── accessors ───────────────────────────────────

def cluster_row(
    sequence_id="AL589148.1",
    cluster_id="AL589148.1_cluster_1",
    start="101",
    end="200",
    average_p="0.9",
    cluster_type="Unknown",
):
    return {
        "sequence_id": sequence_id,
        "cluster_id": cluster_id,
        "start": start,
        "end": end,
        "average_p": average_p,
        "type": cluster_type,
    }


class TestAccessors:
    def test_get_contig(self) -> None:
        assert get_contig(cluster_row()) == "AL589148.1"

    def test_get_contig_missing(self) -> None:
        assert get_contig({}) is None

    def test_get_region_id(self) -> None:
        assert get_region_id(cluster_row()) == "AL589148.1_cluster_1"

    def test_get_coords_converts_1based_inclusive_to_0based_half_open(self) -> None:
        # Real evidence: cluster_1 start=20274 end=53842, LOCUS bp=33569.
        # 1-based inclusive span+1 == LOCUS bp, so start-1, end unchanged.
        assert get_coords(cluster_row(start="20274", end="53842")) == (20273, 53842)

    def test_get_coords_start_1_stays_0(self) -> None:
        assert get_coords(cluster_row(start="1", end="100")) == (0, 100)

    def test_get_coords_missing_field(self) -> None:
        row = cluster_row()
        del row["end"]
        assert get_coords(row) is None

    def test_get_coords_non_numeric(self) -> None:
        assert get_coords(cluster_row(start="abc")) is None

    def test_get_p_bgc(self) -> None:
        assert get_p_bgc(cluster_row(average_p="0.9240075263403591")) == \
            pytest.approx(0.9240075263403591)

    def test_get_p_bgc_non_numeric(self) -> None:
        assert get_p_bgc(cluster_row(average_p="n/a")) is None

    def test_get_predicted_class_populated(self) -> None:
        assert get_predicted_class(cluster_row(cluster_type="NRPS")) == "NRPS"

    def test_get_predicted_class_unknown_kept_as_is(self) -> None:
        # "Unknown" is a real value GECCO emits — not converted to None.
        assert get_predicted_class(cluster_row(cluster_type="Unknown")) == "Unknown"

    def test_get_predicted_class_blank_is_none(self) -> None:
        assert get_predicted_class(cluster_row(cluster_type="")) is None


# ────────────────────────────── row_to_region ───────────────────────────────

class TestRowToRegion:
    def test_happy_path_applies_coordinate_conversion(self) -> None:
        r = row_to_region(cluster_row(
            cluster_id="AL589148.1_cluster_1",
            start="20274", end="53842",
            average_p="0.9240075263403591", cluster_type="Unknown",
        ))
        assert r is not None
        assert r.region_id == "AL589148.1_cluster_1"
        assert r.contig == "AL589148.1"
        assert r.start == 20273
        assert r.end == 53842
        assert r.p_bgc == pytest.approx(0.9240075263403591)
        assert r.predicted_class == "Unknown"

    def test_missing_region_id_returns_none(self) -> None:
        row = cluster_row()
        row["cluster_id"] = ""
        assert row_to_region(row) is None

    def test_missing_score_returns_none(self) -> None:
        row = cluster_row()
        del row["average_p"]
        assert row_to_region(row) is None

    def test_inverted_coords_returns_none(self) -> None:
        assert row_to_region(cluster_row(start="200", end="100")) is None

    def test_equal_coords_after_conversion_returns_none(self) -> None:
        # start=1 -> 0 after conversion; end=0 would be inverted/degenerate.
        assert row_to_region(cluster_row(start="1", end="0")) is None


class TestRowsToRegions:
    def test_skips_bad_rows_keeps_good_ones(self) -> None:
        rows = [cluster_row(), {"sequence_id": "x"}]  # second row is unparseable
        regions = rows_to_regions(rows)
        assert len(regions) == 1


# ────────────────────────────── resolve_clusters_tsv_path ───────────────────

class TestResolveClustersTsvPath:
    def test_direct_file(self, tmp_path: Path) -> None:
        f = tmp_path / "sequence.clusters.tsv"
        f.write_text("sequence_id\tstart\n")
        assert resolve_clusters_tsv_path(f) == f

    def test_directory_with_single_match(self, tmp_path: Path) -> None:
        (tmp_path / "sequence.clusters.tsv").write_text("sequence_id\n")
        (tmp_path / "sequence.features.tsv").write_text("x\n")
        assert resolve_clusters_tsv_path(tmp_path) == tmp_path / "sequence.clusters.tsv"

    def test_directory_with_no_match_raises(self, tmp_path: Path) -> None:
        (tmp_path / "sequence.features.tsv").write_text("x\n")
        with pytest.raises(FileNotFoundError):
            resolve_clusters_tsv_path(tmp_path)

    def test_directory_with_multiple_matches_raises(self, tmp_path: Path) -> None:
        (tmp_path / "a.clusters.tsv").write_text("x\n")
        (tmp_path / "b.clusters.tsv").write_text("x\n")
        with pytest.raises(ValueError):
            resolve_clusters_tsv_path(tmp_path)

    def test_missing_path_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            resolve_clusters_tsv_path(tmp_path / "does_not_exist")


# ────────────────────────────── convert (I/O, real fixture) ─────────────────

class TestConvertEndToEnd:
    def test_real_fixture_produces_expected_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "gecco_predictions.parquet"
        n = convert(FIXTURE, out)
        assert n == 5

        regions = load_predictions_parquet(out)
        by_id = {r.region_id: r for r in regions}
        assert set(by_id) == {
            "AL589148.1_cluster_1",
            "AL589148.1_cluster_2",
            "AL589148.1_cluster_3",
            "AL589148.1_cluster_4",
            "AL589148.1_cluster_5",
        }

        # Real coords from the TSV (1-based inclusive) converted to 0-based
        # half-open: start-1, end unchanged.
        r1 = by_id["AL589148.1_cluster_1"]
        assert (r1.contig, r1.start, r1.end) == ("AL589148.1", 20273, 53842)
        assert r1.p_bgc == pytest.approx(0.9240075263403591, rel=1e-5)
        assert r1.predicted_class == "Unknown"

        r3 = by_id["AL589148.1_cluster_3"]
        assert (r3.start, r3.end) == (179007, 223790)

    def test_convert_accepts_directory(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "sequence.clusters.tsv").write_text(FIXTURE.read_text())
        (out_dir / "sequence.features.tsv").write_text("unrelated\n")

        out = tmp_path / "predictions.parquet"
        n = convert(out_dir, out)
        assert n == 5
