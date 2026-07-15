"""Tests for scripts.convert_deepbgc_to_parquet.

`tests/fixtures/deepbgc_out.bgc.tsv` is the real, unmodified `out.bgc.tsv`
from a DeepBGC 0.1.0 run on AL589148.1 (5 candidate rows, small enough to
check in verbatim). Row 1 has a populated `product_class` (Polyketide); rows
2-5 have it blank — this is the real distribution, not an edge case we
invented, and the parser must handle both. If DeepBGC's TSV layout changes,
--inspect on a fresh real output is what catches it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from convert_deepbgc_to_parquet import (  # noqa: E402
    convert,
    get_coords,
    get_contig,
    get_p_bgc,
    get_predicted_class,
    get_region_id,
    resolve_bgc_tsv_path,
    row_to_region,
    rows_to_regions,
)
from sharp.io import load_predictions_parquet  # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "deepbgc_out.bgc.tsv"


# ────────────────────────────── accessors ───────────────────────────────────

def bgc_row(
    sequence_id="AL589148.1",
    bgc_candidate_id="AL589148.1_100-200.1",
    nucl_start="100",
    nucl_end="200",
    deepbgc_score="0.6",
    product_class="",
):
    return {
        "sequence_id": sequence_id,
        "bgc_candidate_id": bgc_candidate_id,
        "nucl_start": nucl_start,
        "nucl_end": nucl_end,
        "deepbgc_score": deepbgc_score,
        "product_class": product_class,
    }


class TestAccessors:
    def test_get_contig(self) -> None:
        assert get_contig(bgc_row()) == "AL589148.1"

    def test_get_contig_missing(self) -> None:
        assert get_contig({}) is None

    def test_get_region_id(self) -> None:
        assert get_region_id(bgc_row()) == "AL589148.1_100-200.1"

    def test_get_coords(self) -> None:
        assert get_coords(bgc_row(nucl_start="31460", nucl_end="41750")) == (31460, 41750)

    def test_get_coords_missing_field(self) -> None:
        row = bgc_row()
        del row["nucl_end"]
        assert get_coords(row) is None

    def test_get_coords_non_numeric(self) -> None:
        assert get_coords(bgc_row(nucl_start="abc")) is None

    def test_get_p_bgc(self) -> None:
        assert get_p_bgc(bgc_row(deepbgc_score="0.61532")) == pytest.approx(0.61532)

    def test_get_p_bgc_non_numeric(self) -> None:
        assert get_p_bgc(bgc_row(deepbgc_score="n/a")) is None

    def test_get_predicted_class_populated(self) -> None:
        assert get_predicted_class(bgc_row(product_class="Polyketide")) == "Polyketide"

    def test_get_predicted_class_blank_is_none(self) -> None:
        assert get_predicted_class(bgc_row(product_class="")) is None

    def test_get_predicted_class_absent_key_is_none(self) -> None:
        row = bgc_row()
        del row["product_class"]
        assert get_predicted_class(row) is None


# ────────────────────────────── row_to_region ───────────────────────────────

class TestRowToRegion:
    def test_happy_path(self) -> None:
        r = row_to_region(bgc_row(
            nucl_start="31460", nucl_end="41750",
            deepbgc_score="0.61532", product_class="Polyketide",
        ))
        assert r is not None
        assert r.region_id == "AL589148.1_100-200.1"
        assert r.contig == "AL589148.1"
        assert r.start == 31460
        assert r.end == 41750
        assert r.p_bgc == pytest.approx(0.61532)
        assert r.predicted_class == "Polyketide"

    def test_blank_product_class_becomes_none(self) -> None:
        r = row_to_region(bgc_row(product_class=""))
        assert r is not None
        assert r.predicted_class is None

    def test_missing_region_id_returns_none(self) -> None:
        row = bgc_row()
        row["bgc_candidate_id"] = ""
        assert row_to_region(row) is None

    def test_missing_score_returns_none(self) -> None:
        row = bgc_row()
        del row["deepbgc_score"]
        assert row_to_region(row) is None

    def test_inverted_coords_returns_none(self) -> None:
        assert row_to_region(bgc_row(nucl_start="200", nucl_end="100")) is None

    def test_equal_coords_returns_none(self) -> None:
        assert row_to_region(bgc_row(nucl_start="100", nucl_end="100")) is None


class TestRowsToRegions:
    def test_skips_bad_rows_keeps_good_ones(self) -> None:
        rows = [bgc_row(), {"sequence_id": "x"}]  # second row is unparseable
        regions = rows_to_regions(rows)
        assert len(regions) == 1


# ────────────────────────────── resolve_bgc_tsv_path ────────────────────────

class TestResolveBgcTsvPath:
    def test_direct_file(self, tmp_path: Path) -> None:
        f = tmp_path / "out.bgc.tsv"
        f.write_text("sequence_id\tstart\n")
        assert resolve_bgc_tsv_path(f) == f

    def test_directory_with_single_match(self, tmp_path: Path) -> None:
        (tmp_path / "out.bgc.tsv").write_text("sequence_id\n")
        (tmp_path / "out.pfam.tsv").write_text("x\n")
        assert resolve_bgc_tsv_path(tmp_path) == tmp_path / "out.bgc.tsv"

    def test_directory_with_no_match_raises(self, tmp_path: Path) -> None:
        (tmp_path / "out.pfam.tsv").write_text("x\n")
        with pytest.raises(FileNotFoundError):
            resolve_bgc_tsv_path(tmp_path)

    def test_directory_with_multiple_matches_raises(self, tmp_path: Path) -> None:
        (tmp_path / "a.bgc.tsv").write_text("x\n")
        (tmp_path / "b.bgc.tsv").write_text("x\n")
        with pytest.raises(ValueError):
            resolve_bgc_tsv_path(tmp_path)

    def test_missing_path_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            resolve_bgc_tsv_path(tmp_path / "does_not_exist")


# ────────────────────────────── convert (I/O, real fixture) ─────────────────

class TestConvertEndToEnd:
    def test_real_fixture_produces_expected_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "deepbgc_predictions.parquet"
        n = convert(FIXTURE, out)
        assert n == 5

        regions = load_predictions_parquet(out)
        by_id = {r.region_id: r for r in regions}
        assert set(by_id) == {
            "AL589148.1_31460-41750.1",
            "AL589148.1_60518-60743.1",
            "AL589148.1_63898-66577.1",
            "AL589148.1_159898-165658.1",
            "AL589148.1_230262-324563.1",
        }

        r1 = by_id["AL589148.1_31460-41750.1"]
        assert (r1.contig, r1.start, r1.end) == ("AL589148.1", 31460, 41750)
        assert r1.p_bgc == pytest.approx(0.61532)
        assert r1.predicted_class == "Polyketide"

        # Real distribution: most rows have an empty product_class.
        r2 = by_id["AL589148.1_60518-60743.1"]
        assert r2.predicted_class is None
        assert (r2.start, r2.end) == (60518, 60743)

    def test_convert_accepts_directory(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "out.bgc.tsv").write_text(FIXTURE.read_text())
        (out_dir / "out.pfam.tsv").write_text("unrelated\n")

        out = tmp_path / "predictions.parquet"
        n = convert(out_dir, out)
        assert n == 5
