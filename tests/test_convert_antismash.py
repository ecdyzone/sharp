"""Tests for scripts.convert_antismash_to_parquet.

`tests/fixtures/antismash_sequence.json` is a trimmed but real antiSMASH
8.0.4 summary (from a run on AL589148.1): two real `region` features (one
single-product, one hybrid multi-product) plus one real `CDS` feature, kept
to verify the type filter. If antiSMASH's JSON layout changes, --inspect on a
fresh real output is what catches it — this fixture pins the behavior we've
designed for.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from convert_antismash_to_parquet import (  # noqa: E402
    convert,
    feature_to_region,
    get_location,
    get_products,
    get_region_number,
    record_to_regions,
    resolve_summary_path,
)
from sharp.io import load_predictions_parquet  # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "antismash_sequence.json"


# ────────────────────────────── accessors ───────────────────────────────────

def region_feature(location="[100:200](+)", region_number="1", products=("terpene",)):
    return {
        "type": "region",
        "location": location,
        "qualifiers": {
            "region_number": [region_number],
            "product": list(products),
        },
    }


class TestAccessors:
    def test_get_location_parses_half_open_string(self) -> None:
        assert get_location(region_feature(location="[201195:222794](+)")) == (201195, 222794)

    def test_get_location_ignores_strand(self) -> None:
        assert get_location(region_feature(location="[100:200](-)")) == (100, 200)

    def test_get_location_missing(self) -> None:
        assert get_location({"type": "region"}) is None

    def test_get_region_number(self) -> None:
        assert get_region_number(region_feature(region_number="2")) == 2

    def test_get_region_number_missing(self) -> None:
        assert get_region_number({"qualifiers": {}}) is None

    def test_get_products_multi(self) -> None:
        assert get_products(region_feature(products=("furan", "butyrolactone"))) == \
            ["furan", "butyrolactone"]

    def test_get_products_missing(self) -> None:
        assert get_products({"qualifiers": {}}) == []


# ────────────────────────────── feature_to_region ───────────────────────────

class TestFeatureToRegion:
    def test_single_product(self) -> None:
        r = feature_to_region(region_feature(), "AL589148.1")
        assert r is not None
        assert r.region_id == "AL589148.1.region001"
        assert r.contig == "AL589148.1"
        assert r.start == 100
        assert r.end == 200
        assert r.p_bgc == 1.0
        assert r.predicted_class == "terpene"

    def test_hybrid_product_joined(self) -> None:
        r = feature_to_region(
            region_feature(region_number="2", products=("furan", "butyrolactone")),
            "AL589148.1",
        )
        assert r is not None
        assert r.predicted_class == "furan;butyrolactone"

    def test_region_id_padded_and_matches_gbk_filename_convention(self) -> None:
        r = feature_to_region(region_feature(region_number="12"), "CONTIG1")
        assert r is not None
        assert r.region_id == "CONTIG1.region012"

    def test_missing_location_returns_none(self) -> None:
        f = region_feature()
        del f["location"]
        assert feature_to_region(f, "AL589148.1") is None

    def test_missing_region_number_returns_none(self) -> None:
        f = region_feature()
        f["qualifiers"] = {"product": ["terpene"]}
        assert feature_to_region(f, "AL589148.1") is None

    def test_inverted_coords_returns_none(self) -> None:
        f = region_feature(location="[200:100](+)")
        assert feature_to_region(f, "AL589148.1") is None


# ────────────────────────────── record_to_regions ───────────────────────────

class TestRecordToRegions:
    def test_filters_non_region_features(self) -> None:
        record = {
            "id": "AL589148.1",
            "features": [
                region_feature(),
                {"type": "CDS", "location": "[10:20](+)", "qualifiers": {}},
            ],
        }
        regions = record_to_regions(record)
        assert len(regions) == 1

    def test_no_regions(self) -> None:
        record = {"id": "AL589148.1", "features": [{"type": "CDS"}]}
        assert record_to_regions(record) == []

    def test_missing_contig_id(self) -> None:
        record = {"features": [region_feature()]}
        assert record_to_regions(record) == []


# ────────────────────────────── resolve_summary_path ────────────────────────

class TestResolveSummaryPath:
    def test_direct_file(self, tmp_path: Path) -> None:
        f = tmp_path / "sequence.json"
        f.write_text("{}")
        assert resolve_summary_path(f) == f

    def test_directory_with_single_json(self, tmp_path: Path) -> None:
        (tmp_path / "sequence.json").write_text("{}")
        (tmp_path / "index.html").write_text("<html></html>")
        assert resolve_summary_path(tmp_path) == tmp_path / "sequence.json"

    def test_directory_with_no_json_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            resolve_summary_path(tmp_path)

    def test_directory_with_multiple_json_raises(self, tmp_path: Path) -> None:
        (tmp_path / "a.json").write_text("{}")
        (tmp_path / "b.json").write_text("{}")
        with pytest.raises(ValueError):
            resolve_summary_path(tmp_path)

    def test_missing_path_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            resolve_summary_path(tmp_path / "does_not_exist")


# ────────────────────────────── convert (I/O, real fixture) ─────────────────

class TestConvertEndToEnd:
    def test_real_fixture_produces_expected_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "antismash_predictions.parquet"
        n = convert(FIXTURE, out)
        assert n == 2

        regions = load_predictions_parquet(out)
        by_id = {r.region_id: r for r in regions}
        assert set(by_id) == {"AL589148.1.region001", "AL589148.1.region002"}

        r1 = by_id["AL589148.1.region001"]
        assert (r1.contig, r1.start, r1.end) == ("AL589148.1", 201195, 222794)
        assert r1.predicted_class == "terpene"
        assert r1.p_bgc == 1.0

        r2 = by_id["AL589148.1.region002"]
        assert (r2.contig, r2.start, r2.end) == ("AL589148.1", 226409, 255381)
        assert r2.predicted_class == "furan;butyrolactone"

    def test_convert_accepts_directory(self, tmp_path: Path) -> None:
        # Simulate the antiSMASH output dir shape: one summary JSON among
        # other non-JSON output.
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "sequence.json").write_text(FIXTURE.read_text())
        (out_dir / "index.html").write_text("<html></html>")

        out = tmp_path / "predictions.parquet"
        n = convert(out_dir, out)
        assert n == 2
