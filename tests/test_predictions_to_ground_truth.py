"""Tests for scripts/predictions_to_ground_truth.py.

The shim exists so evaluate.py can measure tool-vs-tool agreement. What must
hold is that intervals survive the round trip untouched — an off-by-one here
would quietly move every agreement number.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from predictions_to_ground_truth import run, to_clusters  # noqa: E402
from sharp.io import (  # noqa: E402
    PredictedRegion,
    load_ground_truth_tsv,
    write_predictions_parquet,
)

REGIONS = [
    PredictedRegion("r1", "NC_003155.5", 100, 4_100, 0.9, "t1pks"),
    PredictedRegion("r2", "NC_003155.5", 9_000, 12_000, 0.4, ""),
    PredictedRegion("r3", "NZ_CP157456.1", 500, 2_500, 0.75, None),
]


@pytest.fixture
def preds(tmp_path: Path) -> Path:
    path = tmp_path / "tool.parquet"
    write_predictions_parquet(path, REGIONS)
    return path


class TestToClusters:
    def test_every_region_becomes_a_cluster(self, preds: Path) -> None:
        assert len(to_clusters(preds)) == 3

    def test_coordinates_pass_through_unchanged(self, preds: Path) -> None:
        # Both types are 0-based half-open, so there is nothing to convert.
        got = {c.cluster_id: (c.start, c.end) for c in to_clusters(preds)}
        assert got["r1"] == (100, 4_100)

    def test_contig_is_preserved(self, preds: Path) -> None:
        got = {c.cluster_id: c.contig for c in to_clusters(preds)}
        assert got["r3"] == "NZ_CP157456.1"

    def test_predicted_class_becomes_cluster_class(self, preds: Path) -> None:
        got = {c.cluster_id: c.cluster_class for c in to_clusters(preds)}
        assert got["r1"] == "t1pks"

    def test_min_p_bgc_filters(self, preds: Path) -> None:
        got = to_clusters(preds, min_p_bgc=0.5)
        assert {c.cluster_id for c in got} == {"r1", "r3"}

    def test_prefix_namespaces_the_ids(self, preds: Path) -> None:
        got = to_clusters(preds, prefix="antismash")
        assert all(c.cluster_id.startswith("antismash:") for c in got)


class TestRun:
    def test_writes_a_tsv_evaluate_can_read(self, preds: Path, tmp_path: Path) -> None:
        out = tmp_path / "as_gt.tsv"
        run(preds, out, 0.0, None)
        assert len(load_ground_truth_tsv(out)) == 3

    def test_round_trip_preserves_intervals(self, preds: Path, tmp_path: Path) -> None:
        out = tmp_path / "as_gt.tsv"
        run(preds, out, 0.0, None)
        got = {c.cluster_id: (c.contig, c.start, c.end)
               for c in load_ground_truth_tsv(out)}
        assert got["r2"] == ("NC_003155.5", 9_000, 12_000)

    def test_empty_reference_exits(self, preds: Path, tmp_path: Path) -> None:
        # An empty reference makes every agreement number meaningless, so it
        # must fail loudly rather than produce a confident 0.0.
        with pytest.raises(SystemExit):
            run(preds, tmp_path / "as_gt.tsv", 1.5, None)

    def test_missing_input_exits(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit):
            run(tmp_path / "nope.parquet", tmp_path / "out.tsv", 0.0, None)
