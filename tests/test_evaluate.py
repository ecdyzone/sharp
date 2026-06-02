"""Tests for sharp.evaluate — orchestration, I/O round-trips, CLI."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from sharp.config import EvaluateConfig
from sharp.evaluate import build_parser, run
from sharp.io import (
    KnownCluster,
    PredictedRegion,
    load_ground_truth_tsv,
    load_predictions_parquet,
    write_ground_truth_tsv,
    write_predictions_parquet,
)


# ────────────────────────────── I/O round-trips ────────────────────────────

class TestPredictionsParquet:
    def test_round_trip(self, tmp_path: Path) -> None:
        # p_bgc is stored as float32, so we pick values that round-trip exactly
        # rather than asserting on float64 inputs that get truncated.
        preds = [
            PredictedRegion("R1", "chr1", 0, 1000, 0.5, "T1PKS"),       # exact in fp32
            PredictedRegion("R2", "chr1", 2000, 3500, 0.625, None),      # exact in fp32
        ]
        path = tmp_path / "p.parquet"
        n = write_predictions_parquet(path, preds)
        assert n == 2

        read_back = load_predictions_parquet(path)
        assert read_back == preds

    def test_round_trip_approximate_for_arbitrary_floats(self, tmp_path: Path) -> None:
        # Sanity check: arbitrary float64 inputs round-trip within float32 precision.
        preds = [PredictedRegion("R1", "chr1", 0, 1000, 0.95, "T1PKS")]
        path = tmp_path / "p.parquet"
        write_predictions_parquet(path, preds)
        read_back = load_predictions_parquet(path)
        assert read_back[0].region_id == "R1"
        assert read_back[0].p_bgc == pytest.approx(0.95, abs=1e-6)


class TestGroundTruthTsv:
    def test_round_trip(self, tmp_path: Path) -> None:
        clusters = [
            KnownCluster("BGC0001", "chr1", 0, 30000, "T1PKS"),
            KnownCluster("BGC0002", "chr2", 5000, 60000, "NRPS"),
        ]
        path = tmp_path / "gt.tsv"
        n = write_ground_truth_tsv(path, clusters)
        assert n == 2

        read_back = load_ground_truth_tsv(path)
        assert read_back == clusters

    def test_optional_class_column(self, tmp_path: Path) -> None:
        # Class column present but empty in some rows → None on read.
        path = tmp_path / "gt.tsv"
        path.write_text(
            "cluster_id\tcontig\tstart\tend\tclass\n"
            "C1\tchr1\t0\t100\tT1PKS\n"
            "C2\tchr1\t200\t300\t\n"
        )
        read_back = load_ground_truth_tsv(path)
        assert read_back[0].cluster_class == "T1PKS"
        assert read_back[1].cluster_class is None

    def test_missing_required_column_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "gt.tsv"
        path.write_text("cluster_id\tcontig\tstart\n" "C1\tchr1\t0\n")
        with pytest.raises(ValueError, match="missing columns"):
            load_ground_truth_tsv(path)

    def test_extra_columns_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "gt.tsv"
        path.write_text(
            "cluster_id\tcontig\tstart\tend\tclass\tnotes\textra\n"
            "C1\tchr1\t0\t100\tT1PKS\thello\twhatever\n"
        )
        read_back = load_ground_truth_tsv(path)
        assert len(read_back) == 1
        assert read_back[0].cluster_id == "C1"


# ────────────────────────────── run() orchestration ────────────────────────

class TestRun:
    def _write_inputs(self, tmp_path: Path,
                      preds: list[PredictedRegion],
                      clusters: list[KnownCluster]) -> tuple[Path, Path, Path]:
        p_path = tmp_path / "predictions.parquet"
        g_path = tmp_path / "gt.tsv"
        o_path = tmp_path / "benchmark.json"
        write_predictions_parquet(p_path, preds)
        write_ground_truth_tsv(g_path, clusters)
        return p_path, g_path, o_path

    def test_end_to_end_perfect_match(self, tmp_path: Path) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)

        run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                           output_path=o, min_overlap_frac=0.5))
        assert o.exists()

        result = json.loads(o.read_text())
        assert result["precision"] == 1.0
        assert result["recall"] == 1.0
        assert result["f1"] == 1.0
        assert result["n_predictions"] == 1
        assert result["n_clusters"] == 1
        assert result["recovered_cluster_ids"] == ["C1"]
        assert result["false_positive_prediction_ids"] == []

    def test_end_to_end_mixed_outcome(self, tmp_path: Path) -> None:
        preds = [
            PredictedRegion("R1", "chr1", 0, 100, 0.9),         # TP for C1
            PredictedRegion("R2", "chr1", 500, 600, 0.7),       # FP
        ]
        clusters = [
            KnownCluster("C1", "chr1", 0, 100),                 # recovered
            KnownCluster("C2", "chr1", 800, 900),               # missed
        ]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                           output_path=o, min_overlap_frac=0.5))

        result = json.loads(o.read_text())
        assert result["precision"] == 0.5
        assert result["recall"] == 0.5
        assert result["recovered_cluster_ids"] == ["C1"]
        assert result["missed_cluster_ids"] == ["C2"]
        assert result["false_positive_prediction_ids"] == ["R2"]

    def test_empty_ground_truth_exits(self, tmp_path: Path) -> None:
        p, g, o = self._write_inputs(tmp_path,
                                     [PredictedRegion("R1", "chr1", 0, 100, 0.9)],
                                     [])
        with pytest.raises(SystemExit) as exc_info:
            run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                               output_path=o, min_overlap_frac=0.5))
        assert exc_info.value.code == 1

    def test_creates_output_parent_dirs(self, tmp_path: Path) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p_path = tmp_path / "predictions.parquet"
        g_path = tmp_path / "gt.tsv"
        write_predictions_parquet(p_path, preds)
        write_ground_truth_tsv(g_path, clusters)

        nested = tmp_path / "a" / "b" / "benchmark.json"
        run(EvaluateConfig(predictions_path=p_path, ground_truth_path=g_path,
                           output_path=nested, min_overlap_frac=0.5))
        assert nested.exists()


# ────────────────────────────── CLI ────────────────────────────────────────

class TestCli:
    def test_required_args(self) -> None:
        parser = build_parser()
        args = parser.parse_args([
            "--predictions", "p.parquet",
            "--ground-truth", "gt.tsv",
            "--output", "out.json",
        ])
        assert args.predictions == Path("p.parquet")
        assert args.ground_truth == Path("gt.tsv")
        assert args.output == Path("out.json")
        assert args.min_overlap_frac == 0.5

    def test_predictions_required(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--ground-truth", "gt.tsv", "--output", "out.json"])

    def test_ground_truth_required(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--predictions", "p.parquet", "--output", "out.json"])

    def test_output_required(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--predictions", "p.parquet", "--ground-truth", "gt.tsv"])
