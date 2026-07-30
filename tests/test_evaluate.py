"""Tests for sharp.evaluate — orchestration, I/O round-trips, CLI.

The `TestRealBaselineOutput` class at the bottom is the important one: it runs
the benchmark over checked-in real converter output and pins the numbers. Every
defect this module was rewritten to fix shows up there.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from sharp.config import EvaluateConfig
from sharp.evaluate import build_parser, main, run
from sharp.io import (
    KnownCluster,
    PredictedRegion,
    load_contigs,
    load_ground_truth_tsv,
    load_predictions_parquet,
    write_ground_truth_tsv,
    write_predictions_parquet,
)

FIXTURES = Path(__file__).parent / "fixtures"


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

    def test_degenerate_rows_skipped(self, tmp_path: Path) -> None:
        # A zero-length cluster can never be matched; it would drag recall down
        # silently, so it is dropped at the door.
        path = tmp_path / "gt.tsv"
        path.write_text(
            "cluster_id\tcontig\tstart\tend\tclass\n"
            "C1\tchr1\t0\t100\tT1PKS\n"
            "C2\tchr1\t200\t200\tNRPS\n"
            "C3\tchr1\t400\t300\tterpene\n"
        )
        read_back = load_ground_truth_tsv(path)
        assert [c.cluster_id for c in read_back] == ["C1"]

    def test_accepts_str_path(self, tmp_path: Path) -> None:
        path = tmp_path / "gt.tsv"
        write_ground_truth_tsv(path, [KnownCluster("C1", "chr1", 0, 100)])
        assert len(load_ground_truth_tsv(str(path))) == 1


class TestLoadContigs:
    def test_one_per_line(self, tmp_path: Path) -> None:
        path = tmp_path / "contigs.txt"
        path.write_text("chr1\nchr2\nchr3\n")
        assert load_contigs(path) == {"chr1", "chr2", "chr3"}

    def test_fai_style_takes_first_field(self, tmp_path: Path) -> None:
        path = tmp_path / "genome.fa.fai"
        path.write_text("chr1\t248956422\t112\t60\t61\nchr2\t242193529\t253404903\t60\t61\n")
        assert load_contigs(path) == {"chr1", "chr2"}

    def test_blank_lines_and_comments_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "contigs.txt"
        path.write_text("# analyzed contigs\nchr1\n\n  \nchr2\n")
        assert load_contigs(path) == {"chr1", "chr2"}


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

        run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=o))
        assert o.exists()

        result = json.loads(o.read_text())
        assert result["detection"]["recall"] == 1.0
        assert result["detection"]["matched_prediction_frac"] == 1.0
        assert result["detection"]["f1"] == 1.0
        assert result["nucleotide"]["precision"] == 1.0
        assert result["scope"]["n_predictions_total"] == 1
        assert result["scope"]["n_clusters_total"] == 1
        assert result["recovered_cluster_ids"] == ["C1"]
        assert result["unmatched_prediction_ids"] == []

    def test_end_to_end_mixed_outcome(self, tmp_path: Path) -> None:
        preds = [
            PredictedRegion("R1", "chr1", 0, 100, 0.9),         # matches C1
            PredictedRegion("R2", "chr1", 500, 600, 0.7),       # unmatched
        ]
        clusters = [
            KnownCluster("C1", "chr1", 0, 100),                 # recovered
            KnownCluster("C2", "chr1", 800, 900),               # missed
        ]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=o))

        result = json.loads(o.read_text())
        assert result["detection"]["matched_prediction_frac"] == 0.5
        assert result["detection"]["recall"] == 0.5
        assert result["recovered_cluster_ids"] == ["C1"]
        assert result["missed_cluster_ids"] == ["C2"]
        assert result["unmatched_prediction_ids"] == ["R2"]

    def test_empty_ground_truth_exits(self, tmp_path: Path) -> None:
        p, g, o = self._write_inputs(tmp_path,
                                     [PredictedRegion("R1", "chr1", 0, 100, 0.9)],
                                     [])
        with pytest.raises(SystemExit) as exc_info:
            run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=o))
        assert exc_info.value.code == 1

    def test_creates_output_parent_dirs(self, tmp_path: Path) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, _ = self._write_inputs(tmp_path, preds, clusters)

        nested = tmp_path / "a" / "b" / "benchmark.json"
        run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=nested))
        assert nested.exists()

    def test_inferred_scope_warns(self, tmp_path: Path, caplog) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)

        with caplog.at_level(logging.WARNING, logger="evaluate"):
            run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=o))
        assert "no --contigs given" in caplog.text
        assert json.loads(o.read_text())["scope"]["source"] == "inferred"

    def test_explicit_contigs_no_warning(self, tmp_path: Path, caplog) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100), KnownCluster("C2", "chr2", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        c = tmp_path / "contigs.txt"
        c.write_text("chr1\nchr2\n")

        with caplog.at_level(logging.WARNING, logger="evaluate"):
            run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                               output_path=o, contigs_path=c))
        assert "no --contigs given" not in caplog.text

        result = json.loads(o.read_text())
        assert result["scope"]["source"] == "explicit"
        assert result["scope"]["n_contigs"] == 2
        # chr2 was analyzed but produced nothing — a real miss, unlike inferred scope.
        assert result["detection"]["recall"] == 0.5

    def test_empty_contigs_file_exits(self, tmp_path: Path) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        c = tmp_path / "contigs.txt"
        c.write_text("# nothing here\n")

        with pytest.raises(SystemExit) as exc_info:
            run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                               output_path=o, contigs_path=c))
        assert exc_info.value.code == 1

    def test_contig_name_mismatch_exits_with_diagnostic(
        self, tmp_path: Path, caplog
    ) -> None:
        # The live hazard: BGC Atlas contigs are assembly-qualified, so a bare
        # accession in the predictions matches nothing. Fail loudly, not silently.
        preds    = [PredictedRegion("R1", "contig00011", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "MGYA00004361_contig00011", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)

        with caplog.at_level(logging.ERROR, logger="evaluate"):
            with pytest.raises(SystemExit) as exc_info:
                run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                                   output_path=o))
        assert exc_info.value.code == 1
        assert "no ground-truth clusters on any contig in scope" in caplog.text
        assert "contig00011" in caplog.text
        assert "MGYA00004361_contig00011" in caplog.text
        assert not o.exists()

    def test_min_p_bgc_filters(self, tmp_path: Path) -> None:
        preds = [
            PredictedRegion("R1", "chr1", 0, 100, 0.9),
            PredictedRegion("R2", "chr1", 500, 600, 0.2),
        ]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        run(EvaluateConfig(predictions_path=p, ground_truth_path=g,
                           output_path=o, min_p_bgc=0.5))

        result = json.loads(o.read_text())
        assert result["scope"]["n_predictions_below_threshold"] == 1
        assert result["scope"]["min_p_bgc"] == 0.5
        assert result["unmatched_prediction_ids"] == []

    def test_asymmetric_criterion_written_to_json(self, tmp_path: Path) -> None:
        preds    = [PredictedRegion("R1", "chr1", 0, 100, 0.9)]
        clusters = [KnownCluster("C1", "chr1", 0, 100)]
        p, g, o = self._write_inputs(tmp_path, preds, clusters)
        run(EvaluateConfig(predictions_path=p, ground_truth_path=g, output_path=o,
                           min_cluster_frac=0.4, min_prediction_frac=0.2,
                           reciprocal_frac=0.6))

        result = json.loads(o.read_text())
        assert result["criterion"]["min_cluster_frac"] == 0.4
        assert result["criterion"]["min_prediction_frac"] == 0.2
        assert result["reciprocal_frac"] == 0.6


# ────────────────────────────── CLI ────────────────────────────────────────

class TestCli:
    def test_required_args_and_defaults(self) -> None:
        args = build_parser().parse_args([
            "--predictions", "p.parquet",
            "--ground-truth", "gt.tsv",
            "--output", "out.json",
        ])
        assert args.predictions == Path("p.parquet")
        assert args.ground_truth == Path("gt.tsv")
        assert args.output == Path("out.json")
        assert args.contigs is None
        assert args.min_cluster_frac == 0.5
        assert args.min_prediction_frac == 0.0
        assert args.reciprocal_frac == 0.5
        assert args.min_p_bgc == 0.0
        assert args.max_listed_ids == 1000

    def test_all_flags_parse(self) -> None:
        args = build_parser().parse_args([
            "--predictions", "p.parquet", "--ground-truth", "gt.tsv",
            "--output", "out.json", "--contigs", "c.txt",
            "--min-cluster-frac", "0.4", "--min-prediction-frac", "0.3",
            "--reciprocal-frac", "0.6", "--min-p-bgc", "0.5",
            "--max-listed-ids", "0",
        ])
        assert args.contigs == Path("c.txt")
        assert args.min_cluster_frac == 0.4
        assert args.min_prediction_frac == 0.3
        assert args.reciprocal_frac == 0.6
        assert args.min_p_bgc == 0.5
        assert args.max_listed_ids == 0

    def test_predictions_required(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["--ground-truth", "gt.tsv", "--output", "out.json"])

    def test_ground_truth_required(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["--predictions", "p.parquet", "--output", "out.json"])

    def test_output_required(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["--predictions", "p.parquet", "--ground-truth", "gt.tsv"])

    @pytest.mark.parametrize("flag,value", [
        ("--min-cluster-frac", "0"),
        ("--min-cluster-frac", "1.5"),
        ("--min-prediction-frac", "-0.1"),
        ("--min-prediction-frac", "1.1"),
        ("--reciprocal-frac", "0"),
        ("--reciprocal-frac", "2"),
        ("--max-listed-ids", "-1"),
    ])
    def test_out_of_range_values_exit_2(self, flag: str, value: str) -> None:
        with pytest.raises(SystemExit) as exc_info:
            main(["--predictions", "p.parquet", "--ground-truth", "gt.tsv",
                  "--output", "out.json", flag, value])
        assert exc_info.value.code == 2

    def test_min_prediction_frac_zero_is_valid(self, tmp_path: Path) -> None:
        # 0.0 is the default and must not be rejected as out of range.
        p = tmp_path / "p.parquet"
        g = tmp_path / "gt.tsv"
        write_predictions_parquet(p, [PredictedRegion("R1", "chr1", 0, 100, 0.9)])
        write_ground_truth_tsv(g, [KnownCluster("C1", "chr1", 0, 100)])
        main(["--predictions", str(p), "--ground-truth", str(g),
              "--output", str(tmp_path / "o.json"), "--min-prediction-frac", "0"])
        assert (tmp_path / "o.json").exists()


# ────────────────────────────── real baseline output ───────────────────────

class TestRealBaselineOutput:
    """Regression tests over checked-in real converter output.

    All three parquets come from one run of antiSMASH 8.0.4 / DeepBGC / GECCO
    0.10.3 on AL589148.1 (S. coelicolor plasmid SCP1), converted by the three
    scripts in `scripts/`. MiBIG has exactly one coordinate-resolved cluster on
    that contig: BGC0000914 at [231675, 251017).

    All three tools fully contain that cluster. Before this rewrite they scored
    recall 0.002 / 0.000 / 0.002 against the full 430-cluster ground truth.
    """

    GT = FIXTURES / "AL589148_ground_truth.tsv"

    def _run(self, tmp_path: Path, tool: str, **kwargs) -> dict:
        out = tmp_path / f"benchmark_{tool}.json"
        run(EvaluateConfig(
            predictions_path=FIXTURES / f"{tool}_predictions.parquet",
            ground_truth_path=self.GT,
            output_path=out,
            **kwargs,
        ))
        return json.loads(out.read_text())

    @pytest.mark.parametrize("tool", ["antismash", "deepbgc", "gecco"])
    def test_all_three_tools_find_the_cluster(self, tmp_path: Path, tool: str) -> None:
        result = self._run(tmp_path, tool)
        assert result["detection"]["recall"] == 1.0
        assert result["recovered_cluster_ids"] == ["BGC0000914"]
        # Each tool fully contains the cluster.
        assert result["nucleotide"]["recall"] == 1.0
        assert result["boundary"]["median_cluster_coverage"] == pytest.approx(1.0)

    def test_deepbgc_found_it_but_fails_symmetric_match(self, tmp_path: Path) -> None:
        # The case that motivated asymmetric thresholds: a 94 kb region covering
        # 100% of a 19 kb cluster. Detection says found; reciprocal says no.
        result = self._run(tmp_path, "deepbgc")
        assert result["detection"]["recall"] == 1.0
        assert result["reciprocal"]["recall"] == 0.0
        assert result["boundary"]["median_prediction_coverage"] == pytest.approx(
            0.205, abs=5e-4
        )

    @pytest.mark.parametrize("tool", ["antismash", "gecco"])
    def test_tighter_tools_pass_symmetric_match(self, tmp_path: Path, tool: str) -> None:
        result = self._run(tmp_path, tool)
        assert result["reciprocal"]["recall"] == 1.0

    def test_nucleotide_precision_ranks_tools_by_tightness(self, tmp_path: Path) -> None:
        # The real difference between the three: how much extra territory each
        # calls. antiSMASH tightest, GECCO loosest.
        nt = {t: self._run(tmp_path, t)["nucleotide"]["precision"]
              for t in ("antismash", "deepbgc", "gecco")}
        assert nt["antismash"] == pytest.approx(0.382, abs=5e-4)
        assert nt["deepbgc"] == pytest.approx(0.171, abs=5e-4)
        assert nt["gecco"] == pytest.approx(0.128, abs=5e-4)
        assert nt["antismash"] > nt["deepbgc"] > nt["gecco"]

    def test_unmatched_predictions_are_not_called_false(self, tmp_path: Path) -> None:
        # GECCO calls 5 regions; MiBIG knows of 1. The other 4 are unvalidated,
        # not wrong — the schema must not assert otherwise.
        result = self._run(tmp_path, "gecco")
        assert len(result["unmatched_prediction_ids"]) == 4
        assert "false_positive_prediction_ids" not in result
        assert "precision" not in result["detection"]

    def test_scope_is_reported(self, tmp_path: Path) -> None:
        result = self._run(tmp_path, "antismash")
        assert result["scope"]["source"] == "inferred"
        assert result["scope"]["n_contigs"] == 1
        assert result["scope"]["n_clusters_in_scope"] == 1

    def test_min_p_bgc_drops_low_scoring_gecco_calls(self, tmp_path: Path) -> None:
        # GECCO's average_p ranges 0.924-0.992 across its 5 calls.
        result = self._run(tmp_path, "gecco", min_p_bgc=0.93)
        assert result["scope"]["n_predictions_below_threshold"] == 2
        assert result["detection"]["recall"] == 1.0
