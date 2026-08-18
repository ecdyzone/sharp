"""Tests for scripts/merge_predictions.py.

Builds per-genome output trees from the same checked-in real fixtures the
individual converters are tested against, so the merge is exercised over real
tool output without running any tool.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
_FIXTURES = Path(__file__).resolve().parent / "fixtures"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from merge_predictions import (  # noqa: E402
    TOOLS,
    genome_dirs,
    load_contigs,
    merge,
    report_missing,
    run,
)
from sharp.io import load_predictions_parquet  # noqa: E402


def make_tree(root: Path, fixture: str, filename: str, genomes: list[str]) -> Path:
    """One output subdirectory per genome, each holding a copy of a real fixture."""
    root.mkdir(parents=True, exist_ok=True)
    for g in genomes:
        d = root / g
        d.mkdir()
        shutil.copy(_FIXTURES / fixture, d / filename.format(g=g))
    return root


@pytest.fixture
def deepbgc_tree(tmp_path: Path) -> Path:
    return make_tree(tmp_path / "dbgc", "deepbgc_out.bgc.tsv",
                     "{g}.bgc.tsv", ["GEN_A", "GEN_B"])


class TestGenomeDirs:
    def test_lists_subdirectories_sorted(self, deepbgc_tree: Path) -> None:
        assert [d.name for d in genome_dirs(deepbgc_tree)] == ["GEN_A", "GEN_B"]

    def test_ignores_loose_files(self, deepbgc_tree: Path) -> None:
        (deepbgc_tree / "notes.txt").write_text("stray")
        assert [d.name for d in genome_dirs(deepbgc_tree)] == ["GEN_A", "GEN_B"]


class TestMerge:
    def test_concatenates_every_genome(self, deepbgc_tree: Path) -> None:
        regions, counts, failures = merge(deepbgc_tree, "deepbgc")
        # The real fixture holds 5 rows, so two genomes give 10.
        assert counts == {"GEN_A": 5, "GEN_B": 5}
        assert len(regions) == 10
        assert failures == []

    def test_unreadable_genome_reported_not_raised(self, deepbgc_tree: Path) -> None:
        (deepbgc_tree / "GEN_C").mkdir()  # ran, produced nothing
        regions, counts, failures = merge(deepbgc_tree, "deepbgc")
        assert failures == ["GEN_C"]
        # The readable genomes still merge — one bad run must not lose the rest.
        assert len(regions) == 10

    def test_coordinates_preserved_from_converter(self, deepbgc_tree: Path) -> None:
        regions, _, _ = merge(deepbgc_tree, "deepbgc")
        first = min(regions, key=lambda r: r.start)
        # Straight from the real fixture's first row — the merge must not
        # re-derive or shift coordinates.
        assert (first.start, first.end) == (31460, 41750)

    def test_antismash_tree(self, tmp_path: Path) -> None:
        tree = make_tree(tmp_path / "as", "antismash_sequence.json",
                         "{g}.json", ["G1"])
        regions, counts, failures = merge(tree, "antismash")
        assert failures == []
        assert counts["G1"] == len(regions) > 0

    def test_gecco_tree(self, tmp_path: Path) -> None:
        tree = make_tree(tmp_path / "gc", "gecco_sequence.clusters.tsv",
                         "{g}.clusters.tsv", ["G1"])
        regions, counts, failures = merge(tree, "gecco")
        assert failures == []
        assert counts["G1"] == len(regions) > 0

    def test_every_advertised_tool_is_dispatchable(self, tmp_path: Path) -> None:
        # Guards against a --tool choice with no implementation behind it.
        empty = tmp_path / "empty"
        empty.mkdir()
        for tool in TOOLS:
            regions, counts, failures = merge(empty, tool)
            assert (regions, counts, failures) == ([], {}, [])


class TestReportMissing:
    def test_names_genomes_with_no_output(self) -> None:
        assert report_missing({"A": 3}, ["A", "B", "C"]) == ["B", "C"]

    def test_a_genome_that_predicted_nothing_is_not_missing(self) -> None:
        # It ran and found nothing — different from never having run, and only
        # the latter needs reporting.
        assert report_missing({"A": 0}, ["A"]) == []

    def test_empty_scope_reports_nothing(self) -> None:
        assert report_missing({"A": 1}, []) == []


class TestLoadContigs:
    def test_strips_and_skips_blanks(self, tmp_path: Path) -> None:
        p = tmp_path / "contigs.txt"
        p.write_text("AL645882.2\n\n  CP004370.1  \n\n")
        assert load_contigs(p) == ["AL645882.2", "CP004370.1"]


class TestRun:
    def test_writes_a_parquet_evaluate_can_read(
        self, deepbgc_tree: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "merged.parquet"
        run(deepbgc_tree, "deepbgc", out, None, do_inspect=False)
        regions = load_predictions_parquet(out)
        assert len(regions) == 10

    def test_inspect_writes_nothing(
        self, deepbgc_tree: Path, tmp_path: Path, capsys
    ) -> None:
        out = tmp_path / "merged.parquet"
        run(deepbgc_tree, "deepbgc", out, None, do_inspect=True)
        assert not out.exists()
        assert "MERGE PREVIEW" in capsys.readouterr().out

    def test_output_required_unless_inspecting(self, deepbgc_tree: Path) -> None:
        with pytest.raises(SystemExit):
            run(deepbgc_tree, "deepbgc", None, None, do_inspect=False)

    def test_missing_input_dir_exits(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit):
            run(tmp_path / "nope", "deepbgc", tmp_path / "o.parquet", None,
                do_inspect=False)

    def test_missing_genomes_surfaced_via_contigs(
        self, deepbgc_tree: Path, tmp_path: Path, caplog
    ) -> None:
        scope = tmp_path / "contigs.txt"
        scope.write_text("GEN_A\nGEN_B\nGEN_MISSING\n")
        out = tmp_path / "merged.parquet"
        with caplog.at_level("WARNING"):
            run(deepbgc_tree, "deepbgc", out, scope, do_inspect=False)
        assert "GEN_MISSING" in caplog.text
