"""Tests for scripts/summarize_predictions.py.

The interesting cases are the ones that make raw counts misleading: a tool that
splits what another merges, a threshold that moves the count, and a contig name
that does not join.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from summarize_predictions import (  # noqa: E402
    bp_called_per_contig,
    check_contigs,
    parse_spec,
    run,
    summarize,
    write_assembly_tsv,
    write_tool_tsv,
)
from sharp.io import PredictedRegion, write_predictions_parquet  # noqa: E402

# Two assemblies: one with a chromosome + plasmid, one with a chromosome.
INDEX = {
    "NC_003155.5": ("GCF_000009765.2", 1_000_000),
    "NC_004719.1": ("GCF_000009765.2", 100_000),
    "NZ_CP157456.1": ("GCF_041549525.1", 900_000),
}


def region(contig: str, start: int, end: int, p: float = 1.0) -> PredictedRegion:
    return PredictedRegion(f"{contig}_{start}", contig, start, end, p, "x")


@pytest.fixture
def manifest(tmp_path: Path) -> Path:
    path = tmp_path / "genomes.tsv"
    path.write_text(
        "assembly\tcontig\tlength\tdescription\n"
        + "".join(f"{a}\t{c}\t{n}\tdesc\n" for c, (a, n) in INDEX.items())
    )
    return path


class TestParseSpec:
    def test_splits_name_and_path(self) -> None:
        assert parse_spec("antismash=a/b.parquet") == ("antismash", Path("a/b.parquet"))

    def test_path_may_contain_equals(self) -> None:
        assert parse_spec("t=a=b.parquet")[1] == Path("a=b.parquet")

    @pytest.mark.parametrize("bad", ["nopath", "=path.parquet", "name="])
    def test_malformed_spec_exits(self, bad: str) -> None:
        with pytest.raises(SystemExit):
            parse_spec(bad)


class TestBpCalledPerContig:
    def test_disjoint_regions_add_up(self) -> None:
        got = bp_called_per_contig([region("NC_003155.5", 0, 100),
                                    region("NC_003155.5", 200, 350)])
        assert got == {"NC_003155.5": 250}

    def test_overlapping_regions_are_counted_once(self) -> None:
        # This is what makes a splitter comparable to a merger.
        got = bp_called_per_contig([region("NC_003155.5", 0, 100),
                                    region("NC_003155.5", 50, 150)])
        assert got == {"NC_003155.5": 150}

    def test_contigs_stay_separate(self) -> None:
        got = bp_called_per_contig([region("NC_003155.5", 0, 100),
                                    region("NZ_CP157456.1", 0, 100)])
        assert got == {"NC_003155.5": 100, "NZ_CP157456.1": 100}

    def test_no_regions_is_empty(self) -> None:
        assert bp_called_per_contig([]) == {}


class TestSummarize:
    def test_counts_and_bp_over_the_whole_database(self) -> None:
        overall, _ = summarize("t", [region("NC_003155.5", 0, 40_000)], INDEX, 0.0)
        assert overall.n_regions == 1
        assert overall.bp_called == 40_000
        assert overall.total_bp == 2_000_000

    def test_denominator_includes_genomes_with_no_calls(self) -> None:
        # A tool that finds nothing in a genome must still be charged for it,
        # or calls-per-Mb rewards silence.
        overall, _ = summarize("t", [region("NC_003155.5", 0, 40_000)], INDEX, 0.0)
        assert overall.regions_per_mb == pytest.approx(1 / 2.0)

    def test_threshold_filters_regions(self) -> None:
        regions = [region("NC_003155.5", 0, 100, 0.9),
                   region("NC_003155.5", 500, 600, 0.3)]
        assert summarize("t", regions, INDEX, 0.5)[0].n_regions == 1

    def test_a_splitter_and_a_merger_agree_on_bp_not_on_count(self) -> None:
        merger = [region("NC_003155.5", 0, 40_000)]
        splitter = [region("NC_003155.5", 0, 20_000),
                    region("NC_003155.5", 20_000, 40_000)]
        a, _ = summarize("merger", merger, INDEX, 0.0)
        b, _ = summarize("splitter", splitter, INDEX, 0.0)
        assert a.bp_called == b.bp_called
        assert a.n_regions != b.n_regions

    def test_per_assembly_rows_cover_every_assembly(self) -> None:
        _, rows = summarize("t", [region("NC_003155.5", 0, 100)], INDEX, 0.0)
        assert {r.assembly for r in rows} == {"GCF_000009765.2", "GCF_041549525.1"}

    def test_assembly_bp_sums_its_contigs(self) -> None:
        _, rows = summarize("t", [], INDEX, 0.0)
        by = {r.assembly: r.genome_bp for r in rows}
        assert by == {"GCF_000009765.2": 1_100_000, "GCF_041549525.1": 900_000}

    def test_plasmid_calls_roll_up_to_their_assembly(self) -> None:
        _, rows = summarize("t", [region("NC_004719.1", 0, 5_000)], INDEX, 0.0)
        by = {r.assembly: r.n_regions for r in rows}
        assert by["GCF_000009765.2"] == 1

    def test_n_assemblies_called_counts_only_hits(self) -> None:
        overall, _ = summarize("t", [region("NC_003155.5", 0, 100)], INDEX, 0.0)
        assert overall.n_assemblies_called == 1

    def test_median_region_bp(self) -> None:
        regions = [region("NC_003155.5", 0, 100), region("NC_003155.5", 500, 800),
                   region("NC_003155.5", 1000, 1500)]
        assert summarize("t", regions, INDEX, 0.0)[0].median_region_bp == 300

    def test_empty_predictions_do_not_divide_by_zero(self) -> None:
        overall, _ = summarize("t", [], INDEX, 0.0)
        assert overall.n_regions == 0 and overall.frac_bp_called == 0.0


class TestCheckContigs:
    def test_known_contigs_pass_through(self) -> None:
        regions = [region("NC_003155.5", 0, 100)]
        assert check_contigs("t", regions, INDEX) == regions

    def test_unknown_contig_is_dropped(self) -> None:
        regions = [region("NC_003155.5", 0, 100), region("SOMETHING_ELSE", 0, 100)]
        assert len(check_contigs("t", regions, INDEX)) == 1

    def test_version_suffix_mismatch_is_caught(self, caplog) -> None:
        # The failure this exists for: a tool that reports NC_003155 rather
        # than NC_003155.5 joins to nothing and looks like it found nothing.
        got = check_contigs("t", [region("NC_003155", 0, 100)], INDEX)
        assert got == []
        assert "NOTHING joined" in caplog.text


class TestWriters:
    def test_tool_tsv_has_a_header_and_a_row(self, tmp_path: Path) -> None:
        overall, _ = summarize("t", [region("NC_003155.5", 0, 100)], INDEX, 0.0)
        path = tmp_path / "summary_by_tool.tsv"
        write_tool_tsv(path, [overall])
        lines = path.read_text().splitlines()
        assert lines[0].startswith("tool\tthreshold")
        assert lines[1].split("\t")[0] == "t"

    def test_assembly_tsv_row_per_assembly(self, tmp_path: Path) -> None:
        _, rows = summarize("t", [], INDEX, 0.0)
        path = tmp_path / "summary_by_assembly.tsv"
        assert write_assembly_tsv(path, rows) == 2


class TestRun:
    @pytest.fixture
    def preds(self, tmp_path: Path) -> Path:
        path = tmp_path / "tool.parquet"
        write_predictions_parquet(path, [
            region("NC_003155.5", 0, 40_000, 0.9),
            region("NZ_CP157456.1", 0, 10_000, 0.2),
        ])
        return path

    def test_writes_both_tables(self, manifest: Path, preds: Path,
                                tmp_path: Path) -> None:
        out = tmp_path / "cmp"
        run(manifest, [f"tool={preds}"], [0.0], out)
        assert (out / "summary_by_tool.tsv").is_file()
        assert (out / "summary_by_assembly.tsv").is_file()

    def test_one_row_per_tool_and_threshold(self, manifest: Path, preds: Path,
                                            tmp_path: Path) -> None:
        out = tmp_path / "cmp"
        run(manifest, [f"a={preds}", f"b={preds}"], [0.0, 0.5], out)
        rows = (out / "summary_by_tool.tsv").read_text().splitlines()[1:]
        assert len(rows) == 4

    def test_threshold_changes_the_count(self, manifest: Path, preds: Path,
                                         tmp_path: Path) -> None:
        out = tmp_path / "cmp"
        run(manifest, [f"a={preds}"], [0.0, 0.5], out)
        rows = [r.split("\t") for r in
                (out / "summary_by_tool.tsv").read_text().splitlines()[1:]]
        assert rows[0][2] == "2" and rows[1][2] == "1"

    def test_missing_predictions_file_exits(self, manifest: Path,
                                            tmp_path: Path) -> None:
        with pytest.raises(SystemExit):
            run(manifest, [f"a={tmp_path / 'nope.parquet'}"], [0.0], None)

    def test_no_output_dir_still_prints(self, manifest: Path, preds: Path,
                                        capsys: pytest.CaptureFixture) -> None:
        run(manifest, [f"a={preds}"], [0.0], None)
        assert "RAW COMPARISON" in capsys.readouterr().out
