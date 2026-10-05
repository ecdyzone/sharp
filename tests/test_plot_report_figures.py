"""Tests for scripts/plot_report_figures.py.

Figures are not compared pixel by pixel; what is tested is that every figure
renders from tables build_report_tables.py actually writes, that a placeholder
tool is shown as pending (and stops being pending once numbers are typed in),
and that a missing table skips its figures instead of failing the run.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd
import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import build_report_tables  # noqa: E402
from plot_report_figures import (  # noqa: E402
    FIGURES,
    TOOL_COLORS,
    is_pending,
    main,
    run,
    styles_for,
)
from sharp.io import (  # noqa: E402
    KnownCluster,
    PredictedRegion,
    write_ground_truth_tsv,
    write_predictions_parquet,
)

GT = [
    KnownCluster("BGC1", "c1", 1_000, 11_000, "PKS"),
    KnownCluster("BGC2", "c1", 50_000, 60_000, "NRPS"),
    KnownCluster("BGC3", "c2", 5_000, 15_000, "terpene"),
]
ANTISMASH = [
    PredictedRegion("a1", "c1", 0, 12_000, 1.0, "T1PKS"),
    PredictedRegion("a2", "c2", 0, 20_000, 1.0, "terpene"),
]
DEEPBGC = [
    PredictedRegion("d1", "c1", 2_000, 10_000, 0.9, "Polyketide"),
    PredictedRegion("d2", "c1", 49_000, 61_000, 0.6, None),
    PredictedRegion("d3", "c3", 100, 900, 0.7, "RiPP"),
]


@pytest.fixture(scope="module")
def tables(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Both table sets, written by the real build_report_tables.py."""
    d = tmp_path_factory.mktemp("report")
    gt, ctg = d / "gt.tsv", d / "contigs.txt"
    write_ground_truth_tsv(gt, GT)
    ctg.write_text("c1\nc2\nc3\n")
    a, b = d / "a.parquet", d / "d.parquet"
    write_predictions_parquet(a, ANTISMASH)
    write_predictions_parquet(b, DEEPBGC)
    manifest = d / "genomes.tsv"
    manifest.write_text("assembly\tcontig\tlength\tdescription\n"
                        "GCF_1\tc1\t1000000\tx\nGCF_1\tc2\t100000\tx\n"
                        "GCF_2\tc3\t900000\tx\n")
    assemblies = d / "assemblies.tsv"
    assemblies.write_text("assembly\tn_contigs\ttotal_bp\torganism\tfasta\n"
                          "GCF_1\t2\t1100000\tStreptomyces sp. A\tx\n"
                          "GCF_2\t1\t900000\tNocardia sp. B\tx\n")
    out = d / "tables"
    build_report_tables.main([
        "mibig", "--ground-truth", str(gt), "--contigs", str(ctg),
        "--predictions", f"antiSMASH={a}", f"DeepBGC={b}",
        "--thresholds", "0.5", "0.8", "--output-dir", str(out)])
    build_report_tables.main([
        "raw", "--manifest", str(manifest), "--assemblies", str(assemblies),
        "--predictions", f"antiSMASH={a}", f"DeepBGC={b}",
        "--thresholds", "DeepBGC=0.5,0.8", "--output-dir", str(out)])
    return out


def test_every_figure_renders_in_every_format(tables: Path, tmp_path: Path) -> None:
    written = run(tables, tmp_path, ["png", "svg"], dpi=40)
    assert {p.name for p in written} == {
        f"{name}.{fmt}" for name in FIGURES for fmt in ("png", "svg")}
    assert all(p.stat().st_size > 0 for p in written)


def test_missing_table_skips_only_its_figure(tables: Path, tmp_path: Path) -> None:
    partial = tmp_path / "partial"
    partial.mkdir()
    for f in tables.glob("mibig_*.tsv"):
        (partial / f.name).write_bytes(f.read_bytes())
    written = run(partial, tmp_path / "out", ["png"], dpi=40)
    assert {p.stem for p in written} == {n for n in FIGURES if n.startswith("mibig")}


def test_no_tables_at_all_exits(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        main(["--tables-dir", str(tmp_path), "--output-dir", str(tmp_path / "o")])


def test_only_restricts_the_set(tables: Path, tmp_path: Path) -> None:
    main(["--tables-dir", str(tables), "--output-dir", str(tmp_path),
          "--only", "raw_by_class", "--formats", "png", "--dpi", "40"])
    assert [p.name for p in tmp_path.iterdir()] == ["raw_by_class.png"]


class TestPending:
    def test_placeholder_is_pending(self, tables: Path) -> None:
        df = pd.read_csv(tables / "mibig_metrics.tsv", sep="\t")
        styles = {s.label: s for s in styles_for(df, "tool", ["detection_recall"])}
        assert styles["SHARP"].pending
        assert styles["SHARP"].legend == "SHARP (pending)"
        assert not styles["antiSMASH"].pending

    def test_hand_entered_numbers_end_pending(self, tables: Path, tmp_path: Path) -> None:
        # Someone types S(H)ARP's numbers in but leaves status=placeholder.
        with (tables / "mibig_metrics.tsv").open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for r in rows:
            if r["tool"] == "SHARP":
                r["detection_recall"] = "0.7"
        df = pd.DataFrame(rows)
        assert not is_pending(df[df["tool"] == "SHARP"], ["detection_recall"])

    def test_a_measured_zero_is_not_pending(self) -> None:
        df = pd.DataFrame({"status": ["measured"], "v": [0.0]})
        assert not is_pending(df, ["v"])


class TestColours:
    def test_colour_follows_the_tool_not_its_position(self) -> None:
        a = pd.DataFrame({"tool": ["DeepBGC", "antiSMASH"], "v": [1, 1]})
        b = pd.DataFrame({"tool": ["antiSMASH"], "v": [1]})
        colours_a = {s.label: s.color for s in styles_for(a, "tool", ["v"])}
        colours_b = {s.label: s.color for s in styles_for(b, "tool", ["v"])}
        assert colours_a["antiSMASH"] == colours_b["antiSMASH"] == TOOL_COLORS["antiSMASH"]
        assert colours_a["DeepBGC"] == TOOL_COLORS["DeepBGC"]

    def test_second_cutoff_of_a_tool_is_a_hatched_variant(self) -> None:
        df = pd.DataFrame({"label": ["DeepBGC ≥0.5", "DeepBGC ≥0.8"],
                           "tool": ["DeepBGC", "DeepBGC"], "v": [1, 1]})
        first, second = styles_for(df, "label", ["v"], tool_col="tool")
        assert first.color == TOOL_COLORS["DeepBGC"] and first.hatch is None
        assert second.color != first.color and second.hatch

    def test_unknown_tool_gets_a_colour(self) -> None:
        df = pd.DataFrame({"tool": ["GECCO"], "v": [1]})
        (s,) = styles_for(df, "tool", ["v"])
        assert s.color not in TOOL_COLORS.values()
