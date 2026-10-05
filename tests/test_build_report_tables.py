"""Tests for scripts/build_report_tables.py.

What matters here is that a figure can never disagree with the result file it
illustrates: the tables must reproduce evaluate_predictions / summarize exactly,
per-class counts must add up to the headline, and a placeholder tool must be
visibly a placeholder rather than a tool that found nothing.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from build_report_tables import (  # noqa: E402
    CLASS_ORDER,
    MEASURED,
    PLACEHOLDER,
    MibigSettings,
    Series,
    agreement_rows,
    gt_by_class,
    genus_of,
    load_tools,
    main,
    make_series,
    mibig_class_rows,
    normalize_class,
    parse_thresholds,
    series_label,
    sweep_rows,
)
from summarize_predictions import summarize  # noqa: E402
from sharp.io import (  # noqa: E402
    KnownCluster,
    PredictedRegion,
    write_ground_truth_tsv,
    write_predictions_parquet,
)
from sharp.metrics import evaluate_predictions  # noqa: E402


def region(contig: str, start: int, end: int, p: float = 1.0,
           cls: str | None = None) -> PredictedRegion:
    return PredictedRegion(f"{contig}_{start}", contig, start, end, p, cls)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ── a small benchmark: 2 contigs, 5 clusters over 4 classes ─────────────────

GT = [
    KnownCluster("BGC1", "c1", 1_000, 11_000, "PKS"),
    KnownCluster("BGC2", "c1", 50_000, 60_000, "NRPS"),
    KnownCluster("BGC3", "c1", 100_000, 110_000, "ribosomal"),
    KnownCluster("BGC4", "c2", 5_000, 15_000, "NRPS/PKS"),
    KnownCluster("BGC5", "c2", 40_000, 50_000, "terpene"),
    KnownCluster("BGC6", "c9", 0, 10_000, "PKS"),          # out of scope
]
SCOPE = {"c1", "c2"}

# Unscored tool: finds BGC1, BGC2, BGC4 (wide), misses the rest, plus 1 extra.
ANTISMASH = [
    region("c1", 0, 12_000, cls="T1PKS"),
    region("c1", 48_000, 62_000, cls="NRPS;NRPS-like"),
    region("c2", 0, 40_000, cls="T1PKS;NRPS"),
    region("c2", 200_000, 220_000, cls="terpene"),
]
# Scored tool: finds BGC1 (p .9), BGC3 (p .6), BGC5 (p .55); one extra at .7.
DEEPBGC = [
    region("c1", 2_000, 10_000, p=0.9, cls="Polyketide"),
    region("c1", 101_000, 109_000, p=0.6, cls="RiPP"),
    region("c2", 41_000, 49_000, p=0.55, cls=None),
    region("c2", 300_000, 305_000, p=0.7, cls="NRP-Polyketide"),
]


class TestNormalizeClass:
    @pytest.mark.parametrize("raw,expected", [
        ("T1PKS", "PKS"),
        ("NRPS-like", "NRPS"),                  # antiSMASH name containing "-"
        ("lanthipeptide-class-iii", "RiPP"),
        ("terpene-precursor", "Terpene"),
        ("NI-siderophore", "Other"),
        ("NRPS;NRPS-like", "NRPS"),             # two names, one category
        ("T1PKS;NRPS", "Hybrid"),
        ("terpene;butyrolactone", "Hybrid"),
        ("Polyketide", "PKS"),                  # DeepBGC
        ("NRP-Polyketide", "Hybrid"),
        ("Alkaloid", "Other"),
        ("Terpene", "Terpene"),
        ("ribosomal", "RiPP"),                  # MiBIG
        ("NRPS/PKS", "Hybrid"),
        ("saccharide", "Saccharide"),
        (None, "Unclassified"),
        ("", "Unclassified"),
        (" ", "Unclassified"),
    ])
    def test_known_vocabularies(self, raw: str | None, expected: str) -> None:
        assert normalize_class(raw) == expected

    def test_unknown_name_is_other_and_counted(self) -> None:
        unknown: Counter = Counter()
        assert normalize_class("brand-new-rule", unknown) == "Other"
        assert normalize_class("T1PKS;brand-new-rule", unknown) == "Hybrid"
        assert unknown == Counter({"brand-new-rule": 2})

    def test_every_result_is_a_display_class(self) -> None:
        for raw in ["T1PKS", "Polyketide", "ribosomal", None, "x-y-z"]:
            assert normalize_class(raw) in CLASS_ORDER


class TestGenusOf:
    @pytest.mark.parametrize("organism,expected", [
        ("Streptomyces avermitilis MA-4680 = NBRC 14893", "Streptomyces"),
        ("Candidatus Nanopelagicus limnes isolate MMS-21-122",
         "Candidatus Nanopelagicus"),
        # NCBI brackets a genus the species no longer belongs to — keep it apart.
        ("[Mycobacterium] stephanolepidis", "[Mycobacterium]"),
        ("", "Unknown"),
    ])
    def test_genus(self, organism: str, expected: str) -> None:
        assert genus_of(organism) == expected


class TestThresholdsAndSeries:
    def test_parse_thresholds(self) -> None:
        assert parse_thresholds(["DeepBGC=0.8,0.5"]) == {"DeepBGC": [0.5, 0.8]}

    @pytest.mark.parametrize("bad", ["DeepBGC", "DeepBGC=", "=0.5", "D=x"])
    def test_malformed_thresholds_exit(self, bad: str) -> None:
        with pytest.raises(SystemExit):
            parse_thresholds([bad])

    def test_label_names_the_cutoff_only_when_it_is_a_choice(self) -> None:
        assert series_label("DeepBGC", 0.0, 1) == "DeepBGC"
        assert series_label("DeepBGC", 0.8, 1) == "DeepBGC ≥0.8"
        assert series_label("DeepBGC", 0.5, 2) == "DeepBGC ≥0.5"

    def test_scored_tool_splits_and_filters(self) -> None:
        s = make_series("DeepBGC", DEEPBGC, [0.5, 0.8])
        assert [x.label for x in s] == ["DeepBGC ≥0.5", "DeepBGC ≥0.8"]
        assert [len(x.regions) for x in s] == [4, 1]

    def test_unscored_tool_is_one_series(self) -> None:
        s = make_series("antiSMASH", ANTISMASH, [0.5, 0.8])
        assert [(x.label, len(x.regions)) for x in s] == [("antiSMASH", 4)]


class TestLoadTools:
    def test_placeholder_dropped_once_it_has_data(self, tmp_path: Path) -> None:
        path = tmp_path / "s.parquet"
        write_predictions_parquet(path, DEEPBGC)
        tools, pending = load_tools([f"SHARP={path}"], ["SHARP", "Other"])
        assert [n for n, _ in tools] == ["SHARP"]
        assert pending == ["Other"]

    def test_duplicate_names_exit(self, tmp_path: Path) -> None:
        path = tmp_path / "s.parquet"
        write_predictions_parquet(path, DEEPBGC)
        with pytest.raises(SystemExit):
            load_tools([f"A={path}", f"A={path}"], [])


class TestMibigTables:
    s = MibigSettings()

    def test_class_counts_sum_to_headline(self) -> None:
        in_scope = [c for c in GT if c.contig in SCOPE]
        classes = gt_by_class(in_scope)
        for preds in (ANTISMASH, DEEPBGC):
            rows = mibig_class_rows("t", preds, classes, SCOPE, self.s)
            headline = evaluate_predictions(preds, GT, scope=SCOPE).detection
            assert sum(r["n_recovered"] for r in rows) == headline.n_recovered
            assert sum(r["n_clusters"] for r in rows) == headline.n_clusters

    def test_classes_follow_display_order(self) -> None:
        classes = gt_by_class([c for c in GT if c.contig in SCOPE])
        assert list(classes) == ["PKS", "NRPS", "RiPP", "Terpene", "Hybrid"]

    def test_sweep_is_monotone_for_a_scored_tool(self) -> None:
        rows = sweep_rows("DeepBGC", DEEPBGC, GT, SCOPE, self.s,
                          [0.5, 0.6, 0.8, 0.95])
        recalls = [r["detection_recall"] for r in rows]
        assert recalls == sorted(recalls, reverse=True)
        assert [r["n_predictions"] for r in rows] == [4, 3, 1, 0]

    def test_sweep_is_flat_for_an_unscored_tool(self) -> None:
        rows = sweep_rows("antiSMASH", ANTISMASH, GT, SCOPE, self.s, [0.5, 0.9])
        assert rows[0]["detection_recall"] == rows[1]["detection_recall"]

    def test_end_to_end(self, tmp_path: Path) -> None:
        gt, ctg = tmp_path / "gt.tsv", tmp_path / "contigs.txt"
        write_ground_truth_tsv(gt, GT)
        ctg.write_text("c1\nc2\n")
        a, d = tmp_path / "a.parquet", tmp_path / "d.parquet"
        write_predictions_parquet(a, ANTISMASH)
        write_predictions_parquet(d, DEEPBGC)
        out = tmp_path / "tables"
        main(["mibig", "--ground-truth", str(gt), "--contigs", str(ctg),
              "--predictions", f"antiSMASH={a}", f"DeepBGC={d}",
              "--output-dir", str(out)])

        metrics = {r["tool"]: r for r in read_tsv(out / "mibig_metrics.tsv")}
        assert list(metrics) == ["antiSMASH", "DeepBGC", "SHARP"]
        expected = evaluate_predictions(ANTISMASH, GT, scope=SCOPE)
        assert int(metrics["antiSMASH"]["n_recovered"]) == expected.detection.n_recovered
        assert float(metrics["antiSMASH"]["detection_recall"]) == pytest.approx(
            expected.detection.recall, rel=1e-5)
        assert int(metrics["antiSMASH"]["n_clusters"]) == 5   # BGC6 out of scope

        sharp = metrics["SHARP"]
        assert sharp["status"] == PLACEHOLDER
        assert float(sharp["detection_recall"]) == 0
        assert int(sharp["n_clusters"]) == 5               # GT is not the tool's

        by_class = read_tsv(out / "mibig_by_class.tsv")
        assert {r["status"] for r in by_class if r["tool"] == "SHARP"} == {PLACEHOLDER}
        assert (out / "mibig_threshold_sweep.tsv").is_file()


# ── raw comparison: two assemblies, three contigs ──────────────────────────

INDEX = {
    "c1": ("GCF_1", 1_000_000),
    "c2": ("GCF_1", 100_000),
    "c3": ("GCF_2", 900_000),
}


class TestRawTables:
    def test_agreement_with_itself_is_total(self) -> None:
        a = Series("A", "A", 0.0, ANTISMASH)
        b = Series("B", "B", 0.0, list(ANTISMASH))     # same calls, new tool
        rows = {(r["reference"], r["query"]): r
                for r in agreement_rows([a, b], set(INDEX))}
        assert rows[("A", "B")]["region_agreement"] == 1.0
        assert rows[("A", "B")]["bp_agreement"] == 1.0

    def test_agreement_is_directional(self) -> None:
        wide = Series("wide", "wide", 0.0, [region("c1", 0, 100_000)])
        narrow = Series("narrow", "narrow", 0.0, [region("c1", 0, 10_000)])
        rows = {(r["reference"], r["query"]): r
                for r in agreement_rows([wide, narrow], set(INDEX))}
        # narrow covers 10% of wide; wide covers all of narrow
        assert rows[("wide", "narrow")]["bp_agreement"] == pytest.approx(0.1)
        assert rows[("narrow", "wide")]["bp_agreement"] == 1.0

    def test_placeholder_agreement_is_blank_not_zero(self) -> None:
        a = Series("A", "A", 0.0, ANTISMASH)
        p = Series("SHARP", "SHARP", 0.0, None)
        rows = agreement_rows([a, p], set(INDEX))
        blanks = [r for r in rows if "SHARP" in (r["reference"], r["query"])]
        assert len(blanks) == 3
        assert all(r["status"] == PLACEHOLDER and r["bp_agreement"] is None
                   for r in blanks)

    def test_end_to_end(self, tmp_path: Path) -> None:
        manifest = tmp_path / "genomes.tsv"
        manifest.write_text("assembly\tcontig\tlength\tdescription\n" + "".join(
            f"{a}\t{c}\t{n}\td\n" for c, (a, n) in INDEX.items()))
        assemblies = tmp_path / "assemblies.tsv"
        assemblies.write_text(
            "assembly\tn_contigs\ttotal_bp\torganism\tfasta\n"
            "GCF_1\t2\t1100000\tStreptomyces coelicolor A3(2)\tx\n"
            "GCF_2\t1\t900000\tNocardia farcinica\tx\n")
        a, d = tmp_path / "a.parquet", tmp_path / "d.parquet"
        write_predictions_parquet(a, ANTISMASH)
        write_predictions_parquet(d, DEEPBGC)
        out = tmp_path / "tables"
        main(["raw", "--manifest", str(manifest), "--assemblies", str(assemblies),
              "--predictions", f"antiSMASH={a}", f"DeepBGC={d}",
              "--thresholds", "DeepBGC=0.5,0.8", "--output-dir", str(out)])

        totals = {r["label"]: r for r in read_tsv(out / "raw_totals.tsv")}
        assert list(totals) == ["antiSMASH", "DeepBGC ≥0.5", "DeepBGC ≥0.8", "SHARP"]
        overall, _ = summarize("antiSMASH", ANTISMASH, INDEX, 0.0)
        assert int(totals["antiSMASH"]["n_regions"]) == overall.n_regions
        assert int(totals["antiSMASH"]["bp_called"]) == overall.bp_called
        assert int(totals["DeepBGC ≥0.8"]["n_regions"]) == 1
        assert totals["SHARP"]["status"] == PLACEHOLDER
        assert int(totals["SHARP"]["total_bp"]) == 2_000_000

        by_class = read_tsv(out / "raw_by_class.tsv")
        deep = {r["class"]: int(r["n_regions"]) for r in by_class
                if r["label"] == "DeepBGC ≥0.5"}
        assert deep == {"PKS": 1, "NRPS": 0, "RiPP": 1, "Terpene": 0,
                        "Saccharide": 0, "Other": 0, "Hybrid": 1,
                        "Unclassified": 1}

        genus = [r for r in read_tsv(out / "raw_by_genus.tsv")
                 if r["label"] == "antiSMASH"]
        assert {r["genus"] for r in genus} == {"Streptomyces", "Nocardia"}
        strep = next(r for r in genus if r["genus"] == "Streptomyces")
        assert float(strep["mean_regions_per_genome"]) == 4.0   # all 4 on GCF_1

        agreement = read_tsv(out / "raw_agreement.tsv")
        assert len(agreement) == 16                         # 4 series squared
        assert {r["status"] for r in agreement if r["query"] == "SHARP"} == {PLACEHOLDER}
        assert {r["status"] for r in agreement
                if "SHARP" not in (r["reference"], r["query"])} == {MEASURED}
