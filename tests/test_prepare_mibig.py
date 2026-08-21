"""Tests for scripts.prepare_mibig_ground_truth.

We construct synthetic MiBIG-shaped JSON for both the 4.0 (flat) and 3.x
(nested-under-"cluster") layouts, so the parser's defensive fallbacks are
actually exercised. If the real 4.0 schema differs from what's encoded here,
`--inspect` on the real data is what catches it — these tests pin the
behavior we've designed for.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from prepare_mibig_ground_truth import (  # noqa: E402
    build_ground_truth,
    entry_to_clusters,
    get_classes,
    get_cluster_id,
    get_loci,
    get_locus_coords,
    get_taxonomy_name,
    rejection_reason,
)
from sharp.io import load_ground_truth_tsv  # noqa: E402


# ────────────────────────────── fixtures ───────────────────────────────────

def entry_40() -> dict:
    """A MiBIG 4.0-style (flat) entry."""
    return {
        "accession": "BGC0000001",
        "biosynthesis": {"classes": [{"class": "PKS"}, {"class": "NRPS"}]},
        "taxonomy": {"name": "Streptomyces coelicolor A3(2)"},
        "loci": [
            {"accession": "AB000001.1", "location": {"from": 1, "to": 41001}},
        ],
    }


def entry_3x() -> dict:
    """A MiBIG 3.x-style (nested) entry."""
    return {
        "cluster": {
            "mibig_accession": "BGC0000999",
            "biosyn_class": ["Polyketide"],
            "organism_name": "Amycolatopsis mediterranei",
            "loci": {
                "accession": "CP000010.1",
                "start_coord": 100,
                "end_coord": 5100,
            },
        }
    }


def entry_multi_locus() -> dict:
    return {
        "accession": "BGC0002000",
        "biosynthesis": {"classes": [{"class": "Terpene"}]},
        "taxonomy": {"name": "Streptomyces avermitilis"},
        "loci": [
            {"accession": "C1", "location": {"from": 1, "to": 1000}},
            {"accession": "C2", "location": {"from": 2000, "to": 3000}},
        ],
    }


# ────────────────────────────── accessors ──────────────────────────────────

class TestAccessors:
    def test_cluster_id_40(self) -> None:
        assert get_cluster_id(entry_40(), "fallback") == "BGC0000001"

    def test_cluster_id_3x(self) -> None:
        assert get_cluster_id(entry_3x(), "fallback") == "BGC0000999"

    def test_cluster_id_fallback(self) -> None:
        assert get_cluster_id({}, "BGC_from_filename") == "BGC_from_filename"

    def test_classes_40(self) -> None:
        assert get_classes(entry_40()) == ["PKS", "NRPS"]

    def test_classes_3x(self) -> None:
        assert get_classes(entry_3x()) == ["Polyketide"]

    def test_classes_missing(self) -> None:
        assert get_classes({}) == []

    def test_taxonomy_40(self) -> None:
        assert get_taxonomy_name(entry_40()) == "Streptomyces coelicolor A3(2)"

    def test_taxonomy_3x(self) -> None:
        assert get_taxonomy_name(entry_3x()) == "Amycolatopsis mediterranei"

    def test_loci_list_40(self) -> None:
        assert len(get_loci(entry_40())) == 1

    def test_loci_dict_3x_normalized_to_list(self) -> None:
        # 3.x single-object loci should be wrapped into a list.
        loci = get_loci(entry_3x())
        assert isinstance(loci, list)
        assert len(loci) == 1

    def test_locus_coords_40(self) -> None:
        locus = entry_40()["loci"][0]
        assert get_locus_coords(locus) == ("AB000001.1", 1, 41001)

    def test_locus_coords_3x(self) -> None:
        locus = get_loci(entry_3x())[0]
        assert get_locus_coords(locus) == ("CP000010.1", 100, 5100)

    def test_locus_coords_missing(self) -> None:
        assert get_locus_coords({}) == (None, None, None)


# ────────────────────────────── entry_to_clusters ──────────────────────────

class TestEntryToClusters:
    def test_40_single_locus(self) -> None:
        clusters = entry_to_clusters(entry_40(), "BGC0000001.json")
        assert len(clusters) == 1
        c = clusters[0]
        assert c.cluster_id == "BGC0000001"
        assert c.contig == "AB000001.1"
        # 1-based inclusive [1, 41001] → 0-based half-open [0, 41001)
        assert c.start == 0
        assert c.end == 41001
        assert c.cluster_class == "PKS/NRPS"

    def test_3x_single_locus(self) -> None:
        clusters = entry_to_clusters(entry_3x(), "BGC0000999.json")
        assert len(clusters) == 1
        c = clusters[0]
        assert c.cluster_id == "BGC0000999"
        # [100, 5100] → [99, 5100)
        assert c.start == 99
        assert c.end == 5100
        assert c.cluster_class == "Polyketide"

    def test_multi_locus_gets_suffixed_ids(self) -> None:
        clusters = entry_to_clusters(entry_multi_locus(), "BGC0002000.json")
        assert len(clusters) == 2
        assert clusters[0].cluster_id == "BGC0002000.1"
        assert clusters[1].cluster_id == "BGC0002000.2"

    def test_skips_locus_without_coords(self) -> None:
        entry = {
            "accession": "BGC0003000",
            "loci": [
                {"accession": "C1", "location": {"from": 1, "to": 1000}},
                {"accession": "C2"},   # no location → skipped
            ],
        }
        clusters = entry_to_clusters(entry, "x.json")
        assert len(clusters) == 1
        # only one valid locus, so it should NOT be suffixed... but the entry
        # had 2 loci entries, so multi=True and it IS suffixed. Pin that.
        assert clusters[0].cluster_id == "BGC0003000.1"

    def test_empty_entry_produces_nothing(self) -> None:
        assert entry_to_clusters({}, "x.json") == []

    def test_inverted_coords_skipped(self) -> None:
        entry = {
            "accession": "BGC0004000",
            "loci": [{"accession": "C1", "location": {"from": 5000, "to": 100}}],
        }
        assert entry_to_clusters(entry, "x.json") == []


# ────────────────────────────── build_ground_truth (I/O) ───────────────────

class TestBuildGroundTruth:
    def _write_entries(self, d: Path, entries: dict[str, dict]) -> None:
        for name, entry in entries.items():
            (d / name).write_text(json.dumps(entry))

    def test_end_to_end_writes_valid_tsv(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        self._write_entries(in_dir, {
            "BGC0000001.json": entry_40(),
            "BGC0000999.json": entry_3x(),
        })
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 2

        # The output must be readable by the same loader evaluate.py uses.
        clusters = load_ground_truth_tsv(out)
        ids = {c.cluster_id for c in clusters}
        assert ids == {"BGC0000001", "BGC0000999"}

    def test_genus_filter(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        self._write_entries(in_dir, {
            "strep.json": entry_40(),              # Streptomyces
            "amyco.json": entry_3x(),              # Amycolatopsis
            "strep2.json": entry_multi_locus(),    # Streptomyces
        })
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out, genus="Streptomyces")
        # entry_40 → 1 row, entry_multi_locus → 2 rows, amyco filtered out
        assert n == 3
        clusters = load_ground_truth_tsv(out)
        assert all("BGC0000999" not in c.cluster_id for c in clusters)

    def test_exclude_eukaryotes_drops_fungal_entries(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        fungal = entry_40()
        fungal["accession"] = "BGC0000777"
        fungal["taxonomy"] = {"name": "Aspergillus nidulans FGSC A4",
                              "ncbiTaxId": 227321}
        self._write_entries(in_dir, {
            "strep.json": entry_40(),      # Streptomyces — bacterial, kept
            "amyco.json": entry_3x(),      # Amycolatopsis — bacterial, kept
            "asp.json": fungal,            # Aspergillus — dropped
        })
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out, exclude_eukaryotes=True)
        assert n == 2
        clusters = load_ground_truth_tsv(out)
        assert "BGC0000777" not in {c.cluster_id for c in clusters}

    def test_exclude_eukaryotes_off_by_default(self, tmp_path: Path) -> None:
        # The deny-list must not apply unless explicitly requested.
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        fungal = entry_40()
        fungal["accession"] = "BGC0000777"
        fungal["taxonomy"] = {"name": "Aspergillus nidulans", "ncbiTaxId": 162425}
        self._write_entries(in_dir, {"asp.json": fungal})
        out = tmp_path / "gt.tsv"

        assert build_ground_truth(in_dir, out) == 1

    def test_exclude_eukaryotes_matches_genus_not_substring(
        self, tmp_path: Path
    ) -> None:
        # The deny-list is checked against the first word only. A bacterium
        # whose species epithet happens to contain a listed genus name must
        # survive — this is why the check is not a substring test.
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        tricky = entry_40()
        tricky["accession"] = "BGC0000778"
        tricky["taxonomy"] = {"name": "Streptomyces mucor-like sp. XY",
                              "ncbiTaxId": 1}
        self._write_entries(in_dir, {"tricky.json": tricky})
        out = tmp_path / "gt.tsv"

        assert build_ground_truth(in_dir, out, exclude_eukaryotes=True) == 1

    def test_recursive_search(self, tmp_path: Path) -> None:
        # Files may be nested in subdirectories after tar extraction.
        in_dir = tmp_path / "json"
        nested = in_dir / "subdir"
        nested.mkdir(parents=True)
        (nested / "BGC0000001.json").write_text(json.dumps(entry_40()))
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 1

    def test_malformed_json_skipped(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        (in_dir / "good.json").write_text(json.dumps(entry_40()))
        (in_dir / "bad.json").write_text("{not valid json")
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 1   # the good one still gets through

    def test_empty_dir_returns_zero(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        n = build_ground_truth(in_dir, tmp_path / "gt.tsv")
        assert n == 0


# ══════════════════════════════════════════════════════════════════════════
# Accession / span validation (added 2026-08-18 after the real 4.0 dump was
# found to contain 11 loci no coordinate benchmark can score).
# ══════════════════════════════════════════════════════════════════════════

def _entry(cluster_id: str, contig: str, start: int, end: int) -> dict:
    """Minimal 4.0-shaped entry with one locus at the given coordinates."""
    return {
        "accession": cluster_id,
        "biosynthesis": {"classes": [{"class": "NRPS"}]},
        "taxonomy": {"name": "Streptomyces coelicolor A3(2)"},
        "loci": [{"accession": contig, "location": {"from": start, "to": end}}],
    }


class TestRejectionReason:
    """Pure validation rules — the real offenders from MiBIG 4.0 by name."""

    @pytest.mark.parametrize("contig", [
        "AL645882.2",       # EMBL chromosome
        "AP009493.1",       # DDBJ — starts with 'AP' but is NOT a protein
        "NZ_DS999644.1",    # RefSeq scaffold
        "JJOB01000001.1",   # WGS contig 1 (not the master record)
        "CP114200",         # no version suffix, still usable
    ])
    def test_accepts_real_nucleotide_accessions(self, contig: str) -> None:
        assert rejection_reason(contig, 0, 50_000) is None

    @pytest.mark.parametrize("contig,expected", [
        ("WP_071967254.1", "protein accession"),
        ("NP_123456.1", "protein accession"),
        ("GCA_028752555", "assembly accession"),
        ("GCF_000203835.1", "assembly accession"),
        ("ASM2282770v1", "assembly accession"),
        ("NZ_JOGD01000000", "WGS master accession"),
        ("NZ_LLZK01000000", "WGS master accession"),
    ])
    def test_rejects_unusable_accessions(self, contig: str, expected: str) -> None:
        assert rejection_reason(contig, 0, 50_000) == expected

    def test_rejects_empty_accession(self) -> None:
        assert rejection_reason("   ", 0, 50_000) == "empty accession"

    def test_rejects_degenerate_span(self) -> None:
        reason = rejection_reason("CP021118.1", 965_781, 965_969)  # 188 bp, real
        assert reason is not None and "degenerate span" in reason

    def test_accepts_smallest_real_mibig_locus(self) -> None:
        # BGC0000848 on AP009493.1 is 945 bp — the shortest genuine locus in
        # the Streptomyces set. The threshold must not reject it.
        assert rejection_reason("AP009493.1", 8_273_443, 8_274_388) is None

    def test_wgs_master_rejected_but_contig_one_accepted(self) -> None:
        # The distinction that matters: ...01000000 addresses the whole WGS
        # project, ...01000001 is a real sequence.
        assert rejection_reason("NZ_JOGD01000000", 0, 50_000) is not None
        assert rejection_reason("NZ_JOGD01000001", 0, 50_000) is None


class TestBuildGroundTruthFiltering:
    """The validation and dedup rules as applied by build_ground_truth."""

    @staticmethod
    def _write(in_dir: Path, entries: dict) -> None:
        for name, entry in entries.items():
            (in_dir / name).write_text(json.dumps(entry))

    def test_invalid_accessions_are_dropped(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        self._write(in_dir, {
            "BGC0000001.json": _entry("BGC0000001", "AL645882.2", 1, 50_000),
            "BGC0003020.json": _entry("BGC0003020", "WP_071967254.1", 1, 245),
            "BGC0003030.json": _entry("BGC0003030", "NZ_JOGD01000000", 1, 50_000),
            "BGC0003064.json": _entry("BGC0003064", "GCA_000719695.1", 1, 3),
        })
        out = tmp_path / "gt.tsv"
        assert build_ground_truth(in_dir, out) == 1
        assert {c.cluster_id for c in load_ground_truth_tsv(out)} == {"BGC0000001"}

    def test_duplicate_loci_collapsed_keeping_first_id(self, tmp_path: Path) -> None:
        # Real case: BGC0002850 and BGC0002868 are both OR050662.1:0-58983.
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        self._write(in_dir, {
            "BGC0002850.json": _entry("BGC0002850", "OR050662.1", 1, 58_983),
            "BGC0002868.json": _entry("BGC0002868", "OR050662.1", 1, 58_983),
        })
        out = tmp_path / "gt.tsv"
        assert build_ground_truth(in_dir, out) == 1
        # Sorted filename order means the lower id wins — deterministic.
        assert {c.cluster_id for c in load_ground_truth_tsv(out)} == {"BGC0002850"}

    def test_same_contig_different_interval_is_not_a_duplicate(
        self, tmp_path: Path
    ) -> None:
        in_dir = tmp_path / "json"
        in_dir.mkdir()
        self._write(in_dir, {
            "BGC0000001.json": _entry("BGC0000001", "AL645882.2", 1, 50_000),
            "BGC0000002.json": _entry("BGC0000002", "AL645882.2", 60_000, 90_000),
        })
        out = tmp_path / "gt.tsv"
        assert build_ground_truth(in_dir, out) == 2
