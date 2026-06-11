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
