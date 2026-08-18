"""Tests for scripts/select_benchmark_genomes.py.

Everything here runs offline. The only network-touching function is
`fetch_record_lengths`, which is exercised through the cache path instead —
the selection rules themselves are pure and are what can silently corrupt a
benchmark, so that is where the tests are aimed.

Fixtures use the real accessions and lengths that motivated each rule, so a
failure names the actual record that would break.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from select_benchmark_genomes import (  # noqa: E402
    DEFAULT_MIN_LENGTH,
    RecordInfo,
    caption_of,
    group_into_genomes,
    infer_lower_bounds,
    normalize_ground_truth,
    organism_key,
    primary_accession,
    read_length_cache,
    select_genomes,
    write_benchmark_set,
    write_contigs_file,
    write_length_cache,
)
from sharp.io import KnownCluster  # noqa: E402


def kc(cluster_id: str, contig: str, start: int, end: int) -> KnownCluster:
    return KnownCluster(
        cluster_id=cluster_id, contig=contig, start=start, end=end,
        cluster_class="NRPS",
    )


# The real S. coelicolor twin pair: same sequence, two accessions, clusters
# filed against both.
COELICOLOR = {
    "AL645882": RecordInfo("AL645882.2", 8_667_507,
                           "Streptomyces coelicolor A3(2) complete genome"),
    "NC_003888": RecordInfo("NC_003888.3", 8_667_507,
                            "Streptomyces coelicolor A3(2), complete sequence"),
}


class TestCaptionOf:
    def test_strips_version(self) -> None:
        assert caption_of("AL645882.2") == "AL645882"

    def test_unversioned_passes_through(self) -> None:
        assert caption_of("CP114200") == "CP114200"

    def test_refseq_underscore_preserved(self) -> None:
        assert caption_of("NZ_DS999644.1") == "NZ_DS999644"


class TestOrganismKey:
    def test_refseq_and_genbank_titles_agree(self) -> None:
        # The two S. coelicolor titles differ after the strain, which is
        # exactly why only the first three words are used.
        a = organism_key("Streptomyces coelicolor A3(2) complete genome")
        b = organism_key("Streptomyces coelicolor A3(2), complete sequence")
        assert a == b

    def test_scaffold_wording_ignored(self) -> None:
        a = organism_key("Streptomyces sp. CS113 genomic scaffold scaffold00001")
        b = organism_key("Streptomyces sp. CS113 scaffold00001")
        assert a == b

    def test_different_strains_stay_distinct(self) -> None:
        assert organism_key("Streptomyces sp. CS113 scaffold00001") != \
               organism_key("Streptomyces sp. CS147 scaffold00001")


class TestPrimaryAccession:
    def test_prefers_insdc_over_refseq(self) -> None:
        assert primary_accession(("AL645882.2", "NC_003888.3")) == "AL645882.2"
        assert primary_accession(("NZ_CM001889.1", "CM001889.1")) == "CM001889.1"

    def test_all_refseq_falls_back_to_refseq(self) -> None:
        assert primary_accession(("NZ_DS999644.1",)) == "NZ_DS999644.1"

    def test_deterministic_on_ties(self) -> None:
        assert primary_accession(("CP000002.1", "CP000001.1")) == "CP000001.1"


class TestGroupIntoGenomes:
    def test_twins_merge_into_one_genome(self) -> None:
        clusters = [
            kc("BGC0000038", "AL645882.2", 100, 50_000),
            kc("BGC0000910", "NC_003888.3", 2_943_455, 2_944_875),
        ]
        genomes, unresolved = group_into_genomes(clusters, COELICOLOR)
        assert unresolved == {}
        assert len(genomes) == 1
        g = genomes[0]
        assert g.accession == "AL645882.2"
        assert g.accessions == ("AL645882.2", "NC_003888.3")
        # The whole point: both clusters survive on one genome.
        assert g.n_clusters == 2

    def test_same_length_different_organism_not_merged(self) -> None:
        info = {
            "AAA111": RecordInfo("AAA111.1", 8_000_000, "Streptomyces alpha one x"),
            "BBB222": RecordInfo("BBB222.1", 8_000_000, "Streptomyces beta two y"),
        }
        clusters = [kc("B1", "AAA111.1", 0, 10_000), kc("B2", "BBB222.1", 0, 10_000)]
        genomes, _ = group_into_genomes(clusters, info)
        assert len(genomes) == 2

    def test_same_organism_different_length_not_merged(self) -> None:
        # S. clavuligerus chromosome and its plasmid pSCL4 share an organism
        # but are separate replicons.
        info = {
            "CM000913": RecordInfo("CM000913.1", 6_760_392,
                                   "Streptomyces clavuligerus ATCC 27064 chromosome"),
            "CM000914": RecordInfo("CM000914.1", 1_796_500,
                                   "Streptomyces clavuligerus ATCC 27064 plasmid pSCL4"),
        }
        clusters = [kc("C1", "CM000913.1", 0, 10_000), kc("C2", "CM000914.1", 0, 10_000)]
        genomes, _ = group_into_genomes(clusters, info)
        assert len(genomes) == 2

    def test_unresolved_contig_reported_not_dropped_silently(self) -> None:
        clusters = [kc("X1", "NOSUCH.1", 0, 10_000)]
        genomes, unresolved = group_into_genomes(clusters, COELICOLOR)
        assert genomes == []
        assert "NOSUCH.1" in unresolved["no length from NCBI"]

    def test_ground_truth_without_version_still_matches(self) -> None:
        # Ground truth has bare `CP114200`; NCBI reports `CP114200.1`.
        info = {"CP114200": RecordInfo("CP114200.1", 8_379_354,
                                       "Streptomyces sp. 71268 chromosome")}
        genomes, unresolved = group_into_genomes([kc("G", "CP114200", 0, 9_000)], info)
        assert unresolved == {}
        assert genomes[0].accession == "CP114200.1"


class TestInferLowerBounds:
    def test_bound_is_max_cluster_end(self) -> None:
        clusters = [
            kc("A", "JAFMOF010000002.1", 1_721_958, 1_725_178),
            kc("B", "JAFMOF010000002.1", 10, 5_000),
        ]
        bounds = infer_lower_bounds(clusters, {})
        assert bounds["JAFMOF010000002"].length == 1_725_178
        assert bounds["JAFMOF010000002"].length_known is False

    def test_resolved_accessions_are_left_alone(self) -> None:
        clusters = [kc("A", "AL645882.2", 0, 50_000)]
        assert infer_lower_bounds(clusters, COELICOLOR) == {}

    def test_lower_bound_records_are_never_merged(self) -> None:
        # Two unresolved records that happen to share a bound must not be
        # treated as the same sequence — twin detection needs an exact length.
        clusters = [
            kc("A", "JAAVVL010000001.1", 0, 4_318_177),
            kc("B", "JARAKF010000001.1", 0, 4_318_177),
        ]
        info = infer_lower_bounds(clusters, {})
        genomes, _ = group_into_genomes(clusters, info)
        assert len(genomes) == 2

    def test_short_deposit_stays_below_threshold(self) -> None:
        # The bound can only under-estimate, so a BGC-only deposit can never
        # be promoted into the genome-scale set by this fallback.
        clusters = [kc("A", "MG742725.1", 0, 119_206)]
        info = infer_lower_bounds(clusters, {})
        genomes, _ = group_into_genomes(clusters, info)
        selected, rejected = select_genomes(genomes, top_n=0)
        assert selected == []
        assert len(rejected) == 1


class TestSelectGenomes:
    def test_below_min_length_rejected(self) -> None:
        info = {
            "BIG": RecordInfo("BIG.1", 8_000_000, "Streptomyces big one a"),
            "SMALL": RecordInfo("SMALL.1", 119_206, "Streptomyces small one b"),
        }
        clusters = [kc("A", "BIG.1", 0, 40_000), kc("B", "SMALL.1", 0, 119_206)]
        selected, rejected = self._select(info, clusters, top_n=0)
        assert [g.accession for g in selected] == ["BIG.1"]
        assert [g.accession for g in rejected] == ["SMALL.1"]

    @staticmethod
    def _select(info, clusters, **kw):
        genomes = group_into_genomes(clusters, info)[0]
        return select_genomes(genomes, **kw)

    def test_ranked_by_cluster_count_first(self) -> None:
        info = {
            "MANY": RecordInfo("MANY.1", 6_000_000, "Streptomyces many one a"),
            "HUGE": RecordInfo("HUGE.1", 11_000_000, "Streptomyces huge one b"),
        }
        clusters = [
            kc("A", "MANY.1", 0, 10_000), kc("B", "MANY.1", 20_000, 30_000),
            kc("C", "HUGE.1", 0, 10_000),
        ]
        selected, _ = self._select(info, clusters, top_n=0)
        # Two clusters beats a bigger genome with one.
        assert [g.accession for g in selected] == ["MANY.1", "HUGE.1"]

    def test_length_breaks_cluster_count_ties(self) -> None:
        info = {
            "SHORT": RecordInfo("SHORT.1", 2_000_000, "Streptomyces short one a"),
            "LONG": RecordInfo("LONG.1", 9_000_000, "Streptomyces long one b"),
        }
        clusters = [kc("A", "SHORT.1", 0, 10_000), kc("B", "LONG.1", 0, 10_000)]
        selected, _ = self._select(info, clusters, top_n=0)
        assert [g.accession for g in selected] == ["LONG.1", "SHORT.1"]

    def test_top_n_truncates(self) -> None:
        info = {
            f"G{i}": RecordInfo(f"G{i}.1", 9_000_000 - i, f"Streptomyces g{i} strain")
            for i in range(5)
        }
        clusters = [kc(f"C{i}", f"G{i}.1", 0, 10_000) for i in range(5)]
        selected, _ = self._select(info, clusters, top_n=2)
        assert len(selected) == 2

    def test_top_n_zero_keeps_everything(self) -> None:
        info = {
            f"G{i}": RecordInfo(f"G{i}.1", 9_000_000 - i, f"Streptomyces g{i} strain")
            for i in range(5)
        }
        clusters = [kc(f"C{i}", f"G{i}.1", 0, 10_000) for i in range(5)]
        selected, _ = self._select(info, clusters, top_n=0)
        assert len(selected) == 5

    def test_default_threshold_keeps_smallest_real_replicon(self) -> None:
        # pSCL4 (1.79 Mb) and the S. ambofaciens arm (1.37 Mb) must survive.
        assert 1_367_117 >= DEFAULT_MIN_LENGTH
        assert 1_796_500 >= DEFAULT_MIN_LENGTH


class TestNormalizeGroundTruth:
    def test_contigs_rewritten_to_primary_accession(self) -> None:
        clusters = [
            kc("BGC0000038", "AL645882.2", 100, 50_000),
            kc("BGC0000910", "NC_003888.3", 2_943_455, 2_944_875),
        ]
        genomes, _ = group_into_genomes(clusters, COELICOLOR)
        out = normalize_ground_truth(genomes)
        assert {c.contig for c in out} == {"AL645882.2"}
        # Nothing lost — this is the cluster that would otherwise vanish.
        assert {c.cluster_id for c in out} == {"BGC0000038", "BGC0000910"}

    def test_coordinates_unchanged(self) -> None:
        clusters = [kc("BGC0000910", "NC_003888.3", 2_943_455, 2_944_875)]
        genomes, _ = group_into_genomes(clusters, COELICOLOR)
        out = normalize_ground_truth(genomes)
        assert (out[0].start, out[0].end) == (2_943_455, 2_944_875)

    def test_class_preserved(self) -> None:
        genomes, _ = group_into_genomes(
            [kc("B", "AL645882.2", 0, 10_000)], COELICOLOR)
        assert normalize_ground_truth(genomes)[0].cluster_class == "NRPS"

    def test_sorted_by_contig_then_start(self) -> None:
        info = {
            "AAA": RecordInfo("AAA.1", 8_000_000, "Streptomyces aaa strain x"),
            "BBB": RecordInfo("BBB.1", 8_000_000, "Streptomyces bbb strain y"),
        }
        clusters = [
            kc("C3", "BBB.1", 500, 1_000),
            kc("C2", "AAA.1", 9_000, 10_000),
            kc("C1", "AAA.1", 100, 200),
        ]
        genomes, _ = group_into_genomes(clusters, info)
        out = normalize_ground_truth(genomes)
        assert [c.cluster_id for c in out] == ["C1", "C2", "C3"]


class TestLengthCacheRoundTrip:
    def test_write_then_read_preserves_lookup_key(self, tmp_path: Path) -> None:
        # Regression: the cache is keyed by unversioned caption while storing
        # the versioned accession. Writing one and reading the other back made
        # every rerun miss the cache and re-query NCBI.
        path = tmp_path / "lengths.tsv"
        write_length_cache(path, COELICOLOR)
        back = read_length_cache(path)
        assert set(back) == set(COELICOLOR)
        assert back["AL645882"].accession == "AL645882.2"
        assert back["AL645882"].length == 8_667_507

    def test_missing_file_is_empty_not_an_error(self, tmp_path: Path) -> None:
        assert read_length_cache(tmp_path / "nope.tsv") == {}

    def test_malformed_rows_skipped(self, tmp_path: Path) -> None:
        path = tmp_path / "lengths.tsv"
        path.write_text(
            "caption\taccession\tlength\ttitle\n"
            "GOOD\tGOOD.1\t1000\tStreptomyces good strain\n"
            "BAD\tBAD.1\tnot-a-number\tStreptomyces bad strain\n"
            "SHORT\tonly-two-columns\n"
        )
        assert set(read_length_cache(path)) == {"GOOD"}


class TestOutputFiles:
    @pytest.fixture
    def selected(self):
        clusters = [
            kc("BGC0000038", "AL645882.2", 100, 50_000),
            kc("BGC0000910", "NC_003888.3", 2_943_455, 2_944_875),
        ]
        return group_into_genomes(clusters, COELICOLOR)[0]

    def test_contigs_file_lists_primary_accessions_only(
        self, tmp_path: Path, selected
    ) -> None:
        path = tmp_path / "analyzed_contigs.txt"
        write_contigs_file(path, selected)
        lines = path.read_text().split()
        # Only the accession we actually fetch and analyze — listing the twin
        # too would put a cluster in scope that no prediction can ever match.
        assert lines == ["AL645882.2"]

    def test_benchmark_set_records_merge_provenance(
        self, tmp_path: Path, selected
    ) -> None:
        path = tmp_path / "benchmark_genomes.tsv"
        write_benchmark_set(path, selected)
        rows = path.read_text().strip().split("\n")
        header, row = rows[0].split("\t"), rows[1].split("\t")
        assert header[:6] == ["rank", "accession", "all_accessions", "length",
                              "length_known", "n_clusters"]
        assert row[1] == "AL645882.2"
        assert row[2] == "AL645882.2,NC_003888.3"
        assert row[5] == "2"

    def test_creates_missing_output_dir(self, tmp_path: Path, selected) -> None:
        path = tmp_path / "nested" / "deep" / "out.tsv"
        write_benchmark_set(path, selected)
        assert path.exists()
