"""Tests for scripts/build_genome_manifest.py.

Builds small NCBI-shaped assembly directories on disk rather than mocking the
filesystem: the whole point of the script is that it reads a real layout, and
the layout is the part that can be wrong.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from build_genome_manifest import (  # noqa: E402
    assembly_key,
    find_fasta,
    index_database,
    link_genomes,
    load_contig_index,
    load_genomes_tsv,
    organism_of,
    parse_fasta_index,
    run,
)


def write_fasta(path: Path, records: list[tuple[str, str, int]]) -> None:
    """records: (contig, description, length) — sequence is filler of that length."""
    lines = []
    for contig, desc, length in records:
        lines.append(f">{contig} {desc}")
        seq = "ACGT" * (length // 4) + "A" * (length % 4)
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def db(tmp_path: Path) -> Path:
    """Two assemblies, one multi-contig, shaped like the NCBI dump."""
    root = tmp_path / "db"
    a = root / "GCF_000009765.2_ASM976v2"
    a.mkdir(parents=True)
    write_fasta(a / "GCF_000009765.2_ASM976v2_genomic.fna", [
        ("NC_003155.5", "Streptomyces avermitilis MA-4680, complete sequence", 8000),
        ("NC_004719.1", "Streptomyces avermitilis MA-4680 plasmid SAP1, complete sequence", 400),
    ])
    # Files the script must ignore.
    (a / "GCF_000009765.2_ASM976v2_protein.faa").write_text(">p1\nMA\n")
    (a / "GCF_000009765.2_ASM976v2_genomic.gff").write_text("##gff-version 3\n")

    b = root / "GCF_041549525.1_P_aurescens_B2879"
    b.mkdir(parents=True)
    write_fasta(b / "GCF_041549525.1_P_aurescens_B2879_genomic.fna", [
        ("NZ_CP157456.1", "Paenarthrobacter aurescens strain NRRL B-2879 chromosome, complete genome", 4000),
    ])
    return root


class TestAssemblyKey:
    def test_strips_the_assembly_name_suffix(self) -> None:
        assert assembly_key("GCF_000009765.2_ASM976v2") == "GCF_000009765.2"

    def test_handles_genbank_accessions(self) -> None:
        assert assembly_key("GCA_000203835.1_ASM20383v1") == "GCA_000203835.1"

    def test_non_ncbi_directory_keeps_its_whole_name(self) -> None:
        # Better indexed under an odd key than dropped silently.
        assert assembly_key("my_custom_genome") == "my_custom_genome"


class TestFindFasta:
    def test_prefers_genomic_fna(self, db: Path) -> None:
        got = find_fasta(db / "GCF_000009765.2_ASM976v2")
        assert got is not None and got.name.endswith("_genomic.fna")

    def test_ignores_the_proteome(self, db: Path) -> None:
        got = find_fasta(db / "GCF_000009765.2_ASM976v2")
        assert got is not None and not got.name.endswith(".faa")

    def test_skips_cds_and_rna_companions(self, tmp_path: Path) -> None:
        d = tmp_path / "GCF_1.1_x"
        d.mkdir()
        write_fasta(d / "GCF_1.1_x_cds_from_genomic.fna", [("c", "d", 40)])
        write_fasta(d / "GCF_1.1_x_genomic.fna", [("NZ_1.1", "d", 40)])
        got = find_fasta(d)
        assert got is not None and "cds_from" not in got.name

    def test_returns_none_when_absent(self, tmp_path: Path) -> None:
        d = tmp_path / "empty"
        d.mkdir()
        assert find_fasta(d) is None


class TestParseFastaIndex:
    def test_lengths_exclude_headers_and_newlines(self, db: Path) -> None:
        recs = parse_fasta_index(
            db / "GCF_000009765.2_ASM976v2" / "GCF_000009765.2_ASM976v2_genomic.fna",
            "GCF_000009765.2",
        )
        assert [r.length for r in recs] == [8000, 400]

    def test_contig_is_the_first_header_token(self, db: Path) -> None:
        recs = parse_fasta_index(
            db / "GCF_000009765.2_ASM976v2" / "GCF_000009765.2_ASM976v2_genomic.fna",
            "GCF_000009765.2",
        )
        assert [r.contig for r in recs] == ["NC_003155.5", "NC_004719.1"]

    def test_version_suffix_is_preserved(self, db: Path) -> None:
        # The whole join depends on this: NC_003155 and NC_003155.5 are
        # different strings and only one of them matches tool output.
        recs = parse_fasta_index(
            db / "GCF_000009765.2_ASM976v2" / "GCF_000009765.2_ASM976v2_genomic.fna",
            "GCF_000009765.2",
        )
        assert recs[0].contig.endswith(".5")

    def test_reads_gzipped_fasta(self, tmp_path: Path) -> None:
        plain = tmp_path / "x.fna"
        write_fasta(plain, [("NZ_9.1", "desc", 80)])
        gz = tmp_path / "x.fna.gz"
        gz.write_bytes(gzip.compress(plain.read_bytes()))
        recs = parse_fasta_index(gz, "GCF_9.1")
        assert len(recs) == 1 and recs[0].length == 80


class TestOrganismOf:
    def test_uses_the_longest_contig_not_the_plasmid(self, db: Path) -> None:
        recs = parse_fasta_index(
            db / "GCF_000009765.2_ASM976v2" / "GCF_000009765.2_ASM976v2_genomic.fna",
            "GCF_000009765.2",
        )
        assert "plasmid" not in organism_of(recs).lower()

    def test_strips_trailing_descriptors(self, db: Path) -> None:
        recs = parse_fasta_index(
            db / "GCF_041549525.1_P_aurescens_B2879"
            / "GCF_041549525.1_P_aurescens_B2879_genomic.fna",
            "GCF_041549525.1",
        )
        assert organism_of(recs) == "Paenarthrobacter aurescens strain NRRL B-2879"

    def test_empty_input_is_empty_label(self) -> None:
        assert organism_of([]) == ""


class TestIndexDatabase:
    def test_finds_every_assembly(self, db: Path) -> None:
        assemblies, _, _ = index_database(db)
        assert [a.assembly for a in assemblies] == [
            "GCF_000009765.2", "GCF_041549525.1"]

    def test_total_bp_sums_all_contigs(self, db: Path) -> None:
        assemblies, _, _ = index_database(db)
        assert assemblies[0].total_bp == 8400

    def test_contigs_carry_their_assembly(self, db: Path) -> None:
        _, contigs, _ = index_database(db)
        assert {c.contig: c.assembly for c in contigs} == {
            "NC_003155.5": "GCF_000009765.2",
            "NC_004719.1": "GCF_000009765.2",
            "NZ_CP157456.1": "GCF_041549525.1",
        }

    def test_limit_truncates(self, db: Path) -> None:
        assemblies, _, _ = index_database(db, limit=1)
        assert len(assemblies) == 1

    def test_directory_without_fasta_is_reported_not_raised(self, db: Path) -> None:
        (db / "GCF_999.1_broken").mkdir()
        assemblies, _, problems = index_database(db)
        assert len(assemblies) == 2
        assert any("GCF_999.1_broken" in p for p in problems)

    def test_duplicate_assembly_key_is_reported(self, db: Path) -> None:
        dup = db / "GCF_000009765.2_OTHERNAME"
        dup.mkdir()
        write_fasta(dup / "GCF_000009765.2_OTHERNAME_genomic.fna",
                    [("NZ_X.1", "d", 100)])
        assemblies, _, problems = index_database(db)
        assert len(assemblies) == 2
        assert any("duplicate assembly key" in p for p in problems)


class TestRun:
    def test_writes_all_four_files(self, db: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(db, out, None, None, False)
        for name in ("genomes.tsv", "assemblies.tsv", "assemblies.txt",
                     "analyzed_contigs.txt"):
            assert (out / name).is_file(), name

    def test_assemblies_txt_is_the_array_line_list(self, db: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(db, out, None, None, False)
        lines = (out / "assemblies.txt").read_text().split()
        assert lines == ["GCF_000009765.2", "GCF_041549525.1"]

    def test_analyzed_contigs_is_the_scope_file(self, db: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(db, out, None, None, False)
        lines = (out / "analyzed_contigs.txt").read_text().split()
        assert lines == ["NC_003155.5", "NC_004719.1", "NZ_CP157456.1"]

    def test_manifest_round_trips(self, db: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(db, out, None, None, False)
        assert load_genomes_tsv(out / "genomes.tsv") == {
            "GCF_000009765.2": ["NC_003155.5", "NC_004719.1"],
            "GCF_041549525.1": ["NZ_CP157456.1"],
        }

    def test_contig_index_round_trips(self, db: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        run(db, out, None, None, False)
        assert load_contig_index(out / "genomes.tsv")["NC_004719.1"] == (
            "GCF_000009765.2", 400)

    def test_output_dir_required_unless_inspecting(self, db: Path) -> None:
        with pytest.raises(SystemExit):
            run(db, None, None, None, False)

    def test_inspect_writes_nothing(self, db: Path, tmp_path: Path,
                                    capsys: pytest.CaptureFixture) -> None:
        out = tmp_path / "out"
        run(db, out, None, 5, True)
        assert not out.exists()
        assert "DATABASE LAYOUT" in capsys.readouterr().out

    def test_missing_db_dir_exits(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit):
            run(tmp_path / "nope", tmp_path / "out", None, None, False)


class TestLinkGenomes:
    def test_links_are_named_after_the_assembly(self, db: Path, tmp_path: Path) -> None:
        assemblies, _, _ = index_database(db)
        links = tmp_path / "links"
        n, skipped = link_genomes(links, assemblies)
        assert n == 2 and not skipped
        assert (links / "GCF_000009765.2.fasta").is_symlink()

    def test_link_resolves_to_the_real_fasta(self, db: Path, tmp_path: Path) -> None:
        assemblies, _, _ = index_database(db)
        links = tmp_path / "links"
        link_genomes(links, assemblies)
        target = (links / "GCF_041549525.1.fasta").resolve()
        assert target.name.endswith("_genomic.fna") and target.is_file()

    def test_rerun_replaces_an_existing_link(self, db: Path, tmp_path: Path) -> None:
        assemblies, _, _ = index_database(db)
        links = tmp_path / "links"
        link_genomes(links, assemblies)
        n, _ = link_genomes(links, assemblies)
        assert n == 2

    def test_gzipped_genome_is_refused_not_linked(self, tmp_path: Path) -> None:
        # A .gz link would fail thousands of array tasks in; report instead.
        d = tmp_path / "db" / "GCF_5.1_x"
        d.mkdir(parents=True)
        plain = tmp_path / "tmp.fna"
        write_fasta(plain, [("NZ_5.1", "d", 80)])
        (d / "GCF_5.1_x_genomic.fna.gz").write_bytes(gzip.compress(plain.read_bytes()))
        assemblies, _, _ = index_database(tmp_path / "db")
        n, skipped = link_genomes(tmp_path / "links", assemblies)
        assert n == 0 and len(skipped) == 1 and "gzipped" in skipped[0]
