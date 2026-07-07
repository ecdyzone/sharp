"""Tests for scripts.prepare_bgcatlas_ground_truth.

We synthesize antiSMASH-style GenBank records (structured comment carrying
``Orig. start``/``Orig. end``, one ``region`` feature with ``/product``) and
write them to disk with BioPython, so the real parse path is exercised without
touching the 10 GB ``complete-bgcs`` dump. ``--inspect`` on the real data is the
backstop that catches any drift between these fixtures and the actual schema.
"""
from __future__ import annotations

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from prepare_bgcatlas_ground_truth import (  # noqa: E402
    build_ground_truth,
    get_assembly_id,
    get_cluster_id,
    get_contig,
    get_orig_coords,
    get_region_class,
    record_to_cluster,
)
from sharp.io import load_ground_truth_tsv  # noqa: E402


# ────────────────────────────── fixtures ───────────────────────────────────

def make_record(
    contig: str = "contig00011",
    orig_start: str | None = "36256",
    orig_end: str | None = "61780",
    products: list[str] | None = ("thioamitides", "thiopeptide"),
    region_number: str = "1",
    seq_len: int | None = None,
    with_comment: bool = True,
    with_region: bool = True,
) -> SeqRecord:
    """Build an antiSMASH-region-style SeqRecord. Defaults mirror a real file.

    ``seq_len`` defaults to ``orig_end - orig_start`` (the real invariant); pass
    an explicit value to break it if a test needs to.
    """
    if seq_len is None:
        try:
            seq_len = int(orig_end) - int(orig_start)
        except (TypeError, ValueError):
            seq_len = 500
    rec = SeqRecord(Seq("A" * max(seq_len, 1)), id=contig, name=contig, description="")
    rec.annotations["molecule_type"] = "DNA"
    if with_comment:
        data = {"Version": "7.0.0"}
        if orig_start is not None:
            data["Orig. start"] = orig_start
        if orig_end is not None:
            data["Orig. end"] = orig_end
        rec.annotations["structured_comment"] = {"antiSMASH-Data": data}
    if with_region:
        quals: dict = {"region_number": [region_number]}
        if products is not None:
            quals["product"] = list(products)
        rec.features.append(
            SeqFeature(FeatureLocation(0, max(seq_len, 1)), type="region", qualifiers=quals)
        )
    return rec


def write_gbk(directory: Path, filename: str, record: SeqRecord) -> Path:
    """Write a record to ``directory/filename`` as GenBank; return the path."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / filename
    with path.open("w") as fh:
        SeqIO.write(record, fh, "genbank")
    return path


def roundtrip(record: SeqRecord, tmp_path: Path, filename: str = "x.gbk") -> SeqRecord:
    """Write then re-read a record, so tests see it exactly as the script will."""
    path = write_gbk(tmp_path, filename, record)
    return SeqIO.read(path, "genbank")


# ────────────────────────────── accessors ──────────────────────────────────

class TestAccessors:
    def test_assembly_id(self) -> None:
        assert get_assembly_id("MGYA00004361_contig00011.region001.gbk") == "MGYA00004361"

    def test_assembly_id_contig_with_underscore(self) -> None:
        # Assembly = up to the FIRST underscore; underscores in the contig part
        # (e.g. ckmer87_13914) must not confuse it.
        assert get_assembly_id("MGYA00004368_ckmer87_13914.region001.gbk") == "MGYA00004368"

    def test_assembly_id_no_underscore(self) -> None:
        assert get_assembly_id("noseparator.gbk") is None

    def test_cluster_id_is_stem(self) -> None:
        assert (get_cluster_id("MGYA00004361_contig00011.region001.gbk")
                == "MGYA00004361_contig00011.region001")

    def test_cluster_id_distinguishes_regions(self) -> None:
        # region001 vs region002 on the same contig → different cluster ids.
        a = get_cluster_id("MGYA1_NODE-9.region001.gbk")
        b = get_cluster_id("MGYA1_NODE-9.region002.gbk")
        assert a != b

    def test_contig_is_assembly_qualified(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(contig="contig00011"), tmp_path)
        assert get_contig(rec, "MGYA00004361_contig00011.region001.gbk") == "MGYA00004361_contig00011"

    def test_contig_falls_back_to_recid_without_assembly(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(contig="contig00011"), tmp_path)
        assert get_contig(rec, "noseparator.gbk") == "contig00011"

    def test_orig_coords_parsed_as_is(self, tmp_path: Path) -> None:
        # 0-based half-open — returned verbatim, NOT decremented like MiBIG.
        rec = roundtrip(make_record(orig_start="36256", orig_end="61780"), tmp_path)
        assert get_orig_coords(rec) == (36256, 61780)

    def test_orig_coords_missing_comment(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(with_comment=False), tmp_path)
        assert get_orig_coords(rec) == (None, None)

    def test_region_class_joined(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(products=["thioamitides", "thiopeptide"]), tmp_path)
        assert get_region_class(rec) == "thioamitides/thiopeptide"

    def test_region_class_single(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(products=["NRPS"]), tmp_path)
        assert get_region_class(rec) == "NRPS"

    def test_region_class_none_when_no_region(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(with_region=False), tmp_path)
        assert get_region_class(rec) is None


# ────────────────────────────── record_to_cluster ──────────────────────────

class TestRecordToCluster:
    def test_happy_path(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(), tmp_path)
        c = record_to_cluster(rec, "MGYA00004361_contig00011.region001.gbk")
        assert c is not None
        assert c.cluster_id == "MGYA00004361_contig00011.region001"
        assert c.contig == "MGYA00004361_contig00011"
        assert (c.start, c.end) == (36256, 61780)   # no conversion
        assert c.cluster_class == "thioamitides/thiopeptide"

    def test_missing_coords_returns_none(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(with_comment=False), tmp_path)
        assert record_to_cluster(rec, "MGYA1_c.region001.gbk") is None

    def test_inverted_coords_returns_none(self, tmp_path: Path) -> None:
        rec = roundtrip(make_record(orig_start="5000", orig_end="100", seq_len=500), tmp_path)
        assert record_to_cluster(rec, "MGYA1_c.region001.gbk") is None

    def test_degenerate_coords_returns_none(self, tmp_path: Path) -> None:
        # from == to (the MiBIG-style 0,0 placeholder, should any appear here).
        rec = roundtrip(make_record(orig_start="0", orig_end="0", seq_len=500), tmp_path)
        assert record_to_cluster(rec, "MGYA1_c.region001.gbk") is None

    def test_no_class_still_produces_row(self, tmp_path: Path) -> None:
        # Coordinates are what matter; a missing product should not drop the row.
        rec = roundtrip(make_record(products=None), tmp_path)
        c = record_to_cluster(rec, "MGYA1_c.region001.gbk")
        assert c is not None
        assert c.cluster_class is None


# ────────────────────────────── build_ground_truth (I/O) ───────────────────

class TestBuildGroundTruth:
    def test_end_to_end_writes_valid_tsv(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "gbk"
        write_gbk(in_dir, "MGYA1_contigA.region001.gbk",
                  make_record(contig="contigA", orig_start="100", orig_end="600"))
        write_gbk(in_dir, "MGYA2_contigB.region001.gbk",
                  make_record(contig="contigB", orig_start="2000", orig_end="2500"))
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 2

        # Must round-trip through the loader evaluate.py uses.
        clusters = load_ground_truth_tsv(out)
        assert {c.cluster_id for c in clusters} == {
            "MGYA1_contigA.region001", "MGYA2_contigB.region001",
        }
        assert all(c.start < c.end for c in clusters)

    def test_two_regions_same_contig_are_distinct_rows(self, tmp_path: Path) -> None:
        # A single physical contig can carry region001 and region002 in separate
        # files — both must survive as distinct clusters on the same contig.
        in_dir = tmp_path / "gbk"
        write_gbk(in_dir, "MGYA1_NODE-9.region001.gbk",
                  make_record(contig="NODE-9", orig_start="20217", orig_end="65119",
                              region_number="1"))
        write_gbk(in_dir, "MGYA1_NODE-9.region002.gbk",
                  make_record(contig="NODE-9", orig_start="144865", orig_end="198266",
                              region_number="2"))
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 2
        clusters = load_ground_truth_tsv(out)
        assert len({c.cluster_id for c in clusters}) == 2
        # Same contig, different (non-overlapping) coordinates.
        assert {c.contig for c in clusters} == {"MGYA1_NODE-9"}

    def test_limit_processes_subset(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "gbk"
        for i in range(5):
            write_gbk(in_dir, f"MGYA{i}_c{i}.region001.gbk",
                      make_record(contig=f"c{i}", orig_start="100", orig_end="600"))
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out, limit=2)
        assert n == 2

    def test_recursive_search(self, tmp_path: Path) -> None:
        nested = tmp_path / "gbk" / "sub"
        write_gbk(nested, "MGYA1_c.region001.gbk", make_record())
        out = tmp_path / "gt.tsv"
        n = build_ground_truth(tmp_path / "gbk", out)
        assert n == 1

    def test_skips_file_without_coords(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "gbk"
        write_gbk(in_dir, "MGYA1_good.region001.gbk", make_record())
        write_gbk(in_dir, "MGYA2_bad.region001.gbk", make_record(with_comment=False))
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 1   # the coord-less one is dropped

    def test_malformed_file_skipped(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "gbk"
        in_dir.mkdir()
        write_gbk(in_dir, "MGYA1_good.region001.gbk", make_record())
        (in_dir / "MGYA2_bad.region001.gbk").write_text("this is not genbank\n")
        out = tmp_path / "gt.tsv"

        n = build_ground_truth(in_dir, out)
        assert n == 1   # the good one still gets through

    def test_empty_dir_returns_zero(self, tmp_path: Path) -> None:
        in_dir = tmp_path / "gbk"
        in_dir.mkdir()
        n = build_ground_truth(in_dir, tmp_path / "gt.tsv")
        assert n == 0
