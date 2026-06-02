"""Tests for sharp.io — FASTA parser, FASTA writer, parquet writer."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pytest

from sharp.io import (
    ProteinRecord,
    parse_fasta,
    write_embeddings_parquet,
    write_fasta,
)


# ────────────────────────────── parse_fasta ────────────────────────────────

class TestParseFasta:
    def test_basic_records(self, fasta_path: Path) -> None:
        records = list(parse_fasta(fasta_path))
        assert len(records) == 4
        assert records[0] == ProteinRecord("p1", "R001", "MKTAYIAKQRQI")
        assert records[3].region_id == "R002"

    def test_multiline_sequence(self, tmp_path: Path) -> None:
        path = tmp_path / "multiline.faa"
        path.write_text(
            ">p1 region_id=R001\n"
            "MKTAYIAKQRQI\n"
            "SFVKSHFSRQLE\n"
            "EERLGLIEVQ\n"
        )
        records = list(parse_fasta(path))
        assert len(records) == 1
        assert records[0].sequence == "MKTAYIAKQRQISFVKSHFSRQLEEERLGLIEVQ"

    def test_skips_record_missing_region_id(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        path = tmp_path / "missing.faa"
        path.write_text(
            ">p1 region_id=R001\nMKTAY\n"
            ">p2 product=foo\nACDEF\n"          # no region_id
            ">p3 region_id=R002\nGHIKL\n"
        )
        records = list(parse_fasta(path))
        assert [r.protein_id for r in records] == ["p1", "p3"]
        assert "no region_id" in caplog.text

    def test_skips_empty_sequence(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        path = tmp_path / "empty_seq.faa"
        path.write_text(
            ">p1 region_id=R001\nMKTAY\n"
            ">p2 region_id=R001\n"              # no sequence
            ">p3 region_id=R002\nACDEF\n"
        )
        records = list(parse_fasta(path))
        assert [r.protein_id for r in records] == ["p1", "p3"]
        assert "empty sequence" in caplog.text

    def test_ignores_blank_lines(self, tmp_path: Path) -> None:
        path = tmp_path / "blanks.faa"
        path.write_text(
            "\n\n>p1 region_id=R001\n\nMKTAY\n\n>p2 region_id=R002\nACDEF\n\n"
        )
        records = list(parse_fasta(path))
        assert len(records) == 2

    def test_empty_file(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.faa"
        path.write_text("")
        assert list(parse_fasta(path)) == []

    def test_extra_header_tokens_ignored(self, tmp_path: Path) -> None:
        path = tmp_path / "verbose.faa"
        path.write_text(
            ">p1 region_id=R001 product=hypothetical contig=chr1 strand=+\n"
            "MKTAY\n"
        )
        records = list(parse_fasta(path))
        assert records[0].region_id == "R001"
        assert records[0].sequence == "MKTAY"


# ────────────────────────────── write_fasta ────────────────────────────────

class TestWriteFasta:
    def test_round_trip(
        self, tmp_path: Path, sample_records: list[ProteinRecord]
    ) -> None:
        path = tmp_path / "out.faa"
        n = write_fasta(path, sample_records)
        assert n == len(sample_records)

        round_tripped = list(parse_fasta(path))
        assert round_tripped == sample_records

    def test_line_wrapping(self, tmp_path: Path) -> None:
        path = tmp_path / "wrapped.faa"
        long_seq = "A" * 150
        record = ProteinRecord("p1", "R001", long_seq)
        write_fasta(path, [record], line_width=60)

        content = path.read_text().splitlines()
        # header + 60 + 60 + 30
        assert content == [
            ">p1 region_id=R001",
            "A" * 60,
            "A" * 60,
            "A" * 30,
        ]

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        nested = tmp_path / "a" / "b" / "c" / "out.faa"
        assert not nested.parent.exists()
        write_fasta(nested, [ProteinRecord("p1", "R001", "MKTAY")])
        assert nested.exists()


# ────────────────────────────── write_embeddings_parquet ──────────────────

class TestWriteEmbeddingsParquet:
    def test_schema_and_count(
        self, tmp_path: Path, sample_records: list[ProteinRecord]
    ) -> None:
        path = tmp_path / "emb.parquet"
        dim = 16
        vecs = np.random.default_rng(0).standard_normal((4, dim)).astype(np.float32)

        n = write_embeddings_parquet(path, [(sample_records, vecs)], embedding_dim=dim)
        assert n == 4

        table = pq.read_table(path)
        assert table.column_names == ["protein_id", "region_id", "embedding"]
        assert table.num_rows == 4
        # Embedding column is FixedSizeList → each row has exactly `dim` floats
        first = table.column("embedding")[0].as_py()
        assert len(first) == dim

    def test_values_round_trip(
        self, tmp_path: Path, sample_records: list[ProteinRecord]
    ) -> None:
        path = tmp_path / "emb.parquet"
        dim = 8
        vecs = np.arange(4 * dim, dtype=np.float32).reshape(4, dim)

        write_embeddings_parquet(path, [(sample_records, vecs)], embedding_dim=dim)
        table = pq.read_table(path)

        read_back = np.array(table.column("embedding").to_pylist(), dtype=np.float32)
        np.testing.assert_array_equal(read_back, vecs)

        pids_out = table.column("protein_id").to_pylist()
        assert pids_out == [r.protein_id for r in sample_records]

    def test_multiple_batches_streamed(self, tmp_path: Path) -> None:
        path = tmp_path / "emb.parquet"
        dim = 4

        def gen_batches():
            for i in range(3):
                recs = [ProteinRecord(f"p{i}_{j}", f"R{i:03d}", "MKTAY") for j in range(5)]
                vecs = np.full((5, dim), float(i), dtype=np.float32)
                yield recs, vecs

        n = write_embeddings_parquet(path, gen_batches(), embedding_dim=dim)
        assert n == 15

        table = pq.read_table(path)
        assert table.num_rows == 15
        # Verify the three batches kept their values (0, 1, 2 per group of 5)
        emb = np.array(table.column("embedding").to_pylist(), dtype=np.float32)
        assert (emb[0:5] == 0.0).all()
        assert (emb[5:10] == 1.0).all()
        assert (emb[10:15] == 2.0).all()

    def test_dim_mismatch_raises(
        self, tmp_path: Path, sample_records: list[ProteinRecord]
    ) -> None:
        path = tmp_path / "emb.parquet"
        # Vecs have dim=16 but we tell the writer dim=8 → assertion should fire.
        bad_vecs = np.zeros((4, 16), dtype=np.float32)

        with pytest.raises(AssertionError, match="unexpected vec shape"):
            write_embeddings_parquet(path, [(sample_records, bad_vecs)], embedding_dim=8)
