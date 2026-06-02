"""Tests for sharp.extract_embeddings — batching, orchestration, CLI.

The orchestration test (`test_run`) uses monkeypatching to swap in a stub
Embedder. NOTE: we patch the symbols on `sharp.extract_embeddings`, not on
`sharp.model_management` — because `extract_embeddings` does
`from sharp.model_management import Embedder, ensure_model_available`,
the names are bound in the importing module's namespace and that's where
the patch has to land.
"""
from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest

from sharp.config import EmbeddingConfig
from sharp.extract_embeddings import (
    build_parser,
    length_bucketed_batches,
    run,
)
from sharp.io import ProteinRecord


# ────────────────────────────── length_bucketed_batches ────────────────────

class TestLengthBucketedBatches:
    def _records_with_lengths(self, lengths: list[int]) -> list[ProteinRecord]:
        return [ProteinRecord(f"p{i}", "R001", "A" * L) for i, L in enumerate(lengths)]

    def test_sorted_within_each_batch(self) -> None:
        # Deliberately out-of-order input.
        records = self._records_with_lengths([300, 100, 500, 200, 400, 150, 350])
        batches = list(length_bucketed_batches(records, batch_size=3))

        for batch in batches:
            lens = [r.length for r in batch]
            assert lens == sorted(lens)

    def test_globally_sorted_across_batches(self) -> None:
        # Length-bucketing should also be monotonic across batch boundaries.
        records = self._records_with_lengths([300, 100, 500, 200, 400])
        batches = list(length_bucketed_batches(records, batch_size=2))
        flat = [r.length for batch in batches for r in batch]
        assert flat == sorted(flat)

    def test_count_preserved(self) -> None:
        records = self._records_with_lengths(list(range(50, 73)))  # 23 records
        batches = list(length_bucketed_batches(records, batch_size=5))
        assert sum(len(b) for b in batches) == 23
        assert len(batches) == 5  # 5+5+5+5+3
        assert len(batches[-1]) == 3

    def test_empty_input(self) -> None:
        assert list(length_bucketed_batches([], batch_size=4)) == []

    def test_exact_multiple(self) -> None:
        records = self._records_with_lengths([100] * 12)
        batches = list(length_bucketed_batches(records, batch_size=4))
        assert len(batches) == 3
        assert all(len(b) == 4 for b in batches)


# ────────────────────────────── run() ──────────────────────────────────────

class TestRun:
    def test_orchestration_end_to_end(
        self,
        tmp_path: Path,
        fasta_path: Path,
        stub_embedder,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """run() should: read FASTA, embed every record, write a valid parquet."""
        # Patch on the importing module (see docstring at top).
        monkeypatch.setattr(
            "sharp.extract_embeddings.Embedder",
            lambda *args, **kwargs: stub_embedder,
        )
        monkeypatch.setattr(
            "sharp.extract_embeddings.ensure_model_available",
            lambda *args, **kwargs: "stub",
        )

        output = tmp_path / "embeddings.parquet"
        run(EmbeddingConfig(
            input_path=fasta_path,
            output_path=output,
            batch_size=2,
            log_every=10,
        ))

        assert output.exists()
        table = pq.read_table(output)
        assert table.num_rows == 4
        assert table.column_names == ["protein_id", "region_id", "embedding"]
        assert len(table.column("embedding")[0].as_py()) == stub_embedder.dim

        # All input proteins ended up in the output (regardless of batch order).
        assert set(table.column("protein_id").to_pylist()) == {"p1", "p2", "p3", "p4"}

    def test_exits_on_empty_input(
        self,
        tmp_path: Path,
        stub_embedder,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        empty = tmp_path / "empty.faa"
        empty.write_text("")
        monkeypatch.setattr(
            "sharp.extract_embeddings.Embedder",
            lambda *a, **kw: stub_embedder,
        )
        monkeypatch.setattr(
            "sharp.extract_embeddings.ensure_model_available",
            lambda *a, **kw: "stub",
        )

        with pytest.raises(SystemExit) as exc_info:
            run(EmbeddingConfig(
                input_path=empty,
                output_path=tmp_path / "out.parquet",
            ))
        assert exc_info.value.code == 1

    def test_creates_output_parent_dirs(
        self,
        tmp_path: Path,
        fasta_path: Path,
        stub_embedder,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.setattr(
            "sharp.extract_embeddings.Embedder",
            lambda *a, **kw: stub_embedder,
        )
        monkeypatch.setattr(
            "sharp.extract_embeddings.ensure_model_available",
            lambda *a, **kw: "stub",
        )

        nested = tmp_path / "a" / "b" / "embeddings.parquet"
        assert not nested.parent.exists()
        run(EmbeddingConfig(input_path=fasta_path, output_path=nested))
        assert nested.exists()


# ────────────────────────────── CLI ────────────────────────────────────────

class TestCli:
    def test_required_args(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["--input", "in.faa", "--output", "out.parquet"])
        assert args.input == Path("in.faa")
        assert args.output == Path("out.parquet")

    def test_defaults(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["--input", "in.faa", "--output", "out.parquet"])
        assert args.model == "esm2_t6_8M_UR50D"
        assert args.batch_size == 8
        assert args.max_length == 1024
        assert args.device == "auto"

    def test_input_is_required(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--output", "out.parquet"])

    def test_output_is_required(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--input", "in.faa"])

    def test_rejects_unknown_model(self) -> None:
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                "--input", "in.faa",
                "--output", "out.parquet",
                "--model", "not_a_real_model",
            ])
