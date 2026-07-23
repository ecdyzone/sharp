"""Tests for scripts.parquet_to_tsv.

Generic parquet -> TSV dump used to hand collaborators a plain-text copy of
any pipeline parquet file. Two things are under test: (1) scalar-only
tables (predictions.parquet's shape) round-trip cleanly, and (2) list-typed
columns (embeddings.parquet's `embedding` column) are joined into a single
TSV cell rather than raising or losing data.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from parquet_to_tsv import cell, convert  # noqa: E402

from sharp.io import PredictedRegion, write_predictions_parquet  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestCell:
    def test_scalar_passthrough(self) -> None:
        assert cell("AL589148.1") == "AL589148.1"
        assert cell(0.9) == 0.9
        assert cell(None) is None

    def test_list_joined_with_comma(self) -> None:
        assert cell([0.1, 0.2, 0.3]) == "0.1,0.2,0.3"

    def test_empty_list_joined_to_empty_string(self) -> None:
        assert cell([]) == ""


class TestConvertScalarTable:
    def test_predictions_parquet_round_trips(self, tmp_path: Path) -> None:
        predictions = [
            PredictedRegion(
                region_id="AL589148.1_cluster_1",
                contig="AL589148.1",
                start=100,
                end=200,
                p_bgc=0.9,
                predicted_class="NRPS",
            ),
            PredictedRegion(
                region_id="AL589148.1_cluster_2",
                contig="AL589148.1",
                start=300,
                end=400,
                p_bgc=0.5,
                predicted_class=None,
            ),
        ]
        parquet_path = tmp_path / "predictions.parquet"
        write_predictions_parquet(parquet_path, predictions)

        tsv_path = tmp_path / "predictions.tsv"
        n = convert(parquet_path, tsv_path)
        assert n == 2

        rows = read_tsv(tsv_path)
        assert len(rows) == 2
        assert rows[0]["region_id"] == "AL589148.1_cluster_1"
        assert rows[0]["start"] == "100"
        assert rows[0]["predicted_class"] == "NRPS"
        # A None predicted_class comes back as an empty TSV cell.
        assert rows[1]["predicted_class"] == ""

    def test_creates_output_parent_dirs(self, tmp_path: Path) -> None:
        parquet_path = tmp_path / "predictions.parquet"
        write_predictions_parquet(parquet_path, [])

        tsv_path = tmp_path / "nested" / "dir" / "predictions.tsv"
        n = convert(parquet_path, tsv_path)
        assert n == 0
        assert tsv_path.exists()


class TestConvertListTypedColumn:
    def test_list_column_joined_not_dropped(self, tmp_path: Path) -> None:
        import numpy as np

        from sharp.io import ProteinRecord, write_embeddings_parquet

        records = [ProteinRecord(protein_id="p1", region_id="R001", sequence="MK")]
        vecs = np.array([[0.1, 0.2, 0.3]], dtype=np.float32)

        parquet_path = tmp_path / "embeddings.parquet"
        write_embeddings_parquet(parquet_path, [(records, vecs)], embedding_dim=3)

        tsv_path = tmp_path / "embeddings.tsv"
        n = convert(parquet_path, tsv_path)
        assert n == 1

        rows = read_tsv(tsv_path)
        assert rows[0]["protein_id"] == "p1"
        assert rows[0]["region_id"] == "R001"
        # Embedding vector survives as a comma-joined string, not a Python
        # list repr and not silently truncated.
        values = [float(x) for x in rows[0]["embedding"].split(",")]
        assert len(values) == 3
        assert values[0] == 0.10000000149011612  # float32 -> python float
