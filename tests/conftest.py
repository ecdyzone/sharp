"""Shared fixtures for the test suite."""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

from sharp.io import ProteinRecord


@pytest.fixture
def sample_records() -> list[ProteinRecord]:
    """Four hand-written protein records spanning two regions."""
    return [
        ProteinRecord("p1", "R001", "MKTAYIAKQRQI"),
        ProteinRecord("p2", "R001", "SFVKSHFSRQLE"),
        ProteinRecord("p3", "R002", "ACDEFGHIKLMNPQR"),
        ProteinRecord("p4", "R002", "STVWY"),
    ]


@pytest.fixture
def fasta_path(tmp_path: Path, sample_records: list[ProteinRecord]) -> Path:
    """Write `sample_records` to a FASTA in tmp_path and return its path."""
    path = tmp_path / "proteins.faa"
    lines = []
    for r in sample_records:
        lines.append(f">{r.protein_id} region_id={r.region_id}")
        lines.append(r.sequence)
    path.write_text("\n".join(lines) + "\n")
    return path


class StubEmbedder:
    """Drop-in replacement for the real Embedder. Deterministic output
    per (sequence, dim) so tests can assert reproducibility."""

    def __init__(self, dim: int = 32):
        self.dim = dim

    def embed_batch(self, sequences: Sequence[str]) -> np.ndarray:
        out = np.empty((len(sequences), self.dim), dtype=np.float32)
        for i, seq in enumerate(sequences):
            rng = np.random.default_rng(abs(hash(seq)) % (2**31))
            out[i] = rng.standard_normal(self.dim).astype(np.float32)
        return out


@pytest.fixture
def stub_embedder() -> StubEmbedder:
    return StubEmbedder(dim=32)
