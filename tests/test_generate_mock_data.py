"""Tests for scripts.generate_mock_data — mock protein generation.

We import the script as a module; since `scripts/` is not part of the
`sharp` package, conftest in this folder adds it to sys.path.
"""
from __future__ import annotations

import sys
from pathlib import Path

# scripts/ lives at the project root, not under src/.
_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from generate_mock_data import AMINO_ACIDS, generate_mock_proteins  # noqa: E402


class TestGenerateMockProteins:
    def test_count(self) -> None:
        assert len(generate_mock_proteins(50)) == 50

    def test_deterministic_with_same_seed(self) -> None:
        a = generate_mock_proteins(20, seed=42)
        b = generate_mock_proteins(20, seed=42)
        assert [r.sequence for r in a] == [r.sequence for r in b]
        assert [r.region_id for r in a] == [r.region_id for r in b]

    def test_different_seeds_diverge(self) -> None:
        a = generate_mock_proteins(20, seed=1)
        b = generate_mock_proteins(20, seed=2)
        # Vanishingly unlikely to match by accident at 20 records.
        assert [r.sequence for r in a] != [r.sequence for r in b]

    def test_sequences_use_only_valid_amino_acids(self) -> None:
        valid = set(AMINO_ACIDS)
        for r in generate_mock_proteins(100):
            assert set(r.sequence) <= valid, f"invalid AA in {r.protein_id}"

    def test_lengths_in_expected_range(self) -> None:
        # Generator clips to [50, 800].
        for r in generate_mock_proteins(200):
            assert 50 <= r.length <= 800

    def test_region_distribution(self) -> None:
        # Generator targets ~n/10 distinct regions.
        records = generate_mock_proteins(100, seed=0)
        unique_regions = {r.region_id for r in records}
        assert 1 <= len(unique_regions) <= 10

    def test_unique_protein_ids(self) -> None:
        records = generate_mock_proteins(50)
        pids = [r.protein_id for r in records]
        assert len(pids) == len(set(pids))
