"""Tests for sharp.model_management — registry, device, pooling.

The actual model loading (`Embedder.__init__`, `ensure_model_available`'s
download path) is exercised the first time you run the real pipeline; it's
intentionally not unit-tested because that just tests transformers and
huggingface_hub. We DO test:
  - registry shape and lookup
  - device selection logic
  - the residue-only mean-pool, with synthetic tensors so no real model is needed
"""
from __future__ import annotations

import pytest
import torch

from sharp.model_management import (
    MODEL_REGISTRY,
    ensure_model_available,
    residue_mean_pool,
    select_device,
)


# ────────────────────────────── registry ───────────────────────────────────

class TestRegistry:
    def test_default_model_present(self) -> None:
        assert "esm2_t6_8M_UR50D" in MODEL_REGISTRY
        hub_id, dim = MODEL_REGISTRY["esm2_t6_8M_UR50D"]
        assert hub_id == "facebook/esm2_t6_8M_UR50D"
        assert dim == 320

    def test_dimensions_are_monotonic(self) -> None:
        # Bigger models in the registry should have larger hidden dims.
        dims = [d for _, d in MODEL_REGISTRY.values()]
        assert dims == sorted(dims), "registry should be ordered by capacity"

    def test_ensure_model_available_rejects_unknown(self) -> None:
        with pytest.raises(ValueError, match="unknown model"):
            ensure_model_available("not_a_real_model")


# ────────────────────────────── device ─────────────────────────────────────

class TestSelectDevice:
    def test_explicit_cpu(self) -> None:
        assert select_device("cpu").type == "cpu"

    def test_auto_returns_valid_device(self) -> None:
        # Just verify it returns something in the known set; the actual
        # value depends on the machine running the test.
        device = select_device("auto")
        assert device.type in {"cpu", "cuda", "mps"}

    def test_auto_default(self) -> None:
        # Calling with no args should be equivalent to "auto".
        assert select_device().type == select_device("auto").type


# ────────────────────────────── residue_mean_pool ──────────────────────────

class TestResidueMeanPool:
    """The mean-pool must exclude CLS, EOS, and pad. We construct synthetic
    hidden states where position i holds the value i, so the pooled output
    equals the arithmetic mean of the residue positions — which we can
    verify by hand."""

    def test_excludes_cls_eos_and_pad(self) -> None:
        # Batch of two proteins.
        # Row 0: protein of length 3 → [CLS] r1 r2 r3 [EOS] PAD PAD
        # Row 1: protein of length 5 → [CLS] r1 r2 r3 r4 r5 [EOS]
        attention_mask = torch.tensor([
            [1, 1, 1, 1, 1, 0, 0],
            [1, 1, 1, 1, 1, 1, 1],
        ])
        B, L = attention_mask.shape
        D = 4

        # hidden[b, i, :] = i  →  pooled value should equal mean of residue indices
        hidden = (
            torch.arange(L, dtype=torch.float32)
            .repeat(B, 1)
            .unsqueeze(-1)
            .expand(B, L, D)
            .clone()
        )

        pooled = residue_mean_pool(hidden, attention_mask)

        # Row 0: residues at positions 1, 2, 3 → mean = 2.0
        # Row 1: residues at positions 1, 2, 3, 4, 5 → mean = 3.0
        assert torch.allclose(pooled[0], torch.full((D,), 2.0))
        assert torch.allclose(pooled[1], torch.full((D,), 3.0))

    def test_output_shape(self) -> None:
        B, L, D = 5, 12, 7
        hidden = torch.randn(B, L, D)
        mask = torch.ones(B, L, dtype=torch.long)
        out = residue_mean_pool(hidden, mask)
        assert out.shape == (B, D)

    def test_zero_residue_protein_does_not_div_by_zero(self) -> None:
        # Pathological: a "protein" with only CLS+EOS (no residues).
        # After masking, the denominator would be 0, but clamp(min=1)
        # prevents NaN — the output should be the zero vector.
        attention_mask = torch.tensor([[1, 1, 0, 0, 0]])
        hidden = torch.randn(1, 5, 4)
        out = residue_mean_pool(hidden, attention_mask)
        assert not torch.isnan(out).any()
        assert torch.allclose(out, torch.zeros(1, 4))

    def test_does_not_mutate_attention_mask(self) -> None:
        # Calling code reuses tok.attention_mask after pooling.
        attention_mask = torch.tensor([[1, 1, 1, 1, 0]])
        before = attention_mask.clone()
        hidden = torch.randn(1, 5, 3)
        residue_mean_pool(hidden, attention_mask)
        assert torch.equal(attention_mask, before)
