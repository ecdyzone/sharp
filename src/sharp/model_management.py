"""Model registry, cache management, device selection, and the Embedder class.

Owns the full lifecycle of the embedding model:
  - which models exist (`MODEL_REGISTRY`)
  - is this one cached locally? (`ensure_model_available`)
  - where to run it (`select_device`)
  - tokenize → forward → pool (`Embedder`)

The masking logic inside `Embedder.embed_batch` is tightly coupled to the
tokenizer's CLS/EOS conventions, so the two live together.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

LOG = logging.getLogger(__name__)


# short name → (HuggingFace hub id, hidden dimension)
MODEL_REGISTRY: dict[str, tuple[str, int]] = {
    "esm2_t6_8M_UR50D":    ("facebook/esm2_t6_8M_UR50D",    320),
    "esm2_t12_35M_UR50D":  ("facebook/esm2_t12_35M_UR50D",  480),
    "esm2_t30_150M_UR50D": ("facebook/esm2_t30_150M_UR50D", 640),
    "esm2_t33_650M_UR50D": ("facebook/esm2_t33_650M_UR50D", 1280),
}


# ────────────────────────────── device ─────────────────────────────────────

def select_device(spec: str = "auto") -> torch.device:
    """Resolve a device spec. `auto` prefers CUDA → MPS → CPU."""
    if spec != "auto":
        return torch.device(spec)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# ────────────────────────────── cache ──────────────────────────────────────

def ensure_model_available(
    model_name: str, cache_dir: Path | None = None
) -> str:
    """Ensure model files are present locally. Downloads if not cached.

    Called before any compute-intensive work so a download failure happens
    at a known checkpoint, not mid-batch. Defaults to HuggingFace's standard
    cache (`~/.cache/huggingface/hub/`); pass `cache_dir` to override.

    Returns the resolved hub_id so the caller doesn't need MODEL_REGISTRY.
    """
    if model_name not in MODEL_REGISTRY:
        raise ValueError(
            f"unknown model {model_name!r}; available: {list(MODEL_REGISTRY)}"
        )
    hub_id, _ = MODEL_REGISTRY[model_name]

    # Lazy import: huggingface_hub is pulled in by transformers, but keeping
    # the top of this module focused on the public surface.
    from huggingface_hub import snapshot_download, try_to_load_from_cache

    cached = try_to_load_from_cache(hub_id, "config.json", cache_dir=cache_dir)
    if cached is None:
        LOG.info("model %s not cached — downloading", hub_id)
        snapshot_download(repo_id=hub_id, cache_dir=cache_dir)
        LOG.info("download complete")
    else:
        LOG.info("model %s already cached", hub_id)
    return hub_id


# ────────────────────────────── embedder ───────────────────────────────────

class Embedder:
    """ESM-2 wrapper: tokenize → forward → mean-pool over residues.

    Masks CLS, EOS, and pad tokens so each protein vector reflects only
    residue positions. This matters: for a 50-aa protein, naïvely pooling
    with the attention mask includes ~4% noise from CLS+EOS.
    """

    def __init__(self, model_name: str, device: torch.device, max_length: int):
        hub_id, self.dim = MODEL_REGISTRY[model_name]
        LOG.info("loading %s (%d-dim) on %s", hub_id, self.dim, device)
        self.tokenizer = AutoTokenizer.from_pretrained(hub_id)
        self.model = AutoModel.from_pretrained(hub_id).to(device).eval()
        self.device = device
        self.max_length = max_length

    @torch.inference_mode()
    def embed_batch(self, sequences: list[str]) -> np.ndarray:
        tok = self.tokenizer(
            sequences,
            padding=True,
            truncation=True,
            max_length=self.max_length,
            return_tensors="pt",
        ).to(self.device)

        # (B, L, D) — last hidden state for every token position
        hidden = self.model(**tok).last_hidden_state

        # Residue-only mask: zero CLS (position 0) and EOS (last non-pad position).
        mask = tok.attention_mask.clone()
        mask[:, 0] = 0
        last_idx = tok.attention_mask.sum(dim=1) - 1
        mask[torch.arange(mask.size(0), device=self.device), last_idx] = 0

        m = mask.unsqueeze(-1).to(hidden.dtype)
        pooled = (hidden * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
        return pooled.float().cpu().numpy()
