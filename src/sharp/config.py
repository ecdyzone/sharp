"""Project paths and per-step configuration dataclasses.

Paths follow the Cookiecutter Data Science convention:
  raw/       — immutable external inputs (genomes from NCBI, MiBIG dumps)
  interim/   — intermediate pipeline artifacts (between steps)
  processed/ — final consumer-facing outputs (reports, trained models)
  mock/      — synthetic data for testing
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
MOCK_DIR = DATA_DIR / "mock"


@dataclass(frozen=True)
class EmbeddingConfig:
    """Configuration for the embedding extraction step."""
    input_path: Path
    output_path: Path
    model_name: str = "esm2_t6_8M_UR50D"
    batch_size: int = 8
    max_length: int = 1024
    device: str = "auto"
    log_every: int = 50


@dataclass(frozen=True)
class EvaluateConfig:
    """Configuration for the benchmark / evaluation step."""
    predictions_path: Path
    ground_truth_path: Path
    output_path: Path
    min_overlap_frac: float = 0.5
