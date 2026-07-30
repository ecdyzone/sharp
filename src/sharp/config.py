"""Project paths and per-step configuration dataclasses.

Paths follow the Cookiecutter Data Science convention:
  raw/       — immutable external inputs (genomes from NCBI, MiBIG dumps)
  interim/   — intermediate pipeline artifacts (between steps)
  processed/ — final consumer-facing outputs (reports, trained models)
  mock/      — synthetic data for testing

Machine-specific values (where the data lives, which device to embed on) are read
from the environment, optionally seeded by a `.env` file at the project root. This
exists so the same checkout runs unmodified on a laptop and on a server where the
data sits on a different filesystem. See `.env.example` for the supported keys.

Scope rule: the environment carries *machine identity* only — data roots and
device. Pipeline behaviour (thresholds, model choice, batch sizes) stays in the
dataclass defaults below or on the CLI, so a run is reproducible from its command
line rather than from an unversioned file.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ENV_FILE = Path(os.environ.get("SHARP_ENV_FILE", PROJECT_ROOT / ".env"))


def _load_dotenv(path: Path) -> None:
    """Seed `os.environ` from a `.env` file, without overriding what's already set.

    Real environment variables always win, so `SHARP_DATA_ROOT=... pixi run ...`
    overrides the file for a single run. Deliberately minimal (no export
    keyword, no quoting rules beyond stripping matched quotes) — a dependency on
    python-dotenv is not worth it for the handful of keys we support.
    """
    if not path.is_file():
        return
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, sep, value = line.partition("=")
        if not sep:
            continue
        key = key.strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        # `${VAR}` lets .env derive one path from another (e.g. dirs from the root).
        # Expanded against the environment as built so far, so order matters.
        os.environ.setdefault(key, os.path.expandvars(value))


_load_dotenv(ENV_FILE)


def _env_path(key: str, default: Path) -> Path:
    """Read a path from the environment, falling back to the in-repo default.

    `~` is expanded so `.env` can say `~/scratch/sharp-data` on a server.
    """
    value = os.environ.get(key)
    if not value or not value.strip():
        return default
    return Path(value.strip()).expanduser()


DATA_DIR = _env_path("SHARP_DATA_ROOT", PROJECT_ROOT / "data")
RAW_DIR = _env_path("SHARP_RAW_DIR", DATA_DIR / "raw")
INTERIM_DIR = _env_path("SHARP_INTERIM_DIR", DATA_DIR / "interim")
PROCESSED_DIR = _env_path("SHARP_PROCESSED_DIR", DATA_DIR / "processed")
MOCK_DIR = _env_path("SHARP_MOCK_DIR", DATA_DIR / "mock")

# Device for the embedding step. "auto" preserves the previous behaviour
# (CUDA → MPS → CPU); set SHARP_DEVICE=cpu on a laptop without a usable GPU.
DEFAULT_DEVICE = os.environ.get("SHARP_DEVICE", "auto").strip() or "auto"


@dataclass(frozen=True)
class EmbeddingConfig:
    """Configuration for the embedding extraction step."""
    input_path: Path
    output_path: Path
    model_name: str = "esm2_t6_8M_UR50D"
    batch_size: int = 8
    max_length: int = 1024
    device: str = DEFAULT_DEVICE
    log_every: int = 50


@dataclass(frozen=True)
class EvaluateConfig:
    """Configuration for the benchmark / evaluation step.

    Every knob here changes what a benchmark number *means*, so all of them stay
    on the CLI and out of `.env` (see the scope rule above): the same command
    line must reproduce the same numbers on any machine.

    contigs_path        — contigs the tool was actually run on. None infers the
                          scope from the predictions and warns; see evaluate.py.
    min_cluster_frac    — fraction of a cluster that must be covered to count it
                          as found.
    min_prediction_frac — fraction of a prediction that must be covered for it to
                          count as tightly bounded. 0.0 leaves detection and
                          boundary accuracy as separate questions.
    reciprocal_frac     — threshold for the strict symmetric rule reported
                          alongside the asymmetric one.
    min_p_bgc           — drop predictions scoring below this before scoring.
    max_listed_ids      — cap on each id list in benchmark.json; 0 = unlimited.
    """
    predictions_path: Path
    ground_truth_path: Path
    output_path: Path
    contigs_path: Path | None = None
    min_cluster_frac: float = 0.5
    min_prediction_frac: float = 0.0
    reciprocal_frac: float = 0.5
    min_p_bgc: float = 0.0
    max_listed_ids: int = 1000
