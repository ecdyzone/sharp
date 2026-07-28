"""Tests for environment-driven path/device configuration.

`sharp.config` resolves its constants at import time, so every test here reloads
the module inside a patched environment rather than reading the already-imported
values. `_reload_config` also clears any SHARP_* key inherited from the developer's
own `.env` or shell, so the suite behaves the same on a laptop and on a server.
"""
from __future__ import annotations

import importlib
import os
from pathlib import Path

import pytest

import sharp.config


SHARP_KEYS = [
    "SHARP_ENV_FILE",
    "SHARP_DATA_ROOT",
    "SHARP_RAW_DIR",
    "SHARP_INTERIM_DIR",
    "SHARP_PROCESSED_DIR",
    "SHARP_MOCK_DIR",
    "SHARP_DEVICE",
]


@pytest.fixture
def reload_config(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Reload `sharp.config` with a clean env and no real `.env` in the way."""

    def _reload(env: dict[str, str] | None = None, env_file: Path | None = None):
        for key in SHARP_KEYS:
            monkeypatch.delenv(key, raising=False)
        # Point at a nonexistent .env by default so the developer's real one
        # (if any) can't leak into assertions.
        monkeypatch.setenv("SHARP_ENV_FILE",
                           str(env_file or tmp_path / "absent.env"))
        for key, value in (env or {}).items():
            monkeypatch.setenv(key, value)
        return importlib.reload(sharp.config)

    yield _reload
    # Leave the module in its pristine state for tests that import it normally.
    for key in SHARP_KEYS:
        monkeypatch.delenv(key, raising=False)
    importlib.reload(sharp.config)


# ────────────────────────────── defaults ───────────────────────────────────

def test_defaults_are_in_repo_paths(reload_config):
    cfg = reload_config()
    assert cfg.DATA_DIR == cfg.PROJECT_ROOT / "data"
    assert cfg.RAW_DIR == cfg.DATA_DIR / "raw"
    assert cfg.INTERIM_DIR == cfg.DATA_DIR / "interim"
    assert cfg.PROCESSED_DIR == cfg.DATA_DIR / "processed"
    assert cfg.MOCK_DIR == cfg.DATA_DIR / "mock"
    assert cfg.DEFAULT_DEVICE == "auto"


def test_missing_env_file_is_not_an_error(reload_config, tmp_path):
    cfg = reload_config(env_file=tmp_path / "definitely-not-here.env")
    assert cfg.DATA_DIR == cfg.PROJECT_ROOT / "data"


# ────────────────────────── env var overrides ──────────────────────────────

def test_data_root_moves_all_subdirs(reload_config):
    cfg = reload_config({"SHARP_DATA_ROOT": "/scratch/sharp"})
    assert cfg.DATA_DIR == Path("/scratch/sharp")
    assert cfg.RAW_DIR == Path("/scratch/sharp/raw")
    assert cfg.PROCESSED_DIR == Path("/scratch/sharp/processed")


def test_individual_dir_overrides_root(reload_config):
    cfg = reload_config({
        "SHARP_DATA_ROOT": "/scratch/sharp",
        "SHARP_RAW_DIR": "/mnt/shared/raw",
    })
    assert cfg.RAW_DIR == Path("/mnt/shared/raw")
    # The others still follow the root.
    assert cfg.INTERIM_DIR == Path("/scratch/sharp/interim")


def test_tilde_is_expanded(reload_config):
    cfg = reload_config({"SHARP_DATA_ROOT": "~/sharp-data"})
    assert cfg.DATA_DIR == Path.home() / "sharp-data"
    assert "~" not in str(cfg.RAW_DIR)


def test_blank_value_falls_back_to_default(reload_config):
    cfg = reload_config({"SHARP_DATA_ROOT": "   "})
    assert cfg.DATA_DIR == cfg.PROJECT_ROOT / "data"


def test_device_override(reload_config):
    cfg = reload_config({"SHARP_DEVICE": "cpu"})
    assert cfg.DEFAULT_DEVICE == "cpu"
    assert cfg.EmbeddingConfig(input_path=Path("in.faa"),
                               output_path=Path("out.parquet")).device == "cpu"


# ──────────────────────────── .env file parsing ────────────────────────────

def test_env_file_is_read(reload_config, tmp_path):
    env_file = tmp_path / ".env"
    env_file.write_text("SHARP_DATA_ROOT=/data/from-file\nSHARP_DEVICE=cuda\n")
    cfg = reload_config(env_file=env_file)
    assert cfg.DATA_DIR == Path("/data/from-file")
    assert cfg.DEFAULT_DEVICE == "cuda"


def test_comments_blanks_and_quotes(reload_config, tmp_path):
    env_file = tmp_path / ".env"
    env_file.write_text(
        "# a comment\n"
        "\n"
        '   SHARP_DATA_ROOT = "/data/quoted"   \n'
        "SHARP_DEVICE='mps'\n"
        "not_a_pair\n"
    )
    cfg = reload_config(env_file=env_file)
    assert cfg.DATA_DIR == Path("/data/quoted")
    assert cfg.DEFAULT_DEVICE == "mps"


def test_variable_interpolation(reload_config, tmp_path):
    env_file = tmp_path / ".env"
    env_file.write_text(
        "SHARP_DATA_ROOT=/scratch/sharp\n"
        "SHARP_RAW_DIR=${SHARP_DATA_ROOT}/raw-alt\n"
    )
    cfg = reload_config(env_file=env_file)
    assert cfg.RAW_DIR == Path("/scratch/sharp/raw-alt")


def test_real_env_wins_over_env_file(reload_config, tmp_path):
    """A shell override must beat `.env` — that's what makes one-off runs work."""
    env_file = tmp_path / ".env"
    env_file.write_text("SHARP_DEVICE=cuda\n")
    cfg = reload_config({"SHARP_DEVICE": "cpu"}, env_file=env_file)
    assert cfg.DEFAULT_DEVICE == "cpu"


def test_env_file_does_not_clobber_unrelated_vars(reload_config, tmp_path):
    env_file = tmp_path / ".env"
    env_file.write_text("PATH=/nope\n")
    reload_config(env_file=env_file)
    assert os.environ["PATH"] != "/nope"
