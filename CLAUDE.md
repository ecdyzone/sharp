# SHARP benchmarks — assistant instructions

This project benchmarks the team's SHARP tool (called `sharptool` here),
antiSMASH, DeepBGC and GECCO. It does not implement the team's prediction tool.
Read README.md for commands and docs/ for methodology and scope caveats.

## Installation and environment

The source is intended for `iGEM_analysis/benchmarks` in the team repository.
`bash install.sh [DESTINATION]` copies a clean snapshot to `~/sharp-benchmarks`
by default and installs its locked Pixi environment. Existing destinations are
refused. It copies no private .env, data, results, Git metadata or environments.
It does not install baseline tools or download datasets.

Run commands from the independent installation:

```bash
cd ~/sharp-benchmarks
pixi run pytest
pixi run python scripts/generate_mock_benchmark_data.py \
    --n-clusters 20 --recall-rate 0.7 --n-false-positives 5
pixi run python -m sharp.evaluate \
    --predictions data/mock/predictions.parquet \
    --ground-truth data/mock/ground_truth.tsv \
    --output data/processed/benchmark.json
```

Use Pixi rather than bare Python. The internal package name `sharp` intentionally
remains unchanged; do not install it into the team's SHARP environment. The
`report` environment remains available for existing plotting commands.
For development in the source checkout, use `pixi install --locked`.

Slurm scripts use `${SHARP_BENCHMARKS_DIR:-$HOME/sharp-benchmarks}` to locate the
installation. Their tool/database paths and resource directives are specific to
the original cluster. Baseline setup scripts create separate tool environments.

## Layout and contracts

- `src/sharp/config.py`: data paths, optional .env loading, evaluation config.
- `src/sharp/io.py`: data types and file readers/writers.
- `src/sharp/metrics.py`: pure metric calculations; no I/O or logging.
- `src/sharp/evaluate.py`: benchmark orchestration.
- `scripts/`: dataset preparation/downloads, converters, baseline setup/jobs,
  merge/evaluation wrapper, report tables and figures.
- `tests/`: existing unit/integration tests and real/trimmed output fixtures.
- `notebooks/`: benchmark walkthroughs and exports; run from notebooks/.
- `docs/`: architecture, metrics, benchmark scopes and NCBI mirror notes.
- `data/{raw,interim,processed,mock}/`: local data, excluded from delivery.

Keep coordinates 0-based half-open `[start, end)` after ingestion. Each converter
owns its source format's coordinate conversion. Types stay in io.py; config
classes stay in config.py. Shared embedding/FASTA helpers remain for compatibility
but no embedding or model modules are part of this project.

Always pass the same explicit `--contigs` scope file to every compared tool.
MiBIG is primary curated ground truth. BGC Atlas is secondary/noisy; agreement
with its computational predictions is not proof of correctness. The published
three-tool report requires original inputs that are not bundled here. Never
invent results or present mock outputs as measured biological results.

Environment variables/.env carry machine-specific paths and credentials.
Thresholds and other experimental parameters stay on the CLI. Python paths
resolve on import; config tests reload the module with patched environments.

Keep scoring, converters and output schemas unchanged during delivery work.
Tests mirror modules/scripts; patch imported symbols on the importing module.
Use the existing suite and synthetic smoke run; do not launch expensive baseline
analyses or downloads as installation checks.

## Changes and documentation

Keep changes small. For a new file, update the README directory tree and document
usage in its docstring/header, README and this file. Commit meaningful units of
work using Conventional Commits. Do not change the team repository without the
user's permission. Do not commit credentials, environments or bulk data.
