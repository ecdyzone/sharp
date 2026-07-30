# Projeto

> Link GitHub Pages: <https://ecdyzone.github.io/sharp>  
> Repositório GitHub: <https://github.com/ecdyzone/sharp>  
> Diagrama DAG: <https://ecdyzone.github.io/sharp/docs/sharp_dag.html>  

Clique nos links abaixo para ir às páginas HTML interativas:

- [Diagrama DAG (directed acyclic graph)](docs/sharp_dag.html) - Fluxograma do projeto, evidenciando inputs-processos-outputs.
- [Página Descritiva](docs/sharp_pipeline.html) - Praticamente o mesmo conteúdo do Diagrama DAG, mas apresentado com uma interface menos técnica.

## Setting up

First run:

```bash
git clone <repo>
cd <repo>
```

then install with `pixi` or `conda`

### Using pixi (recommended)

Just run:

```bash
pixi install
```

### Using conda

Option 1: use `pixi.lock`

```bash
conda create --name my-env --file pixi.lock
```

Option 2: use `environment.yml`:

```bash
conda env create -f environment.yml
```

### Machine-specific settings (optional `.env`)

Everything works out of the box with no `.env` — data is read from and written to
`./data` and the embedding step picks its device automatically. Create one only
when this machine differs, e.g. a server that keeps the data on a scratch disk:

```bash
cp .env.example .env
$EDITOR .env
```

```bash
# .env — gitignored, never committed
SHARP_DATA_ROOT=/scratch/$USER/sharp-data   # default: <repo>/data
SHARP_DEVICE=cpu                            # default: auto (CUDA -> MPS -> CPU)
```

Setting `SHARP_DATA_ROOT` moves `raw/`, `interim/`, `processed/` and `mock/` with
it; override one individually with `SHARP_RAW_DIR`, `SHARP_INTERIM_DIR`,
`SHARP_PROCESSED_DIR` or `SHARP_MOCK_DIR`. See `.env.example` for all keys.

A real environment variable beats the file, so a single run can be redirected
without editing anything:

```bash
SHARP_DEVICE=cpu pixi run python -m sharp.extract_embeddings ...
```

`.env` holds **machine identity only** — where the data lives, what hardware is
present, and credentials once a step needs them. Pipeline parameters (thresholds,
model choice, batch sizes) stay on the CLI so a run remains reproducible from the
command that produced it.

## Workflow

> For reproducing with conda/mamba you have to:
>
> - `conda activate <environment-name>`
> - run the commands below without `pixi run`

### Extract Embeddings

```bash
## 1. Generate test data
pixi run python scripts/generate_mock_data.py --n 100

## 2. Run the step against it
pixi run python -m sharp.extract_embeddings \
    --input data/mock/neighborhood_proteins.faa \
    --output data/interim/embeddings.parquet
```

### Benchmarks

```bash
# 1. Generate correlated mock data — clusters and predictions that overlap by construction
pixi run python scripts/generate_mock_benchmark_data.py \
    --n-clusters 20 --recall-rate 0.7 --n-false-positives 5

# 2. Evaluate
pixi run python -m sharp.evaluate \
    --predictions data/mock/predictions.parquet \
    --ground-truth data/mock/ground_truth.tsv \
    --output data/processed/benchmark.json

# Output: detection recall=0.700 (14/20), matched 14/19 predictions
#         — matches the generator's --recall-rate 0.7 and 5 injected extras
```

#### Reading `benchmark.json`

| block | what it answers |
|---|---|
| `scope` | how much of the ground truth was evaluable; `explicit` vs `inferred` |
| `detection` | *did the tool find the BGC?* (fraction of the cluster covered) |
| `reciprocal` | the strict symmetric rule, reported for comparison |
| `nucleotide` | bp-level agreement; `precision` = how much extra territory was called |
| `boundary` | tightness of the calls, plus split/merge diagnostics |

Two things the schema is deliberate about:

- **Recall counts only contigs the tool was run on.** Pass `--contigs` (one name
  per line, or a `.fai`), and give **every tool in a comparison the same file**.
  Omitted, the scope is inferred from the predictions and a warning is logged —
  that is optimistic, because a contig analyzed but not called on drops out of
  the denominator.
- **Unmatched is not false.** Ground truth is incomplete by construction, so a
  prediction with no match is unvalidated rather than wrong. The output reports
  `matched_prediction_frac` (a lower bound on precision) and
  `unmatched_prediction_ids` — there is no region-level `precision` field.

See `docs/ARCHITECTURE.md` → "Metrics — methodological choices" for the full
rationale.

### Competitor baselines (antiSMASH / DeepBGC / GECCO)

S(H)ARP does **not** run these tools — each has incompatible dependencies and
installs into its own isolated pixi env via `scripts/setup_<tool>.sh`. You run
the tool yourself, then convert its output to `predictions.parquet` and evaluate
it exactly like S(H)ARP's own predictions.

All three locations the setup scripts write to are set in `.env` (see
`.env.example`), so a laptop and a server can differ without editing a tracked
file. Defaults apply when `.env` is absent, and each script echoes the paths it
resolved:

| Key | Default | What |
|---|---|---|
| `TOOLS_INSTALL_DIR` | `~/.local/src` | one isolated pixi env per tool |
| `ANTISMASH_DOWNLOADS_DIR` | `${DATABASES:-~/.local/share}/antismash/databases` | ~10GB reference data |
| `DEEPBGC_DOWNLOADS_DIR` | `${DATABASES:-~/.local/share}/deepbgc/data` | ~3GB models + Pfam |

The database dirs follow `$DATABASES` when the machine exports it (as the server
does) and fall back to `~/.local/share` otherwise, so neither machine normally
needs an edit. The `cd ~/.local/src/<tool>` commands below assume the default
`TOOLS_INSTALL_DIR` — if you override it, use the path the setup script prints
when it finishes.

```bash
# 1. Install a baseline into its own env (~/.local/src/<tool>/), one-time
bash scripts/setup_antismash.sh

# 2. Run it yourself, from its own env (or on HPC / in a container)
cd ~/.local/src/antismash && pixi run antismash <genome.gbk> --output-dir <out>
# for non-annotated fasta, the code changes a bit:
# cd ~/.local/src/antismash && pixi run antismash <genome.fasta> --output-dir <out> --genefinding-tool prodigal

# 3. Convert its output to predictions.parquet (runs in the S(H)ARP env)
#    (antiSMASH, DeepBGC, and GECCO converters all written)
#    Inspect first to verify the schema against your actual output:
pixi run python scripts/convert_antismash_to_parquet.py --inspect <out>

pixi run python scripts/convert_antismash_to_parquet.py \
    --input <out> --output data/interim/antismash_predictions.parquet

# 4. List the contigs the tool was run on — the denominator of recall.
#    Build it once from the input genome and reuse it for EVERY tool, or the
#    recall denominators differ and the numbers stop being comparable.
grep '^>' <genome.fasta> | cut -c2- | cut -d' ' -f1 \
    > data/interim/analyzed_contigs.txt

# 5. Evaluate against the same ground truth and scope as S(H)ARP
pixi run python -m sharp.evaluate \
    --predictions data/interim/antismash_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --contigs data/interim/analyzed_contigs.txt \
    --output data/processed/benchmark_antismash.json
```

DeepBGC follows the same shape — the tool runs in its own env, S(H)ARP only parses `<prefix>.bgc.tsv`:

```bash
bash scripts/setup_deepbgc.sh
cd ~/.local/src/deepbgc && pixi run deepbgc pipeline <genome.fasta> --output out

pixi run python scripts/convert_deepbgc_to_parquet.py --inspect out
pixi run python scripts/convert_deepbgc_to_parquet.py \
    --input out --output data/interim/deepbgc_predictions.parquet

pixi run python -m sharp.evaluate \
    --predictions data/interim/deepbgc_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --contigs data/interim/analyzed_contigs.txt \
    --output data/processed/benchmark_deepbgc.json
```

GECCO too — its `start`/`end` are 1-based inclusive (the one baseline tool that
needs a coordinate conversion), which the converter applies automatically:

```bash
bash scripts/setup_gecco.sh
cd ~/.local/src/gecco && pixi run gecco run --genome <genome.fasta> --output-dir out

pixi run python scripts/convert_gecco_to_parquet.py --inspect out
pixi run python scripts/convert_gecco_to_parquet.py \
    --input out --output data/interim/gecco_predictions.parquet

pixi run python -m sharp.evaluate \
    --predictions data/interim/gecco_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --contigs data/interim/analyzed_contigs.txt \
    --output data/processed/benchmark_gecco.json
```

### Converting Parquet to TSV

Any parquet file the pipeline produces (`predictions.parquet`,
`embeddings.parquet`, `kg_features.parquet`, ...) can be dumped to a plain
TSV for inspection or sharing. List-typed columns (e.g. `embeddings.parquet`'s
`embedding` vector) have no TSV representation, so they're joined into a
single comma-separated cell — this is a one-way, informational dump, not a
round-trippable format.

```bash
# Inspect the schema first, especially for anything with list-typed columns
pixi run python scripts/parquet_to_tsv.py --inspect data/interim/embeddings.parquet

pixi run python scripts/parquet_to_tsv.py \
    --input data/interim/predictions.parquet \
    --output data/interim/predictions.tsv
```

### Downloading a Benchmark Genome

Fetches one contig by accession from NCBI nuccore and derives the `--contigs`
scope file alongside it. Defaults to `AL645882.2` (*S. coelicolor* A3(2)), which
carries 15 coordinate-resolved MiBiG clusters — the most of any single contig,
and enough for a recall number that actually varies.

```bash
# Default: S. coelicolor A3(2)
scripts/download_genome.sh

# Or any other nuccore accession
scripts/download_genome.sh CP002993.1
```

Writes `data/raw/<ACCESSION>.fasta` and `data/interim/analyzed_contigs.txt`.
Pass that same contigs file to **every** tool in a comparison — see the
benchmark-scope caveat in `CLAUDE.md`.

### Preparing MiBiG Database

```bash
# All clusters
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/mibig_ground_truth.tsv

# Or focused on your organism of interest
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/streptomyces_ground_truth.tsv \
    --genus Streptomyces
```

### Preparing BGC Atlas Database

Secondary, noisy ground truth (labels are themselves antiSMASH predictions —
report alongside MiBiG, never alone). The dump is 204k antiSMASH `.gbk` files
downloaded by `scripts/download_bgc-atlas.sh` (DVC-managed).

```bash
# Inspect a few real files first (verify the schema on disk)
pixi run python scripts/prepare_bgcatlas_ground_truth.py \
    --inspect data/raw/complete-bgcs

# Build the TSV (streams over all ~204k files)
pixi run python scripts/prepare_bgcatlas_ground_truth.py \
    --input-dir data/raw/complete-bgcs \
    --output data/raw/bgcatlas_ground_truth.tsv

# Develop / test against a small subset without walking 10 GB
pixi run python scripts/prepare_bgcatlas_ground_truth.py \
    --input-dir data/raw/complete-bgcs \
    --output data/interim/bgcatlas_sample.tsv --limit 100
```

## Tests

Run tests with:

```bash
pixi run pytest
```

## Directory Structure

```bash
.
├── .env.example                          # template for machine-specific settings (copy to .env)
├── benchmarks
├── config
├── data
├── docs
│   ├── sharp_dag.html
│   └── sharp_pipeline.html
├── environment.yml
├── LICENSE
├── notebooks
│   ├── benchmarks_part1.py
│   ├── benchmarks_part2.py
│   ├── conversions
│   │   ├── benchmarks_part1.html
│   │   ├── benchmarks_part1.ipynb
│   │   ├── benchmarks_part2.html
│   │   └── benchmarks_part2.ipynb
│   └── inspecting_parquet.py
├── pixi.lock
├── pixi.toml
├── pyproject.toml
├── README.md
├── scripts
│   ├── convert_antismash_to_parquet.py   # antiSMASH JSON -> predictions.parquet (no coord conversion)
│   ├── convert_deepbgc_to_parquet.py     # DeepBGC .bgc.tsv -> predictions.parquet (no coord conversion)
│   ├── convert_gecco_to_parquet.py       # GECCO .clusters.tsv -> predictions.parquet (start-1: 1-based -> 0-based)
│   ├── download_bgc-atlas.sh
│   ├── download_genome.sh                # NCBI accession -> data/raw/<ACC>.fasta + --contigs scope file
│   ├── download_mibig.sh
│   ├── generate_mock_benchmark_data.py
│   ├── generate_mock_data.py
│   ├── parquet_to_tsv.py                 # generic parquet -> TSV dump (any pipeline parquet file)
│   ├── prepare_bgcatlas_ground_truth.py
│   ├── prepare_mibig_ground_truth.py
│   ├── _load_env.sh               # sourced by the setup scripts: loads .env, exports it
│   ├── setup_antismash.sh         # install baseline into its own isolated pixi env
│   ├── setup_deepbgc.sh
│   └── setup_gecco.sh
├── src
│   └── sharp
│       ├── __init__.py
│       ├── config.py
│       ├── evaluate.py
│       ├── extract_embeddings.py
│       ├── io.py
│       ├── metrics.py
│       └── model_management.py
└── tests
    ├── conftest.py
    ├── fixtures
    │   ├── AL589148_ground_truth.tsv        # the one MiBIG cluster on AL589148.1
    │   ├── antismash_predictions.parquet    # converted real output, benchmark regression
    │   ├── antismash_sequence.json          # trimmed real antiSMASH 8.0.4 summary
    │   ├── deepbgc_out.bgc.tsv              # real (unmodified) DeepBGC 0.1.0 output
    │   ├── deepbgc_predictions.parquet      # converted real output, benchmark regression
    │   ├── gecco_predictions.parquet        # converted real output, benchmark regression
    │   └── gecco_sequence.clusters.tsv      # real (unmodified) GECCO 0.10.3 output
    ├── test_config.py
    ├── test_convert_antismash.py
    ├── test_convert_deepbgc.py
    ├── test_convert_gecco.py
    ├── test_evaluate.py
    ├── test_extract_embeddings.py
    ├── test_generate_mock_data.py
    ├── test_io.py
    ├── test_metrics.py
    ├── test_model_management.py
    ├── test_parquet_to_tsv.py
    ├── test_prepare_bgcatlas.py
    └── test_prepare_mibig.py
```

## Currently Working on

- prototyping benchmarks
