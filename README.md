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

# Output: precision=0.737, recall=0.700, F1=0.718 — matches the generator's prediction
```

### Competitor baselines (antiSMASH / DeepBGC / GECCO)

S(H)ARP does **not** run these tools — each has incompatible dependencies and
installs into its own isolated pixi env via `scripts/setup_<tool>.sh`. You run
the tool yourself, then convert its output to `predictions.parquet` and evaluate
it exactly like S(H)ARP's own predictions.

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

# 4. Evaluate against the same ground truth as S(H)ARP
pixi run python -m sharp.evaluate \
    --predictions data/interim/antismash_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
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
    --output data/processed/benchmark_gecco.json
```

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
│   ├── download_mibig.sh
│   ├── generate_mock_benchmark_data.py
│   ├── generate_mock_data.py
│   ├── prepare_bgcatlas_ground_truth.py
│   ├── prepare_mibig_ground_truth.py
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
    │   ├── antismash_sequence.json          # trimmed real antiSMASH 8.0.4 summary
    │   ├── deepbgc_out.bgc.tsv              # real (unmodified) DeepBGC 0.1.0 output
    │   └── gecco_sequence.clusters.tsv      # real (unmodified) GECCO 0.10.3 output
    ├── test_convert_antismash.py
    ├── test_convert_deepbgc.py
    ├── test_convert_gecco.py
    ├── test_evaluate.py
    ├── test_extract_embeddings.py
    ├── test_generate_mock_data.py
    ├── test_io.py
    ├── test_metrics.py
    ├── test_model_management.py
    ├── test_prepare_bgcatlas.py
    └── test_prepare_mibig.py
```

## Currently Working on

- prototyping benchmarks
