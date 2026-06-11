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
│   └── inspecting_parquet.py
├── pixi.lock
├── pixi.toml
├── pyproject.toml
├── README.md
├── scripts
│   ├── download_mibig.sh
│   ├── generate_mock_benchmark_data.py
│   ├── generate_mock_data.py
│   └── prepare_mibig_ground_truth.py
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
    ├── test_evaluate.py
    ├── test_extract_embeddings.py
    ├── test_generate_mock_data.py
    ├── test_io.py
    ├── test_metrics.py
    ├── test_model_management.py
    └── test_prepare_mibig.py
```

## Currently Working on

- prototyping benchmarks
