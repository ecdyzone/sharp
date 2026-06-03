# Projeto

> Link GitHub Pages: <https://ecdyzone.github.io/sharp>  
> Repositório GitHub: <https://github.com/ecdyzone/sharp>  
> Diagrama DAG: <https://ecdyzone.github.io/sharp/docs/sharp_dag.html>  

Clique nos links abaixo para ir às páginas HTML interativas:

- [Diagrama DAG (directed acyclic graph)](docs/sharp_dag.html) - Fluxograma do projeto, evidenciando inputs-processos-outputs.
- [Página Descritiva](docs/sharp_pipeline.html) - Praticamente o mesmo conteúdo do Diagrama DAG, mas apresentado com uma interface menos técnica.

# Setting up

First run:

```bash
git clone <repo>
cd <repo>
```

then install with `pixi` or `conda`

## using pixi

Just run:

```bash
pixi install
```

## using conda

Option 1: use `pixi.lock`

```bash
conda create --name my-env --file pixi.lock
```

Option 2: use `environment.yml`:

```bash
conda env create -f environment.yml
```

# Workflow

```bash
# 1. Generate test data
pixi run python scripts/generate_mock_data.py --n 100

# 2. Run the step against it
pixi run python -m sharp.extract_embeddings \
    --input data/mock/neighborhood_proteins.faa \
    --output data/interim/embeddings.parquet
```

Run tests with:

```bash
pixi run pytest
```

# Directory Structure

```bash
.
├── benchmarks
├── config
├── data
├── environment.yml
├── LICENSE
├── notebooks
│   └── inspecting_parquet.py
├── pixi.lock
├── pixi.toml
├── pyproject.toml
├── README.md
├── scripts
│   └── generate_mock_data.py
├── sharp_dag.html
├── sharp_pipeline.html
├── src
│   └── sharp
│       ├── __init__.py
│       ├── config.py
│       ├── extract_embeddings.py
│       ├── io.py
│       └── model_management.py
└── tests
    ├── conftest.py
    ├── test_extract_embeddings.py
    ├── test_generate_mock_data.py
    ├── test_io.py
    └── test_model_management.py
```

# Currently Working on

NOW

- prototyping embeddings with tests

NEXT

- prototyping benchmarks
