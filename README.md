# Projeto

Clique nos links abaixo para ir às páginas HTML interativas:

- [Diagrama DAG (directed acyclic graph)](sharp_dag.html) - Fluxograma do projeto, evidenciando inputs-processos-outputs.
- [Página Descritiva](sharp_pipeline.html) - Praticamente o mesmo conteúdo do Diagrama DAG, mas apresentado com uma interface menos técnica.

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
