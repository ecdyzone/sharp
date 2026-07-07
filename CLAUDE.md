# CLAUDE.md — S(H)ARP Project Context

> Read this before touching any code. Additional detail in `docs/`.

## What this project is

S(H)ARP predicts Biosynthetic Gene Clusters (BGCs) in *Streptomyces* and related actinomycetes using **SARP transcription factors as anchors**. A SARP (Streptomyces Antibiotic Regulatory Protein) is a regulator with an HTH-BTAD domain (± NB-ARC / TPR / AAA / LuxR) that binds heptameric repeats (afsR-box) in BGC promoters.

Differentiator vs. antiSMASH: we use regulatory context + protein language model embeddings, not just biosynthetic enzyme patterns.

See `docs/PIPELINE.md` for the full biological pipeline. See `docs/ARCHITECTURE.md` for module ownership.

---

## Repo layout

```
project_root/
├── CLAUDE.md                   ← you are here
├── pyproject.toml              ← editable install: `pip install -e .`
├── pixi.toml                   ← environment (use pixi, not conda/pip directly)
├── src/sharp/
│   ├── __init__.py
│   ├── config.py               ← paths + config dataclasses (DONE)
│   ├── io.py                   ← all data types + file I/O (DONE)
│   ├── metrics.py              ← pure metric math (DONE)
│   └── evaluate.py             ← benchmark step orchestration (DONE)
├── scripts/
│   ├── generate_mock_data.py           ← synthetic proteins for embedding step (DONE)
│   ├── generate_mock_benchmark_data.py ← synthetic predictions + GT for benchmark (DONE)
│   └── prepare_mibig_ground_truth.py   ← parse MiBIG 4.0 JSON → ground_truth.tsv (DONE)
├── tests/
│   ├── conftest.py
│   ├── test_io.py
│   ├── test_model_management.py
│   ├── test_extract_embeddings.py
│   ├── test_generate_mock_data.py
│   ├── test_metrics.py
│   ├── test_evaluate.py
│   └── test_prepare_mibig.py
└── data/
    ├── raw/          ← immutable inputs (MiBIG dump, downloaded genomes)
    ├── interim/      ← intermediate pipeline artifacts
    ├── processed/    ← final outputs (model.pkl, report.html, benchmark.json)
    └── mock/         ← synthetic data for testing
```

> **Note:** `extract_embeddings.py` and `model_management.py` were implemented and tested but belong in `src/sharp/` — they may already be there if you've been working in this repo. If missing, see `docs/ARCHITECTURE.md` for their specs.

---

## Git workflow

Commit after every meaningful unit of work. Follow **Conventional Commits**:

```
<type>(<scope>): <short description>

[optional body]
```

Types used in this project:

| Type | When |
|---|---|
| `feat` | new pipeline step, new script, new metric |
| `fix` | bug fix in existing code |
| `test` | adding or fixing tests |
| `refactor` | restructuring without behavior change |
| `docs` | CLAUDE.md, docs/, docstrings |
| `chore` | pixi.toml, pyproject.toml, CI config |
| `data` | scripts that produce or transform data files |

Scopes are module names or script names: `io`, `metrics`, `evaluate`,
`extract-embeddings`, `prepare-mibig`, `run-antismash`, etc.

Examples:
```
feat(metrics): add reciprocal_overlap and BenchmarkResult
test(evaluate): add end-to-end orchestration tests
fix(io): handle missing region_id in FASTA header gracefully
feat(prepare-mibig): add --inspect mode for schema verification
data(prepare-mibig): build Streptomyces ground truth from MiBIG 4.0
refactor(extract-embeddings): extract residue_mean_pool as pure function
docs(claude): add benchmark comparison section and backlog tier 0
chore: add antismash and deepbgc to pixi.toml
```

Keep commits atomic — one logical change per commit. Don't batch unrelated
changes (e.g. don't fix a bug and add a feature in the same commit).

---

## Environment

```bash
pixi run python ...          # always use pixi, not bare python
pixi run pytest              # run tests
pixi run python -m sharp.evaluate --help
```

Package is installed editable: `import sharp.io` works anywhere.

---

## Conventions — read before writing any code

**Coordinates:** 0-based half-open `[start, end)` everywhere. MiBIG uses 1-based inclusive `[from, to]`; `prepare_mibig_ground_truth.py` converts on ingest (`start - 1`, `end` unchanged). Never store 1-based coords in any data type or file.

**Data types live in `io.py`**, not in a separate `types.py`. The rule: a type lives in the module that constructs it. Extract to `types.py` only if a third module needs it without going through io.

**Config dataclasses in `config.py`.** One frozen dataclass per pipeline step (e.g. `EmbeddingConfig`, `EvaluateConfig`). Steps receive a config object, not loose `**kwargs`.

**Each pipeline step = one module** with a `run(cfg: StepConfig) -> None` function and a `build_parser() -> argparse.ArgumentParser` function. Entry point: `python -m sharp.<step>`.

**Side-effect isolation:** `metrics.py` is pure (no I/O, no logging). `io.py` owns disk. Orchestration modules (`evaluate.py`, `extract_embeddings.py`, etc.) call both and log.

**Streaming writes for large files.** Parquet is written batch-by-batch via `pq.ParquetWriter` context manager. Never accumulate all rows in memory.

**Tests mirror `src/sharp/`** one-to-one. Test file for `sharp/foo.py` → `tests/test_foo.py`. Scripts tested in `tests/test_<script_name>.py` with `sys.path` injection (see existing examples).

**Monkeypatching rule:** patch on the *importing* module, not the source. If `extract_embeddings.py` does `from sharp.model_management import Embedder`, patch `sharp.extract_embeddings.Embedder`, not `sharp.model_management.Embedder`.

---

## What is DONE (with tests)

| Module / Script | Responsibility | Tests |
|---|---|---|
| `sharp/config.py` | Paths, `EmbeddingConfig`, `EvaluateConfig` | — |
| `sharp/io.py` | `ProteinRecord`, `PredictedRegion`, `KnownCluster`; FASTA r/w; parquet r/w; TSV r/w; JSON w | `test_io.py` |
| `sharp/model_management.py` | ESM-2 registry, device selection, `residue_mean_pool`, `Embedder`, `ensure_model_available` | `test_model_management.py` |
| `sharp/extract_embeddings.py` | Embedding extraction step: load FASTA → embed → write parquet | `test_extract_embeddings.py` |
| `sharp/metrics.py` | `overlap_bp`, `reciprocal_overlap`, `evaluate_predictions`, `BenchmarkResult` | `test_metrics.py` |
| `sharp/evaluate.py` | Benchmark step: load predictions + GT → compute metrics → write JSON | `test_evaluate.py` |
| `scripts/generate_mock_data.py` | Synthetic proteins → FASTA (for embedding step smoke tests) | `test_generate_mock_data.py` |
| `scripts/generate_mock_benchmark_data.py` | Synthetic predictions + GT with controlled overlap (for benchmark smoke tests) | `test_evaluate.py` (integration) |
| `scripts/prepare_mibig_ground_truth.py` | MiBIG 4.0 JSON dir → `ground_truth.tsv`; handles 3.x fallback; `--inspect` mode | `test_prepare_mibig.py` |

---

## What is NOT YET IMPLEMENTED

Implement these in order. Each is a pipeline step; each gets its own module + config dataclass + tests.

### 1. `sharp/annotate.py` — genome annotation
**Input:** `data/raw/genome.fasta`
**Output:** `data/interim/proteins.faa`, `data/interim/annotated.gbk`, `data/interim/genes.gff`
**Tool:** Bakta (shell out via `subprocess`)
**Config:** `AnnotateConfig(input_path, output_dir, threads, min_contig_length)`
**Notes:** Bakta writes its own output dir. Wrapper should copy/symlink the three output files to canonical interim paths. Validate that all three files exist after run.

### 2. `sharp/detect_sarp.py` — SARP detection by HMM
**Input:** `data/interim/proteins.faa`, `data/raw/sarp_models.hmm`
**Output:** `data/interim/anchors_sarp.tsv` (columns: `protein_id, contig, start, end, strand, score, type`)
**Tool:** hmmscan (shell out)
**Logic:** Parse hmmscan tblout format. Add `type=SARP`. Filter by e-value threshold (default 1e-5).
**Config:** `DetectSarpConfig(proteins_path, hmm_path, output_path, evalue_threshold)`

### 3. `sharp/detect_heptarepeats.py` — motif search in DNA
**Input:** `data/raw/genome.fasta`, FIMO motif file (afsR-box PWM)
**Output:** `data/interim/anchors_heptarepeats.tsv` (same columns as above, `type=heptarepeat`)
**Tool:** FIMO from MEME suite (shell out)
**Logic:** Parse FIMO TSV output. Coords are already 0-based in FIMO output — verify this on real output before assuming.
**Config:** `DetectHeptarepeatsConfig(genome_path, motif_path, output_path, pvalue_threshold)`

### 4. `sharp/merge_anchors.py` — unify anchor tables
**Input:** `anchors_sarp.tsv`, `anchors_heptarepeats.tsv`
**Output:** `data/interim/anchors.tsv`
**Logic:** Concatenate, deduplicate by position, sort by contig+start. Pure function, minimal I/O.

### 5. `sharp/extract_neighborhood.py` — genomic window extraction
**Input:** `data/interim/anchors.tsv`, `data/interim/annotated.gbk`
**Output:** `data/interim/neighborhoods.tsv`, `data/interim/neighborhood_proteins.faa`, `data/interim/neighborhood_dna.fna`
**Logic:**
- Window: ±20 genes **or** ±20 kb, whichever is larger
- Merge anchors within 50 kb of each other into one region (avoids double-counting)
- `neighborhoods.tsv` columns: `region_id, contig, start, end, anchor_ids, n_proteins`
- FASTA headers: `>PROTEIN_ID region_id=R001` (required by `parse_fasta`)
**Config:** `ExtractNeighborhoodConfig(anchors_path, genbank_path, output_dir, window_genes, window_bp, merge_distance)`
**Library:** BioPython `SeqIO` for GenBank parsing.

### 6. `sharp/annotate_domains.py` — Pfam domain annotation
**Input:** `data/interim/neighborhood_proteins.faa`, `data/raw/pfam_models.hmm`
**Output:** `data/interim/domains.tsv` (columns: `protein_id, region_id, domain, e_value, start, end`)
**Tool:** hmmscan (shell out, domtblout format)
**Logic:** Parse domtblout. One row per domain hit per protein. Filter by e-value. Note: `region_id` must be recovered from the FASTA header (use `parse_fasta` then build a `protein_id → region_id` map).
**Config:** `AnnotateDomainsConfig(proteins_path, hmm_path, output_path, evalue_threshold)`

### 7. `sharp/extract_kg_features.py` — knowledge graph context features ⭐ yours
**Input:** `data/interim/neighborhoods.tsv`, `data/interim/domains.tsv`, `data/raw/kg.gpickle`
**Output:** `data/interim/kg_features.parquet` (columns: `region_id, n_similar_clusters, modal_class, has_large_sarp, ...`)
**Logic:** For each region, query the KG for clusters with similar domain architecture. Extract tabular features.
**Config:** `KgFeaturesConfig(neighborhoods_path, domains_path, kg_path, output_path)`
**Note:** KG is built by a separate one-time script (`scripts/build_kg.py`) — see `docs/PIPELINE.md`.

### 8. `sharp/train.py` — classifier training ⭐ yours
**Input:** `data/interim/embeddings.parquet`, `data/interim/domains.tsv`, `data/interim/kg_features.parquet`, `data/raw/mibig_ground_truth.tsv`
**Output:** `data/processed/model.pkl`, `data/processed/metrics.json`, `data/processed/feature_importance.tsv`
**Logic:**
- Aggregate proteins → regions: mean-pool embeddings by `region_id`; one-hot domains by `region_id`; KG features already per-region
- Label regions: positive if overlaps a MiBIG cluster (use `reciprocal_overlap` from `metrics.py`), negative otherwise
- Train LightGBM with k-fold CV (k=5)
- Serialize with `joblib.dump`
**Config:** `TrainConfig(embeddings_path, domains_path, kg_features_path, ground_truth_path, output_dir, n_folds, min_overlap_frac)`

### 9. `sharp/predict.py` — inference
**Input:** `data/processed/model.pkl`, same features as train
**Output:** `data/interim/predictions.parquet` (columns: `region_id, contig, start, end, p_bgc, predicted_class`)
**Logic:** Load model, run feature pipeline (same aggregation as train), predict. Output format must match `PredictedRegion` schema in `io.py`.
**Config:** `PredictConfig(...)`

### 10. `sharp/filter.py` — heuristic post-filter
**Input:** `data/interim/predictions.parquet`, `data/interim/domains.tsv`
**Output:** `data/interim/filtered_predictions.parquet`
**Logic (rules to start with):**
- Drop regions where `p_bgc < threshold` (default 0.5)
- Drop regions containing only ribosomal protein domains
- Drop regions with fewer than 3 proteins
**Config:** `FilterConfig(predictions_path, domains_path, output_path, p_bgc_threshold)`

### 11. `sharp/generate_report.py` — HTML report
**Input:** `data/interim/filtered_predictions.parquet`, `data/interim/neighborhoods.tsv`, `data/interim/domains.tsv`
**Output:** `data/processed/report.html`
**Tool:** Jinja2
**Logic:** One section per predicted BGC: region coordinates, class, p_bgc score, domain architecture diagram (SVG or simple HTML table).

---

## Benchmark comparison — competitor baselines

**Priority: high.** The team wants S(H)ARP benchmarked against antiSMASH and
DeepBGC (at minimum). Check recent literature for others.

The architecture is already correct: any tool's output can be converted to
`predictions.parquet` and passed through `evaluate.py` unchanged. Each tool
gets one conversion script in `scripts/`.

### Ground truth sources

| Source | Reliability | Use as GT |
|---|---|---|
| MiBIG 4.0 | ✅ Manually curated | Primary — always use |
| BGC Atlas | ⚠️ Computationally predicted, no manual curation | Secondary — noisier, interpret separately |

BGC Atlas results should be reported with a caveat in any paper/presentation:
benchmark numbers on BGC Atlas are optimistic by nature (the positive labels are
themselves predictions, so agreement with them doesn't prove correctness).

### Baseline conversion scripts (not yet written)

**`scripts/prepare_bgcatlas_ground_truth.py`**
Similar to `prepare_mibig_ground_truth.py`. BGC Atlas distributes as TSV or GFF.
Inspect the actual format before writing the parser. Output: `data/raw/bgcatlas_ground_truth.tsv`
(same schema as `mibig_ground_truth.tsv`).

**`scripts/run_antismash_baseline.py`**
Runs antiSMASH on a genome and converts its output to `predictions.parquet`.
antiSMASH output format: per-region GenBank files + a JSON summary (`regions.js`
or `<genome>.json`). The JSON summary is the cleanest parse target.
Key fields to extract: `records[].features[]` where `type == "region"`,
with `qualifiers.region_number`, `qualifiers.product`, and coordinates from `location`.
Output: `data/interim/antismash_predictions.parquet` (same schema as `PredictedRegion`).
Set `p_bgc=1.0` for all antiSMASH predictions (it doesn't produce scores).

**`scripts/run_deepbgc_baseline.py`**
Runs DeepBGC and converts its output. DeepBGC produces a `.bgc.tsv` with columns
including `sequence_id`, `start`, `end`, `deepbgc_score`, `product_class`.
Output: `data/interim/deepbgc_predictions.parquet`.
Use `deepbgc_score` as `p_bgc`.

### Running a full comparison

```bash
# Build ground truth
python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/mibig_ground_truth.tsv --genus Streptomyces

# Run baselines on the same genome
python scripts/run_antismash_baseline.py \
    --genome data/raw/genome.fasta \
    --output data/interim/antismash_predictions.parquet

python scripts/run_deepbgc_baseline.py \
    --genome data/raw/genome.fasta \
    --output data/interim/deepbgc_predictions.parquet

# Evaluate all three against the same ground truth
python -m sharp.evaluate \
    --predictions data/interim/antismash_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --output data/processed/benchmark_antismash.json

python -m sharp.evaluate \
    --predictions data/interim/deepbgc_predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --output data/processed/benchmark_deepbgc.json

python -m sharp.evaluate \
    --predictions data/interim/predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --output data/processed/benchmark_sharp.json
```

---

## Deliberate omissions (do NOT add unless asked)

These were scoped out of the MVP intentionally. Add only when the feature is explicitly needed.

| Feature | Where it belongs | When to add |
|---|---|---|
| Evo nucleotide embeddings | `extract_embeddings.py` | After ESM-2 baseline is benchmarked |
| ESM-IF / Foldseek structural embeddings | `extract_embeddings.py` | After Evo |
| GNN embeddings from KG | `extract_kg_features.py` | After tabular KG features are validated |
| Multi-class classification | `train.py` | After binary classifier AUROC > 0.85 |
| Ensemble of modality-specific models | `train.py` | After multi-class |
| Asymmetric overlap thresholds | `metrics.py` | If team decides one threshold is insufficient |
| AUROC in benchmark | `metrics.py` | Once `predict.py` scores ALL candidates (not just positives) |
| Per-class benchmark breakdown | `evaluate.py` | When team asks "why is NRPS recall low?" |
| DeepBGC / antiSMASH per-class breakdown | `evaluate.py` extension | When team asks "which BGC class is each tool best at?" |
| Resumable embedding extraction | `extract_embeddings.py` | When datasets exceed ~100k proteins |
| fp16/bf16 inference | `model_management.py` | When running on GPU cluster |
| `BaseStep` abstraction | new `pipeline.py` | When 3+ steps need to share boilerplate |
| Logging to file / JSON structured logs | `config.py` | When deploying beyond laptop |
| Docker / Snakemake / Nextflow | new | When moving to HPC |

---

## Key domain facts (don't get these wrong)

- **SARP = Streptomyces Antibiotic Regulatory Protein.** HTH-BTAD domain is obligatory. Larger SARPs carry NB-ARC, TPR, AAA, or LuxR domains additionally.
- **afsR-box** = the heptameric DNA repeat that SARPs bind. FIMO searches for these.
- **BGC** = Biosynthetic Gene Cluster. Classes: T1PKS, T2PKS, NRPS, terpene, RiPP, etc.
- **MiBIG** = ground truth database. We use v4.0. Coordinates are 1-based inclusive.
- **`neighborhood_dna.fna`** is generated but not consumed by any MVP step. Reserved for Evo (nucleotide language model).

---

## Running the benchmark (current state)

```bash
# Smoke test with synthetic data (no real genome needed)
pixi run python scripts/generate_mock_benchmark_data.py \
    --n-clusters 20 --recall-rate 0.7 --n-false-positives 5
pixi run python -m sharp.evaluate \
    --predictions data/mock/predictions.parquet \
    --ground-truth data/mock/ground_truth.tsv \
    --output data/processed/benchmark.json

# Verify MiBIG 4.0 JSON schema (do once after download)
pixi run python scripts/prepare_mibig_ground_truth.py \
    --inspect data/raw/mibig_json_4.0

# Build real ground truth
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/mibig_ground_truth.tsv \
    --genus Streptomyces

# Full competitor comparison (once baseline scripts are written)
# See "Benchmark comparison" section above for full command sequence
```
