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
├── .env.example                ← template for machine-specific settings (copy → .env, gitignored)
├── pyproject.toml              ← editable install: `pip install -e .`
├── pixi.toml                   ← environment (use pixi, not conda/pip directly)
├── src/sharp/
│   ├── __init__.py
│   ├── config.py               ← paths + `.env` loading + config dataclasses (DONE)
│   ├── io.py                   ← all data types + file I/O (DONE)
│   ├── metrics.py              ← pure metric math (DONE)
│   └── evaluate.py             ← benchmark step orchestration (DONE)
├── scripts/
│   ├── generate_mock_data.py           ← synthetic proteins for embedding step (DONE)
│   ├── generate_mock_benchmark_data.py ← synthetic predictions + GT for benchmark (DONE)
│   └── prepare_mibig_ground_truth.py   ← parse MiBIG 4.0 JSON → ground_truth.tsv (DONE)
├── tests/
│   ├── conftest.py
│   ├── test_config.py
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

**`.env` is for machine identity, not pipeline behaviour.** `config.py` reads an optional gitignored `.env` at the project root (template: `.env.example`) to resolve the data-dir constants and `DEFAULT_DEVICE`. Precedence: real env var → `.env` → in-repo default.

| Belongs in `.env` | Never in `.env` |
|---|---|
| `SHARP_DATA_ROOT` (+ `SHARP_{RAW,INTERIM,PROCESSED,MOCK}_DIR`) | Thresholds (`min_cluster_frac`, `min_p_bgc`, e-value cutoffs) |
| `SHARP_DEVICE` — hardware physically present | `model_name`, `batch_size`, `max_length` |
| Credentials, once a step actually needs one | Anything that changes what a run *means* |

Reason: a benchmark number must be reproducible from its command line. If a threshold could come from an unversioned file, the same command could produce different numbers on the laptop and the server. When adding a new step, put its knobs in the config dataclass and on the CLI — reach for `.env` only if the value describes the machine rather than the experiment.

Paths are resolved once at import, so tests that need different values reload the module under a patched environment (see `tests/test_config.py`). Don't add `python-dotenv` — the loader in `config.py` is intentionally ~15 lines.

**Each pipeline step = one module** with a `run(cfg: StepConfig) -> None` function and a `build_parser() -> argparse.ArgumentParser` function. Entry point: `python -m sharp.<step>`.

**When you add a new file (script or module), do these two things in the same change:**
1. **Update the directory-structure tree** in `README.md` (`## Directory Structure`) so it stays accurate.
2. **Document its usage** in both `README.md` (a runnable command block, like the "Preparing MiBiG / BGC Atlas Database" sections) and `CLAUDE.md`. The single source of truth is the file's own module docstring — mirror its `Usage:` block; don't invent new invocations. Keep the three (docstring, README, CLAUDE.md) in sync.

**Side-effect isolation:** `metrics.py` is pure (no I/O, no logging). `io.py` owns disk. Orchestration modules (`evaluate.py`, `extract_embeddings.py`, etc.) call both and log.

**Streaming writes for large files.** Parquet is written batch-by-batch via `pq.ParquetWriter` context manager. Never accumulate all rows in memory.

**Tests mirror `src/sharp/`** one-to-one. Test file for `sharp/foo.py` → `tests/test_foo.py`. Scripts tested in `tests/test_<script_name>.py` with `sys.path` injection (see existing examples).

**Monkeypatching rule:** patch on the *importing* module, not the source. If `extract_embeddings.py` does `from sharp.model_management import Embedder`, patch `sharp.extract_embeddings.Embedder`, not `sharp.model_management.Embedder`.

---

## What is DONE (with tests)

| Module / Script | Responsibility | Tests |
|---|---|---|
| `sharp/config.py` | Paths (env-overridable), `.env` loading, `DEFAULT_DEVICE`, `EmbeddingConfig`, `EvaluateConfig` | `test_config.py` |
| `sharp/io.py` | `ProteinRecord`, `PredictedRegion`, `KnownCluster`; FASTA r/w; parquet r/w; TSV r/w; JSON w | `test_io.py` |
| `sharp/model_management.py` | ESM-2 registry, device selection, `residue_mean_pool`, `Embedder`, `ensure_model_available` | `test_model_management.py` |
| `sharp/extract_embeddings.py` | Embedding extraction step: load FASTA → embed → write parquet | `test_extract_embeddings.py` |
| `sharp/metrics.py` | `overlap_bp`, `merge_intervals`, `covered_bp`, `matches`, `reciprocal_overlap`, `MatchCriterion`, `evaluate_predictions`, `BenchmarkResult` | `test_metrics.py` |
| `sharp/evaluate.py` | Benchmark step: load predictions + GT (+ optional `--contigs` scope) → compute metrics → write JSON | `test_evaluate.py` |
| `scripts/generate_mock_data.py` | Synthetic proteins → FASTA (for embedding step smoke tests) | `test_generate_mock_data.py` |
| `scripts/generate_mock_benchmark_data.py` | Synthetic predictions + GT with controlled overlap (for benchmark smoke tests) | `test_evaluate.py` (integration) |
| `scripts/prepare_mibig_ground_truth.py` | MiBIG 4.0 JSON dir → `ground_truth.tsv`; handles 3.x fallback; `--inspect` mode; **rejects unusable accessions** (protein `WP_`/`NP_`, assembly `GCA_`/`ASMnnnvn`, WGS master `...01000000`) and spans <500 bp, and **collapses duplicate loci** filed under two cluster ids | `test_prepare_mibig.py` |
| `scripts/prepare_bgcatlas_ground_truth.py` | BGC Atlas `.gbk` dump → `bgcatlas_ground_truth.tsv` (secondary/noisy GT) | `test_prepare_bgcatlas.py` |
| `scripts/convert_antismash_to_parquet.py` | antiSMASH `sequence.json` → `predictions.parquet`; no coord conversion; `--inspect` mode | `test_convert_antismash.py` |
| `scripts/convert_deepbgc_to_parquet.py` | DeepBGC `.bgc.tsv` → `predictions.parquet`; no coord conversion; `--inspect` mode | `test_convert_deepbgc.py` |
| `scripts/convert_gecco_to_parquet.py` | GECCO `.clusters.tsv` → `predictions.parquet`; `start-1` coord conversion; `--inspect` mode | `test_convert_gecco.py` |
| `scripts/run_antismash.sbatch` | Slurm job: antiSMASH on the benchmark genome. CPU-only like DeepBGC. It takes `--cpus` and hands it to its own module scheduler, but that scheduler parallelises very little, so a wide allocation is wasted. **Sizing measured** (`seff` on array job 45315 index 1, which ran this script's own default genome `AL645882.2`: 5.40% CPU efficiency of 16 cores, 1.62 GB peak) → 4 cores / 4G, matching `run_antismash_array.sbatch`; walltime stays 12h since `$1` may be an arbitrary genome. Paths are explicit at the top of the file | shell, no test |
| `scripts/run_deepbgc.sbatch` | Slurm job: DeepBGC on the benchmark genome. CPU-only on a GPU-free node (the `python=3.7` env predates CUDA-capable TF). Sized from `seff` on job 42995: 0.86 cores, 1.71 GB, 32 min → 2 cores / 8G / 2h. hmmscan dominates the runtime but DeepBGC does not thread it, so the pipeline is serial. Paths are explicit at the top of the file | shell, no test |
| `scripts/run_antismash_array.sbatch` | Slurm **job array**: antiSMASH over the whole benchmark set, one task per genome (index → line of `analyzed_contigs.txt`). Resumable (skips genomes with existing output), per-genome output dir. Sizing **measured** (`seff` on job 45315, indices 1–2, 2026-08-19): ~3 min wall, 5.4%/7.8% CPU efficiency of 16 cores (~1.3 cores of real parallelism), 1.6 GB peak → 4 cores / 4G / 1h. antiSMASH ignores most of `--cpus`, the same way DeepBGC did (sized 8, measured 0.86) | shell, no test |
| `scripts/run_deepbgc_array.sbatch` | Slurm **job array**: DeepBGC over the benchmark set, same shape. Per-task sizing carried from the measured single-genome run (2 cores / 8G / 2h); DeepBGC is single-threaded so the per-task shape does not change with the array | shell, no test |
| `scripts/merge_predictions.py` | Per-genome baseline outputs → one `predictions.parquet`. Reuses each `convert_<tool>_to_parquet.py` (no reparsing — format assumptions stay in the tested converters); `--tool antismash\|deepbgc\|gecco`; `--contigs` reports genomes that produced no output at all, which would otherwise be scored as "tool found nothing"; `--inspect` mode | `test_merge_predictions.py` |
| `scripts/download_benchmark_genomes.sh` | Batch downloader: `benchmark_genomes.tsv` → `data/raw/genomes/<ACC>.fasta`, one per genome. Resumable (skips valid existing files, re-fetches truncated ones), retries with backoff, and verifies each FASTA header equals the expected accession — a mismatch would silently score that genome zero. Fetches by nucleotide accession, **not** from the cluster's NCBI mirror (that mirror is *assembly*-indexed while the GT is *nucleotide*-keyed, and no index file bridges the two — see `docs/NCBI_MIRROR.md`) | shell, no test |
| `scripts/_fetch_nuccore.sh` | Sourced helper (not executed): `fetch_nuccore` (efetch + retry + "is this really FASTA" validation) and `fasta_contig_ids`. Shared by both downloaders so the rules live in one place | shell, no test |
| `scripts/download_genome.sh` | NCBI nuccore accession → `data/raw/<ACC>.fasta` + `data/interim/analyzed_contigs.txt`; defaults to `AL645882.2` | shell, no test |
| `scripts/select_benchmark_genomes.py` | Ground truth → benchmark genome set: drops BGC-only deposits (`--min-length`), merges RefSeq/GenBank twins, ranks, and emits `benchmark_genomes.tsv` + `analyzed_contigs.txt` + `benchmark_ground_truth.tsv` (contigs normalized onto the primary accession). NCBI lengths cached in `data/interim/record_lengths.tsv`; `--inspect`/`--offline` modes | `test_select_benchmark_genomes.py` |
| `scripts/parquet_to_tsv.py` | Generic dump: any pipeline parquet file → TSV; list-typed columns (e.g. `embeddings.parquet`'s `embedding` vector) comma-joined per cell; `--inspect` mode | `test_parquet_to_tsv.py` |

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

**Benchmark scope caveat (verified 2026-07-29).** Recall is measured only over
ground-truth clusters on contigs the tool was actually run on — `evaluate.py`
takes `--contigs` for this. Ground truth spans a database while a run spans one
assembly, so without scoping, recall is capped by the ratio between them: on a
real run of all three baselines against `streptomyces_ground_truth.tsv` (430
clusters / 363 contigs), the analyzed contig `AL589148.1` carried exactly **one**
coordinate-resolved cluster, and every tool scored recall ≤ 0.002. **Pass the same
`--contigs` file to every tool in a comparison** — omitted, the scope is inferred
from the predictions, which is optimistic (a contig analyzed but not called on
drops out of the denominator) and logs a warning.

BGC Atlas results should be reported with a caveat in any paper/presentation:
benchmark numbers on BGC Atlas are optimistic by nature (the positive labels are
themselves predictions, so agreement with them doesn't prove correctness).

**MiBIG 4.0 coordinate-coverage caveat (verified 2026-07-07).** ~45% of all MiBIG
4.0 entries — and **478 of 905 (53%) of *Streptomyces* entries** — store their
locus as `location: {from: 0, to: 0}`, i.e. the compound is characterized but the
genomic coordinates are unknown. `prepare_mibig_ground_truth.py` correctly drops
these (a coordinate-based benchmark can't score a cluster with no interval; the
drop count is logged as "N entries had no locus with usable coordinates"). The
resulting *Streptomyces* ground truth is **~430 loci from 427 clusters, not ~900**.
Two consequences to report in any paper/presentation:
- The recall **denominator is ~half** of MiBIG's *Streptomyces* content by design.
- The dropped half is **not random** — it skews toward older, compound-first
  submissions (dropped IDs cluster in the low `BGC00000xx` range), so the benchmark
  over-represents well-characterized PKS/NRPS clusters. This affects *every* tool
  (S(H)ARP, antiSMASH, DeepBGC) equally, so it doesn't bias the *comparison* — but
  it does mean absolute recall numbers are "recall over coordinate-resolved MiBIG,"
  not "recall over all known *Streptomyces* BGCs."

### Baseline integration — converters, not wrappers

**S(H)ARP never invokes the baseline tools.** antiSMASH, DeepBGC, and GECCO each
install into their own isolated pixi env under `~/.local/src/<tool>/` (via
`scripts/setup_<tool>.sh`) — they have mutually incompatible dependencies and
must stay isolated. You run each tool yourself (its own env, or HPC, or a
container); S(H)ARP only parses the *output files* it leaves behind.

**Setup scripts and `.env`.** All three `setup_<tool>.sh` scripts source
`scripts/_load_env.sh` (sourced, not executed) to read their write locations from
`.env`: `TOOLS_INSTALL_DIR` (where the per-tool pixi envs go, formerly a hardcoded
`~/.local/src` repeated in all three) and, for antiSMASH/DeepBGC,
`ANTISMASH_DOWNLOADS_DIR` / `DEEPBGC_DOWNLOADS_DIR` for the ~10GB / ~3GB reference
databases. The helper uses `set -a` so the values are exported — DeepBGC has no
path flag and reads `DEEPBGC_DOWNLOADS_DIR` from the environment itself; each
script also applies a `: "${VAR:=default}"` fallback so it works with no `.env`.
Note the README's `cd ~/.local/src/<tool>` examples assume the default
`TOOLS_INSTALL_DIR`; the scripts' closing `echo`s interpolate the real value.
These two keys are the exception to the "no shell expansion" rule for `.env`: they
use `${DATABASES:-$HOME/.local/share}`, which bash resolves when sourcing but
`config.py`'s minimal parser does not (`os.path.expandvars` has no `:-` support,
so Python sees the literal string — inert, since no Python step reads them).
`DATABASES` is intentionally *not* defined in `.env.example`: the server exports
it, a laptop falls back to `~/.local/share`, so neither machine needs an edit.
Keep shell-only keys in that block and `SHARP_*` keys — read by Python — above it.

So each baseline gets one **converter** script (not a subprocess wrapper):

```
scripts/convert_<tool>_to_parquet.py --input <tool output> --output <predictions.parquet>
```

Runs entirely in the S(H)ARP env, no external binary, no tool-path config.
Each converter isolates every tool-format assumption (column names, coordinate
base) in one clearly-marked block and provides an `--inspect` mode that prints a
real output file's structure — verify the schema against actual output before
trusting the parser (same pattern as `prepare_mibig_ground_truth.py`).

**Coordinate base is tool-specific — verified per tool 2026-07-15 against a real
run (`antismash 8.0.4`, `deepbgc`, `gecco 0.10.3` on the same input FASTA,
`AL589148.1`).**

| Tool | `p_bgc` source | Coordinate base | Conversion |
|---|---|---|---|
| antiSMASH | none → set `1.0` | 0-based half-open (verified) | none |
| DeepBGC | `deepbgc_score` | 0-based half-open (verified) | none — refutes old hypothesis |
| GECCO | `average_p` | 1-based inclusive (verified) | `start - 1`, `end` unchanged — confirms old hypothesis |

Evidence (span = `end - start` from the TSV/JSON row, cross-checked against the
matching region/cluster `.gbk` LOCUS bp length, across every row in each output —
not just one — since a single row can't distinguish the two conventions if
mis-signed. If `span == LOCUS bp`, the source is 0-based half-open. If
`span == LOCUS bp - 1`, the source is 1-based inclusive (the `.gbk` extraction
naturally has `LOCUS bp = end - start + 1` bases for an inclusive interval):
- **antiSMASH**: `sequence.json` region `location` string is `"[201195:222794](+)"`
  (not plain ints — needs regex parsing). Both regions checked: span exactly equals
  the matching `region00N.gbk` LOCUS bp (`21599`/`21599`, `28972`/`28972`). 0-based
  half-open, no conversion. Same value also available as `Orig. start`/`Orig. end`
  in the region `.gbk` COMMENT block.
- **DeepBGC**: `out.bgc.tsv` columns are `nucl_start`/`nucl_end` (not `start`/`end`).
  All 5 rows checked against `out.bgc.gbk`'s 5 LOCUS records: span exactly equals
  LOCUS bp every time (`10290/10290`, `225/225`, `2679/2679`, `5760/5760`,
  `94301/94301`). 0-based half-open, no conversion needed.
- **GECCO**: `sequence.clusters.tsv` columns are `start`/`end`. All 5 clusters
  checked against their `_cluster_N.gbk` LOCUS bp: span is `LOCUS bp - 1` every
  time (e.g. cluster_1 span `33568` vs LOCUS `33569`; cluster_3 span `44782` vs
  LOCUS `44783`). Consistent off-by-one across every row confirms 1-based
  inclusive, not 0-based half-open — convert with `start - 1`, `end` as-is, same
  as the MiBIG ingest pattern.

**Caution for implementation:** a coordinate check on a single row can look
consistent with either convention if you only compare one direction of the
off-by-one; always check the full span-vs-LOCUS relationship (`==` vs `== -1`)
across multiple rows before trusting the parser, as done above.

Tests parse a small, checked-in, real (trimmed) output fixture per tool — no tool
execution in the suite (that would break env isolation). Same approach as
`test_prepare_mibig.py`.

**`scripts/prepare_bgcatlas_ground_truth.py`** ✅ done (2026-07-07)
Parses the BGC Atlas `complete-bgcs` dump — 204,661 antiSMASH-produced `.gbk`
files, one region per file (downloaded by `scripts/download_bgc-atlas.sh`, DVC-managed
under `data/raw/complete-bgcs/`). Output: `data/raw/bgcatlas_ground_truth.tsv`
(same schema as `mibig_ground_truth.tsv`). Verified schema facts:
- Genomic coords are the antiSMASH `Orig. start`/`Orig. end` structured-comment
  fields (NOT the region-local LOCUS coords), and are **already 0-based half-open**
  (`end - start == len(seq)` across thousands of files) — so, unlike MiBIG, **no
  coordinate conversion is applied**.
- `cluster_id` = filename stem (unique; includes `.regionNNN`, so region001 and
  region002 on one contig stay distinct). `contig` = `<MGYA assembly>_<rec.id>`
  (assembly-qualified, because `rec.id` alone repeats across assemblies).
- `--limit N` for dev/tests (walks a deterministic subset instead of all 10 GB);
  `--inspect DIR` to re-verify the schema. Tests: `tests/test_prepare_bgcatlas.py`.
Secondary/noisy GT — report alongside MiBIG with the optimism caveat above.

**`scripts/convert_antismash_to_parquet.py`** ✅ written (2026-07-15, verified
against a real `antismash 8.0.4` run; tests: `tests/test_convert_antismash.py`,
fixture: `tests/fixtures/antismash_sequence.json` — trimmed real summary JSON)
Parses `sequence.json` → `data/interim/antismash_predictions.parquet`. Iterate
`data['records'][*]['features']` where `type == "region"`. Per region:
- `contig` = `record['id']`
- `start`/`end` = regex-parsed from the `location` string `"[start:end](strand)"`
  (it is a string, not plain ints) — 0-based half-open, no conversion
- `region_id` = `f"{contig}.region{region_number:03d}"` (mirrors the
  `<contig>.region00N.gbk` filename, so a region row can always be traced back to
  its source `.gbk`; `qualifiers.region_number` alone resets per contig and is not
  globally unique)
- `predicted_class` = `";".join(qualifiers['product'])` (hybrid regions can list
  multiple products, e.g. `["furan", "butyrolactone"]` — join rather than truncate)
- `p_bgc` = `1.0` (antiSMASH is rule-based, no score)
Multi-contig genomes: loop all of `data['records']`, not just the first.

**`scripts/convert_deepbgc_to_parquet.py`** ✅ written (2026-07-15, verified
against a real DeepBGC 0.1.0 run; tests: `tests/test_convert_deepbgc.py`,
fixture: `tests/fixtures/deepbgc_out.bgc.tsv` — real, unmodified output)
Parses `out.bgc.tsv` → `data/interim/deepbgc_predictions.parquet`. Real header
(28 columns) confirms `product_class` exists as documented, but it sits at column
18, after several unrelated columns (`detector_version`, `num_proteins`,
`product_activity`, per-activity probabilities) — don't assume column order.
Per row:
- `contig` = `sequence_id`
- `start`/`end` = `nucl_start`/`nucl_end` (**not** `start`/`end` — different column
  names than GECCO despite the same convention) — 0-based half-open, no conversion
- `region_id` = `bgc_candidate_id` (already unique, e.g. `AL589148.1_31460-41750.1`)
- `p_bgc` = `deepbgc_score`
- `predicted_class` = `product_class` — **frequently empty string** in practice
  (4 of 5 rows in the verification run); downstream consumers must handle a blank
  class, not assume it's always populated.

**`scripts/convert_gecco_to_parquet.py`** ✅ written (2026-07-15, verified
against a real `gecco 0.10.3` run; tests: `tests/test_convert_gecco.py`,
fixture: `tests/fixtures/gecco_sequence.clusters.tsv` — real, unmodified output)
Parses `sequence.clusters.tsv` → `data/interim/gecco_predictions.parquet`. Real
header confirms `sequence_id`, `cluster_id`, `start`, `end`, `average_p`, `max_p`,
`type`, plus per-class probability columns (`nrp_probability`,
`polyketide_probability`, etc.) and `proteins`/`domains` list columns. Per row:
- `contig` = `sequence_id`
- `start`/`end` — **1-based inclusive, requires conversion**: `start - 1`, `end`
  unchanged (verified across all 5 rows against `.gbk` LOCUS lengths — see
  coordinate table above; this is the one tool where CLAUDE.md's original
  hypothesis held)
- `region_id` = `cluster_id` (already unique, e.g. `AL589148.1_cluster_1`)
- `p_bgc` = `average_p`
- `predicted_class` = `type` — **was `"Unknown"` for every row** in the
  verification run (5/5); the per-class probability columns are more informative
  when this happens but v1 keeps `type` as-is (simple, matches `PredictedRegion`
  schema) rather than deriving an argmax class, which is a modeling decision to
  revisit later if `"Unknown"` turns out to dominate real runs too.

### Running a full comparison

**The scaled run (50 genomes, 113 clusters) is the default now** — a
single-genome run caps the recall denominator at 16 clusters. The
single-genome flow below still works for a smoke test; see README
"Benchmarking" → "Scaling up: the 50-genome benchmark" for the full
array-job sequence.

```bash
# 1. Select the genome set from the ground truth. This is not just "sort by
#    cluster count": ~58% of MiBIG records are BGC-only deposits (the record IS
#    the cluster, so every tool scores ~1.0 by construction) and one physical
#    sequence can carry several accessions (NC_003888.3 and AL645882.2 are the
#    same S. coelicolor chromosome, with 15 clusters filed under one and 1 under
#    the other). Emits the set, the scope file, and a contig-normalized GT.
python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/streptomyces_ground_truth.tsv \
    --output-dir data/interim/benchmark_set

# 2. Fetch them (resumable, ~420 Mb). By nucleotide accession, so the FASTA
#    header IS the name the ground truth uses.
scripts/download_benchmark_genomes.sh

# 3. Run each baseline as a job array, then merge (see README), and evaluate
#    against benchmark_ground_truth.tsv + analyzed_contigs.txt from step 1 —
#    NOT the raw MiBIG ground truth, whose contig names are not normalized.

# ── single-genome smoke test ────────────────────────────────────────────────
# Fetch one genome and its --contigs scope file in one step.
scripts/download_genome.sh

# Build ground truth
python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/mibig_ground_truth.tsv --genus Streptomyces

# Run each baseline yourself in its own env (see scripts/setup_<tool>.sh), then
# convert its output — S(H)ARP never invokes the tools:
python scripts/convert_antismash_to_parquet.py \
    --input <antismash output dir/json> \
    --output data/interim/antismash_predictions.parquet

python scripts/convert_deepbgc_to_parquet.py \
    --input <deepbgc .bgc.tsv> \
    --output data/interim/deepbgc_predictions.parquet

python scripts/convert_gecco_to_parquet.py \
    --input <gecco .clusters.tsv> \
    --output data/interim/gecco_predictions.parquet

# Evaluate all against the same ground truth AND the same scope.
# --contigs lists the contigs the tools were run on (one per line, or a .fai);
# every tool must get the same file or the recall denominators differ.
# download_genome.sh already wrote analyzed_contigs.txt. For a genome obtained
# some other way, derive it the same way:
#   grep '^>' <genome.fasta> | cut -c2- | cut -d' ' -f1 > data/interim/analyzed_contigs.txt

for tool in antismash deepbgc gecco; do
    python -m sharp.evaluate \
        --predictions data/interim/${tool}_predictions.parquet \
        --ground-truth data/raw/mibig_ground_truth.tsv \
        --contigs data/interim/analyzed_contigs.txt \
        --output data/processed/benchmark_${tool}.json
done

python -m sharp.evaluate \
    --predictions data/interim/predictions.parquet \
    --ground-truth data/raw/mibig_ground_truth.tsv \
    --contigs data/interim/analyzed_contigs.txt \
    --output data/processed/benchmark_sharp.json
```

Reading the output (`benchmark.json`, see `docs/ARCHITECTURE.md` → "Metrics"):

| block | what it answers |
|---|---|
| `scope` | how much of the ground truth was evaluable, and whether scope was `explicit` or `inferred` |
| `detection` | *did the tool find the BGC?* — `min_cluster_frac` only |
| `reciprocal` | the strict symmetric rule, for comparison |
| `nucleotide` | bp-level agreement; `precision` says how much extra territory was called |
| `boundary` | `median_prediction_coverage` (tightness), split/merge diagnostics |

`matched_prediction_frac` is a **lower bound** on precision, not precision — the
ground truth is incomplete, so an unmatched prediction is unvalidated rather than
wrong. There is deliberately no region-level `precision` or `false_positive` field.

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
# → detection recall=0.700 (14/20), matched 14/19 predictions — matches the
#   generator's --recall-rate 0.7 and 5 injected unmatched predictions

# Verify MiBIG 4.0 JSON schema (do once after download)
pixi run python scripts/prepare_mibig_ground_truth.py \
    --inspect data/raw/mibig_json_4.0

# Build real ground truth
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/mibig_ground_truth.tsv \
    --genus Streptomyces

# BGC Atlas secondary ground truth (noisy — report alongside MiBIG, never alone)
pixi run python scripts/prepare_bgcatlas_ground_truth.py \
    --inspect data/raw/complete-bgcs          # verify schema first
pixi run python scripts/prepare_bgcatlas_ground_truth.py \
    --input-dir data/raw/complete-bgcs \
    --output data/raw/bgcatlas_ground_truth.tsv
#   add --limit N to build against a small subset for dev/tests

# Full competitor comparison (once baseline scripts are written)
# See "Benchmark comparison" section above for full command sequence
```
