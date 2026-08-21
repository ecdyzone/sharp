# S(H)ARP

S(H)ARP predicts Biosynthetic Gene Clusters (BGCs) in *Streptomyces* and related
actinomycetes using **SARP transcription factors as anchors**, combining regulatory
context with protein language model embeddings rather than biosynthetic enzyme
patterns alone.

> GitHub Pages: <https://ecdyzone.github.io/sharp>  
> GitHub repository: <https://github.com/ecdyzone/sharp>  
> DAG diagram: <https://ecdyzone.github.io/sharp/docs/sharp_dag.html>  

Click the links below for the interactive HTML pages:

- [DAG diagram (directed acyclic graph)](docs/sharp_dag.html) - Project flowchart, showing inputs-processes-outputs.
- [Descriptive page](docs/sharp_pipeline.html) - Practically the same content as the DAG diagram, but presented with a less technical interface.

## Table of Contents

- [Quickstart](#quickstart)
- [Setting up](#setting-up)
- [Pipeline steps](#pipeline-steps)
- [Ground truth](#ground-truth)
- [Benchmarking](#benchmarking)
  - [Run once, slice many](#run-once-slice-many--the-benchmarking-approach)
- [Utilities](#utilities)
- [Tests](#tests)
- [Directory Structure](#directory-structure)
- [Currently Working on](#currently-working-on)

## Quickstart

Just cloned the repo? This runs end to end with no data download:

```bash
pixi install                                   # set up the environment
pixi run pytest                                # confirm it works

# Benchmark against synthetic data — clusters and predictions that overlap
# by construction, so the expected number is known in advance.
pixi run python scripts/generate_mock_benchmark_data.py \
    --n-clusters 20 --recall-rate 0.7 --n-false-positives 5
pixi run python -m sharp.evaluate \
    --predictions data/mock/predictions.parquet \
    --ground-truth data/mock/ground_truth.tsv \
    --output data/processed/benchmark.json
# → detection recall=0.700 (14/20), matched 14/19 predictions
```

Where to go next:

| If you want to... | Read |
|---|---|
| build the ground truth from MiBiG | [Ground truth](#ground-truth) |
| compare against antiSMASH / DeepBGC / GECCO | [Benchmarking](#benchmarking) |
| run a S(H)ARP pipeline step | [Pipeline steps](#pipeline-steps) |
| move data onto a scratch disk | [Machine-specific settings](#machine-specific-settings-optional-env) |

## Setting up

First run:

```bash
git clone https://github.com/ecdyzone/sharp
cd sharp
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


## Pipeline steps

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

## Ground truth

Every benchmark number is measured against one of these two tables. **MiBiG is the
primary ground truth** — manually curated, and what any reported result should be
based on. **BGC Atlas is secondary and noisy**: its labels are themselves antiSMASH
predictions, so agreement with them does not prove correctness. Report it alongside
MiBiG, never alone.

Build MiBiG first — [Selecting the genome set](#1-selecting-the-genome-set) consumes
its output.

> **Coverage caveat.** ~53% of MiBiG 4.0 *Streptomyces* entries store no genomic
> coordinates, so they cannot be scored and are dropped. The resulting ground truth
> is around 400 loci, not 900, and the dropped half skews toward older compound-first
> submissions. This affects every tool equally, so it does not bias the comparison —
> but absolute recall is "recall over coordinate-resolved MiBiG". See `CLAUDE.md`
> for the full analysis.

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

# Bacteria only — drops fungal/plant/animal entries (22% of the
# coordinate-resolved set: 327 fungal + 26 plant/animal)
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/bacterial_ground_truth.tsv \
    --exclude-eukaryotes
```

MiBiG's taxonomy block carries only `{name, ncbiTaxId}` — no lineage — so
"keep only bacteria" cannot be expressed as a field test. `--exclude-eukaryotes`
applies a curated genus deny-list (`EUKARYOTIC_GENERA`, ~100 genera) against the
first word of the taxonomy name, and logs which genera it skipped. It is a
deny-list rather than an allow-list of bacteria so that an unlisted eukaryote is
*kept* (and visible in the log) rather than a novel bacterium vanishing silently.

**Not every entry is scoreable.** Of MiBiG 4.0's 3,013 entries, **1,363 (45%)
store `location: {from: 0, to: 0}`** — the compound is characterized but the
genomic locus is unknown, and a coordinate-based benchmark has nothing to score
against. These are dropped on ingest and the count is logged. What remains:

| scope | clusters | accessions |
|---|---|---|
| coordinate-resolved, all genera | 1,634 | 1,420 |
| bacteria only (`--exclude-eukaryotes`) | 1,280 | 1,112 |
| *Streptomyces* only (`--genus`) | 414 | 352 |

Absolute recall figures are therefore "recall over coordinate-resolved MiBiG",
not "recall over all known BGCs". The dropped half skews toward older,
compound-first submissions — but it is dropped identically for every tool, so
the *comparison* stays unbiased.

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


## Benchmarking

Read [How evaluation works](#how-evaluation-works) first — the scope rule it
describes governs every command in this section.

### How evaluation works

Every tool, including S(H)ARP itself, is scored the same way: its output is
converted to a `predictions.parquet` and passed to `sharp.evaluate`, which writes a
`benchmark.json`.

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

### Smoke test with mock data

Needs no genome and no ground truth — the generator builds clusters and
predictions that overlap by construction, so the expected recall is known before
the run:

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

### Competitor baselines (antiSMASH / DeepBGC / GECCO)

S(H)ARP does **not** run these tools — each has incompatible dependencies and
installs into its own isolated pixi env via `scripts/setup_<tool>.sh`. You run
the tool yourself, then convert its output to `predictions.parquet` and evaluate
it exactly like S(H)ARP's own predictions.

#### Installing a baseline

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

#### Running and converting one tool

antiSMASH is the worked example; DeepBGC and GECCO differ only where noted below.

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

**DeepBGC** — same shape; S(H)ARP parses `<prefix>.bgc.tsv` instead:

```bash
bash scripts/setup_deepbgc.sh
cd ~/.local/src/deepbgc && pixi run deepbgc pipeline <genome.fasta> --output out

pixi run python scripts/convert_deepbgc_to_parquet.py --inspect out
pixi run python scripts/convert_deepbgc_to_parquet.py \
    --input out --output data/interim/deepbgc_predictions.parquet
```

**GECCO** — same shape, but its `start`/`end` are 1-based inclusive (the one
baseline tool needing a coordinate conversion), which the converter applies
automatically:

```bash
bash scripts/setup_gecco.sh
cd ~/.local/src/gecco && pixi run gecco run --genome <genome.fasta> --output-dir out

pixi run python scripts/convert_gecco_to_parquet.py --inspect out
pixi run python scripts/convert_gecco_to_parquet.py \
    --input out --output data/interim/gecco_predictions.parquet
```

Both then evaluate exactly like step 5 above, with the same `--ground-truth` and
the same `--contigs` file, writing `benchmark_deepbgc.json` / `benchmark_gecco.json`.

#### Running on Slurm

`scripts/run_antismash.sbatch` and `scripts/run_deepbgc.sbatch` wrap the tool
invocation with explicit paths and a preflight check. Paths are set at the top of
each file — edit them for your machine:

```bash
sbatch scripts/run_antismash.sbatch                    # benchmark genome
sbatch scripts/run_antismash.sbatch /path/to/other.fasta

sbatch scripts/run_deepbgc.sbatch                      # benchmark genome
sbatch scripts/run_deepbgc.sbatch /path/to/other.fasta

squeue -u $USER
tail -f deepbgc-scoe-<jobid>.out
```

Both request CPU only, and no GPU. DeepBGC is Prodigal → `hmmscan` vs Pfam (both
CPU-bound, and together the bulk of the runtime) → a small Keras classifier. Only
that last stage could use a GPU and it is a rounding error in the total, so a GPU
allocation would sit idle — leave the GPU nodes for the ESM-2 embedding step.
antiSMASH has no GPU code path at all.

DeepBGC's sizing (2 cores / 8G / 2h) comes from `seff` on a real run of this
script over `AL645882.2`: **0.86 cores, 1.71 GB, 32m24s**. Despite `hmmscan`
dominating the runtime, DeepBGC does not thread it, so the pipeline is serial —
an earlier 8-core allocation ran at 10.73% CPU efficiency.

antiSMASH's sizing (4 cores / 4G) is measured the same way. It accepts `--cpus`
and hands it to its own module scheduler, but that scheduler parallelises very
little in practice: `seff` on array job 45315 index 1 — which ran the same
`AL645882.2` genome this script defaults to — reported 5.40% CPU efficiency of
16 cores and 1.62 GB peak, i.e. ~1.3 cores of real parallelism. Both
`run_antismash.sbatch` and `run_antismash_array.sbatch` are sized from that.
The single-genome script keeps a loose 12h walltime because it accepts an
arbitrary genome as `$1`, where the array runs a known set. Measure your
own runs rather than trusting these numbers on a different genome:

```bash
seff <jobid>
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

### Run once, slice many — the benchmarking approach

**This is how every benchmark is produced.** Run each baseline over the broadest
genome set you are willing to pay for, *once*, then derive as many benchmark
numbers as you like by re-scoping — no recomputation.

This works because the expensive step is keyed by genome, not by experiment:

```
~/projects/<tool>/out_benchmark/<ACCESSION>/     ← one shared pool, keyed by accession
        │                                          array tasks skip genomes already done
        ├─ merge_predictions.py --contigs A ──▶ <tool>_predictions_A.parquet ──▶ benchmark_A_<tool>.json
        ├─ merge_predictions.py --contigs B ──▶ <tool>_predictions_B.parquet ──▶ benchmark_B_<tool>.json
        └─ merge_predictions.py --contigs C ──▶ ...
```

`--contigs` filters **both** the predictions and the ground-truth denominator
(`metrics.py`), so a scope file plus its matching ground truth fully define a
benchmark. Producing a new one costs minutes.

**Two rules that make this safe:**

1. **Never partition the output pool per experiment.** `OUTROOT` is keyed by
   accession precisely so two experiments sharing a genome share one directory —
   that is what makes re-slicing free. Keep the *derived* artifacts separate
   instead, named per scope:
   `benchmark_set_<name>/`, `<tool>_predictions_<name>.parquet`,
   `benchmark_<name>_<tool>.json`.
2. **`--contigs` is mandatory, not optional.** Because the pool holds every
   genome ever run, `merge_predictions.py` without `--contigs` would sweep in
   genomes from unrelated experiments. It is what scopes the pool down to one
   benchmark.

A scope is always a **pair** of files from one `select_benchmark_genomes.py
--output-dir`: `analyzed_contigs.txt` *and* its `benchmark_ground_truth.tsv`.
They must travel together — the ground truth is contig-normalized for that
selection (RefSeq/GenBank twins collapsed onto a primary accession), so pairing
a scope file with the raw MiBiG ground truth silently drops clusters filed under
a twin.

See [TODO.md](TODO.md) for the scenarios queued against the shared pool
(*Streptomyces*, bacteria-only, all genera, BGC-only deposits).

### Worked example: the 50-genome benchmark

The four steps below run in order — each consumes what the previous one wrote.
A single genome caps the recall denominator at ~15 clusters; this set is 50
genomes / 113 clusters. Substitute any other scope by changing `--output-dir`
and the `--ground-truth` it is built from.

#### 1. Selecting the genome set

The ground truth cannot be ranked naively:

- **~58% of MiBiG records are not genomes.** They are BGC-only deposits where
  the record *is* the cluster, so every tool scores detection recall ~1.0 on
  them by construction. `--min-length` drops them.
- **The same sequence appears under several accessions.** `NC_003888.3` and
  `AL645882.2` are both the 8,667,507 bp *S. coelicolor* chromosome, with 15
  clusters filed under one and 1 under the other. They are merged, and the
  ground truth is rewritten onto the primary accession so that 16th cluster is
  not lost.

```bash
# Inspect the selection first: ranking, merges, what got filtered
pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/streptomyces_ground_truth.tsv --inspect

# Write the set (default: top 50 genome-scale records)
pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/streptomyces_ground_truth.tsv \
    --output-dir data/interim/benchmark_set

# Or take every genome-scale record
pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/streptomyces_ground_truth.tsv \
    --output-dir data/interim/benchmark_set --top-n 0
```

Writes three files into `--output-dir`:

| file | use |
|---|---|
| `benchmark_genomes.tsv` | the set, with merge provenance and cluster ids |
| `analyzed_contigs.txt` | the `--contigs` scope file for `sharp.evaluate` |
| `benchmark_ground_truth.tsv` | ground truth for the selection, contigs normalized |

Pass `benchmark_ground_truth.tsv` and `analyzed_contigs.txt` — not the raw
MiBiG ground truth — to every tool in the comparison.

Record lengths are fetched once from NCBI esummary and cached in
`data/interim/record_lengths.tsv`; `--offline` reuses the cache and never
queries. Accessions NCBI cannot resolve (some WGS contigs) fall back to a
length lower bound derived from ground-truth coordinates, and are reported.

#### 2. Downloading the genome set

Fetches every genome named by `benchmark_genomes.tsv`, one FASTA per genome.
Resumable — an existing valid FASTA is skipped, so re-running after an
interruption or a rate-limit failure costs nothing.

```bash
# Uses data/interim/benchmark_set/benchmark_genomes.tsv by default
scripts/download_benchmark_genomes.sh

# Or point at another set / output directory
scripts/download_benchmark_genomes.sh path/to/set.tsv data/raw/genomes
```

Writes `data/raw/genomes/<ACCESSION>.fasta` (~420 Mb for the default 50).
Each download is verified to start with `>` and its header checked against the
expected accession — a mismatch there would make the genome score zero against
the ground truth for reasons invisible in the metrics, so it is reported loudly.

Genomes are fetched by **nucleotide accession** rather than read from the
cluster's local NCBI mirror. That mirror is indexed by *assembly* (GCA/GCF)
while the ground truth is keyed by *nucleotide record* (`AL645882.2`), and no
index file maps between the two — the mapping lives only inside each
per-assembly `*_assembly_report.txt`, so using it means building and maintaining
a reverse index over ~12.8k files. Fetching by accession makes the FASTA header
*be* the name the ground truth uses, and pins the benchmark to
`accession.version` rather than to whenever the mirror was last synced. See
[docs/NCBI_MIRROR.md](docs/NCBI_MIRROR.md) for the full inspection and the
conditions under which the mirror becomes worth it.

#### 3. Running the baselines as job arrays

Each baseline runs as a Slurm job array — one task per genome, so failures are
isolated and resubmitting retries only what failed. Both array scripts skip
genomes that already have output, so a partially-failed array can be resubmitted
wholesale.

```bash
CONTIGS=data/interim/benchmark_set/analyzed_contigs.txt
N=$(wc -l < $CONTIGS)

# Both scripts take the contigs file as $1 — that is how one shared output
# pool serves several scopes. Submit DeepBGC first; it dominates the runtime.

# DeepBGC: single-threaded, 2 cores / 8G / 2h per task, 8 concurrent.
sbatch --array=1-${N}%8 scripts/run_deepbgc_array.sbatch $CONTIGS

# antiSMASH: 4 cores / 4G / 1h per task, 8 concurrent.
sbatch --array=1-${N}%8 scripts/run_antismash_array.sbatch $CONTIGS
```

Output lands in `~/projects/<tool>/out_benchmark/<ACCESSION>/` — the shared
pool. Genomes already present are skipped, so widening the scope later only
pays for the genomes that are new. Raise `--time` in the antiSMASH script
before running scopes that include fungal genomes; the 1h sizing was measured
on *Streptomyces* chromosomes.

Both per-task sizings are **measured**, not assumed. antiSMASH comes from `seff`
on job 45315 (indices 1–2): ~3 min wall, 5.4%/7.8% CPU efficiency of 16 cores —
about 1.3 cores of real parallelism — and 1.6 GB peak, hence 4 cores / 4G / 1h.
DeepBGC carries the single-genome measurement (0.86 cores, 1.71 GB, 32m24s)
unchanged, since it is single-threaded and the per-task shape does not vary with
the array.

#### 4. Merging and evaluating

**The quick way — `scripts/run_benchmark.sh`.** Steps 1-3 stay manual (they are
slow, need the network, and the array is something you want to watch). Step 4 is
cheap and repeatable, so it is wrapped:

```bash
# Both tools; derives every path from the one scope name
scripts/run_benchmark.sh benchmark_set_strep

# One tool
scripts/run_benchmark.sh benchmark_set_strep antismash

# Sweep a score threshold — args after `--` go to sharp.evaluate.
# The parquet is reused, so each point costs seconds.
scripts/run_benchmark.sh benchmark_set_strep deepbgc -- --min-p-bgc 0.5

# Show the resolved commands without running them
scripts/run_benchmark.sh benchmark_set_strep --dry-run
```

It reads `data/interim/<scope>/{analyzed_contigs.txt,benchmark_ground_truth.tsv}`
and the pool at `$POOL_ROOT/<tool>/out_benchmark` (default `~/projects`), then
writes `data/interim/<tool>_predictions_<scope>.parquet` and
`data/processed/benchmark_<scope>_<tool>.json`, printing recall per tool at the
end.

Deriving both scope files from one name is the point: pairing a scope file with
the *wrong* ground truth yields a plausible-looking `benchmark.json` scored
against the wrong denominator, with nothing in the output to flag it. It also
refuses to overwrite an existing result (`--force`), warns when the pool holds
fewer genomes than the scope lists (the array is probably still running), and
rebuilds the parquet only when the pool has changed since it was written
(`--remerge` forces).

The two steps it wraps, if you want to run them by hand:

```bash
SCOPE=benchmark_set            # the scope name; one per experiment

for tool in antismash deepbgc; do
  pixi run python scripts/merge_predictions.py --tool ${tool} \
      --input-dir ~/projects/${tool}/out_benchmark \
      --contigs data/interim/${SCOPE}/analyzed_contigs.txt \
      --output data/interim/${tool}_predictions_${SCOPE}.parquet
done
```

Pass `--contigs` here too: a genome whose run never completed produces no
predictions but still counts in the recall denominator, so it is scored as
though the tool looked and found nothing. `merge_predictions.py` names those
genomes rather than letting them disappear into a lower recall number. Use
`--inspect` to see per-genome region counts before writing.

Finally, evaluate each tool against the **normalized** ground truth and the
shared scope file:

```bash
for tool in antismash deepbgc; do
    pixi run python -m sharp.evaluate \
        --predictions data/interim/${tool}_predictions_${SCOPE}.parquet \
        --ground-truth data/interim/${SCOPE}/benchmark_ground_truth.tsv \
        --contigs data/interim/${SCOPE}/analyzed_contigs.txt \
        --output data/processed/benchmark_${SCOPE}_${tool}.json
done
```

To produce another benchmark from the same pool, change `SCOPE` and re-run
these last two blocks — the arrays do not run again.

```bash
```

## Utilities

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

### Downloading a single genome

For a one-genome smoke test, rather than the 50-genome set above.
`download_genome.sh` fetches one contig by accession from NCBI nuccore and
derives the `--contigs` scope file alongside it. Defaults to `AL645882.2`
(*S. coelicolor* A3(2)), which carries 15 coordinate-resolved MiBiG clusters —
the most of any single contig, and enough for a recall number that actually
varies.

```bash
# Default: S. coelicolor A3(2)
scripts/download_genome.sh

# Or any other nuccore accession
scripts/download_genome.sh CP002993.1
```

Writes `data/raw/<ACCESSION>.fasta` and `data/interim/analyzed_contigs.txt`.
Pass that same contigs file to **every** tool in a comparison — see
[How evaluation works](#how-evaluation-works).

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
│   ├── NCBI_MIRROR.md                    # why the benchmark downloads instead of reading the cluster mirror
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
│   ├── download_benchmark_genomes.sh     # benchmark set TSV -> data/raw/genomes/<ACC>.fasta (resumable)
│   ├── download_genome.sh                # NCBI accession -> data/raw/<ACC>.fasta + --contigs scope file
│   ├── download_mibig.sh
│   ├── generate_mock_benchmark_data.py
│   ├── merge_predictions.py              # per-genome tool outputs -> one predictions.parquet
│   ├── run_benchmark.sh                  # score one scope against the shared pool (merge + evaluate)
│   ├── generate_mock_data.py
│   ├── parquet_to_tsv.py                 # generic parquet -> TSV dump (any pipeline parquet file)
│   ├── prepare_bgcatlas_ground_truth.py
│   ├── prepare_mibig_ground_truth.py
│   ├── select_benchmark_genomes.py       # ground truth -> benchmark genome set + scope + normalized GT
│   ├── _fetch_nuccore.sh          # sourced by the downloaders: efetch + "is this really FASTA" check
│   ├── _load_env.sh               # sourced by the setup scripts: loads .env, exports it
│   ├── setup_antismash.sh         # install baseline into its own isolated pixi env
│   ├── setup_deepbgc.sh
│   ├── setup_gecco.sh
│   ├── run_antismash.sbatch              # Slurm job: antiSMASH on the benchmark genome (CPU, n01)
│   ├── run_antismash_array.sbatch        # same, as a Slurm job array over the benchmark set
│   ├── run_deepbgc.sbatch                # Slurm job: DeepBGC on the benchmark genome (CPU, n01)
│   └── run_deepbgc_array.sbatch          # same, as a Slurm job array over the benchmark set
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
    ├── test_merge_predictions.py
    ├── test_metrics.py
    ├── test_model_management.py
    ├── test_parquet_to_tsv.py
    ├── test_prepare_bgcatlas.py
    ├── test_prepare_mibig.py
    └── test_select_benchmark_genomes.py
```

## Currently Working on

- Finishing benchmarks
