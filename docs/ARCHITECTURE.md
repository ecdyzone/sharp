# Architecture — S(H)ARP

## Module ownership

Every module has exactly one stated responsibility. If you're unsure where
something belongs, this table is the tie-breaker.

| Module | Owns | Does NOT own |
|---|---|---|
| `config.py` | `PROJECT_ROOT`, data dir constants, `DEFAULT_DEVICE`, `.env` loading, frozen config dataclasses | Pipeline logic, I/O beyond reading `.env`, behavioural defaults that vary at runtime |
| `io.py` | Data types (`ProteinRecord`, `PredictedRegion`, `KnownCluster`); all disk read/write functions | Metric math, model logic, domain knowledge |
| `model_management.py` | `MODEL_REGISTRY`, `select_device`, `ensure_model_available`, `residue_mean_pool`, `Embedder` | Batching strategy, logging, file I/O |
| `extract_embeddings.py` | Embedding step orchestration: load FASTA → batch → embed → write | Model loading, device selection (delegates to model_management), I/O (delegates to io) |
| `metrics.py` | Pure metric functions: `overlap_bp`, `reciprocal_overlap`, `evaluate_predictions`, `BenchmarkResult` | Any I/O, logging, or coordination |
| `evaluate.py` | Benchmark step orchestration: load → evaluate → write | Metric math (delegates to metrics), I/O (delegates to io) |

### Machine-specific configuration (`.env`)

`config.py` resolves its data-dir constants and `DEFAULT_DEVICE` from the
environment at import time, seeded by an optional `.env` at the project root
(gitignored; `.env.example` is the committed template). This exists so the same
checkout runs unmodified on a laptop and on a server where the data lives on a
different filesystem.

The boundary is **machine identity vs. pipeline behaviour**:

| Goes in `.env` | Stays in dataclass defaults / CLI |
|---|---|
| `SHARP_DATA_ROOT` and the four dir overrides | `min_cluster_frac`, `min_p_bgc`, `batch_size`, `max_length` |
| `SHARP_DEVICE` (hardware that's physically present) | `model_name` (a modelling choice) |
| Credentials, once a step needs one | Anything that changes what a run *means* |

Rationale: a benchmark number has to be reproducible from its command line. If a
threshold could come from an unversioned file, two machines could report
different numbers for the same command — so behavioural knobs never move into the
environment. Precedence is real env var → `.env` → in-repo default, which keeps
`SHARP_DEVICE=cpu pixi run ...` working as a one-off override.

Path constants are read once at import. Tests that need different values reload
the module under a patched environment (see `tests/test_config.py`).

### Where new things go

- New pipeline step → new module (`sharp/<step>.py`) + new config dataclass in `config.py`
- New data type → `io.py` (if constructed there), or `types.py` if 3+ modules need it
- New metric → `metrics.py` (pure function only)
- New external tool wrapper → either inside the step module (if it's specific) or a new `tools.py` (if reused)
- New reusable utility with no clear home → resist the urge; wait until it's needed twice

---

## Data type contracts

### `ProteinRecord`
```python
protein_id: str    # unique per genome
region_id: str     # from FASTA header: ">ID region_id=R001"
sequence: str      # uppercase amino acids, no gaps
```

### `PredictedRegion`
```python
region_id: str
contig: str        # NCBI accession or contig name from the genome FASTA
start: int         # 0-based half-open
end: int           # 0-based half-open
p_bgc: float       # model probability [0, 1]
predicted_class: str | None
```

### `KnownCluster`
```python
cluster_id: str    # e.g. "BGC0000001" or "BGC0000001.2" for multi-locus
contig: str        # NCBI accession (from MiBIG JSON loci[].accession)
start: int         # 0-based half-open (converted from MiBIG 1-based inclusive on ingest)
end: int           # 0-based half-open
cluster_class: str | None   # e.g. "PKS/NRPS"
```

**Coordinate invariant:** `start < end` always. Functions that produce coordinates
must assert or check this. `overlap_bp` and `reciprocal_overlap` handle degenerate
cases defensively (return 0) but no data type should store invalid coordinates.

---

## File formats

| File | Format | Schema |
|---|---|---|
| `proteins.faa`, `neighborhood_proteins.faa`, `neighborhood_dna.fna` | FASTA | Header: `>PROTEIN_ID region_id=R001 [optional fields]` |
| `anchors_sarp.tsv`, `anchors_heptarepeats.tsv`, `anchors.tsv` | TSV | `protein_id, contig, start, end, strand, score, type` |
| `neighborhoods.tsv` | TSV | `region_id, contig, start, end, anchor_ids, n_proteins` |
| `domains.tsv` | TSV | `protein_id, region_id, domain, e_value, start, end` |
| `embeddings.parquet` | Parquet (zstd) | `protein_id: str, region_id: str, embedding: list<float32>[D]` |
| `kg_features.parquet` | Parquet (zstd) | `region_id: str, n_similar_clusters: int, modal_class: str, ...` |
| `predictions.parquet` | Parquet (zstd) | `region_id, contig, start, end, p_bgc: float32, predicted_class: str` |
| `ground_truth.tsv` | TSV | `cluster_id, contig, start, end, class` |
| `benchmark.json` | JSON | `BenchmarkResult` dataclass (see `metrics.py`) |
| `model.pkl` | joblib | Serialized LightGBM + metadata |

---

## Metrics — methodological choices

### Scope: recall counts only contigs that were analyzed

Ground truth spans a database; a run spans one assembly. Recall is measured over
ground-truth clusters on contigs in **scope** — the contigs the tool was actually
run on — because a cluster on a contig the tool never saw cannot be found and
must not count as missed.

Pass `--contigs` (one name per line, or a `.fai`) and give **every tool in a
comparison the same file**. Without it, scope is inferred from the contigs
present in the predictions and a warning is logged: inferred scope is optimistic,
since a contig that was analyzed but produced no prediction silently leaves the
denominator, which flatters a tool that calls less.

*Why this is not a corner case:* on a real run of all three baselines against
`streptomyces_ground_truth.tsv` (430 clusters / 363 contigs), the analyzed contig
`AL589148.1` carried exactly **one** coordinate-resolved cluster. Dividing by all
430 capped recall at 0.002 for every tool.

### Match definition: asymmetric by default

`MatchCriterion` has two knobs, and they answer different questions:

| knob | question | default |
|---|---|---|
| `min_cluster_frac` | *Did the tool find this BGC?* — fraction of the **cluster** covered | `0.5` |
| `min_prediction_frac` | *Did it bound the BGC tightly?* — fraction of the **prediction** covered | `0.0` (off) |

Setting both to the same value reproduces the symmetric reciprocal rule, which is
still computed and reported alongside under `reciprocal_frac`
(`metrics.py:reciprocal_overlap` survives for labelling in `train.py`).

This **reverses** the earlier decision to defer asymmetric thresholds as "two
knobs instead of one; add later if needed". The trigger: DeepBGC called a 94 kb
region that covered **100%** of a 19 kb true cluster, and the symmetric rule
scored it recall `0.000` — reporting "not found" when only the bounds were loose.
Fusing detection and boundary accuracy into one pass/fail hides which of the two
actually failed, and penalizes tools that call wide regions on an axis that has
nothing to do with whether they found anything.

Boundary accuracy is not discarded — it moves to where it can be read directly:
`nucleotide.precision` and `boundary.median_prediction_coverage`. On the same
real run these order the tools as antiSMASH `0.382` > DeepBGC `0.171` >
GECCO `0.128`, i.e. by how much extra territory each calls.

Alternatives still rejected:
- **Any overlap**: too lenient — a 100 kb prediction touching a 10 kb cluster by 1 bp "matches"
- **Jaccard as the match rule**: equivalent ranking, less intuitive threshold (it *is* reported, as a bp-level summary)

### Recovery semantics: set-based, per unique cluster

A cluster is **recovered** if ≥1 prediction matches it (not all predictions).
Recall = |recovered clusters in scope| / |clusters in scope|.

Two overlapping predictions that together cover a cluster do NOT count as
recovery unless at least one of them individually passes the threshold. **Finding
a BGC means finding it as a unit** — this is unchanged. Such clusters are counted
in `boundary.n_clusters_recovered_by_union_only` so a tool that habitually splits
clusters is visible rather than silently penalized.

### Unmatched is not false — there is no region-level "precision"

MiBIG is deliberately incomplete: ~53% of *Streptomyces* entries are dropped for
missing coordinates (see CLAUDE.md), so **absence from the ground truth carries no
information**. A prediction with no match is *unvalidated*, not wrong.

So the output reports `matched_prediction_frac` and `unmatched_prediction_ids` —
never `precision` or `false_positive_prediction_ids` at region level. Treat
`matched_prediction_frac` as a lower bound on precision. Nucleotide-level
`precision` keeps its name because there it is a plain bp ratio, not a claim
about correctness.

*Concretely:* GECCO called 5 regions on `AL589148.1` where MiBIG knows of 1. The
old schema reported `precision=0.200`, asserting the other four were wrong.

### The criterion is recorded in every `benchmark.json`

A benchmark number without its thresholds is uninterpretable. Every run writes
`criterion` (both fracs), `reciprocal_frac`, and the full `scope` block —
including `n_clusters_total` beside `n_clusters_in_scope`, so how much of the
ground truth was excluded is always visible — plus `min_p_bgc`, the score cutoff
applied before scoring.

---

## External tool invocation pattern

Steps that shell out to bioinformatics tools (Bakta, hmmscan, FIMO) should:

```python
import subprocess
import logging

LOG = logging.getLogger(__name__)

def run_tool(cmd: list[str], step_name: str) -> None:
    LOG.info("running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{step_name} failed (exit {result.returncode}):\n{result.stderr}"
        )
    if result.stderr:
        LOG.debug("%s stderr:\n%s", step_name, result.stderr)
```

Always validate that the expected output file(s) exist after the call.
Never swallow non-zero exit codes silently.

---

## Testing conventions

- **Unit tests**: pure functions (metrics, accessors, parsers) — no I/O, deterministic
- **Integration tests**: use `tmp_path` (pytest fixture) for all disk I/O; never touch `data/`
- **Stub pattern for external models**: patch at the importing module's namespace (not the source module). The `StubEmbedder` in `conftest.py` is the reference implementation.
- **Property tests over edge cases**: test empty inputs, missing fields, inverted coordinates, degenerate intervals — not just the happy path
- **No mocking of `subprocess`**: shell-out wrappers are tested with real tool invocations in a separate `tests/integration/` dir (not run by default; marked `@pytest.mark.slow`)

### Marks

```python
# In pyproject.toml:
[tool.pytest.ini_options]
markers = [
    "slow: requires external tools (Bakta, hmmscan, FIMO) or GPU",
]
```

Run slow tests: `pytest -m slow`
Run fast tests only: `pytest -m "not slow"` (default in CI)

---

## One-time setup scripts (not pipeline steps)

These run once to build pre-requisites. They are in `scripts/` and do NOT
have corresponding `sharp/` modules.

| Script | Status | Output |
|---|---|---|
| `download_mibig.sh` | ✅ written | `data/raw/mibig_json_4.0/`, `mibig_gbk_4.0/`, `mibig_prot_seqs_4.0.fasta` |
| `prepare_mibig_ground_truth.py` | ✅ written | `data/raw/mibig_ground_truth.tsv` |
| `download_bgc-atlas.sh` | ✅ written | `data/raw/complete-bgcs/` (204k antiSMASH `.gbk`) |
| `prepare_bgcatlas_ground_truth.py` | ✅ written | `data/raw/bgcatlas_ground_truth.tsv` (secondary/noisy GT) |
| `build_sarp_hmm.py` | 🔲 not written | `data/raw/sarp_models.hmm` |
| `build_kg.py` | 🔲 not written | `data/raw/kg.gpickle` |

`build_sarp_hmm.py`: collect SARP sequences from UniProt/literature → align (MUSCLE or MAFFT) → `hmmbuild`.
`build_kg.py`: parse `mibig_json_4.0/` + BGC Atlas → build NetworkX graph → serialize.

---

## Competitor baselines (not pipeline steps)

antiSMASH, DeepBGC, and GECCO have mutually incompatible dependencies, so each
installs into its own isolated pixi env under `~/.local/src/<tool>/`. S(H)ARP
never invokes them: the user runs the tool, and S(H)ARP only **parses the output
files** into `predictions.parquet` — one converter script per tool, no subprocess.

| Script | Status | Role |
|---|---|---|
| `setup_antismash.sh` | ✅ written | install antiSMASH into its own pixi env |
| `setup_deepbgc.sh` | ✅ written | install DeepBGC into its own pixi env |
| `setup_gecco.sh` | ✅ written | install GECCO into its own pixi env |
| `convert_antismash_to_parquet.py` | ✅ written | antiSMASH JSON → `antismash_predictions.parquet` |
| `convert_deepbgc_to_parquet.py` | ✅ written | DeepBGC `.bgc.tsv` → `deepbgc_predictions.parquet` |
| `convert_gecco_to_parquet.py` | ✅ written | GECCO `.clusters.tsv` → `gecco_predictions.parquet` |

Each converter isolates the tool's column names + coordinate base in one place and
offers an `--inspect` mode (like `prepare_mibig_ground_truth.py`). Coordinate base
verified 2026-07-15 against a real run of all three tools on the same input FASTA
(evidence: span vs. matching `.gbk` LOCUS bp length, checked across every output
row, not just one): antiSMASH and DeepBGC are **both** already 0-based half-open,
no conversion; GECCO is 1-based inclusive, needs `start-1`. This confirms the old
hypothesis for GECCO but refutes it for DeepBGC (was assumed 1-based).
See `CLAUDE.md` → "Baseline integration" for the full spec and evidence.
