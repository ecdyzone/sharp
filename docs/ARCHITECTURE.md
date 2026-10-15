# Benchmark architecture

The benchmark package keeps four responsibilities separate:

| Module | Responsibility |
|---|---|
| `config.py` | Data paths, optional `.env` loading and evaluation configuration |
| `io.py` | Data types and file readers/writers |
| `metrics.py` | Pure interval matching and metric calculations |
| `evaluate.py` | Load inputs, apply scope, evaluate and write JSON |

Scripts normalize external tool outputs into predictions, prepare ground truth,
and build report tables and figures. Baseline tools run in separate environments.
The benchmark package is installed independently of the team's SHARP package.
Coordinates are always 0-based half-open `[start, end)` after ingestion.

## Metrics — methodological choices

> This section is the *rationale*. For a field-by-field walkthrough of an
> actual `benchmark.json`, aimed at readers who have not read the code, see
> [METRICS.md](METRICS.md).

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
