# TODO

## Benchmarking approach — run once, slice many

**Decided 2026-08-21. This supersedes the per-experiment run model.**

The baselines run over the broadest genome set we are willing to pay for,
*once*. Every benchmark after that is a re-scope, not a re-run: the array output
pool is keyed by accession (`~/projects/<tool>/out_benchmark/<ACCESSION>/`), and
`--contigs` filters both the predictions and the ground-truth denominator
(`metrics.py:288-289`). See README → "Run once, slice many".

Rules:
- **Never partition the output pool per experiment** — shared accession keys are
  what make re-slicing free and let a widened scope pay only for new genomes.
- **`--contigs` is mandatory on `merge_predictions.py`**, since the pool holds
  every genome ever run.
- A scope is a **pair**: `analyzed_contigs.txt` + its `benchmark_ground_truth.tsv`
  from the same `select_benchmark_genomes.py --output-dir`. Never pair a scope
  file with the raw MiBiG GT — twins-merged clusters vanish silently.
- Name derived artifacts per scope: `benchmark_set_<name>/`,
  `<tool>_predictions_<name>.parquet`, `benchmark_<name>_<tool>.json`.

### Scenarios queued against the shared pool

Numbers measured 2026-08-21 against `data/raw/mibig_json_4.0`. "Genomes" is
after `--min-length` drops BGC-only deposits; `--min-length 0` keeps them.

| # | Scope | Clusters | Genomes | Ground truth | Status |
|---|---|---|---|---|---|
| 0 | 50-genome *Streptomyces* (current) | 113 | 50 | `benchmark_set/` | ✅ done, see `sharp-davinci-copy/data/processed-50genomes/` |
| 1 | **Full *Streptomyces*** | **156** | **93** | `streptomyces_ground_truth.tsv` | selected → `data/interim/benchmark_set_strep/`, not run |
| 2 | **Bacteria only** | ~1,280 GT | TBD | `bacterial_ground_truth.tsv` ✅ built | not selected |
| 3 | All genera (entire MiBiG) | ~1,634 GT | TBD | `mibig_ground_truth.tsv` | not selected |
| 4 | BGC-only deposits (`--min-length 0` minus #1) | ~245 (*Strep*) | ~241 | any of the above | not selected |
| 5 | Per-genus slices (*Amycolatopsis*, *Micromonospora*, …) | — | — | `--genus <name>` | idea only |

- [ ] **#1 Full *Streptomyces* — run next.** Selected: 93 genomes / 156
      clusters. **Note the correction:** an earlier estimate of ~352 genomes /
      ~414 clusters was wrong — it counted all coordinate-resolved accessions
      without `--min-length`. 263 of 352 *Streptomyces* accessions are BGC-only
      deposits (median 67 kb) carrying 245 clusters, dropped by design. 156 is
      the real *Streptomyces* ceiling, and the current run at 113 is already 72%
      of it. Clean superset of #0: all 50 genomes retained, 43 added.
      Download ~0.66 GB. Full commands in `../TODO-next-experiments.md`.
- [ ] **#2 Bacteria only.** `bacterial_ground_truth.tsv` is built (1,280
      clusters / 1,112 accessions via `--exclude-eukaryotes`). Selection needs
      ~1,100 uncached NCBI esummary lookups — run under tmux, it is resumable.
      Recommended over #3: removes the fungal walltime trap, and DeepBGC is a
      bacterial model so scoring it on *Aspergillus* measures
      domain-of-applicability rather than detection quality.
- [ ] **#3 All genera.** Only if the team leader specifically means this by
      "entire MiBiG". Adds 353 eukaryotic clusters (327 fungal + 26
      plant/animal) over #2. **Raise antiSMASH `--time` to 12h first** — the 1h
      sizing was measured on *Streptomyces* chromosomes and fungal genomes are
      much slower.
- [ ] **#4 BGC-only deposits — report separately, never merged into a headline
      table.** The record *is* the cluster (median cluster covers 68% of its
      record; 97 of 241 are ≥90%), so detection recall approaches 1.0 by
      construction for every tool. Merging them would push recall to ~0.95+ for
      both tools and dilute the antiSMASH-vs-DeepBGC difference into noise, and
      would invert precision (little non-cluster territory to be wrong about).
      Worth running *because they are cheap* (~67 kb vs ~8 Mb, ~1% the compute)
      and having them in the pool means the sanity table costs no second
      campaign. Open refinement: `--min-length` is blunt; a coverage-based
      filter (drop records where the cluster is >50% of the sequence) would
      recover ~90 clusters of real signal from the 95 records currently under
      50% coverage. Needs a stated threshold in the methods.
- [ ] **Disk check before any large scope.** The pool grows monotonically;
      antiSMASH writes HTML + region `.gbk`s + JSON per genome. Measure with
      `du -sh ~/projects/antismash/out_benchmark` on the existing 50 and
      extrapolate before submitting ~1,100.

## Benchmarks — real data

- [x] Run the full comparison with real data — done for the SCP1 smoke test
      (`AL589148.1`, all three tools). This exposed three defects in the
      benchmark core, now fixed (see `docs/BACKLOG.md` Tier 0). Still to do at
      scale: a genome with more than one coordinate-resolved MiBIG cluster on it.
- [x] Re-run the comparison on a genome with real recall signal — done on
      *S. coelicolor* A3(2) (`AL645882.2`, 15 coordinate-resolved MiBIG
      clusters), run on the davinci server via `scripts/run_antismash.sbatch`
      and `scripts/run_deepbgc.sbatch`. Both used the same explicit `--contigs`
      scope (`source: "explicit"`, `n_contigs: 1`, 15 clusters in scope out of
      1675 in the GT file). Results below.
- [~] **GECCO — paused.** Not part of the `AL645882.2` comparison and no
      `scripts/run_gecco.sbatch` exists. If it comes back: mirror
      `run_deepbgc.sbatch`, size it from `seff` after the first run, and convert
      with `scripts/convert_gecco_to_parquet.py` (GECCO is the one tool needing
      `start - 1`).
- [ ] S(H)ARP itself can't be benchmarked yet — `predict.py` and the rest of the
      pipeline (`annotate.py` → `train.py`) aren't implemented yet (see CLAUDE.md
      "What is NOT YET IMPLEMENTED").
- [x] **Scale the benchmark past one genome** — done 2026-08-18. The
      `AL645882.2` run had a denominator of 15 clusters; the selected set is 50
      genomes / 113 clusters. Tooling: `select_benchmark_genomes.py` (filters
      BGC-only deposits, merges RefSeq/GenBank twins, emits a contig-normalized
      GT), `download_benchmark_genomes.sh`, `run_{antismash,deepbgc}_array.sbatch`,
      `merge_predictions.py`. **Not yet run** — needs the davinci server.
- [ ] **Run it** on davinci. Full sequence:
      ```bash
      scripts/download_benchmark_genomes.sh          # ~420 Mb, resumable
      N=$(wc -l < data/interim/benchmark_set/analyzed_contigs.txt)   # 50

      # Both sizings are now measured, so submit the full arrays directly.
      # (antiSMASH indices 1-2 already ran under job 45315; the script is
      # resumable, so re-running them just skips.)
      sbatch --array=1-${N}%8 scripts/run_antismash_array.sbatch
      sbatch --array=1-${N}%8 scripts/run_deepbgc_array.sbatch

      pixi run python scripts/merge_predictions.py --tool antismash \
          --input-dir ~/projects/antismash/out_benchmark \
          --contigs data/interim/benchmark_set/analyzed_contigs.txt \
          --output data/interim/antismash_predictions.parquet

      pixi run python -m sharp.evaluate \
          --predictions data/interim/antismash_predictions.parquet \
          --ground-truth data/interim/benchmark_set/benchmark_ground_truth.tsv \
          --contigs data/interim/benchmark_set/analyzed_contigs.txt \
          --output data/processed/benchmark_antismash.json
      ```
      **Use `benchmark_ground_truth.tsv`, not `streptomyces_ground_truth.tsv`** —
      the benchmark set merges RefSeq/GenBank twins onto one primary accession,
      and only the normalized GT has contigs renamed to match. Passing the raw GT
      silently drops every cluster filed under a non-primary twin (this is how
      *S. coelicolor*'s 16th cluster goes missing). Pass the same `--contigs` to
      every tool or the recall denominators differ.
- [x] **antiSMASH array sizing measured** (2026-08-19, `seff` on job 45315,
      indices 1-2): ~3 min wall per genome, CPU efficiency 5.4%/7.8% of 16
      cores, 1.6 GB peak. antiSMASH's own module scheduler parallelises far
      less than `--cpus` suggests — the same lesson DeepBGC taught (sized 8
      cores, measured 0.86). `run_antismash_array.sbatch` is now 4 cores / 4G /
      1h, throttle `%8`.
- [ ] Write up the numbers below with the MiBIG coordinate-coverage,
      benchmark-scope, and BGC Atlas optimism caveats (CLAUDE.md → "Benchmark
      comparison"). Report `detection` and `reciprocal` recall together.
- [x] Decide whether `min_prediction_frac` should be non-zero for the headline
      table — **the two rankings do disagree**, so keep 0.0 and always report
      both. See "Open question" below for what still needs deciding.

> **Superseded numbers.** The ground truth was corrected on 2026-08-18
> (unusable accessions and duplicate loci dropped at ingest; RefSeq/GenBank
> twins merged), which raises *S. coelicolor* from 15 to **16** clusters — the
> 16th was filed under `NC_003888.3`, the RefSeq copy of the same sequence.
> The table below is from before that fix and will be regenerated by the
> scaled run.

### Results — `AL645882.2` (S. coelicolor A3(2)), 15 clusters in scope

Ground truth: `data/raw/mibig_ground_truth.tsv` (all genera, 1675 clusters;
scoping to the one contig leaves the same 15 *Streptomyces* clusters the
genus-filtered GT would). Full write-up with caveats:
`../sharp-davinci-copy/data/processed/AL645882.2.md`; raw JSON alongside it.

| | antiSMASH | DeepBGC |
|---|---|---|
| predictions in scope | 29 | 167 |
| **detection** recall | **1.000** (15/15) | 0.733 (11/15) |
| **reciprocal** recall | 0.267 (4/15) | **0.467** (7/15) |
| matched prediction frac | 0.483 (14/29) | 0.060 (10/167) |
| nucleotide recall | 0.986 | 0.897 |
| nucleotide precision | 0.251 | 0.142 |
| predicted bp | 1,128,277 | 1,806,162 |
| median prediction coverage | 0.317 | 0.702 |
| clusters recovered by union only | 0 | 0 |

Reading it:

- **The two recall criteria rank the tools in opposite directions.** antiSMASH
  finds every cluster (detection 1.000) but its regions are wide — median
  prediction coverage 0.317 means the typical matched region is ~2/3 territory
  that isn't the cluster, so only 4/15 survive the symmetric 50% rule. DeepBGC
  misses 4 clusters outright but its matched calls are tighter (0.702), so it
  wins on reciprocal. Neither number alone is a fair headline; this is exactly
  the boundary-tightness story flagged after the SCP1 run, now with n=15 instead
  of n=1.
- **DeepBGC's precision problem is volume, not just width.** 167 calls covering
  1.8 Mb of an 8.7 Mb chromosome (~21%), 157 of them unmatched, many only a few
  hundred bp. `matched_prediction_frac` 0.060 is a lower bound on precision, not
  precision — but the gap to antiSMASH's 0.483 is large enough that the ranking
  is unlikely to be an artefact of incomplete ground truth.
- **Both tools recover clusters with single predictions.** `n_clusters_recovered_by_union_only`
  is 0 for both, and each has exactly one merged prediction (a single call
  spanning ≥2 clusters), so split/merge pathology is not driving the numbers.
- DeepBGC's 4 misses: `BGC0000551`, `BGC0000660`, `BGC0000940`, `BGC0001181`.
  antiSMASH missed none. Worth a look at what those four have in common before
  the write-up.

### Open question

- Headline table: keep `min_prediction_frac = 0.0` (detection) as the primary
  recall and report `reciprocal` beside it, or promote reciprocal? Now that the
  rankings are known to disagree, this is a presentation decision for the team,
  not a metrics bug. Recommendation: report both columns side by side and never
  a single "recall" number.
