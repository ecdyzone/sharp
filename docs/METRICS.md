# Reading a `benchmark.json`

What each field of a benchmark result means. Written for reading a number, not
producing one — to produce one see [the README's benchmarking
section](../README.md#benchmarking).

Examples throughout come from one real result,
`benchmark_pool_bact_antismash.mzamaral.json` — antiSMASH over the whole
bacterial pool, 1,280 clusters on 1,087 contigs.

## Comparing two tools — the short version

Compare these six fields. Everything else is context for them.

| Read | Field | antiSMASH on `pool_bact` |
|---|---|---|
| Did it find the BGCs? | `detection.recall` | `0.810` |
| Were the calls tight too? | `reciprocal.recall` | `0.461` |
| How many calls landed on a known cluster? | `detection.matched_prediction_frac` | `0.127` |
| How much called territory is supported? | `nucleotide.precision` | `0.103` |
| How much of a typical call is real cluster? | `boundary.median_prediction_coverage` | `0.552` |
| Does it split or fuse clusters? | `boundary.n_clusters_recovered_by_union_only` / `n_merged_predictions` | `2` / `26` |

Read as one story: antiSMASH found 81% of the clusters, but a typical call is
only 55% cluster and 45% flanking sequence, which is why the strict rule drops
recall to 46%. It called 338 Mb where 84 Mb is known.

**Before comparing anything, check three things.**

1. `scope.source` is `"explicit"` in both files. `"inferred"` is optimistic and
   not comparable.
2. `scope.n_clusters_in_scope` matches. Different denominators, different
   benchmark.
3. `criterion` and `min_p_bgc` match. The thresholds change the numbers.

`provenance.scope` names the genome set both files should share. Look it up in
[BENCHMARK_SCOPES.md](BENCHMARK_SCOPES.md).

## Three traps

- **`matched_prediction_frac` is not precision.** Ground truth is incomplete on
  purpose, so an unmatched prediction is unvalidated, not wrong. It is a *lower
  bound* on precision, and it is why no region-level `precision` field exists.
  For a precision claim use `nucleotide.precision`, which is a plain bp ratio.
- **Recall counts only in-scope clusters.** A cluster on a contig the tool never
  saw is excluded, not missed. That is what `--contigs` is for.
- **antiSMASH has no score**, so all its predictions carry `p_bgc = 1.0`. A
  `--min-p-bgc` sweep moves DeepBGC and GECCO only.

---

## Field reference

Coordinates are 0-based half-open, so `end - start` is a length in bp. Overlap is
counted only between items on the same contig.

### `scope` — what was measured

| Field | Meaning |
|---|---|
| `source` | `explicit` if `--contigs` was passed, `inferred` if reconstructed from the predictions |
| `n_contigs` | contigs the tool was run on; everything else is excluded |
| `n_clusters_in_scope` | the recall denominator |
| `n_clusters_total` | clusters in the ground-truth file before scoping |
| `n_predictions_in_scope` / `n_predictions_total` | same split on the prediction side |
| `n_predictions_below_threshold` | predictions dropped for scoring too low |
| `min_p_bgc` | the score cutoff that dropped them; `0.0` keeps everything |

In the example, all 1,280 clusters are in scope because the ground truth was
built for this pool, while 389 of 8,332 predictions sit on contigs outside it.
A prediction pool wider than the scope is normal: baselines run once over
everything and each experiment re-slices.

### `criterion` and `reciprocal_frac` — the match rule

| Field | Question | Default |
|---|---|---|
| `min_cluster_frac` | did the tool cover enough of the **cluster**? | `0.5` |
| `min_prediction_frac` | is the **prediction** mostly cluster? | `0.0`, off |
| `reciprocal_frac` | threshold for the strict symmetric rule, reported separately | `0.5` |

By default a prediction matches when it covers half the cluster, and its own size
is not held against it. Boundary quality is measured separately rather than
folded into the pass/fail. See [ARCHITECTURE.md → Metrics](ARCHITECTURE.md#metrics--methodological-choices)
for why.

### `detection` — the headline block

| Field | Meaning |
|---|---|
| `recall` | `n_recovered / n_clusters` |
| `n_recovered` | clusters found by at least one prediction, counted once each |
| `n_clusters` | same as `scope.n_clusters_in_scope` |
| `matched_prediction_frac` | `n_matched_predictions / n_predictions`; a lower bound on precision |
| `f1` | harmonic mean of the two fractions, so also a lower bound |

A cluster counts as recovered only if a **single** prediction clears the
threshold. Clusters covered only by several predictions together land in
`boundary.n_clusters_recovered_by_union_only` instead.

### `reciprocal` — same counts, strict rule

Identical fields, recomputed symmetrically at `reciprocal_frac`. Its value is the
comparison with `detection`: equal means the calls were tight, much lower means
the tool found the clusters but called far wider regions. The example halves,
from `0.810` to `0.461`.

### `nucleotide` — base pairs, not regions

Compares the merged union of all predictions against the merged union of all
clusters, so it is unaffected by region boundaries and by clusters split across
calls.

| Field | Meaning |
|---|---|
| `gt_bp` / `predicted_bp` / `intersect_bp` | bp of ground truth, of predictions, and of both |
| `recall` | `intersect_bp / gt_bp` |
| `precision` | `intersect_bp / predicted_bp` — the honest precision number |
| `jaccard` | `intersect / (gt + predicted - intersect)`, a symmetric summary |

### `boundary` — tightness, splits, merges

| Field | Meaning |
|---|---|
| `median_prediction_coverage` | median fraction of a matched prediction that is cluster; `1.0` is a perfect call |
| `median_cluster_coverage` | median fraction of a recovered cluster covered by its best prediction |
| `n_clusters_recovered_by_union_only` | clusters no single call covered enough, but several together would; **not** in `detection.n_recovered` |
| `n_merged_predictions` | predictions that matched two or more clusters at once |

The last two are opposite failure modes, splitting and fusing, and neither is
visible in recall.

### The id lists

`recovered_cluster_ids` and `missed_cluster_ids` partition the in-scope ground
truth as MiBIG accessions. `matched_prediction_ids` and
`unmatched_prediction_ids` partition the in-scope predictions, using each tool's
own region ids, which trace back to its output files:

| Tool | Id shape | Source |
|---|---|---|
| antiSMASH | `AL645882.2.region004` | the matching `.region004.gbk` |
| DeepBGC | `AL589148.1_31460-41750.1` | a row of `out.bgc.tsv` |
| GECCO | `AL589148.1_cluster_1` | the matching `_cluster_1.gbk` |

`ids_truncated: true` means at least one list was cut at 1,000 entries
(`--max-listed-ids`, `0` for unlimited). The counts stay complete; only the lists
are short. It is `true` in the example.

### `provenance` — who produced the file

Added by `scripts/run_benchmark.sh`, absent if `sharp.evaluate` was called
directly. Records `user`, `host`, `written_at`, `scope`, `label`, `tool`, `pool`,
`pool_n_genomes`, `predictions`, `predictions_reused`, `git_commit` and
`evaluate_args`. Several people benchmark inside one server clone, so this is how
a result stays readable months later. `pool_n_genomes` says how many genomes had
finished when the merge ran, which is the field to check if a number looks stale.

---

## Not in the schema, on purpose

Region-level `precision` and any `false_positive` list, for the reason above.
AUROC requires scores for all candidates, not only reported positive regions.
Per-class breakdown is not built into the evaluation JSON; report tables provide
a separate per-class recovery view.

One caveat applies to every absolute number: denominators are *coordinate-resolved*
MiBIG. About half of MiBIG's entries have no coordinates and are dropped, skewed
toward older, well-characterized PKS/NRPS submissions. It hits every tool equally,
so comparisons stay fair, but quote absolute recall with that qualifier.

Rationale for the design: [ARCHITECTURE.md → Metrics](ARCHITECTURE.md#metrics--methodological-choices).
Implementation: `src/sharp/metrics.py`.
