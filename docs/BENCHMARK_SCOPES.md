# Benchmark pools and scopes

Which genomes the baselines are run over, which slices of that we score, and why
each one exists.

Read [README → Run once, slice many](../README.md#run-once-slice-many--the-benchmarking-approach)
first for the mechanism. This document is the *catalogue*: the pool we committed
to, the scopes carved out of it, and the commands that build both.

---

## The one distinction that matters

A **pool** and a **scope** are different things, and only one of them is
expensive.

|  | pool | scope |
|---|---|---|
| what it is | genomes the tools were actually run on | a slice we score |
| where it lives | `$POOL_ROOT/<tool>/out_benchmark/<ACCESSION>/` | `data/interim/<scope>/` (28 KB) |
| cost to make | hours to days of cluster time | seconds |
| how many | **one** | as many as are useful |
| chosen | once, deliberately, as wide as affordable | freely, forever after |

Every scope below is a subset of the one pool. Adding a scope re-reads output
that already exists — no download, no Slurm, no waiting. That is why the pool is
built wider than any single experiment needs: the marginal genome is cheap now
and impossible to add cheaply later.

The corollary, from `TODO.md`: **never partition the pool per experiment.** The
pool is keyed by accession precisely so that two scopes sharing a genome share
one directory.

---

## The pool: all coordinate-resolved bacterial MiBiG

**1,280 clusters across 1,112 accessions** (`data/raw/bacterial_ground_truth.tsv`,
built with `--exclude-eukaryotes`).

Three decisions are baked into that.

**Bacteria, not all genera.** Eukaryotes add 353 clusters (327 fungal, 26
plant/animal) and two problems: fungal genomes blow through the antiSMASH array's
1-hour walltime, which was measured on *Streptomyces* chromosomes; and DeepBGC is
a bacterial model, so scoring it on *Aspergillus* measures domain-of-applicability
rather than detection quality. If "entire MiBiG" is ever asked for literally,
rebuild from `mibig_ground_truth.tsv` (1,634 clusters) and **raise
`run_antismash_array.sbatch --time` to 12h first**.

**Coordinate-resolved only, which is about half of MiBiG.** 3,013 entries → 1,634
with usable coordinates → 1,280 bacterial. The 1,363 dropped entries store
`location: {from: 0, to: 0}` — the compound is characterized but the locus is
unknown, so no amount of compute can score them. This is a property of the
database, not of our filtering.

**BGC-only deposits are in the pool, and out of every headline scope.** ~800 of
the 1,112 accessions are short records where the deposit *is* the cluster (median
~67 kb). They must not appear in a headline number — detection recall approaches
1.0 by construction for every tool, which would compress the antiSMASH-vs-DeepBGC
difference into noise. But they are ~1% of the compute, and having them in the
pool means the scenario-#4 sanity table and any future coverage-based refinement
cost no second campaign. So: **pool them, filter them at scope time** with
`--min-length` (the default, 1 Mb, drops them).

### What it costs

Only ~300 of the 1,112 accessions are genome-scale; the rest are ~67 kb deposits
that cost almost nothing:

| | genome-scale (~300) | BGC-only (~800) | total | wall at `%8` |
|---|---|---|---|---|
| antiSMASH | ~15 core-h | ~10 core-h | ~25 core-h | **~3 h** |
| DeepBGC | ~160 core-h | ~27 core-h | ~190 core-h | **~1 day** |

Extrapolated from measured `seff` numbers (antiSMASH ~3 min/genome, DeepBGC ~32
min/genome on ~8 Mb *Streptomyces* chromosomes). FASTA download is ~2-3 GB.
**Measure pool disk on the existing set before submitting** — it grows
monotonically and antiSMASH writes HTML, region `.gbk`s and JSON per genome.

---

## The scopes

Genus counts below are over the bacterial ground truth, genus taken as the first
word of the MiBiG taxonomy name. "Clusters" is before `--min-length`; the
selected count after it is roughly 35-40% of that, by analogy with
*Streptomyces* (414 → 156).

### Tier 1 — the ones to build

| scope | ground truth | clusters | accessions | answers |
|---|---|---|---|---|
| `benchmark_set_strep` | `streptomyces_ground_truth.tsv` | 414 → **156** | 352 → **93** | the headline number |
| `benchmark_set_actino` | `actinomycete_ground_truth.tsv` | **575** | 489 | where S(H)ARP actually applies |
| `benchmark_set_bact` | `bacterial_ground_truth.tsv` | **1,280** | 1,112 | "the entire (bacterial) MiBiG" |
| `benchmark_set_smoke` | derived from `_strep` | ~36 | 3 | is my setup working? |

**`benchmark_set_strep`** — already selected, at 156 clusters / 93 genomes. Note
the correction to an earlier estimate of ~352 genomes / ~414 clusters: that
counted every coordinate-resolved accession without `--min-length`. 263 of the
352 are BGC-only deposits. **156 is the real *Streptomyces* ceiling**, so the
published 50-genome run at 113 clusters is already 72% of it.

**`benchmark_set_actino`** — 34 actinomycete genera, 575 clusters / 489
accessions. *Streptomyces* is 407 of those, so this adds ~168 clusters beyond the
headline set. **This is the scope that matches the biology.** SARPs are an
actinomycete protein family, not a *Streptomyces* one, so this — not
`_strep` — is S(H)ARP's true domain of applicability, and it is the concrete
answer to "is the filtering only for actinobacteria?": yes, deliberately, because
outside the actinomycetes S(H)ARP has no anchor to work from.

The genera, by cluster count:

```
Streptomyces 407   Amycolatopsis 28   Micromonospora 22   Salinispora 16
Nocardia 12        Actinomadura 11    Saccharopolyspora 8 Kitasatospora 8
Actinoplanes 5     Mycobacterium 5    Planomonospora 4    Streptacidiphilus 4
Thermobifida 4     Nocardiopsis 4     Rhodococcus 4       Kutzneria 4
Frankia 3          Nonomuraea 2       Catenulispora 2     Corynebacterium 2
Microbispora 2     Verrucosispora 2   Actinokineospora 2  Planobispora 2
Pseudonocardia 2   Streptomonospora 2 Dactylosporangium 1 Actinosynnema 1
Clavibacter 1      Rathayibacter 1    Saccharothrix 1     Streptosporangium 1
Lentzea 1          Actinoallomurus 1
```

**`benchmark_set_bact`** — the whole pool minus BGC-only deposits. Expect
~300-450 genomes / ~470-600 clusters after `--min-length`; the range is wide
because the only measured shrink rate (27% of accessions, 37% of clusters
retained) comes from *Streptomyces*, where BGC-only deposits are unusually
common. *Pseudomonas*, *Bacillus* and *Burkholderia* entries are more often
complete genomes, so the true retention is likely higher.

**`benchmark_set_smoke`** — three genomes, not for science. A coworker should be
able to confirm their setup works in two minutes rather than discovering a typo
after merging 300 genomes. Taking the top three of `_strep` (they rank by cluster
count) gives ~36 clusters, so the smoke test still exercises real matching rather
than scoring an empty denominator.

### Tier 2 — negative controls

| scope | clusters | accessions |
|---|---|---|
| `benchmark_set_pseudomonas` | 61 | 49 |
| `benchmark_set_burkholderia` | 34 | 29 |
| `benchmark_set_bacillus` | 32 | 27 |

These are where S(H)ARP **should not** win: no SARPs, no anchor. Being able to
show antiSMASH ahead of S(H)ARP outside the actinomycetes, on purpose, is a
stronger claim than only reporting wins — it demonstrates the anchor is doing the
work rather than the classifier memorizing MiBiG. *Pseudomonas* at 61 clusters is
the only one large enough to stand alone; below ~25 clusters a recall estimate is
too noisy to interpret.

### Tier 3 — free slices, no new compute

**Per-class.** Same `analyzed_contigs.txt`, ground truth filtered on the `class`
column. Over the bacterial GT: NRPS 302, ribosomal 281, PKS 235, NRPS/PKS
hybrids 125, other 166, terpene 69, saccharide 20. This answers "which BGC class
is each tool best at?", which `CLAUDE.md` lists as a deliberate omission from
`evaluate.py` — but as a *scope* it needs no change to `evaluate.py` at all.

**BGC-only deposits** (scenario #4 in `TODO.md`). The ceiling table: recall
approaches 1.0 for every tool by construction. Report it separately, labelled,
and never merge it into a headline number.

---

## Building all of it

Everything below runs on the server, in the clone that owns the pool. Steps 1-4
are the pool owner's job and run once; step 5 onward is what everyone else uses.

### 1. Ground truths

```bash
# The pool's ground truth — all bacteria (already built)
pixi run python scripts/prepare_mibig_ground_truth.py \
    --input-dir data/raw/mibig_json_4.0 \
    --output data/raw/bacterial_ground_truth.tsv \
    --exclude-eukaryotes

# Single-genus scopes (already built for Streptomyces)
for g in Streptomyces Pseudomonas Burkholderia Bacillus; do
    pixi run python scripts/prepare_mibig_ground_truth.py \
        --input-dir data/raw/mibig_json_4.0 \
        --output "data/raw/$(echo "$g" | tr '[:upper:]' '[:lower:]')_ground_truth.tsv" \
        --genus "$g"
done
```

The actinomycete ground truth needs one run per genus, because `--genus` takes a
single substring (`prepare_mibig_ground_truth.py`). Concatenating the results is
a workaround until it accepts a list:

```bash
ACTINO="Streptomyces Amycolatopsis Micromonospora Salinispora Nocardia
Actinomadura Saccharopolyspora Kitasatospora Actinoplanes Mycobacterium
Planomonospora Streptacidiphilus Thermobifida Nocardiopsis Rhodococcus Kutzneria
Frankia Nonomuraea Catenulispora Corynebacterium Microbispora Verrucosispora
Actinokineospora Planobispora Pseudonocardia Streptomonospora Dactylosporangium
Actinosynnema Clavibacter Rathayibacter Saccharothrix Streptosporangium Lentzea
Actinoallomurus"

TMP="$(mktemp -d)"
for g in $ACTINO; do
    pixi run python scripts/prepare_mibig_ground_truth.py \
        --input-dir data/raw/mibig_json_4.0 \
        --output "${TMP}/${g}.tsv" --genus "$g"
done
head -1 "${TMP}/Streptomyces.tsv"           > data/raw/actinomycete_ground_truth.tsv
tail -q -n +2 "${TMP}"/*.tsv | sort -u     >> data/raw/actinomycete_ground_truth.tsv
rm -rf "${TMP}"

wc -l data/raw/actinomycete_ground_truth.tsv    # expect ~576 (575 clusters + header)
```

### 2. The pool list

`--min-length 0` keeps the BGC-only deposits, `--top-n 0` keeps everything. This
is the only selection run that needs the network — ~780 uncached NCBI `esummary`
lookups, so use tmux. It is resumable, and the lengths land in
`data/interim/record_lengths.tsv`, which every later scope reuses.

```bash
pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/bacterial_ground_truth.tsv \
    --output-dir data/interim/pool_bact \
    --top-n 0 --min-length 0
```

Deliberately named `pool_bact`, not `benchmark_set_*`: `run_benchmark.sh` lists
`benchmark_set*/` as available scopes, and the pool is not a scope — scoring it
directly would mix BGC-only deposits into the number.

### 3. Download the genomes — once, for every scope

```bash
scripts/download_benchmark_genomes.sh data/interim/pool_bact/benchmark_genomes.tsv
```

**This one download serves every scope.** It writes `data/raw/genomes/<ACC>.fasta`,
keyed by accession, exactly like the output pool — so a scope needs no download of
its own. ~2-3 GB, resumable (a file is trusted only if it is non-empty and starts
with `>`, so a truncated fetch is retried), and it verifies that each FASTA header
equals the expected accession. Use tmux.

### 4. Run the baselines over the pool

The pool is larger than Slurm will accept as a single array. `MaxArraySize`
(commonly 1001) caps the highest legal array *index* at 1000 — it is not a
concurrency limit, so letting a running array drain does not make index 1001
legal; `--array=1001-1087` is rejected at parse time with `Invalid job array
specification`. Both scripts therefore take an **index offset as `$2`** and read
line `SLURM_ARRAY_TASK_ID + OFFSET`, so successive windows reuse indices
`1..1000` over the same list rather than cutting it into per-chunk files that
would each have to be regenerated and frozen separately.

```bash
mkdir -p logs        # Slurm writes task logs here and will NOT create it
CONTIGS=data/interim/pool_bact/analyzed_contigs.txt
N=$(wc -l < "$CONTIGS")
MAX=1000             # MaxArraySize - 1; scontrol show config | grep -i maxarraysize

# submit() <script> <throttle> — one chained submission per window of $MAX
submit() {
    local prev="" off n dep
    for (( off = 0; off < N; off += MAX )); do
        n=$(( N - off < MAX ? N - off : MAX ))
        dep=${prev:+--dependency=afterany:$prev}
        prev=$(sbatch --parsable $dep --array=1-${n}%$2 "$1" "$CONTIGS" "$off")
        echo "$(basename "$1") offset $off -> job $prev"
    done
}

submit scripts/run_deepbgc_array.sbatch   8    # submit first, ~1 day
submit scripts/run_antismash_array.sbatch 8    # ~3 h
```

For the 1,087-genome pool that is two submissions per tool: `--array=1-1000` at
offset 0, then `--array=1-87` at offset 1000, chained with `--dependency` so the
second starts as the first drains. The offset defaults to 0, so a pool under the
cap still submits with a plain `--array=1-${N}%8 ... "$CONTIGS"`.

**Do not edit `analyzed_contigs.txt` while an array is in flight.** Each task
reads it with `sed -n "${LINE}p"` when that task *starts* — only the script body
is snapshotted at submit time. Renumbering the lines mid-run silently skips some
genomes and double-runs others. This is why dropping unusable accessions from
the scope has to wait until the pool has finished.

Both arrays are resumable: a task whose result file already exists
(`<ACC>/<ACC>.json` for antiSMASH, `<ACC>/<ACC>.bgc.tsv` for DeepBGC) exits
immediately. Resubmit the whole array to fill gaps; already-done genomes cost
seconds. This is also what makes widening a scope cheap later — only new genomes
run.

Watch for one failure mode: a task killed on walltime leaves a partial output
directory with no result file, so it is retried — but antiSMASH may refuse to
write into a non-empty output directory. Check that on the first failed task
rather than across 800.

### 5. Carve the scopes

All offline — the lengths are cached from step 2.

```bash
# Tier 1
pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/streptomyces_ground_truth.tsv \
    --output-dir data/interim/benchmark_set_strep  --top-n 0 --offline

pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/actinomycete_ground_truth.tsv \
    --output-dir data/interim/benchmark_set_actino --top-n 0 --offline

pixi run python scripts/select_benchmark_genomes.py \
    --ground-truth data/raw/bacterial_ground_truth.tsv \
    --output-dir data/interim/benchmark_set_bact   --top-n 0 --offline

# Tier 2 — negative controls
for g in pseudomonas burkholderia bacillus; do
    pixi run python scripts/select_benchmark_genomes.py \
        --ground-truth "data/raw/${g}_ground_truth.tsv" \
        --output-dir "data/interim/benchmark_set_${g}" --top-n 0 --offline
done
```

The smoke scope is three genomes lifted out of `_strep`, with its ground truth
restricted to match:

```bash
SRC=data/interim/benchmark_set_strep
DST=data/interim/benchmark_set_smoke
mkdir -p "$DST"
head -3 "$SRC/analyzed_contigs.txt" > "$DST/analyzed_contigs.txt"
awk -F'\t' '{sub(/\r$/, "")} NR==FNR {c[$1]; next} FNR==1 || ($2 in c)' \
    "$DST/analyzed_contigs.txt" "$SRC/benchmark_ground_truth.tsv" \
    > "$DST/benchmark_ground_truth.tsv"
```

Per-class scopes are the same trick on the other axis — same contigs, ground
truth filtered by class:

```bash
SRC=data/interim/benchmark_set_actino
for cls in NRPS PKS ribosomal terpene; do
    DST="data/interim/benchmark_set_actino_${cls}"
    mkdir -p "$DST"
    cp "$SRC/analyzed_contigs.txt" "$DST/"
    awk -F'\t' -v c="$cls" '{sub(/\r$/, "")} FNR == 1 || $5 == c' \
        "$SRC/benchmark_ground_truth.tsv" > "$DST/benchmark_ground_truth.tsv"
done
```

`$5 == c` keeps pure classes only; hybrids (`NRPS/PKS`) fall out of every
single-class scope. Use `$5 ~ c` instead to count a hybrid toward each of its
components — state which you used, since it changes the denominator.

> **`{sub(/\r$/, "")}` is not decoration.** Every `*_ground_truth.tsv` is written
> with **CRLF** line endings (`io.py` uses `csv.writer`, whose default
> `lineterminator` is `\r\n`). Python readers do not care, so `sharp.evaluate` is
> unaffected — but in `awk`, `cut` or `grep` the *last* column carries a trailing
> `\r`, so `$5 == "NRPS"` matches **nothing**. Without the strip these scopes come
> out with a header and zero clusters, and a benchmark scored against an empty
> denominator looks like a finished run, not an error. Verify any hand-built scope
> with `wc -l` before scoring it.

**Check every scope against the pool** before trusting a number. Anything printed
here is a genome the scope expects but the pool never ran, which scores as "the
tool found nothing":

```bash
for s in data/interim/benchmark_set_*/; do
    missing=$(comm -23 <(sort "$s/analyzed_contigs.txt") \
                       <(sort data/interim/pool_bact/analyzed_contigs.txt) | wc -l)
    echo "$(basename "$s"): ${missing} genome(s) not in the pool"
done
```

### 6. Score them

```bash
for s in benchmark_set_strep benchmark_set_actino benchmark_set_bact; do
    scripts/run_benchmark.sh "$s" --label "$USER"
done
```

See [README → Benchmarking in a shared clone](../README.md#benchmarking-in-a-shared-clone).

---

## Caveats that belong in any write-up

- **Recall is over coordinate-resolved MiBiG, not all known BGCs.** The
  denominator is ~half of MiBiG by design, and the dropped half is not random —
  it skews toward older, compound-first submissions. This hits every tool
  equally, so the *comparison* is unbiased; the absolute numbers are not
  "recall over all known BGCs".
- **`matched_prediction_frac` is a lower bound on precision, not precision.** The
  ground truth is incomplete, so an unmatched prediction is unvalidated rather
  than wrong.
- **Report `detection` and `reciprocal` recall together.** They ranked antiSMASH
  and DeepBGC in opposite directions on the *S. coelicolor* run; neither is a
  fair single headline.
- **A scope is a pair of files.** `analyzed_contigs.txt` and
  `benchmark_ground_truth.tsv` from one `--output-dir`. The ground truth is
  contig-normalized for that selection (RefSeq/GenBank twins collapsed onto a
  primary accession), so pairing a scope file with the raw MiBiG ground truth
  silently drops clusters filed under a twin. `run_benchmark.sh` derives both
  from one scope name so this cannot happen by hand.
- **BGC Atlas is a secondary ground truth.** Its labels are themselves
  predictions, so agreement with it does not prove correctness. Never report it
  alone.
