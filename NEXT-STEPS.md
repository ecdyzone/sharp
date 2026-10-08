# NEXT STEPS

Working notes, **2026-10-08**. Scratch file — not a doc. Anything here that
outlives the week belongs in `TODO.md`, `docs/BENCHMARK_SCOPES.md` or `CLAUDE.md`
instead; this is only "where we stopped".

**State: the bacterial pool is computed and scored.** antiSMASH and DeepBGC have
run over all of `pool_bact`, the merge-time version fix has landed, and
`scripts/run_benchmark.sh pool_bact --label $USER --remerge --force` produces
numbers. Two checks remain before those numbers are reportable (1 and 2 below).

---

## 1. Confirm the 13 unusable accessions are out of the scope

8 accessions never downloaded (3 protein accessions, 5 WGS master records) and 5
`GPC_`/`GPS_` identifiers resolve to unrelated sequences. If they are still in
`analyzed_contigs.txt` and `benchmark_ground_truth.tsv`, they sit in the recall
denominator as clusters no tool was given a chance to find.

```bash
wc -l data/interim/pool_bact/analyzed_contigs.txt     # 1074 = done, 1087 = not
```

If 1087 — safe now, no array is in flight:

```bash
cd data/interim/pool_bact
cat > /tmp/drop.txt <<'EOF'
EGF94505
RVT50611
SEN48742
BMMK00000000.1
MDEQ00000000.2
AJJQ00000000.1
SUMB00000000.1
JAPMUZ000000000.1
GPC_000001832
GPC_000011789
GPC_000011790
GPS_020388193
GPS_020388126
EOF
sed -i 's/\r$//' analyzed_contigs.txt
grep -vxFf /tmp/drop.txt analyzed_contigs.txt > t && mv t analyzed_contigs.txt
awk -F'\t' 'NR==FNR{d[$1];next} FNR==1 || !($2 in d)' \
    /tmp/drop.txt benchmark_ground_truth.tsv > t && mv t benchmark_ground_truth.tsv
wc -l analyzed_contigs.txt          # expect 1074
cd -
scripts/run_benchmark.sh pool_bact --label $USER --force   # parquet is reused
```

The permanent fix — rejecting these at ground-truth build time in
`prepare_mibig_ground_truth.py` — is in `TODO.md`.

---

## 2. Decide the headline scope: `pool_bact` or `benchmark_set_bact`

`pool_bact` keeps the BGC-only deposits (the record *is* the cluster, so every
tool scores ~1.0 on them by construction), which inflates detection recall.
The README recommends scoring through `benchmark_set_bact`, built with
`select_benchmark_genomes.py` at its default `--min-length`. That needs no new
runs — selection and scoring only — but it does not exist in this checkout yet;
recipe in `docs/BENCHMARK_SCOPES.md`.

Either build it, or report `pool_bact` with the caveat stated next to the number.

---

## 3. Then, optional

- **DeepBGC `--min-p-bgc` sweep.** Precision so far is at threshold 0.0. Each
  point costs seconds:
  `scripts/run_benchmark.sh pool_bact deepbgc --label p50 -- --min-p-bgc 0.5`.
- **Report figures** against the final scope: `build_report_tables.py mibig`,
  then `plot_report_figures.py` (commands in `CLAUDE.md` → *Running the
  benchmark*).
- **The extra 130 genomes** (`data/interim/extra_130/genomes.txt`, from
  `~/eduardo/benchmarks/missing_mibig_actino`). Both tools have run on all 130
  into the shared pool, so they are in the merged parquets, but **no scope lists
  them, so no benchmark JSON includes them**. Only 11 MiBiG clusters fall on
  them and all 11 are already in `pool_bact` under another accession/version
  (twins, e.g. `NC_003888.3` ↔ `AL645882.2`) — so they add nothing to the MiBiG
  benchmark; **do not** add them to `pool_bact`'s scope, it would only add
  unmatched calls. Their use is the no-ground-truth comparison
  (`summarize_predictions.py` / `build_report_tables.py raw`), which needs a
  small `genomes.tsv` manifest for them — not built yet.

---

## What landed since 2026-09-01

- `71f9ba8 fix(convert-deepbgc)` — lift csv's 128 KiB field cap; long DeepBGC
  candidates made the pool merge fail with `field larger than field limit`.
- `c1ed69b fix(run-deepbgc)` — Prodigal meta mode for inputs under 20 kb.
  BGC-only deposits produced no proteins and no `.bgc.tsv`, which failed the
  task ("DeepBGC finished but … is missing or empty").
- `d1d9b37`, `da9f153 feat(sbatch)` — larger per-task configs, `max50` partition.
- `4ec2f3d fix(merge-predictions)` — rename versioned prediction contigs onto
  unversioned scope entries (`X.1` → `X`). Recovered 49 antiSMASH / 48 DeepBGC
  contigs on `pool_bact` that were scoring a silent zero. The root cause in
  `select_benchmark_genomes.py` is still open in `TODO.md`.
- Pool runs: the second array window (last 87 of the pool list, offset
  `N - 87`), the 31 failed DeepBGC tasks from job 50668 resubmitted by index,
  and the 130 extra genomes.

## Still open from before

- **Wire the `.sbatch` scripts to `POOL_ROOT`.** They still hardcode
  `$HOME/projects/<tool>` while `run_benchmark.sh` reads `POOL_ROOT` from `.env`.
- **`--genus` takes only one substring.** The 34-genus actinomycete scope needs a
  loop; there is no builder for class-filtered scopes either, only the awk recipe
  in `BENCHMARK_SCOPES.md`.
- **GECCO paused** — no `run_gecco_array.sbatch`. If it returns, its converter
  needs the same csv field-cap fix as DeepBGC's (`proteins`/`domains` columns).
- **No S(H)ARP row yet** — `predict.py` and everything upstream of it are
  unimplemented, so every table stays baselines-only.
- **Pre-existing test failure**, unrelated:
  `test_extract_embeddings.py::TestRun::test_orchestration_end_to_end` (old
  NVIDIA driver → `UserWarning` → error under `filterwarnings = ["error"]`).
  498 passed / 1 failed as of `4ec2f3d`.
