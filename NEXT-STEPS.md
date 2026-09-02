# NEXT STEPS

Working notes, **2026-09-01**. Scratch file — not a doc. Anything here that
outlives the week belongs in `TODO.md`, `docs/BENCHMARK_SCOPES.md` or `CLAUDE.md`
instead; this is only "where we stopped".

**State: the bacterial pool is being computed.** The genomes are downloaded, the
arrays are (re)submitted, and the next real work is the merge-time fix below,
which must land before any number can be trusted.

---

## 1. Blocking: version-suffix normalization in `merge_predictions.py`

**Nothing can be scored until this is done.** Full diagnosis in `TODO.md` →
*Known defects*; the short version:

51 of the 1,087 pool accessions are unversioned in the pool list, so the FASTA
NCBI returned carries a version the ground truth lacks (`JADBID010000001` vs
`>JADBID010000001.1`). Every tool copies that header into its output, so those
predictions no longer match `--contigs` — they get filtered out while their
clusters stay in the recall denominator. **51 genomes score a silent zero.**

The fix is in `merge_predictions.py`, not in the FASTAs: when a prediction's
contig is absent from the scope but its caption (`accession.split(".")[0]`) is
present, map it to the caption. Symmetric across tools, and it also absorbs a
future NCBI version bump.

Do **not** fix it by rewriting FASTA headers — the two pools would then disagree
with each other and need two different normalizations.

One thing to confirm first, once antiSMASH tasks have landed. This assumes
antiSMASH preserves the header verbatim, but GenBank LOCUS names cap at 16
characters and `JADBID010000001.1` is 17:

```bash
python3 -c "
import json, os
d = 'JADBID010000001'
j = os.path.expanduser(f'~/projects/antismash/out_benchmark/{d}/{d}.json')
print([r['id'] for r in json.load(open(j))['records']])
"
```

Expected `['JADBID010000001.1']`. Anything else and the normalization rule
changes.

---

## 2. Drop the 13 unusable accessions from the scope

8 accessions never downloaded (3 protein accessions, 5 WGS master records) and 5
`GPC_`/`GPS_` identifiers resolve to unrelated sequences. They must leave
`analyzed_contigs.txt` **and** `benchmark_ground_truth.tsv`, or they sit in the
recall denominator as clusters no tool was given a chance to find.

**Only while no array is in flight** — each task reads the contigs file with
`sed -n "${LINE}p"` when it *starts*, so editing it mid-run renumbers the lines
under pending tasks.

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
```

The permanent fix — teaching `prepare_mibig_ground_truth.py` to reject these at
ground-truth build time — is in `TODO.md`.

---

## 3. Check on the arrays

```bash
squeue -u $USER
sacct -X --format=JobID%18,JobName%18,State,Elapsed,MaxRSS -S today | tail -20
seff <jobid>_<index>              # note the array index suffix
```

Expect ~8 failed tasks if the drop in step 2 had not happened before submission
— those are the missing FASTAs, and they are harmless.

Three things that would have gone wrong overnight, worth ruling out first:

- **`logs/` missing.** The scripts now write there and Slurm does not create it;
  every task dies instantly if it is absent. `ls logs/ | wc -l`.
- **Disk.** Never verified. `du -sh ~/projects/deepbgc/out_benchmark`,
  `df -h ~/projects`.
- **Walltime.** antiSMASH 4h, DeepBGC 6h, both raised this session. A timeout
  leaves a partial output dir with no result file, which is correctly retried —
  but antiSMASH may refuse to write into a non-empty output dir, so `rm -rf` the
  directory before retrying.

Then merge and score:

```bash
pixi run python scripts/merge_predictions.py --tool antismash \
    --input-dir ~/projects/antismash/out_benchmark \
    --contigs data/interim/pool_bact/analyzed_contigs.txt \
    --output data/interim/antismash_predictions_pool_bact.parquet
scripts/run_benchmark.sh pool_bact
```

---

## 4. Then: carve the scopes

The pool exists to be sliced. `docs/BENCHMARK_SCOPES.md` is the catalogue — the
genus scopes, the per-class slices, the negative controls and the exact commands
that build each. Scoring a scope is seconds; only the pool was expensive.

Note the CRLF trap when hand-building a class scope (`TODO.md` → *Known
defects*): every `*_ground_truth.tsv` is CRLF, so an `awk` filter on the last
column silently matches nothing. The recipes in `BENCHMARK_SCOPES.md` already
strip it with `{sub(/\r$/, "")}`.

---

## What landed this session

Seven commits, `86c6951..bb0f736` — history was rewritten at the end (rebase),
so the server needs `git fetch origin && git reset --hard origin/main`, not a
plain `git pull`. Nothing running is affected: Slurm snapshots batch scripts at
submit time.

- `fix(arrays)` — Slurm task logs to `logs/`, and `logs/`/`*.out`/`*.err`
  gitignored. At ~1,100 tasks that was ~2,200 files in the repo root.
- `feat(arrays)` — **index offset as `$2`**. Slurm's `MaxArraySize` caps the
  highest legal array *index* at 1000, so a 1,087-genome pool cannot be submitted
  as `--array=1001-1087`; successive windows now reuse indices `1..1000` over the
  one list. Also documents that the contigs file is read per task at task start.
- `chore(arrays)` — walltimes 1h→4h (antiSMASH), 2h→6h (DeepBGC).
- `docs(todo)` — the CRLF defect, plus the two accession defects from this
  session.
- The `.env.example` change was folded back into `docs(env)` where it belongs.

Backup branch `backup/pre-rebase` still points at the pre-rewrite tip; delete it
once the force-push has settled.

## Still open from before

- **Sweep `--min-p-bgc` for DeepBGC.** Every precision figure so far is at
  threshold 0.0 — a ranked list compared against rule-based antiSMASH, which has
  no score. Recall is unaffected and reportable now. The wrapper reuses the
  parquet, so each point costs seconds; use `--label` to keep the curve.
- **Wire the `.sbatch` scripts to `POOL_ROOT`.** They still hardcode
  `$HOME/projects/<tool>` while `run_benchmark.sh` reads `POOL_ROOT` from `.env`.
  Override one and the wrapper reads a pool the arrays never wrote to.
- **`--genus` takes only one substring.** The 34-genus actinomycete scope needs a
  loop; there is no builder for class-filtered scopes either, only the awk recipe
  in `BENCHMARK_SCOPES.md`.
- **GECCO paused** — no `run_gecco_array.sbatch`. Mirror the DeepBGC one if it
  returns; GECCO is the tool needing `start - 1`.
- **No S(H)ARP row yet** — `predict.py` and everything upstream of it are
  unimplemented, so every table stays baselines-only.
- **Pre-existing test failure**, unrelated to any of this:
  `test_extract_embeddings.py::TestRun::test_orchestration_end_to_end`. An old
  NVIDIA driver makes `torch.cuda.is_available()` emit a `UserWarning`, which
  `filterwarnings = ["error"]` (`pyproject.toml:18`) promotes to a failure.
  364 passed / 1 failed, before and after.
