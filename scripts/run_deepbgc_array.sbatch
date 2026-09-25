#!/usr/bin/env bash
#
# Run DeepBGC across the whole benchmark genome set as a Slurm job array, one
# task per genome, on the davinci cluster.
#
# Submit (the array range must match the number of genomes):
#   N=$(wc -l < data/interim/benchmark_set/analyzed_contigs.txt)
#   sbatch --array=1-${N}%8 scripts/run_deepbgc_array.sbatch
#
# Or against another set / list:
#   sbatch --array=1-50%8 scripts/run_deepbgc_array.sbatch /path/to/contigs.txt
#
# Args: $1 = genome list (default below), $2 = index offset (default 0),
#       $3 = genome directory (default below),
#       $4 = output pool root (default below).
#
# ── Running over a genome database instead of the MiBiG benchmark set ───────
# The list is one *key* per line and the genome is `$3/<key>.fasta`, so the same
# script covers both layouts. For a database of NCBI assembly dumps (one
# directory per assembly), `scripts/build_genome_manifest.py --link-dir` stages
# exactly that shape — `<link-dir>/<assembly>.fasta` symlinks — and emits the
# matching `assemblies.txt`:
#
#   sbatch --array=1-1000%24 scripts/run_deepbgc_array.sbatch \
#       "$SHARP/data/interim/actino_db/assemblies.txt" 0 \
#       "$SHARP/data/raw/actino_db" "$DEEPBGC_DIR/out_actino"
#
# The output pool is then keyed by assembly rather than by contig accession.
# That is fine for every downstream step except `merge_predictions.py`'s
# "produced no output" report, which needs `--manifest genomes.tsv` to map an
# assembly back to its contigs.
#
# ── Submitting a pool larger than MaxArraySize ──────────────────────────────
# Slurm's MaxArraySize (commonly 1001) caps the highest legal array *index* at
# 1000. It is not a concurrency limit, so waiting for a running array to drain
# does not raise it -- `--array=1001-1087` is rejected at parse time with
#
#   sbatch: error: Batch job submission failed: Invalid job array specification
#
# Submit in windows over the same list, moving the offset (arg 2) instead of
# cutting the list up. Chunk 2 covers lines 1001+ using indices 1..87 again:
#
#   scontrol show config | grep -i maxarraysize      # confirm the cap
#   mkdir -p logs                     # Slurm will NOT create it; tasks die if absent
#   CONTIGS=data/interim/pool_bact/analyzed_contigs.txt
#   N=$(wc -l < "$CONTIGS")
#   MAX=1000                          # MaxArraySize - 1
#
#   prev=""
#   for (( off = 0; off < N; off += MAX )); do
#       n=$(( N - off < MAX ? N - off : MAX ))
#       dep=${prev:+--dependency=afterany:$prev}
#       prev=$(sbatch --parsable $dep --array=1-${n}%16 \
#           scripts/run_deepbgc_array.sbatch "$CONTIGS" "$off")
#       echo "chunk at offset $off -> job $prev"
#   done
#
# Then convert every genome's output in one pass, from the S(H)ARP env:
#   pixi run python scripts/merge_predictions.py --tool deepbgc \
#       --input-dir ~/projects/deepbgc/out_benchmark \
#       --output data/interim/deepbgc_predictions.parquet
#
# ── Why an array and not one long job ───────────────────────────────────────
# 50 genomes at ~30 min each is ~25 hours serially, and one failure loses the
# lot. As an array each genome is an independent task: failures are isolated and
# retried by resubmitting just that index, and the scheduler backfills small
# 2-core tasks far more readily than one huge reservation.
#
# `%8` throttles to 8 concurrent tasks. It is a courtesy limit, not a technical
# one -- raise it if the partition is idle, lower it if you are sharing.
#
# ── Sizing ──────────────────────────────────────────────────────────────────
# Per task, carried over from the measured single-genome run (`seff` on job
# 42995, AL645882.2, 8.7 Mb): 0.86 cores, 1.71 GB, 32 min wall. DeepBGC is
# single-threaded end to end -- it shells out to hmmscan without threading it --
# so the per-task shape does not change with the array; only the count does.
#
# The largest genome in the Streptomyces set is ~11.9 Mb (CP002047.1) against the
# 8.7 Mb measured, so ~1.4x the work.
#
# Walltime is 6h against a measured 32 min. That is deliberate slack: the
# measurement is a Streptomyces chromosome, the bacterial pool reaches ~13 Mb
# (myxobacteria), DeepBGC's hmmscan stage is serial, and a contended node
# stretches wall-clock without changing CPU time. A timeout costs the task plus a
# retry; an unused reservation on a pinned node costs nothing. The one real cost
# is that a hung task holds a %-slot for 6h -- lower it if you see tasks hanging.
# Re-measure with `seff <jobid>_<index>` -- note the array index suffix.

#SBATCH --job-name=deepbgc-bench
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail

# ── Explicit paths ──────────────────────────────────────────────────────────
DEEPBGC_DIR="$HOME/projects/deepbgc"
SHARP_DIR="$HOME/projects/igem2026/git-clones/sharp"

# `deepbgc download` and `deepbgc pipeline` read this from the environment --
# there is no CLI flag for it. Set here so the job does not depend on the tool's
# pixi.toml being intact, or on DATABASES being exported to batch jobs.
export DEEPBGC_DOWNLOADS_DIR="$HOME/projects/databases/deepbgc/data"

# Read fresh by every task, at the moment that task starts -- NOT snapshotted
# at submit time like the script body is. Editing this file while the array is
# in flight shifts the line numbering under the tasks that have not started
# yet, silently skipping some genomes and running others twice. Point it at a
# frozen pool list and leave that file alone until the array drains.
CONTIGS="${1:-$SHARP_DIR/data/interim/benchmark_set/analyzed_contigs.txt}"

# Slurm caps the highest array index at MaxArraySize - 1 (commonly 1000), so a
# pool larger than that needs more than one submission. Only the starting point
# differs between them, so take it as an offset rather than slicing $CONTIGS
# into per-chunk files -- one list stays the single source of truth, and there
# is no derived copy to keep in sync or freeze separately.
OFFSET="${2:-0}"

# Where `<key>.fasta` lives. Overridable so this script can run over a genome
# database staged as symlinks (see the header) without a second copy of itself.
GENOME_DIR="${3:-$SHARP_DIR/data/raw/genomes}"
# The output pool root. Overridable so a separate campaign gets its own
# pool: merge_predictions.py converts EVERY directory under --input-dir
# (--contigs only drives the missing-genome report, it does not filter
# the parquet), so two genome sets sharing one pool end up in each
# other's predictions file.
OUTROOT="${4:-$DEEPBGC_DIR/out_benchmark}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: not running as a job array. Submit with --array, e.g." >&2
    echo "  sbatch --array=1-\$(wc -l < $CONTIGS)%8 \$0" >&2
    exit 1
fi

# One key per line (contig accession, or assembly for a genome database);
# this task takes the line at its index plus the offset.
LINE=$(( SLURM_ARRAY_TASK_ID + OFFSET ))
GENOME_KEY="$(sed -n "${LINE}p" "$CONTIGS")"
if [[ -z "$GENOME_KEY" ]]; then
    echo "ERROR: no genome key at line ${LINE} of $CONTIGS" >&2
    echo "       ($(wc -l < "$CONTIGS") lines; array index ${SLURM_ARRAY_TASK_ID}, offset ${OFFSET})" >&2
    exit 1
fi

GENOME="$GENOME_DIR/${GENOME_KEY}.fasta"
# DeepBGC names its output files after the basename of --output, so a per-genome
# directory gives <ACC>/<ACC>.bgc.tsv and keeps the runs from colliding.
OUTDIR="$OUTROOT/${GENOME_KEY}"
RESULT="$OUTDIR/${GENOME_KEY}.bgc.tsv"

# ── Preflight ───────────────────────────────────────────────────────────────
for path in "$DEEPBGC_DIR" "$DEEPBGC_DOWNLOADS_DIR" "$GENOME"; do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: does not exist on $(hostname): $path" >&2
        exit 1
    fi
done

echo "Host:      $(hostname)"
echo "Array idx: ${SLURM_ARRAY_TASK_ID} (+${OFFSET} = line ${LINE})"
echo "Genome key: $GENOME_KEY"
echo "Genome:    $GENOME"
echo "Output:    $OUTDIR"
echo "Pool:      $OUTROOT"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-?}"
echo

# Resumable: a completed genome is skipped, so a partially-failed array can be
# resubmitted wholesale instead of picking out indices by hand.
if [[ -s "$RESULT" ]]; then
    echo "Already done ($RESULT exists) — skipping."
    exit 0
fi

mkdir -p "$OUTROOT"
cd "$DEEPBGC_DIR"

pixi run deepbgc pipeline "$GENOME" --output "$OUTDIR"

if [[ ! -s "$RESULT" ]]; then
    echo "ERROR: DeepBGC finished but $RESULT is missing or empty" >&2
    exit 1
fi

echo
echo "Done: $RESULT"
