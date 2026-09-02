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
# The largest genome in the set is ~11.9 Mb (CP002047.1) against the 8.7 Mb that
# was measured, so ~1.4x the work. 2 hours keeps the same ~3.7x headroom the
# single-genome script deliberately carries. Re-measure with
# `seff <jobid>_<index>` -- note the array index suffix.

#SBATCH --job-name=deepbgc-bench
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
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

CONTIGS="${1:-$SHARP_DIR/data/interim/benchmark_set/analyzed_contigs.txt}"
GENOME_DIR="$SHARP_DIR/data/raw/genomes"
OUTROOT="$DEEPBGC_DIR/out_benchmark"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: not running as a job array. Submit with --array, e.g." >&2
    echo "  sbatch --array=1-\$(wc -l < $CONTIGS)%8 \$0" >&2
    exit 1
fi

# One accession per line; this task takes the line matching its array index.
ACCESSION="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIGS")"
if [[ -z "$ACCESSION" ]]; then
    echo "ERROR: no accession at line ${SLURM_ARRAY_TASK_ID} of $CONTIGS" >&2
    exit 1
fi

GENOME="$GENOME_DIR/${ACCESSION}.fasta"
# DeepBGC names its output files after the basename of --output, so a per-genome
# directory gives <ACC>/<ACC>.bgc.tsv and keeps the runs from colliding.
OUTDIR="$OUTROOT/${ACCESSION}"
RESULT="$OUTDIR/${ACCESSION}.bgc.tsv"

# ── Preflight ───────────────────────────────────────────────────────────────
for path in "$DEEPBGC_DIR" "$DEEPBGC_DOWNLOADS_DIR" "$GENOME"; do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: does not exist on $(hostname): $path" >&2
        exit 1
    fi
done

echo "Host:      $(hostname)"
echo "Array idx: ${SLURM_ARRAY_TASK_ID}"
echo "Accession: $ACCESSION"
echo "Genome:    $GENOME"
echo "Output:    $OUTDIR"
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
