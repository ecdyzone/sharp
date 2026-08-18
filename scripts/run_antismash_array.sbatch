#!/usr/bin/env bash
#
# Run antiSMASH across the whole benchmark genome set as a Slurm job array, one
# task per genome, on the davinci cluster.
#
# Submit (the array range must match the number of genomes):
#   N=$(wc -l < data/interim/benchmark_set/analyzed_contigs.txt)
#   sbatch --array=1-${N}%4 scripts/run_antismash_array.sbatch
#
# Or against another set / list:
#   sbatch --array=1-50%4 scripts/run_antismash_array.sbatch /path/to/contigs.txt
#
# Then convert every genome's output in one pass, from the S(H)ARP env:
#   pixi run python scripts/merge_predictions.py --tool antismash \
#       --input-dir ~/projects/antismash/out_benchmark \
#       --output data/interim/antismash_predictions.parquet
#
# ── Why an array and not one long job ───────────────────────────────────────
# Same reasoning as run_deepbgc_array.sbatch: independent tasks isolate
# failures, resubmitting retries only what failed, and the scheduler backfills
# many small tasks more readily than one long reservation.
#
# `%4` rather than DeepBGC's `%8` because each task asks for 16 cores; four
# concurrent tasks is already 64 cores. Tune to the partition.
#
# ── Sizing ──────────────────────────────────────────────────────────────────
# STILL UNMEASURED, as in the single-genome script. Run `seff <jobid>_<index>`
# on the first few completed tasks and cut hard if the CPU efficiency comes back
# low -- DeepBGC was sized at 8 cores on the reasoning that hmmscan would
# parallelise and measured 0.86.
#
# An array makes that measurement cheap: submit a couple of indices first
#   sbatch --array=1-2 scripts/run_antismash_array.sbatch
# check `seff`, then submit the rest with the corrected --cpus-per-task.
#
# 12 hours per task carries over from the single-genome script. The largest
# genome here is ~11.9 Mb against the 8.7 Mb benchmarked, so ~1.4x the work.

#SBATCH --job-name=antismash-bench
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

set -euo pipefail

# ── Explicit paths ──────────────────────────────────────────────────────────
ANTISMASH_DIR="$HOME/projects/antismash"
SHARP_DIR="$HOME/projects/igem2026/git-clones/sharp"
ANTISMASH_DB="$HOME/projects/databases/antismash/databases"

CONTIGS="${1:-$SHARP_DIR/data/interim/benchmark_set/analyzed_contigs.txt}"
GENOME_DIR="$SHARP_DIR/data/raw/genomes"
OUTROOT="$ANTISMASH_DIR/out_benchmark"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: not running as a job array. Submit with --array, e.g." >&2
    echo "  sbatch --array=1-\$(wc -l < $CONTIGS)%4 \$0" >&2
    exit 1
fi

# One accession per line; this task takes the line matching its array index.
ACCESSION="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIGS")"
if [[ -z "$ACCESSION" ]]; then
    echo "ERROR: no accession at line ${SLURM_ARRAY_TASK_ID} of $CONTIGS" >&2
    exit 1
fi

GENOME="$GENOME_DIR/${ACCESSION}.fasta"
OUTDIR="$OUTROOT/${ACCESSION}"
# antiSMASH names its summary JSON after the input FASTA basename.
RESULT="$OUTDIR/${ACCESSION}.json"

# ── Preflight ───────────────────────────────────────────────────────────────
for path in "$ANTISMASH_DIR" "$ANTISMASH_DB" "$GENOME"; do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: does not exist on $(hostname): $path" >&2
        exit 1
    fi
done

echo "Host:      $(hostname)"
echo "Array idx: ${SLURM_ARRAY_TASK_ID}"
echo "Accession: $ACCESSION"
echo "Genome:    $GENOME"
echo "Databases: $ANTISMASH_DB"
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
cd "$ANTISMASH_DIR"

# --genefinding-tool prodigal is required for a plain FASTA input: without gene
# annotations antiSMASH has no CDS features to run its detection rules over.
# (A GenBank input already carries them, and then this flag must be omitted.)
pixi run antismash "$GENOME" \
    --output-dir "$OUTDIR" \
    --genefinding-tool prodigal \
    --databases "$ANTISMASH_DB" \
    --cpus "${SLURM_CPUS_PER_TASK:-16}"

if [[ ! -s "$RESULT" ]]; then
    echo "ERROR: antiSMASH finished but $RESULT is missing or empty" >&2
    echo "Contents of $OUTDIR:" >&2
    ls -1 "$OUTDIR" >&2 || true
    exit 1
fi

echo
echo "Done: $RESULT"
