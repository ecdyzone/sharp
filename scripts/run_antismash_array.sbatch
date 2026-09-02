#!/usr/bin/env bash
#
# Run antiSMASH across the whole benchmark genome set as a Slurm job array, one
# task per genome, on the davinci cluster.
#
# Submit (the array range must match the number of genomes):
#   N=$(wc -l < data/interim/benchmark_set/analyzed_contigs.txt)
#   sbatch --array=1-${N}%8 scripts/run_antismash_array.sbatch
#
# Or against another set / list:
#   sbatch --array=1-50%8 scripts/run_antismash_array.sbatch /path/to/contigs.txt
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
# `%8` throttles to 8 concurrent tasks, matching DeepBGC now that the measured
# sizing below brings each task down to 4 cores. Tune to the partition.
#
# ── Sizing ──────────────────────────────────────────────────────────────────
# MEASURED on the davinci cluster, `seff` on array job 45315, indices 1 and 2
# (AL645882.2, 8.7 Mb, and the index-2 genome):
#
#   index 1: 2:56 wall, CPU efficiency 5.40% of 16 cores, 1.62 GB peak
#   index 2: 3:05 wall, CPU efficiency 7.84% of 16 cores, 1.28 GB peak
#
# So antiSMASH averages ~1.3 cores of real parallelism and ~1.7 GB, not the 16
# cores / 16 GB originally reserved on the untested assumption that --cpus would
# be used. Same story as DeepBGC, which was sized at 8 cores and measured 0.86:
# the module scheduler parallelises far less than the flag suggests.
#
# Sized to 4 cores / 4G -- ~3x the measured CPU and ~2.4x the measured memory,
# so a heavier genome has room without hoarding a node. The largest genome in
# the set is ~11.9 Mb (CP002047.1) against the 8.7 Mb measured, ~1.4x the work,
# which 1 hour covers with an order of magnitude to spare.
#
# `%8` at 4 cores per task is 32 cores concurrently -- the same courtesy-limit
# reasoning as DeepBGC, not a technical ceiling. Raise it if the partition is
# idle.

#SBATCH --job-name=antismash-bench
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

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
    --cpus "${SLURM_CPUS_PER_TASK:-4}"

if [[ ! -s "$RESULT" ]]; then
    echo "ERROR: antiSMASH finished but $RESULT is missing or empty" >&2
    echo "Contents of $OUTDIR:" >&2
    ls -1 "$OUTDIR" >&2 || true
    exit 1
fi

echo
echo "Done: $RESULT"
