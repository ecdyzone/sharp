#!/usr/bin/env bash
#
# Run DeepBGC on the benchmark genome, on the davinci cluster.
#
# Submit from anywhere:
#   sbatch ~/projects/igem2026/git-clones/sharp/scripts/run_deepbgc.sbatch
#
# Optionally override the genome (default: the S. coelicolor benchmark contig):
#   sbatch scripts/run_deepbgc.sbatch /path/to/other.fasta
#
# Then convert the output for benchmarking, from the S(H)ARP env:
#   pixi run python scripts/convert_deepbgc_to_parquet.py \
#       --input ~/projects/deepbgc/out/out.bgc.tsv \
#       --output data/interim/deepbgc_predictions.parquet
#
# ── Why CPU-only, on n01 ────────────────────────────────────────────────────
# DeepBGC is Prodigal (single-threaded, CPU) -> hmmscan vs Pfam (CPU-bound,
# dominates runtime on an 8.7 Mb genome) -> a small Keras classifier scoring a
# few thousand Pfam-domain vectors. Only the last stage could use a GPU, and it
# is a rounding error in the total, so a GPU allocation would sit idle. The env
# is pinned to python=3.7 (see scripts/setup_deepbgc.sh), which predates any
# CUDA-enabled TensorFlow that could address an L40S or A30 anyway.
#
# n01 is the only GPU-free node, so running here leaves gn01-03 for work that
# is actually GPU-bound (the ESM-2 embedding step).

#SBATCH --job-name=deepbgc-scoe
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# ── Explicit paths ──────────────────────────────────────────────────────────
DEEPBGC_DIR="$HOME/projects/deepbgc"
SHARP_DIR="$HOME/projects/igem2026/git-clones/sharp"

# `deepbgc download` and `deepbgc pipeline` read this from the environment --
# there is no CLI flag for it. setup_deepbgc.sh also bakes it into the tool's
# pixi.toml [activation.env], but set it here too so this job does not depend on
# that file being intact, or on DATABASES being exported to batch jobs.
export DEEPBGC_DOWNLOADS_DIR="$HOME/projects/databases/deepbgc/data"

# First argument overrides the genome; otherwise the benchmark contig.
GENOME="${1:-$SHARP_DIR/data/raw/AL645882.2.fasta}"

# DeepBGC output goes next to the tool, per the run layout.
OUTDIR="$DEEPBGC_DIR/out"

# ── Preflight ───────────────────────────────────────────────────────────────
# Check before the expensive stages: a missing path here fails in seconds
# rather than after Prodigal and a partial Pfam scan.
for path in "$DEEPBGC_DIR" "$DEEPBGC_DOWNLOADS_DIR" "$GENOME"; do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: does not exist on $(hostname): $path" >&2
        exit 1
    fi
done

echo "Host:     $(hostname)"
echo "Genome:   $GENOME"
echo "Models:   $DEEPBGC_DOWNLOADS_DIR"
echo "Output:   $OUTDIR"
echo "CPUs:     ${SLURM_CPUS_PER_TASK:-?}"
echo

cd "$DEEPBGC_DIR"

# Confirm the models actually resolve from this node before running.
pixi run deepbgc info

pixi run deepbgc pipeline "$GENOME" --output "$OUTDIR"

echo
echo "Done. Convert with:"
echo "  cd $SHARP_DIR && pixi run python scripts/convert_deepbgc_to_parquet.py \\"
echo "      --input $OUTDIR/out.bgc.tsv \\"
echo "      --output data/interim/deepbgc_predictions.parquet"
