#!/usr/bin/env bash
#
# Run antiSMASH on the benchmark genome, on the davinci cluster.
#
# Submit from anywhere:
#   sbatch ~/projects/igem2026/git-clones/sharp/scripts/run_antismash.sbatch
#
# Optionally override the genome (default: the S. coelicolor benchmark contig):
#   sbatch scripts/run_antismash.sbatch /path/to/other.fasta
#
# Then convert the output for benchmarking, from the S(H)ARP env:
#   pixi run python scripts/convert_antismash_to_parquet.py \
#       --input ~/projects/antismash/out \
#       --output data/interim/antismash_predictions.parquet
#
# ── Why CPU-only, and why more cores than DeepBGC ───────────────────────────
# antiSMASH is Prodigal gene calling -> HMM detection rules -> per-region
# analysis modules (NRPS/PKS domain prediction, ClusterBlast comparisons against
# the reference DB, etc.). All of it is CPU-bound HMMER and DIAMOND work with no
# GPU code path at all, so as with DeepBGC there is nothing for an L40S or A30
# to do here -- leave the GPU nodes for the ESM-2 embedding step.
#
# Unlike DeepBGC, antiSMASH takes an explicit `--cpus` and parallelises its
# analysis modules across it, so a wider allocation genuinely helps. 16 is a
# reasonable point on the curve for a single 8.7 Mb chromosome; the ClusterBlast
# stage is also I/O-heavy against the 10 GB database, so returns flatten well
# before the node's 64 threads.

#SBATCH --job-name=antismash-scoe
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# ── Explicit paths ──────────────────────────────────────────────────────────
ANTISMASH_DIR="$HOME/projects/antismash"
SHARP_DIR="$HOME/projects/igem2026/git-clones/sharp"
ANTISMASH_DB="$HOME/projects/databases/antismash/databases"

# First argument overrides the genome; otherwise the benchmark contig.
GENOME="${1:-$SHARP_DIR/data/raw/AL645882.2.fasta}"

# antiSMASH output goes next to the tool, matching the DeepBGC run layout.
OUTDIR="$ANTISMASH_DIR/out"

# ── Preflight ───────────────────────────────────────────────────────────────
# Check before the expensive stages: a missing path here fails in seconds
# rather than part-way through a multi-hour analysis.
for path in "$ANTISMASH_DIR" "$ANTISMASH_DB" "$GENOME"; do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: does not exist on $(hostname): $path" >&2
        exit 1
    fi
done

echo "Host:      $(hostname)"
echo "Genome:    $GENOME"
echo "Databases: $ANTISMASH_DB"
echo "Output:    $OUTDIR"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-?}"
echo

cd "$ANTISMASH_DIR"

# --genefinding-tool prodigal is required for a plain FASTA input: without gene
# annotations antiSMASH has no CDS features to run its detection rules over.
# (A GenBank input already carries them, and then this flag must be omitted.)
pixi run antismash "$GENOME" \
    --output-dir "$OUTDIR" \
    --genefinding-tool prodigal \
    --databases "$ANTISMASH_DB" \
    --cpus "${SLURM_CPUS_PER_TASK:-16}"

echo
echo "Done. Convert with (the converter globs the output dir for the summary"
echo ".json, which antiSMASH names after the input FASTA):"
echo "  cd $SHARP_DIR && pixi run python scripts/convert_antismash_to_parquet.py \\"
echo "      --input $OUTDIR \\"
echo "      --output data/interim/antismash_predictions.parquet"
