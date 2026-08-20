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
# ── Why CPU-only ────────────────────────────────────────────────────────────
# antiSMASH is Prodigal gene calling -> HMM detection rules -> per-region
# analysis modules (NRPS/PKS domain prediction, ClusterBlast comparisons against
# the reference DB, etc.). All of it is CPU-bound HMMER and DIAMOND work with no
# GPU code path at all, so as with DeepBGC there is nothing for an L40S or A30
# to do here -- leave the GPU nodes for the ESM-2 embedding step.
#
# ── Sizing, measured ────────────────────────────────────────────────────────
# MEASURED on the davinci cluster. This script previously reserved 16 cores /
# 16G on the untested argument that antiSMASH takes an explicit `--cpus` and
# hands it to its own module scheduler. That argument turned out to be wrong.
#
# `seff` on array job 45315 index 1 ran this script's own default genome
# (AL645882.2, 8.7 Mb): 2:56 wall, CPU efficiency 5.40% of 16 cores, 1.62 GB
# peak. Index 2 agreed at 7.84% and 1.28 GB. So antiSMASH delivers ~1.3 cores
# of real parallelism regardless of what `--cpus` is set to -- the same story as
# DeepBGC, sized at 8 cores and measured at 0.86 (job 42995, 10.73%).
#
# Sized to 4 cores / 4G to match run_antismash_array.sbatch: ~3x the measured
# CPU and ~2.4x the measured memory, so an override genome heavier than the
# 8.7 Mb default has room without hoarding a node.
#
# Time stays at 12h rather than the array's 1h. The array runs a known set whose
# largest member is ~11.9 Mb, while this script takes an arbitrary genome as $1
# -- so the walltime is deliberately loose here. Measure your own runs with
# `seff <jobid>` if you feed it something much larger.
#
# Memory is 4G rather than the 64G first written here. ClusterBlast and HMMER
# stream against the on-disk databases instead of loading the 10 GB into RAM,
# which the 1.62 GB measurement confirms.

#SBATCH --job-name=antismash-scoe
#SBATCH --partition=basic
#SBATCH --nodelist=n01
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=12:00:00
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
    --cpus "${SLURM_CPUS_PER_TASK:-4}"

echo
echo "Done. Convert with (the converter globs the output dir for the summary"
echo ".json, which antiSMASH names after the input FASTA):"
echo "  cd $SHARP_DIR && pixi run python scripts/convert_antismash_to_parquet.py \\"
echo "      --input $OUTDIR \\"
echo "      --output data/interim/antismash_predictions.parquet"
