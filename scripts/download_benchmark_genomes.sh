#!/usr/bin/env bash

# Download every genome in a benchmark set produced by
# scripts/select_benchmark_genomes.py.
#
# Reads the `accession` column of benchmark_genomes.tsv and fetches one FASTA
# per genome into data/raw/genomes/. Resumable: a genome that already has a
# valid FASTA on disk is skipped, so re-running after an interruption or a
# rate-limit failure costs nothing.
#
# Usage:
#   scripts/download_benchmark_genomes.sh
#   scripts/download_benchmark_genomes.sh path/to/benchmark_genomes.tsv
#   scripts/download_benchmark_genomes.sh path/to/set.tsv path/to/outdir
#
# Outputs:
#   data/raw/genomes/<ACCESSION>.fasta   (one per selected genome)
#
# ── Why download rather than read a mirrored NCBI tree ──────────────────────
# A local NCBI FTP mirror is indexed by RefSeq *assembly* (GCF_*), while the
# ground truth is keyed by *nucleotide* accession, and 44 of the 50 selected
# accessions are GenBank-style. Going through a mirror therefore needs a
# GenBank->RefSeq contig-name translation, and any mistake there is silent: the
# contig names in the FASTA stop matching the ground truth and every genome
# scores zero.
#
# Fetching by nucleotide accession sidesteps that entirely — the FASTA header
# *is* the accession the ground truth names, which this script verifies. It also
# pins the benchmark to accession.version rather than to whenever a mirror was
# last rsynced. At ~420 Mb for 50 genomes the download is a one-off cost.
#
# Revisit the mirror if the benchmark ever scales to thousands of genomes, where
# downloading stops being reasonable.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SET_TSV="${1:-${REPO_ROOT}/data/interim/benchmark_set/benchmark_genomes.tsv}"
OUT_DIR="${2:-${REPO_ROOT}/data/raw/genomes}"

if [[ ! -f "${SET_TSV}" ]]; then
    echo "ERROR: benchmark set not found: ${SET_TSV}" >&2
    echo "Build it first:" >&2
    echo "  pixi run python scripts/select_benchmark_genomes.py \\" >&2
    echo "      --ground-truth data/raw/streptomyces_ground_truth.tsv \\" >&2
    echo "      --output-dir data/interim/benchmark_set" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# shellcheck source=scripts/_load_env.sh
source "${SCRIPT_DIR}/_load_env.sh"
# shellcheck source=scripts/_fetch_nuccore.sh
source "${SCRIPT_DIR}/_fetch_nuccore.sh"

# Without an API key E-utilities allows 3 requests/sec; with one, 10. Stay
# under whichever applies rather than relying on the retry path to absorb it.
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    SLEEP_BETWEEN=0.15
else
    SLEEP_BETWEEN=0.4
fi

# Column 2 is `accession` — the primary accession of each merged genome, and the
# name the normalized ground truth and the --contigs scope file both use.
mapfile -t ACCESSIONS < <(tail -n +2 "${SET_TSV}" | cut -f2 | grep -v '^$')

TOTAL=${#ACCESSIONS[@]}
if (( TOTAL == 0 )); then
    echo "ERROR: no accessions in ${SET_TSV}" >&2
    exit 1
fi

echo "Benchmark set: ${SET_TSV}"
echo "Genomes:       ${TOTAL}"
echo "Output:        ${OUT_DIR}"
echo

n_skipped=0
n_fetched=0
failed=()
mismatched=()

for i in "${!ACCESSIONS[@]}"; do
    acc="${ACCESSIONS[$i]}"
    fasta="${OUT_DIR}/${acc}.fasta"
    label="[$((i + 1))/${TOTAL}] ${acc}"

    # Resume: trust an existing file only if it still looks like FASTA, so a
    # truncated download from an interrupted run is retried rather than kept.
    if [[ -s "${fasta}" ]] && [[ "$(head -c 1 "${fasta}")" == ">" ]]; then
        echo "${label} — already present, skipping"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    echo -n "${label} — downloading... "
    if fetch_nuccore "${acc}" "${fasta}"; then
        bp=$(grep -v '^>' "${fasta}" | tr -d '\n' | wc -c)
        echo "ok (${bp} bp)"
        n_fetched=$((n_fetched + 1))
    else
        echo "FAILED"
        failed+=("${acc}")
        continue
    fi

    # The check that matters: the contig id the tools will report must equal the
    # accession the ground truth names. If NCBI answers with a different header
    # the benchmark would score zero on this genome for reasons invisible in the
    # metrics, so surface it here.
    got=$(fasta_contig_ids "${fasta}" | head -1)
    if [[ "${got}" != "${acc}" ]]; then
        echo "  WARNING: header is '${got}', expected '${acc}'" >&2
        mismatched+=("${acc} (got ${got})")
    fi

    sleep "${SLEEP_BETWEEN}"
done

echo
echo "Fetched ${n_fetched}, skipped ${n_skipped}, failed ${#failed[@]}."

if (( ${#mismatched[@]} > 0 )); then
    echo
    echo "Contig-name mismatches — these will NOT match the ground truth:" >&2
    printf '  %s\n' "${mismatched[@]}" >&2
fi

if (( ${#failed[@]} > 0 )); then
    echo
    echo "Failed accessions (re-run to retry just these):" >&2
    printf '  %s\n' "${failed[@]}" >&2
    exit 1
fi

echo
echo "All ${TOTAL} genomes present in ${OUT_DIR}"
echo "Total: $(du -sh "${OUT_DIR}" | cut -f1)"
