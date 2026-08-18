#!/usr/bin/env bash
#
# Shared helper: fetch one nuccore record as FASTA and prove it is one.
# `source` this, do not execute it. Sourced by download_genome.sh and
# download_benchmark_genomes.sh so the retry and validation rules live in one
# place.
#
# Provides:
#   fetch_nuccore <ACCESSION> <OUTPUT_PATH> [MAX_ATTEMPTS]
#     Returns 0 on success, 1 on failure. Leaves no file behind on failure.
#
# Expects NCBI_API_KEY to already be in the environment if it is set at all
# (source _load_env.sh first). It is optional and raises the E-utilities rate
# limit from 3 to 10 requests/sec.

EFETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

fetch_nuccore() {
    local accession="$1"
    local output="$2"
    local max_attempts="${3:-3}"

    local query="db=nuccore&id=${accession}&rettype=fasta&retmode=text"
    if [[ -n "${NCBI_API_KEY:-}" ]]; then
        query="${query}&api_key=${NCBI_API_KEY}"
    fi

    local attempt
    for (( attempt = 1; attempt <= max_attempts; attempt++ )); do
        # -q keeps batch output readable; the caller reports progress itself.
        if wget -q -O "${output}" "${EFETCH_URL}?${query}"; then
            # efetch answers HTTP 200 with an error body for a bad accession, so
            # wget's exit status proves nothing. A real FASTA starts with '>';
            # anything else is an error page we must not leave on disk
            # pretending to be a genome.
            if [[ -s "${output}" ]] && [[ "$(head -c 1 "${output}")" == ">" ]]; then
                return 0
            fi
        fi
        if (( attempt < max_attempts )); then
            # NCBI throttles by returning a short error body rather than a 429,
            # so back off before retrying instead of hammering.
            sleep $(( attempt * 3 ))
        fi
    done

    echo "ERROR: ${accession} did not return FASTA after ${max_attempts} attempt(s). NCBI said:" >&2
    head -c 300 "${output}" 2>/dev/null >&2 || true
    echo >&2
    rm -f "${output}"
    return 1
}

# The contig id a tool will report for this FASTA: the first whitespace-delimited
# token of the header, which is what every baseline uses as its contig name and
# therefore what the --contigs scope file and the ground truth must agree with.
fasta_contig_ids() {
    grep '^>' "$1" | cut -c2- | cut -d' ' -f1
}
