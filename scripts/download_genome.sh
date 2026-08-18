#!/usr/bin/env bash

# Download a single genome contig by accession from NCBI, for benchmark runs.
#
# Default: AL645882.2 — Streptomyces coelicolor A3(2) chromosome. It carries 15
# coordinate-resolved MiBIG clusters in data/raw/streptomyces_ground_truth.tsv,
# the most of any single contig, which is what makes it a usable recall target
# (see TODO.md: AL589148.1 has exactly one, so every recall number from it is
# 0/1 or 1/1).
#
# Also writes the `--contigs` scope file that every tool in a comparison must
# share. Omitting --contigs makes recall optimistic and lets the three baselines
# drift apart — see CLAUDE.md "Benchmark scope caveat".
#
# Usage:
#   scripts/download_genome.sh                 # AL645882.2, the default
#   scripts/download_genome.sh CP002993.1      # any other nuccore accession
#
# Outputs:
#   data/raw/<ACCESSION>.fasta
#   data/interim/analyzed_contigs.txt
#
# ── Why efetch and not ncbi-datasets-cli ────────────────────────────────────
# This fetches one nuccore record, which is exactly the unit the benchmark is
# scoped to: one contig, one --contigs line, one recall denominator. It needs no
# new pixi dependency.
#
# The alternative is `ncbi-datasets-cli` (bioconda) fetching a whole assembly:
#
#   datasets download genome accession GCF_000203835.1 --include genome
#
# That returns the chromosome plus the SCP1/SCP2 plasmids, and you would then
# subset to the contig you want. Prefer it once a later pipeline step
# (annotate.py onward) actually needs the complete assembly rather than the one
# benchmark contig; for now it is a dependency and an extra subsetting step for
# no benefit.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/../data/raw"
INTERIM_DIR="${SCRIPT_DIR}/../data/interim"

# Accession to fetch. Override by passing one as the first argument.
ACCESSION="${1:-AL645882.2}"

FASTA="${RAW_DIR}/${ACCESSION}.fasta"
CONTIGS="${INTERIM_DIR}/analyzed_contigs.txt"

mkdir -p "${RAW_DIR}" "${INTERIM_DIR}"

# NCBI_API_KEY is optional — it raises the E-utilities rate limit from 3 to 10
# requests/sec. One download does not need it, but .env reserves the key, so
# honour it when set. `source`d rather than parsed so .env's shell syntax works.
# shellcheck source=scripts/_load_env.sh
source "${SCRIPT_DIR}/_load_env.sh"

# Fetch and validation live in _fetch_nuccore.sh, shared with the batch
# downloader so both apply the same retry and "is this really FASTA" rules.
# shellcheck source=scripts/_fetch_nuccore.sh
source "${SCRIPT_DIR}/_fetch_nuccore.sh"

echo "Downloading ${ACCESSION} from NCBI nuccore..."
fetch_nuccore "${ACCESSION}" "${FASTA}" || exit 1

# The --contigs scope file: one contig id per line, matching how the tools name
# contigs in their output (first whitespace-delimited token of the FASTA header).
fasta_contig_ids "${FASTA}" > "${CONTIGS}"

echo
echo "Genome:  ${FASTA}"
echo "Contigs: ${CONTIGS}"
sed 's/^/  /' "${CONTIGS}"
echo
echo "Pass the same --contigs file to every tool in a comparison."
