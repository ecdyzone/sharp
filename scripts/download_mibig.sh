#!/usr/bin/env bash

# FILES DOWNLOADED FROM:
# https://mibig.secondarymetabolites.org/download
#
# Version 4.0 (November 15, 2024)
# All entries in JSON format.
# All entries in GBK format.
# All genes from MIBiG entries in FASTA format (amino acid sequences).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/../data/raw"

mkdir -p "${RAW_DIR}"

# download files
wget -O "${RAW_DIR}/mibig_json_4.0.tar.gz" \
  https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz

wget -O "${RAW_DIR}/mibig_gbk_4.0.tar.gz" \
  https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz

wget -O "${RAW_DIR}/mibig_prot_seqs_4.0.fasta" \
  https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_4.0.fasta

# extract archives (they contain top-level directories already)
tar -xzf "${RAW_DIR}/mibig_json_4.0.tar.gz" -C "${RAW_DIR}"
tar -xzf "${RAW_DIR}/mibig_gbk_4.0.tar.gz" -C "${RAW_DIR}"

# remove archives
rm -f "${RAW_DIR}/mibig_json_4.0.tar.gz"
rm -f "${RAW_DIR}/mibig_gbk_4.0.tar.gz"
