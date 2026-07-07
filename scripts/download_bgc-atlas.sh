#!/usr/bin/env bash

# FILES DOWNLOADED FROM:
# https://bgc-atlas.cs.uni-tuebingen.de/downloads
#
# complete-bgcs.tar.gz (3.5 GB)
# GBK files for complete BGCs only (partial/fragmented BGCs excluded).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/../data/raw"

mkdir -p "${RAW_DIR}"

# download
wget -O "${RAW_DIR}/complete-bgcs.tar.gz" \
  https://bgc-atlas.cs.uni-tuebingen.de/downloads/complete-bgcs.tar.gz

# extract
tar -xzf "${RAW_DIR}/complete-bgcs.tar.gz" -C "${RAW_DIR}"

# remove archive
rm -f "${RAW_DIR}/complete-bgcs.tar.gz"
