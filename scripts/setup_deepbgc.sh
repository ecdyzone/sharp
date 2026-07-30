#!/usr/bin/env bash
set -euo pipefail

# Install and downloads dirs come from `.env` (see .env.example); fall back if unset.
# `deepbgc download` has no path flag — it reads DEEPBGC_DOWNLOADS_DIR from the
# environment, which the loader exports for us.
source "$(dirname "${BASH_SOURCE[0]}")/_load_env.sh"
: "${TOOLS_INSTALL_DIR:=$HOME/.local/src}"
: "${DEEPBGC_DOWNLOADS_DIR:=$HOME/.local/share/deepbgc/data}"
export DEEPBGC_DOWNLOADS_DIR

INSTALL_DIR="$TOOLS_INSTALL_DIR"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Create project
mkdir -p deepbgc
cd deepbgc

# Initialize Pixi
pixi init

# Add channels
pixi workspace channel add conda-forge
pixi workspace channel add bioconda

# Add dependencies
pixi add python=3.7 hmmer prodigal pip

# Install DeepBGC from PyPI
pixi add --pypi deepbgc

# don't know why but solved some warnings
pixi add "protobuf=3.20.*"
pixi add --pypi "deepbgc[hmm]"

# add env var to pixi runs
cat >> pixi.toml <<EOF

[activation.env]
DEEPBGC_DOWNLOADS_DIR = "$DEEPBGC_DOWNLOADS_DIR"
EOF

# Solve and install environment
pixi install

# Verify installation
pixi run python --version
pixi run deepbgc --help



# Download DeepBGC models/data
# Before you can use DeepBGC, download trained models and Pfam database.
# Destination is DEEPBGC_DOWNLOADS_DIR, resolved at the top of this script.
echo "Downloading DeepBGC models and Pfam database (~3GB) to $DEEPBGC_DOWNLOADS_DIR ..."
pixi run deepbgc download # downloads almost 3GB
# You can display downloaded dependencies and models using:
pixi run deepbgc info

echo "DeepBGC installation complete at $INSTALL_DIR/deepbgc"
echo ""
echo "Run DeepBGC from its own env, e.g.:"
echo "  cd $INSTALL_DIR/deepbgc && pixi run deepbgc pipeline <genome.fasta> --output <out>"
echo "Then convert its output for benchmarking:"
echo "  pixi run python scripts/convert_deepbgc_to_parquet.py --input <out> --output data/interim/deepbgc_predictions.parquet"
