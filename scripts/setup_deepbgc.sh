#!/usr/bin/env bash
set -euo pipefail

# Install to ~/.local/src
INSTALL_DIR="${HOME}/.local/src"
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

# Solve and install environment
pixi install

# Verify installation
pixi run python --version
pixi run deepbgc --help

# Download DeepBGC models/data
# Before you can use DeepBGC, download trained models and Pfam database:
pixi run deepbgc download # downloads almost 3GB
# You can display downloaded dependencies and models using:
pixi run deepbgc info

echo "DeepBGC installation complete at $INSTALL_DIR/deepbgc"
echo ""
echo "Run DeepBGC from its own env, e.g.:"
echo "  cd $INSTALL_DIR/deepbgc && pixi run deepbgc pipeline <genome.fasta> --output <out>"
echo "Then convert its output for benchmarking:"
echo "  pixi run python scripts/convert_deepbgc_to_parquet.py --input <out> --output data/interim/deepbgc_predictions.parquet"
