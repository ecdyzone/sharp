#!/usr/bin/env bash
set -euo pipefail

# Install and database dirs come from `.env` (see .env.example); fall back if unset.
source "$(dirname "${BASH_SOURCE[0]}")/_load_env.sh"
: "${TOOLS_INSTALL_DIR:=$HOME/.local/src}"
: "${ANTISMASH_DOWNLOADS_DIR:=$HOME/.local/share/antismash/databases}"

INSTALL_DIR="$TOOLS_INSTALL_DIR"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

mkdir -p antismash
cd antismash

pixi init

pixi workspace channel add conda-forge
pixi workspace channel add bioconda

pixi add antismash

pixi install

echo "antiSMASH installation complete at $INSTALL_DIR/antismash"
echo ""
echo "Downloading antiSMASH databases (10GB) to $ANTISMASH_DOWNLOADS_DIR ..."
pixi run download-antismash-databases --database-dir "$ANTISMASH_DOWNLOADS_DIR"

echo ""
echo "Run antiSMASH from its own env, e.g.:"
echo "  cd $INSTALL_DIR/antismash && pixi run antismash <genome.gbk> --output-dir <out>"
echo "  or"
echo "  cd $INSTALL_DIR/antismash && pixi run antismash <genome.fasta> --output-dir <out> --genefinding-tool prodigal"
echo "Then convert its output for benchmarking:"
echo "  pixi run python scripts/convert_antismash_to_parquet.py --input <out> --output data/interim/antismash_predictions.parquet"
