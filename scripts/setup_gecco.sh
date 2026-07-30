#!/usr/bin/env bash
set -euo pipefail

# Install dir comes from `.env` (see .env.example); fall back if unset.
# GECCO ships its models in the conda package, so there is no database download.
source "$(dirname "${BASH_SOURCE[0]}")/_load_env.sh"
: "${TOOLS_INSTALL_DIR:=$HOME/.local/src}"

INSTALL_DIR="$TOOLS_INSTALL_DIR"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

mkdir -p gecco
cd gecco

pixi init
pixi workspace channel add conda-forge
pixi workspace channel add bioconda

pixi add gecco

pixi install

echo "GECCO installation complete at $INSTALL_DIR/gecco"
echo ""
echo "Run GECCO from its own env, e.g.:"
echo "  cd $INSTALL_DIR/gecco && pixi run gecco run --genome <genome.fasta> --output-dir <out>"
echo "Then convert its output for benchmarking:"
echo "  pixi run python scripts/convert_gecco_to_parquet.py --input <out> --output data/interim/gecco_predictions.parquet"
