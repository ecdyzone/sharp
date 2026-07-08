#!/usr/bin/env bash
set -euo pipefail

# Install to ~/.local/src
INSTALL_DIR="${HOME}/.local/src"
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
