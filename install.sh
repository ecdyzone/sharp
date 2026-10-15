#!/usr/bin/env bash
# Install an independent benchmark snapshot.
# Usage:
#   bash install.sh [DESTINATION]
# Default destination: ~/sharp-benchmarks
set -euo pipefail

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo 'Usage: bash install.sh [DESTINATION]'
    echo 'Copies benchmark source and installs its Pixi environment.'
    echo 'Default: ~/sharp-benchmarks. Existing destinations are refused.'
    exit 0
fi
if [[ $# -gt 1 || ( $# -eq 1 && -z "$1" ) ]]; then
    echo 'Usage: bash install.sh [DESTINATION]' >&2
    exit 1
fi

for command in pixi tar; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: $command is required; install it before running this script." >&2
        exit 1
    fi
done

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
DESTINATION="${1:-$HOME/sharp-benchmarks}"
if [[ -e "$DESTINATION" || -L "$DESTINATION" ]]; then
    echo "ERROR: destination already exists: $DESTINATION" >&2
    echo 'Choose a new destination; existing data and results are never overwritten.' >&2
    exit 1
fi

mkdir -p "$(dirname "$DESTINATION")"
DESTINATION="$(cd "$(dirname "$DESTINATION")" && pwd -P)/$(basename "$DESTINATION")"
case "$DESTINATION/" in
    "$SOURCE_DIR/"*)
        echo 'ERROR: choose a destination outside the benchmark source folder.' >&2
        exit 1
        ;;
esac

# An explicit list keeps private configuration, data, Git metadata and local
# environments out of the delivery. Exclude caches inside source directories too.
mkdir "$DESTINATION"
tar -C "$SOURCE_DIR" \
    --exclude='__pycache__' --exclude='*.pyc' --exclude='*.pyo' \
    --exclude='.ipynb_checkpoints' --exclude='*.egg-info' \
    -cf - README.md CLAUDE.md LICENSE .gitignore .env.example \
    pyproject.toml pixi.toml pixi.lock environment.yml install.sh \
    src scripts tests docs notebooks |
    tar -C "$DESTINATION" -xf -
mkdir -p "$DESTINATION/data/"{raw,interim,processed,mock}
cp "$SOURCE_DIR/data/.gitignore" "$DESTINATION/data/.gitignore"

echo "Installing benchmark environment in $DESTINATION"
if ! pixi install --manifest-path "$DESTINATION/pixi.toml" --locked; then
    echo "Installation failed. The copied source remains at: $DESTINATION" >&2
    echo 'After resolving the error, retry inside that folder with: pixi install --locked' >&2
    exit 1
fi

printf '\nInstalled. Run independently from this folder:\n  cd %q\n' "$DESTINATION"
echo '  pixi run pytest'
echo '  pixi run python scripts/generate_mock_benchmark_data.py --n-clusters 20 --recall-rate 0.7 --n-false-positives 5'
echo '  pixi run python -m sharp.evaluate --predictions data/mock/predictions.parquet --ground-truth data/mock/ground_truth.tsv --output data/processed/benchmark.json'
