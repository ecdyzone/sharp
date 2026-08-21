#!/usr/bin/env bash
# Score one benchmark scope against the shared baseline output pool.
#
# This is step 4 of a benchmark run — the cheap, repeatable part. Steps 1-3
# (select the genome set, download it, submit the sbatch arrays) stay manual;
# they are slow, need the network, and the array is something you want to watch
# and resize yourself.
#
# Why this exists: under "run once, slice many" the tool output pool
# (~/projects/<tool>/out_benchmark/<ACCESSION>/) holds every genome ever run,
# and a benchmark is defined by a *pair* of files — a scope file and its
# matching contig-normalized ground truth. Pairing a scope with the wrong
# ground truth produces a plausible-looking benchmark.json scored against the
# wrong denominator, with nothing in the output to flag it. Deriving both from
# one scope name makes that mismatch structurally impossible.
#
# Usage:
#     scripts/run_benchmark.sh <scope> [tool ...] [-- evaluate-args ...]
#
#     # Both tools, default paths
#     scripts/run_benchmark.sh benchmark_set_strep
#
#     # One tool
#     scripts/run_benchmark.sh benchmark_set_strep antismash
#
#     # Sweep a score threshold (extra args go to sharp.evaluate). The parquet
#     # is reused, so a sweep costs seconds per point.
#     scripts/run_benchmark.sh benchmark_set_strep deepbgc -- --min-p-bgc 0.5
#
#     # See the resolved commands without running them
#     scripts/run_benchmark.sh benchmark_set_strep --dry-run
#
# Reads (per scope, both required — they must come from one
# select_benchmark_genomes.py --output-dir):
#     data/interim/<scope>/analyzed_contigs.txt
#     data/interim/<scope>/benchmark_ground_truth.tsv
#
# Reads (per tool, the shared pool):
#     $POOL_ROOT/<tool>/out_benchmark/          (default: ~/projects)
#
# Writes:
#     data/interim/<tool>_predictions_<scope>.parquet
#     data/processed/benchmark_<scope>_<tool>.json
#
# Options:
#     --force      overwrite an existing benchmark_<scope>_<tool>.json
#     --remerge    rebuild the parquet even if it is newer than the pool
#     --dry-run    print the commands and exit
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/_load_env.sh"

: "${POOL_ROOT:=$HOME/projects}"
DEFAULT_TOOLS=(antismash deepbgc)

FORCE=0
REMERGE=0
DRY_RUN=0
SCOPE=""
TOOLS=()
EVAL_ARGS=()

# ── argument parsing ────────────────────────────────────────────────────────
# First positional is the scope; any further positionals are tool names.
# Everything after `--` is forwarded verbatim to sharp.evaluate, which is how
# threshold sweeps work without this script knowing about evaluate's flags.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)    FORCE=1; shift ;;
        --remerge)  REMERGE=1; shift ;;
        --dry-run)  DRY_RUN=1; shift ;;
        -h|--help)  sed -n '2,47p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        --)         shift; EVAL_ARGS=("$@"); break ;;
        -*)         echo "ERROR: unknown option $1 (evaluate args go after --)" >&2; exit 2 ;;
        *)          if [[ -z "$SCOPE" ]]; then SCOPE="$1"; else TOOLS+=("$1"); fi; shift ;;
    esac
done

if [[ -z "$SCOPE" ]]; then
    echo "ERROR: no scope given." >&2
    echo "Usage: scripts/run_benchmark.sh <scope> [tool ...] [-- evaluate-args ...]" >&2
    echo "Available scopes:" >&2
    for d in "${REPO_ROOT}"/data/interim/benchmark_set*/; do
        [[ -d "$d" ]] && echo "  $(basename "$d")" >&2
    done
    exit 2
fi
[[ ${#TOOLS[@]} -eq 0 ]] && TOOLS=("${DEFAULT_TOOLS[@]}")

# Accept either a bare scope name or a path to the scope directory.
SCOPE="$(basename "$SCOPE")"
SCOPE_DIR="${REPO_ROOT}/data/interim/${SCOPE}"
CONTIGS="${SCOPE_DIR}/analyzed_contigs.txt"
GT="${SCOPE_DIR}/benchmark_ground_truth.tsv"

# ── preflight ───────────────────────────────────────────────────────────────
# Everything is checked before any work starts, so a bad invocation fails in a
# second rather than after a long merge.
fail=0
if [[ ! -d "$SCOPE_DIR" ]]; then
    echo "ERROR: no scope directory ${SCOPE_DIR}" >&2
    echo "       build one with: select_benchmark_genomes.py --output-dir data/interim/${SCOPE}" >&2
    exit 2
fi
for f in "$CONTIGS" "$GT"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty $f" >&2
        echo "       a scope needs BOTH analyzed_contigs.txt and benchmark_ground_truth.tsv," >&2
        echo "       from the same select_benchmark_genomes.py --output-dir" >&2
        fail=1
    fi
done
[[ $fail -eq 1 ]] && exit 2

N_CONTIGS=$(grep -cve '^[[:space:]]*$' "$CONTIGS" || true)
N_CLUSTERS=$(( $(wc -l < "$GT") - 1 ))
echo "scope ${SCOPE}: ${N_CONTIGS} contigs, ${N_CLUSTERS} ground-truth clusters"

mkdir -p "${REPO_ROOT}/data/interim" "${REPO_ROOT}/data/processed"

run() {
    if [[ $DRY_RUN -eq 1 ]]; then printf '  +'; printf ' %q' "$@"; printf '\n'; return 0; fi
    "$@"
}

# ── per tool ────────────────────────────────────────────────────────────────
failed=()
for tool in "${TOOLS[@]}"; do
    echo
    echo "── ${tool} ─────────────────────────────────────────────────────────"
    POOL="${POOL_ROOT}/${tool}/out_benchmark"
    PARQUET="${REPO_ROOT}/data/interim/${tool}_predictions_${SCOPE}.parquet"
    JSON="${REPO_ROOT}/data/processed/benchmark_${SCOPE}_${tool}.json"

    if [[ ! -d "$POOL" ]]; then
        echo "ERROR: no output pool at ${POOL}" >&2
        echo "       set POOL_ROOT, or run the array first:" >&2
        echo "       sbatch --array=1-${N_CONTIGS}%8 scripts/run_${tool}_array.sbatch ${CONTIGS}" >&2
        failed+=("$tool"); continue
    fi

    # A pool smaller than the scope usually means the array is still running.
    # Proceeding would score the missing genomes as "the tool found nothing",
    # which is indistinguishable in the metrics from a real miss.
    n_pool=$(find "$POOL" -mindepth 1 -maxdepth 1 -type d | wc -l)
    if [[ "$n_pool" -lt "$N_CONTIGS" ]]; then
        echo "WARNING: pool holds ${n_pool} genome dir(s) but the scope lists ${N_CONTIGS}." >&2
        echo "         If the array is still running, wait — genomes with no output are" >&2
        echo "         scored as if the tool found nothing there." >&2
    fi

    if [[ -f "$JSON" && $FORCE -eq 0 && $DRY_RUN -eq 0 ]]; then
        echo "ERROR: ${JSON} already exists (modified $(date -r "$JSON" '+%Y-%m-%d %H:%M'))." >&2
        echo "       re-run with --force to overwrite." >&2
        failed+=("$tool"); continue
    fi

    # Merging reparses every genome in the pool, so reuse a parquet that is
    # already newer than the newest thing in the pool. This is what makes a
    # --min-p-bgc sweep cheap: only evaluate re-runs.
    newest_pool="$(find "$POOL" -mindepth 1 -maxdepth 1 -type d -newer "$PARQUET" -print -quit 2>/dev/null || true)"
    if [[ -f "$PARQUET" && $REMERGE -eq 0 && -z "$newest_pool" ]]; then
        echo "reusing $(basename "$PARQUET") (newer than the pool; --remerge to rebuild)"
    else
        if ! run pixi run --manifest-path "${REPO_ROOT}/pixi.toml" \
                python "${SCRIPT_DIR}/merge_predictions.py" \
                --tool "$tool" \
                --input-dir "$POOL" \
                --contigs "$CONTIGS" \
                --output "$PARQUET"; then
            echo "ERROR: merge failed for ${tool}" >&2
            failed+=("$tool"); continue
        fi
    fi

    if ! run pixi run --manifest-path "${REPO_ROOT}/pixi.toml" \
            python -m sharp.evaluate \
            --predictions "$PARQUET" \
            --ground-truth "$GT" \
            --contigs "$CONTIGS" \
            --output "$JSON" \
            ${EVAL_ARGS[@]+"${EVAL_ARGS[@]}"}; then
        echo "ERROR: evaluate failed for ${tool}" >&2
        failed+=("$tool"); continue
    fi
done

# ── summary ─────────────────────────────────────────────────────────────────
# Read the numbers back so the result is visible without opening any file.
if [[ $DRY_RUN -eq 0 ]]; then
    echo
    echo "── ${SCOPE} ────────────────────────────────────────────────────────"
    for tool in "${TOOLS[@]}"; do
        JSON="${REPO_ROOT}/data/processed/benchmark_${SCOPE}_${tool}.json"
        [[ -f "$JSON" ]] || continue
        pixi run --manifest-path "${REPO_ROOT}/pixi.toml" python - "$JSON" "$tool" <<'PY' || true
import json, sys
d = json.load(open(sys.argv[1]))
det, rec, sc = d["detection"], d["reciprocal"], d["scope"]
print(f'  {sys.argv[2]:<10s} detection {det["recall"]:.3f} ({det["n_recovered"]}/{det["n_clusters"]})'
      f'   reciprocal {rec["recall"]:.3f}'
      f'   {sc["n_predictions_in_scope"]} predictions')
PY
    done
    echo "  → data/processed/benchmark_${SCOPE}_<tool>.json"
fi

if [[ ${#failed[@]} -gt 0 ]]; then
    echo >&2
    echo "FAILED: ${failed[*]}" >&2
    exit 1
fi
