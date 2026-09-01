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
# ── Several people share one clone ──────────────────────────────────────────
# On the server the pool, the genomes and the scope files all live in one
# checkout, and everyone benchmarks inside it. Three consequences are handled
# here:
#
#   * Results are per-person. `benchmark_<scope>_<tool>.json` is keyed by scope
#     and tool alone, so two people scoring the same scope — or one person
#     sweeping --min-p-bgc — write the same path. Use `--label <name>` to claim
#     your own file. Overwriting someone else's still needs --force, and the
#     refusal names its owner.
#   * The merged parquet is shared on purpose. It is a derived cache keyed by
#     scope and tool, and reusing it is what makes a sweep cost seconds, so it
#     deliberately carries no label. It is written to a temporary path and
#     renamed into place, so a concurrent reader never sees a half-written file.
#   * Every result records who produced it. A `provenance` block is stamped
#     into the JSON (user, host, time, git commit, pool and its genome count,
#     evaluate args) — without it a shared directory of results is unreadable
#     a month later.
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
#     # Your own copy of a result, in a clone other people also benchmark in
#     scripts/run_benchmark.sh benchmark_set_strep --label alice
#
#     # Sweep a score threshold (extra args go to sharp.evaluate). The parquet
#     # is reused, so a sweep costs seconds per point. Label each point, or
#     # every threshold overwrites the same file.
#     for t in 0.0 0.3 0.5 0.7 0.9; do
#         scripts/run_benchmark.sh benchmark_set_strep deepbgc \
#             --label "p${t}" -- --min-p-bgc "$t"
#     done
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
#     data/interim/<tool>_predictions_<scope>.parquet      (shared cache)
#     data/processed/benchmark_<scope>_<tool>[.<label>].json
#
# Options:
#     --label NAME open results as benchmark_<scope>_<tool>.NAME.json, so two
#                  people (or two thresholds) do not overwrite each other
#     --force      overwrite an existing benchmark JSON, including someone else's
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
LABEL=""
TOOLS=()
EVAL_ARGS=()

# Print the header comment block as the help text. Scanning for the first
# non-comment line rather than a fixed line range keeps --help correct when the
# block above is edited.
usage() {
    awk 'NR == 1 { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' \
        "${BASH_SOURCE[0]}"
}

# ── argument parsing ────────────────────────────────────────────────────────
# First positional is the scope; any further positionals are tool names.
# Everything after `--` is forwarded verbatim to sharp.evaluate, which is how
# threshold sweeps work without this script knowing about evaluate's flags.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)    FORCE=1; shift ;;
        --remerge)  REMERGE=1; shift ;;
        --dry-run)  DRY_RUN=1; shift ;;
        --label)    LABEL="${2:-}"; shift 2 ;;
        --label=*)  LABEL="${1#--label=}"; shift ;;
        -h|--help)  usage; exit 0 ;;
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

# The label becomes part of a filename, so keep it to characters that survive a
# shell glob and a scp without quoting.
if [[ -n "$LABEL" && ! "$LABEL" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "ERROR: --label must match [A-Za-z0-9][A-Za-z0-9._-]* (got '${LABEL}')" >&2
    exit 2
fi

# Accept either a bare scope name or a path to the scope directory.
SCOPE="$(basename "$SCOPE")"
SCOPE_DIR="${REPO_ROOT}/data/interim/${SCOPE}"
CONTIGS="${SCOPE_DIR}/analyzed_contigs.txt"
GT="${SCOPE_DIR}/benchmark_ground_truth.tsv"

INTERIM_DIR="${REPO_ROOT}/data/interim"
PROCESSED_DIR="${REPO_ROOT}/data/processed"

# benchmark_<scope>_<tool>.json, or ..._<tool>.<label>.json when labelled. The
# dot keeps the label visually separate from the underscore-joined parts, since
# scope names contain underscores themselves.
json_path() {
    local tool="$1"
    if [[ -n "$LABEL" ]]; then
        printf '%s/benchmark_%s_%s.%s.json' "$PROCESSED_DIR" "$SCOPE" "$tool" "$LABEL"
    else
        printf '%s/benchmark_%s_%s.json' "$PROCESSED_DIR" "$SCOPE" "$tool"
    fi
}

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

# Both output directories are written to, and in a shared clone they are owned
# by whoever ran the first benchmark. Check now: failing here costs a second,
# failing after the merge costs minutes.
for d in "$INTERIM_DIR" "$PROCESSED_DIR"; do
    if [[ ! -d "$d" ]] && ! mkdir -p "$d" 2>/dev/null; then
        echo "ERROR: cannot create ${d}" >&2
        fail=1
        continue
    fi
    if [[ ! -w "$d" ]]; then
        echo "ERROR: ${d} is not writable by $(id -un)." >&2
        echo "       Everyone benchmarks inside this clone, so both output dirs must be" >&2
        echo "       group-writable. Its owner should run, once:" >&2
        echo "         chmod g+rwxs ${d}     # setgid: new files keep the group" >&2
        echo "       and each person should set 'umask 002' so their own results stay" >&2
        echo "       readable and writable by the rest of the group." >&2
        fail=1
    fi
done
[[ $fail -eq 1 ]] && exit 2

N_CONTIGS=$(grep -cve '^[[:space:]]*$' "$CONTIGS" || true)
N_CLUSTERS=$(( $(wc -l < "$GT") - 1 ))
echo "scope ${SCOPE}: ${N_CONTIGS} contigs, ${N_CLUSTERS} ground-truth clusters"
[[ -n "$LABEL" ]] && echo "label: ${LABEL}"

# Recorded in every result so a JSON found later can be traced to the code that
# produced it. A dirty tree is flagged rather than silently reported as clean.
GIT_COMMIT="$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
if ! git -C "$REPO_ROOT" diff --quiet HEAD -- 2>/dev/null; then
    GIT_COMMIT="${GIT_COMMIT}-dirty"
fi

# Merges write to <parquet>.tmp.<pid> and rename, so an interrupted run leaves
# nothing behind for the next person to trip over.
cleanup() { rm -f "${INTERIM_DIR}"/*.tmp."$$" 2>/dev/null || true; }
trap cleanup EXIT

run() {
    if [[ $DRY_RUN -eq 1 ]]; then printf '  +'; printf ' %q' "$@"; printf '\n'; return 0; fi
    "$@"
}

# Append a provenance block to a finished benchmark JSON. Kept here rather than
# in sharp.evaluate because every field is a wrapper concept: evaluate.py knows
# nothing about the pool, the label, or the shared clone.
stamp_provenance() {
    [[ $DRY_RUN -eq 1 ]] && return 0
    SHARP_PROV_JSON="$1" \
    SHARP_PROV_TOOL="$2" \
    SHARP_PROV_POOL="$3" \
    SHARP_PROV_NPOOL="$4" \
    SHARP_PROV_PARQUET="$5" \
    SHARP_PROV_REUSED="$6" \
    SHARP_PROV_SCOPE="$SCOPE" \
    SHARP_PROV_LABEL="$LABEL" \
    SHARP_PROV_COMMIT="$GIT_COMMIT" \
    SHARP_PROV_ARGS="${EVAL_ARGS[*]:-}" \
    pixi run --manifest-path "${REPO_ROOT}/pixi.toml" python - <<'PY'
import datetime, getpass, json, os, pathlib, socket

path = pathlib.Path(os.environ["SHARP_PROV_JSON"])
doc = json.loads(path.read_text())
doc["provenance"] = {
    "user": getpass.getuser(),
    "host": socket.gethostname(),
    "written_at": datetime.datetime.now(datetime.timezone.utc)
    .replace(microsecond=0)
    .isoformat(),
    "scope": os.environ["SHARP_PROV_SCOPE"],
    "label": os.environ["SHARP_PROV_LABEL"] or None,
    "tool": os.environ["SHARP_PROV_TOOL"],
    "pool": os.environ["SHARP_PROV_POOL"],
    "pool_n_genomes": int(os.environ["SHARP_PROV_NPOOL"]),
    "predictions": os.environ["SHARP_PROV_PARQUET"],
    "predictions_reused": os.environ["SHARP_PROV_REUSED"] == "1",
    "git_commit": os.environ["SHARP_PROV_COMMIT"],
    "evaluate_args": os.environ["SHARP_PROV_ARGS"],
}
path.write_text(json.dumps(doc, indent=2) + "\n")
PY
}

# ── per tool ────────────────────────────────────────────────────────────────
failed=()
succeeded=()
for tool in "${TOOLS[@]}"; do
    echo
    echo "── ${tool} ─────────────────────────────────────────────────────────"
    POOL="${POOL_ROOT}/${tool}/out_benchmark"
    PARQUET="${INTERIM_DIR}/${tool}_predictions_${SCOPE}.parquet"
    TMP_PARQUET="${PARQUET}.tmp.$$"
    JSON="$(json_path "$tool")"

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

    # In a shared clone the existing file may well be someone else's, so name
    # its owner and offer --label before offering --force.
    if [[ -f "$JSON" && $FORCE -eq 0 && $DRY_RUN -eq 0 ]]; then
        owner="$(stat -c '%U' "$JSON" 2>/dev/null || echo 'unknown')"
        echo "ERROR: $(basename "$JSON") already exists" >&2
        echo "       (owner ${owner}, modified $(date -r "$JSON" '+%Y-%m-%d %H:%M'))." >&2
        if [[ "$owner" != "$(id -un)" ]]; then
            echo "       That is not your file. Prefer --label <name> to write your own" >&2
            echo "       copy; --force overwrites theirs." >&2
        else
            echo "       re-run with --force to overwrite, or --label <name> to keep both." >&2
        fi
        failed+=("$tool"); continue
    fi

    # Merging reparses every genome in the pool, so reuse a parquet that is
    # already newer than the newest thing in the pool. This is what makes a
    # --min-p-bgc sweep cheap: only evaluate re-runs. The parquet is shared
    # across everyone benchmarking this scope, hence no label in its name.
    newest_pool="$(find "$POOL" -mindepth 1 -maxdepth 1 -type d -newer "$PARQUET" -print -quit 2>/dev/null || true)"
    if [[ -f "$PARQUET" && $REMERGE -eq 0 && -z "$newest_pool" ]]; then
        echo "reusing $(basename "$PARQUET") (newer than the pool; --remerge to rebuild)"
        reused=1
    else
        reused=0
        # Merge to a temporary path and rename, so a concurrent reader either
        # sees the old parquet or the new one, never a partial write.
        if ! run pixi run --manifest-path "${REPO_ROOT}/pixi.toml" \
                python "${SCRIPT_DIR}/merge_predictions.py" \
                --tool "$tool" \
                --input-dir "$POOL" \
                --contigs "$CONTIGS" \
                --output "$TMP_PARQUET"; then
            rm -f "$TMP_PARQUET"
            echo "ERROR: merge failed for ${tool}" >&2
            failed+=("$tool"); continue
        fi
        if ! run mv -f "$TMP_PARQUET" "$PARQUET"; then
            rm -f "$TMP_PARQUET"
            echo "ERROR: could not install $(basename "$PARQUET")" >&2
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

    # Non-fatal: the numbers are already on disk and correct, only the audit
    # trail is missing.
    if ! stamp_provenance "$JSON" "$tool" "$POOL" "$n_pool" "$PARQUET" "$reused"; then
        echo "WARNING: could not stamp provenance into $(basename "$JSON")" >&2
    fi

    succeeded+=("$tool")
done

# ── summary ─────────────────────────────────────────────────────────────────
# Read the numbers back so the result is visible without opening any file.
# Only tools this run actually scored are listed: a tool that was refused
# because someone else's result is already there still has a readable JSON on
# disk, and printing it here would present their numbers as yours.
if [[ $DRY_RUN -eq 0 && ${#succeeded[@]} -gt 0 ]]; then
    echo
    echo "── ${SCOPE} ────────────────────────────────────────────────────────"
    for tool in "${succeeded[@]}"; do
        JSON="$(json_path "$tool")"
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
    if [[ -n "$LABEL" ]]; then
        echo "  → data/processed/benchmark_${SCOPE}_<tool>.${LABEL}.json"
    else
        echo "  → data/processed/benchmark_${SCOPE}_<tool>.json"
    fi
fi

if [[ ${#failed[@]} -gt 0 ]]; then
    echo >&2
    echo "FAILED: ${failed[*]}" >&2
    exit 1
fi
