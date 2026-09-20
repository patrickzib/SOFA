#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

usage() {
    cat <<'USAGE'
Usage: run_core_scaling_experiment.sh [OPTIONS] [RUN_SUITE_OPTIONS]

Runs the 16/32/64-query-core comparison for TRIE/SPARTAN, SOFA/SFA+iSAX,
and MESSI/SAX+iSAX. Index construction uses 64 workers by default.

Options:
  --index-threads N|auto   Index workers (default: 64)
  --core-list LIST         Query workers, comma-separated (default: 16,32,64)
  --experiment-root PATH   Separate logs/results root (default: ./results/core_scaling)
  -h, --help               Show this help

All other options are passed to run_suite.sh. In particular, use
--datasets to select a smaller dataset set and --dry-run to inspect commands.
USAGE
}

INDEX_THREADS=64
CORE_LIST=16,32,64
EXPERIMENT_ROOT=${MESSI_CORE_SCALING_ROOT:-"$PWD/results/core_scaling"}
PASSTHROUGH=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --index-threads)
            [[ $# -ge 2 ]] || { printf 'Error: --index-threads requires a value\n' >&2; exit 2; }
            INDEX_THREADS=$2
            shift 2
            ;;
        --core-list)
            [[ $# -ge 2 ]] || { printf 'Error: --core-list requires a value\n' >&2; exit 2; }
            CORE_LIST=$2
            shift 2
            ;;
        --experiment-root)
            [[ $# -ge 2 ]] || { printf 'Error: --experiment-root requires a path\n' >&2; exit 2; }
            EXPERIMENT_ROOT=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            PASSTHROUGH+=("$1")
            shift
            ;;
    esac
done

run_system() {
    local name=$1 index_type=$2 methods=$3
    shift 3
    local results_root="$EXPERIMENT_ROOT/results/$name"
    local log_root="$EXPERIMENT_ROOT/logs/$name"

    printf '\n=== %s: %s (%s) ===\n' "$name" "$index_type" "$methods" >&2
    MESSI_RESULTS_ROOT="$results_root" \
    MESSI_LOG_ROOT="$log_root" \
        "$SCRIPT_DIR/run_suite.sh" standard \
            --threads "$CORE_LIST" \
            --index-threads "$INDEX_THREADS" \
            --index-type "$index_type" \
            --methods "$methods" \
            --rerun-existing \
            "${PASSTHROUGH[@]}" "$@"
}

# Queue count follows each query-core count because run_suite defaults it to
# the current --threads value.
run_system trie trie "spartan-depth,spartan-width" \
    --trie-mbr-dims 128 --n-segments 64 --trie-split-dims 64 --trie-fanout 8 \
    --leaf-size 20000 --min-leaf-size 20000 --trie-leaf-ivf 16 \
    --trie-leaf-ivf-min-size 4096 --trie-record-mbr-suffix-bound \
    --trie-streaming-leaf-scan --trie-residual-record-only \
    --trie-residual-order symbolic-first --trie-leaf-ivf-radial-bound

run_system sofa isax "sfa-depth,sfa-width" --enable-sofa-v2
run_system messi isax "sax"
