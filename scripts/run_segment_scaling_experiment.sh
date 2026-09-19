#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

usage() {
    cat <<'USAGE'
Usage: run_segment_scaling_experiment.sh [OPTIONS] [RUN_SUITE_OPTIONS]

Compares MESSI/SAX+iSAX, SOFA/SFA+iSAX, and SPARTAN/trie while increasing
the symbolic dimensions. Query and index construction use 64 workers by
default. SOFA-v2 bounds are deliberately excluded.

Options:
  --segment-list LIST       Symbolic dimensions (default: 16,32,64)
  --threads N               Query workers (default: 64)
  --index-threads N|auto    Index workers (default: 64)
  --experiment-root PATH    Logs/results root (default: ./results/segment_scaling)
  -h, --help                Show this help

All other options are passed to run_suite.sh. Use --datasets to select a
smaller dataset set and --dry-run to inspect commands.
USAGE
}

die() { printf 'Error: %s\n' "$*" >&2; exit 2; }

SEGMENT_LIST=16,32,64
QUERY_THREADS=64
INDEX_THREADS=64
EXPERIMENT_ROOT=${MESSI_SEGMENT_SCALING_ROOT:-"$PWD/results/segment_scaling"}
PASSTHROUGH=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --segment-list)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SEGMENT_LIST=$2
            shift 2
            ;;
        --threads)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            QUERY_THREADS=$2
            shift 2
            ;;
        --index-threads)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            INDEX_THREADS=$2
            shift 2
            ;;
        --experiment-root)
            [[ $# -ge 2 ]] || die "$1 requires a path"
            EXPERIMENT_ROOT=$2
            shift 2
            ;;
        --enable-sofa-v2|--isax-node-mbr|--isax-record-mbr-suffix-bound|--isax-record-lb-table|--no-simd)
            die "$1 is controlled by this experiment"
            ;;
        --isax-n-segments|--n-segments|--trie-record-lb-dims|--trie-split-dims|--trie-mbr-dims|--methods|--index-type)
            die "$1 is controlled by this experiment"
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

[[ $QUERY_THREADS =~ ^[1-9][0-9]*$ ]] || die '--threads must be a positive integer'
[[ $INDEX_THREADS == auto || $INDEX_THREADS =~ ^[1-9][0-9]*$ ]] || \
    die '--index-threads must be a positive integer or auto'

IFS=',' read -r -a SEGMENTS <<< "$SEGMENT_LIST"
[[ ${#SEGMENTS[@]} -gt 0 ]] || die '--segment-list must not be empty'
for segments in "${SEGMENTS[@]}"; do
    [[ $segments =~ ^[1-9][0-9]*$ ]] || die "invalid segment count '$segments'"
    (( segments == 16 || segments == 32 || segments == 64 )) || \
        die '--segment-list supports 16, 32, and 64'
done

run_system() {
    local segments=$1 name=$2 index_type=$3 methods=$4
    shift 4
    local results_root="$EXPERIMENT_ROOT/results/$name/segments-$segments"
    local log_root="$EXPERIMENT_ROOT/logs/$name/segments-$segments"

    printf '\n=== segments=%s system=%s layout=%s methods=%s ===\n' \
        "$segments" "$name" "$index_type" "$methods" >&2
    MESSI_RESULTS_ROOT="$results_root" \
    MESSI_LOG_ROOT="$log_root" \
        "$SCRIPT_DIR/run_suite.sh" standard \
            --threads "$QUERY_THREADS" \
            --index-threads "$INDEX_THREADS" \
            --index-type "$index_type" \
            --methods "$methods" \
            --rerun-existing \
            "${PASSTHROUGH[@]}" "$@"
}

for segments in "${SEGMENTS[@]}"; do
    run_system "$segments" messi isax sax \
        --isax-n-segments "$segments"

    run_system "$segments" sofa isax sfa-depth,sfa-width \
        --isax-n-segments "$segments"

    run_system "$segments" trie trie spartan-depth,spartan-width \
        --trie-mbr-dims 128 --n-segments "$segments" \
        --trie-split-dims "$segments" --trie-fanout 8 \
        --leaf-size 20000 --min-leaf-size 20000 --trie-leaf-ivf 16 \
        --trie-leaf-ivf-min-size 4096 --trie-record-mbr-suffix-bound \
        --trie-streaming-leaf-scan --trie-residual-record-only \
        --trie-residual-order symbolic-first --trie-leaf-ivf-radial-bound
done
