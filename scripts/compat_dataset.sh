#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
[[ $# -ge 2 ]] || { printf 'Internal error: dataset and profile are required\n' >&2; exit 2; }
DATASET=$1
PROFILE=$2
shift 2

THREADS=${1:-}
QUEUE_NUMBER=${2:-}
ARGS=("$DATASET" "$PROFILE")
[[ -n $THREADS ]] && ARGS+=(--threads "$THREADS")
[[ -n $QUEUE_NUMBER ]] && ARGS+=(--queue-number "$QUEUE_NUMBER")

if [[ $DATASET == seisbench ]]; then
    [[ -n ${3:-} ]] && ARGS+=(--dataset-file "$3")
    [[ -n ${4:-} ]] && ARGS+=(--query-file "$4")
    [[ -n ${5:-} ]] && ARGS+=(--dataset-size "$5")
    case "$PROFILE" in
        knn) [[ -n ${6:-} ]] && ARGS+=(--k "$6") ;;
        sampling) [[ -n ${6:-} ]] && ARGS+=(--sample-factor "$6") ;;
    esac
else
    case "$PROFILE" in
        standard) [[ -n ${3:-} ]] && ARGS+=(--query-file "$3") ;;
        knn) [[ -n ${3:-} ]] && ARGS+=(--k "$3") ;;
        sampling) [[ -n ${3:-} ]] && ARGS+=(--sample-factor "$3") ;;
    esac
fi

exec "$SCRIPT_DIR/run_dataset.sh" "${ARGS[@]}"
