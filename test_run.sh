#!/usr/bin/env bash
set -euo pipefail

FILE_PATH=data_head/astro_head.bin
QUERIES_PATH=data_queries/astro_queries.bin

TS_SIZE=256
COEFF_NUMBER=64
QUERY_SIZE=100

DATASET_SIZE=10000
SAMPLE_SIZE=10000

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  printf 'Usage: %s [threads] [isax|trie]\n' "$0"
  exit 0
fi

THREADS=${1:-auto}
INDEX_TYPE=${2:-trie}

case "$INDEX_TYPE" in
  isax|trie) ;;
  *)
    printf 'Usage: %s [threads] [isax|trie]\n' "$0" >&2
    exit 2
    ;;
esac

args=(
  --dataset "$FILE_PATH"
  --dataset-size "$DATASET_SIZE"
  --queries "$QUERIES_PATH"
  --queries-size "$QUERY_SIZE"
  --timeseries-size "$TS_SIZE"
  --sample-size "$SAMPLE_SIZE"
  --threads "$THREADS"
  --index-type "$INDEX_TYPE"
  --function-type 5
  --histogram-type 2
  --sfa-n-coefficients "$COEFF_NUMBER"
  --leaf-size 100
)

if [[ $INDEX_TYPE == isax ]]; then
  args+=(--dynamic-root-split-variance)
else
  # Exercise the combined per-record bound: the 16 symbolic dimensions plus
  # the disjoint leaf-MBR suffix over the remaining trie dimensions.
  args+=(--trie-record-mbr-suffix-bound)
fi

"${MESSI_BIN:-./bin/MESSI}" "${args[@]}"
