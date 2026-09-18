#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

usage() {
    cat <<'USAGE'
Usage: run_paper_pruning_experiment.sh [RUN_SUITE_OPTIONS]

Profiles the five lower bounds described in the S3-Trie paper on the eight
datasets used by the historical pruning-curve experiment. The executable must
have been configured with --enable-trie-pruning-trace.

Defaults:
  datasets: deep1b,sift1b,ethc,iquique,lendb,neic,obs,obst2024
  workers/queues: 64
  D=128 (capped by series length), B=64, split dimensions=64, fanout=8
  leaf capacity=20000, IVF=16 for leaves with at least 4096 series
  uniform sample=min(1000000,dataset size), 100 queries, z-normalized input

Additional arguments are passed to run_suite.sh after these defaults. Set
MESSI_BINARY, MESSI_DATA_ROOT, MESSI_QUERY_ROOT, MESSI_SEISBENCH_ROOT,
MESSI_SEISBENCH_QUERY_ROOT, MESSI_LOG_ROOT, and MESSI_RESULTS_ROOT as needed.
USAGE
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    usage
    exit 0
fi

exec "$SCRIPT_DIR/run_suite.sh" standard \
    --datasets deep1b,sift1b,ethc,iquique,lendb,neic,obs,obst2024 \
    --threads 64 --queue-number 64 --methods spartan-depth \
    --query-size 100 --sample-size 1000000 --sample-type 2 --sampling-seed 1 \
    --apply-z-norm --leaf-size 20000 --min-leaf-size 20000 \
    --trie-mbr-dims 128 --n-segments 64 --trie-split-dims 64 --trie-fanout 8 \
    --trie-leaf-ivf 16 --trie-leaf-ivf-min-size 4096 \
    --no-trie-leaf-ivf-raw-ball-bound --trie-leaf-ivf-radial-bound \
    --trie-record-mbr-suffix-bound --trie-streaming-leaf-scan \
    --trie-residual-record-only --trie-residual-order symbolic-first \
    --trie-pruning-curve "$@"
