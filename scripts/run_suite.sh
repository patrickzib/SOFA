#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=lib/datasets.sh
source "$SCRIPT_DIR/lib/datasets.sh"

usage() {
    cat <<'USAGE'
Usage: run_suite.sh SUITE [OPTIONS]

Suites:
  standard, high-frequency, knn, sampling
  generated-queries, hard-queries, noise-workloads

Options:
  --threads LIST          Comma-separated CPU/queue counts (default: available physical CPU cores)
  --k-values LIST         K values for knn (default: 20,50)
  --sample-factors LIST   Factors for sampling (default: 0.15,...,0.5)
  --datasets LIST         Limit regular suites to dataset IDs
  --methods LIST          Comma-separated methods to run
  --index-type TYPE       Index layout: trie (default) or isax
  --enable-sofa-v2        Enable all opt-in iSAX SOFA v2 bounds and variance root splitting
  --isax-node-mbr         Enable iSAX node MBR bounds
  --isax-record-mbr-suffix-bound
                          Enable learned iSAX record-MBR suffix bounds
  --isax-mbr-dims N       Extended iSAX MBR dimensions (default: 32)
  --isax-record-lb-table  Use query-local iSAX record lower-bound tables
  --queue-number N        Fixed queue count (default: same as each --threads value)
  --numa MODE             auto (default), none, or number of NUMA nodes
  --dataset-file PATH     Override the dataset path for every selected run
  --query-file PATH       Override the query path for regular-suite runs
  --dataset-size N        Override dataset records; accepts count suffixes
  --query-size N          Override queries per run; accepts count suffixes
  --leaf-size N           Maximum records per leaf; accepts count suffixes
  --min-leaf-size N       iSAX query-leaf threshold; trie does not enforce it
  --sample-size N         Override binning sample size; accepts count suffixes
  --sample-type 1|2|3     Binning sampling: first values, uniform (default), or random
  --sampling-seed N       Sampling seed (default: 1)
  --no-simd               Disable SIMD for every run
  --trie-mbr-dims N       Trie MBR dimensions (default: 128; capped by series length)
  --n-segments N          Trie record-prefix lower-bound dimensions (default: 64; range: 16--64)
  --trie-split-dims N     Trie split-candidate dimensions (default: min(32, MBR dimensions))
  --trie-record-mbr-suffix-bound
                          Add leaf-MBR contributions outside record-prefix dimensions (default for trie)
  --no-trie-record-mbr-suffix-bound
                          Disable trie record-MBR suffix pruning
  --trie-streaming-leaf-scan
                          Refine each passing record immediately (default)
  --no-trie-streaming-leaf-scan
                          Use the record lower-bound heap instead
  --trie-leaf-ivf K        Build K flat IVF MBR groups inside large trie leaves (default: 16)
  --no-trie-leaf-ivf       Disable flat leaf IVF groups
  --no-trie-leaf-ivf-raw-ball-bound
                          Disable certified raw centroid/radius cluster pruning
  --trie-leaf-ivf-radial-bound
                          Prune IVF records by centroid radius without reordering
  --trie-fanout 2|4|8      Trie symbolic split fanout (default: 8)
  --trie-dynamic-alphabet Use one global variance-weighted alphabet allocation
  --trie-min-fanout N     Minimum dynamic trie fanout (default: 2)
  --trie-max-fanout N     Maximum dynamic trie fanout (default: 16)
  --trie-alphabet-budget-bits N
                          Average dynamic alphabet budget in bits (default: 3)
  --trie-query-parallel   Parallelize each trie query
  --trie-query-batch      Batch independent trie queries
  --query-report-interval N
                          Print query progress every N queries (0 disables it)
  --profile-query-phases  Measure traversal, lower-bound, and exact work
  --dynamic-root-split-variance
                          Enable variance-assigned iSAX root bits
  --no-dynamic-root-split-variance
                          Disable variance-assigned root bits for learned iSAX methods
  --tight-bound           Enable iSAX tight-bound pruning (default for iSAX)
  --binary PATH           MESSI executable
  --data-root PATH        Main dataset root
  --query-root PATH       Main query root
  --seisbench-root PATH   SeisBench dataset root
  --seisbench-query-root PATH
                          SeisBench query root
  --rerun-existing        Run workloads even when an archive directory exists
  --dry-run               Print commands without running or archiving
  -h, --help              Show this help

Path-related environment variables and runner options are documented in
scripts/README.md.
USAGE
}

die() { printf 'Error: %s\n' "$*" >&2; exit 2; }
split_csv() { local value=$1; IFS=',' read -r -a SPLIT_RESULT <<< "$value"; }
dataset_label() {
    case "$1" in
        bigann) printf 'BIGANN' ;; sald) printf 'SALD' ;; sift1b) printf 'SIFT1b' ;;
        deep1b) printf 'DEEP1b' ;; scedc) printf 'SCEDC' ;; astro) printf 'ASTRO' ;;
        ethc) printf 'ETHC' ;; isc_ehb_depthphases) printf 'ISC_EHB_DepthPhases' ;;
        lendb) printf 'LenDB' ;; iquique) printf 'Iquique' ;; neic) printf 'NEIC' ;;
        obs) printf 'OBS' ;; obst2024) printf 'OBST2024' ;; pnw) printf 'PNW' ;;
        meier2019jgr) printf 'Meier2019JGR' ;; stead) printf 'STEAD' ;; txed) printf 'TXED' ;;
        *) printf '%s' "$1" ;;
    esac
}

[[ $# -ge 1 ]] || { usage >&2; exit 2; }
SUITE=$1
shift

THREADS_CSV=$(physical_core_count) || die 'unable to detect physical CPU cores; pass --threads N'
K_VALUES_CSV=20,50
SAMPLE_FACTORS_CSV=0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5
DATASETS_CSV=
METHODS_OVERRIDE=
INDEX_TYPE=trie
QUEUE_NUMBER=
NUMA_MODE=auto
DATASET_FILE=
QUERY_FILE=
DATASET_SIZE=
QUERY_SIZE=
LEAF_SIZE=
MIN_LEAF_SIZE=
SAMPLE_SIZE=
SAMPLE_TYPE=2
SAMPLING_SEED=
NO_SIMD=false
TRIE_FANOUT=8
TRIE_MBR_DIMS=
TRIE_RECORD_LB_DIMS=64
TRIE_RECORD_LB_DIMS_SPECIFIED=false
TRIE_SPLIT_DIMS=
TRIE_RECORD_MBR_SUFFIX_BOUND=
TRIE_STREAMING_LEAF_SCAN=true
TRIE_STREAMING_LEAF_SCAN_SPECIFIED=false
TRIE_LEAF_IVF=16
TRIE_LEAF_IVF_SPECIFIED=false
TRIE_LEAF_IVF_RAW_BALL_BOUND=true
TRIE_LEAF_IVF_RADIAL_BOUND=false
TRIE_DYNAMIC_ALPHABET=false
TRIE_MIN_FANOUT=2
TRIE_MAX_FANOUT=16
TRIE_ALPHABET_BUDGET_BITS=3
TRIE_QUERY_PARALLEL=false
TRIE_QUERY_BATCH=false
QUERY_REPORT_INTERVAL=
PROFILE_QUERY_PHASES=false
DYNAMIC_ROOT_SPLIT_VARIANCE=
DYNAMIC_ROOT_SPLIT_VARIANCE_SPECIFIED=false
ENABLE_SOFA_V2=false
ISAX_NODE_MBR=false
ISAX_RECORD_MBR_SUFFIX_BOUND=false
ISAX_MBR_DIMS=
ISAX_RECORD_LB_TABLE=false
TIGHT_BOUND=true
MESSI_EXECUTABLE=
DATA_ROOT=
QUERY_ROOT=
SEISBENCH_ROOT=
SEISBENCH_QUERY_ROOT=
RERUN_EXISTING=false
DRY_RUN=false
RESULTS_ROOT=${MESSI_RESULTS_ROOT:-"$HOME/MESSI_SFA_logs"}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) [[ $# -ge 2 ]] || die "$1 requires a value"; THREADS_CSV=$2; shift 2 ;;
        --k-values) [[ $# -ge 2 ]] || die "$1 requires a value"; K_VALUES_CSV=$2; shift 2 ;;
        --sample-factors) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_FACTORS_CSV=$2; shift 2 ;;
        --datasets) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASETS_CSV=$2; shift 2 ;;
        --methods) [[ $# -ge 2 ]] || die "$1 requires a value"; METHODS_OVERRIDE=$2; shift 2 ;;
        --index-type) [[ $# -ge 2 ]] || die "$1 requires a value"; INDEX_TYPE=$2; shift 2 ;;
        --enable-sofa-v2) ENABLE_SOFA_V2=true; shift ;;
        --isax-node-mbr) ISAX_NODE_MBR=true; shift ;;
        --isax-record-mbr-suffix-bound) ISAX_RECORD_MBR_SUFFIX_BOUND=true; shift ;;
        --isax-mbr-dims|--isax-mbr-dimensions) [[ $# -ge 2 ]] || die "$1 requires a value"; ISAX_MBR_DIMS=$2; shift 2 ;;
        --isax-record-lb-table) ISAX_RECORD_LB_TABLE=true; shift ;;
        --queue-number) [[ $# -ge 2 ]] || die "$1 requires a value"; QUEUE_NUMBER=$2; shift 2 ;;
        --numa) [[ $# -ge 2 ]] || die "$1 requires a value"; NUMA_MODE=$2; shift 2 ;;
        --dataset-file) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_FILE=$2; shift 2 ;;
        --query-file) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_FILE=$2; shift 2 ;;
        --dataset-size) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_SIZE=$2; shift 2 ;;
        --query-size) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_SIZE=$2; shift 2 ;;
        --leaf-size) [[ $# -ge 2 ]] || die "$1 requires a value"; LEAF_SIZE=$2; shift 2 ;;
        --min-leaf-size) [[ $# -ge 2 ]] || die "$1 requires a value"; MIN_LEAF_SIZE=$2; shift 2 ;;
        --sample-size) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_SIZE=$2; shift 2 ;;
        --sample-type) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_TYPE=$2; shift 2 ;;
        --sampling-seed) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLING_SEED=$2; shift 2 ;;
        --no-simd) NO_SIMD=true; shift ;;
        --trie-mbr-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MBR_DIMS=$2; shift 2 ;;
        --n-segments|--trie-record-lb-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_RECORD_LB_DIMS=$2; TRIE_RECORD_LB_DIMS_SPECIFIED=true; shift 2 ;;
        --trie-split-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_SPLIT_DIMS=$2; shift 2 ;;
        --trie-record-mbr-suffix-bound) TRIE_RECORD_MBR_SUFFIX_BOUND=true; shift ;;
        --no-trie-record-mbr-suffix-bound) TRIE_RECORD_MBR_SUFFIX_BOUND=false; shift ;;
        --trie-streaming-leaf-scan) TRIE_STREAMING_LEAF_SCAN=true; TRIE_STREAMING_LEAF_SCAN_SPECIFIED=true; shift ;;
        --no-trie-streaming-leaf-scan) TRIE_STREAMING_LEAF_SCAN=false; TRIE_STREAMING_LEAF_SCAN_SPECIFIED=true; shift ;;
        --trie-leaf-ivf) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_LEAF_IVF=$2; TRIE_LEAF_IVF_SPECIFIED=true; shift 2 ;;
        --no-trie-leaf-ivf) TRIE_LEAF_IVF=0; TRIE_LEAF_IVF_SPECIFIED=true; shift ;;
        --no-trie-leaf-ivf-raw-ball-bound) TRIE_LEAF_IVF_RAW_BALL_BOUND=false; shift ;;
        --trie-leaf-ivf-radial-bound) TRIE_LEAF_IVF_RADIAL_BOUND=true; shift ;;
        --trie-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_FANOUT=$2; shift 2 ;;
        --trie-dynamic-alphabet) TRIE_DYNAMIC_ALPHABET=true; shift ;;
        --trie-min-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MIN_FANOUT=$2; shift 2 ;;
        --trie-max-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MAX_FANOUT=$2; shift 2 ;;
        --trie-alphabet-budget-bits) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_ALPHABET_BUDGET_BITS=$2; shift 2 ;;
        --trie-query-parallel) TRIE_QUERY_PARALLEL=true; shift ;;
        --trie-query-batch) TRIE_QUERY_BATCH=true; shift ;;
        --query-report-interval) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_REPORT_INTERVAL=$2; shift 2 ;;
        --profile-query-phases) PROFILE_QUERY_PHASES=true; shift ;;
        --dynamic-root-split-variance) DYNAMIC_ROOT_SPLIT_VARIANCE=true; DYNAMIC_ROOT_SPLIT_VARIANCE_SPECIFIED=true; shift ;;
        --no-dynamic-root-split-variance) DYNAMIC_ROOT_SPLIT_VARIANCE=false; DYNAMIC_ROOT_SPLIT_VARIANCE_SPECIFIED=true; shift ;;
        --tight-bound) TIGHT_BOUND=true; shift ;;
        --no-tight-bound) TIGHT_BOUND=false; shift ;;
        --binary) [[ $# -ge 2 ]] || die "$1 requires a value"; MESSI_EXECUTABLE=$2; shift 2 ;;
        --data-root) [[ $# -ge 2 ]] || die "$1 requires a value"; DATA_ROOT=$2; shift 2 ;;
        --query-root) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_ROOT=$2; shift 2 ;;
        --seisbench-root) [[ $# -ge 2 ]] || die "$1 requires a value"; SEISBENCH_ROOT=$2; shift 2 ;;
        --seisbench-query-root) [[ $# -ge 2 ]] || die "$1 requires a value"; SEISBENCH_QUERY_ROOT=$2; shift 2 ;;
        --rerun-existing) RERUN_EXISTING=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option '$1'" ;;
    esac
done

case "$SUITE" in
    standard|high-frequency|knn|sampling|generated-queries|hard-queries|noise-workloads) ;;
    *) die "unknown suite '$SUITE'" ;;
esac
[[ $INDEX_TYPE == isax || $INDEX_TYPE == trie ]] || die '--index-type must be isax or trie'
[[ $SAMPLE_TYPE == 1 || $SAMPLE_TYPE == 2 || $SAMPLE_TYPE == 3 ]] || \
    die '--sample-type must be 1 (first values), 2 (uniform), or 3 (random)'
[[ $TRIE_FANOUT == 2 || $TRIE_FANOUT == 4 || $TRIE_FANOUT == 8 ]] || \
    die '--trie-fanout must be 2, 4, or 8'
[[ $INDEX_TYPE == trie || $TRIE_FANOUT == 8 ]] || \
    die '--trie-fanout requires --index-type trie'
[[ -z $TRIE_MBR_DIMS || $INDEX_TYPE == trie ]] || \
    die '--trie-mbr-dims requires --index-type trie'
[[ $TRIE_RECORD_LB_DIMS_SPECIFIED == false || $INDEX_TYPE == trie ]] || \
    die '--n-segments requires --index-type trie'
[[ -z $TRIE_RECORD_MBR_SUFFIX_BOUND || $INDEX_TYPE == trie ]] || \
    die '--trie-record-mbr-suffix-bound requires --index-type trie'
[[ $TRIE_STREAMING_LEAF_SCAN_SPECIFIED == false || $INDEX_TYPE == trie ]] || \
    die 'trie streaming leaf-scan options require --index-type trie'
[[ $TRIE_QUERY_PARALLEL == false || $INDEX_TYPE == trie ]] || \
    die '--trie-query-parallel requires --index-type trie'
[[ $TRIE_QUERY_BATCH == false || $INDEX_TYPE == trie ]] || \
    die '--trie-query-batch requires --index-type trie'
[[ $TRIE_QUERY_PARALLEL == false || $TRIE_QUERY_BATCH == false ]] || \
    die 'choose at most one of --trie-query-parallel and --trie-query-batch'
[[ $DYNAMIC_ROOT_SPLIT_VARIANCE != true || $INDEX_TYPE == isax ]] || \
    die '--dynamic-root-split-variance requires --index-type isax'
[[ $TRIE_LEAF_IVF_SPECIFIED == false || $INDEX_TYPE == trie ]] || die '--trie-leaf-ivf requires --index-type trie'
[[ $TRIE_LEAF_IVF_RAW_BALL_BOUND == true || $INDEX_TYPE == trie ]] || \
    die '--no-trie-leaf-ivf-raw-ball-bound requires --index-type trie'
[[ $TRIE_LEAF_IVF_RADIAL_BOUND == false || $INDEX_TYPE == trie ]] || \
    die '--trie-leaf-ivf-radial-bound requires --index-type trie'
[[ $TRIE_DYNAMIC_ALPHABET == false || $INDEX_TYPE == trie ]] || \
    die '--trie-dynamic-alphabet requires --index-type trie'
if [[ $INDEX_TYPE == trie ]]; then
    TRIE_RECORD_MBR_SUFFIX_BOUND=${TRIE_RECORD_MBR_SUFFIX_BOUND:-true}
    [[ $TRIE_LEAF_IVF_RADIAL_BOUND == false || $TRIE_LEAF_IVF != 0 ]] || \
        die '--trie-leaf-ivf-radial-bound requires --trie-leaf-ivf'
fi
if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
    [[ $TRIE_FANOUT == 8 ]] || die '--trie-fanout cannot be combined with --trie-dynamic-alphabet'
    [[ $TRIE_MIN_FANOUT =~ ^(2|4|8|16|32|64|128|256)$ ]] || die '--trie-min-fanout must be a power of two between 2 and 256'
    [[ $TRIE_MAX_FANOUT =~ ^(2|4|8|16|32|64|128|256)$ ]] || die '--trie-max-fanout must be a power of two between 2 and 256'
    [[ $TRIE_ALPHABET_BUDGET_BITS =~ ^[1-8]$ ]] || die '--trie-alphabet-budget-bits must be between 1 and 8'
fi

split_csv "$THREADS_CSV"; THREADS=("${SPLIT_RESULT[@]}")
split_csv "$K_VALUES_CSV"; K_VALUES=("${SPLIT_RESULT[@]}")
split_csv "$SAMPLE_FACTORS_CSV"; SAMPLE_FACTORS=("${SPLIT_RESULT[@]}")

archive_exists() {
    local label=$1
    [[ -d "$RESULTS_ROOT/$label" ]]
}

if [[ -n $DATASETS_CSV ]]; then
    split_csv "$DATASETS_CSV"; DATASETS=("${SPLIT_RESULT[@]}")
else
    mapfile -t DATASETS < <({ active_datasets; seisbench_datasets; })
fi

run_one() {
    local dataset=$1 profile=$2 threads=$3 result_value=$4
    shift 4
    if [[ $DRY_RUN == false && $RERUN_EXISTING == false && $profile != high-frequency ]] && \
       archive_exists "$(dataset_label "$dataset")"; then
        printf 'Skipping dataset=%s profile=%s run=%s: archive directory already exists at %s/%s\n' \
            "$dataset" "$profile" "$result_value" "$RESULTS_ROOT" \
            "$(dataset_label "$dataset")" >&2
        return 0
    fi
    local queue_number=${QUEUE_NUMBER:-$threads}
    local -a command=("$SCRIPT_DIR/run_dataset.sh" "$dataset" "$profile" --threads "$threads" \
        --queue-number "$queue_number" --numa "$NUMA_MODE" --index-type "$INDEX_TYPE")
    [[ -n $DATASET_FILE ]] && command+=(--dataset-file "$DATASET_FILE")
    [[ -n $QUERY_FILE ]] && command+=(--query-file "$QUERY_FILE")
    [[ -n $DATASET_SIZE ]] && command+=(--dataset-size "$DATASET_SIZE")
    [[ -n $QUERY_SIZE ]] && command+=(--query-size "$QUERY_SIZE")
    [[ -n $SAMPLE_SIZE ]] && command+=(--sample-size "$SAMPLE_SIZE")
    command+=(--sample-type "$SAMPLE_TYPE")
    [[ -n $SAMPLING_SEED ]] && command+=(--sampling-seed "$SAMPLING_SEED")
    [[ -n $MESSI_EXECUTABLE ]] && command+=(--binary "$MESSI_EXECUTABLE")
    [[ -n $DATA_ROOT ]] && command+=(--data-root "$DATA_ROOT")
    [[ -n $QUERY_ROOT ]] && command+=(--query-root "$QUERY_ROOT")
    [[ -n $SEISBENCH_ROOT ]] && command+=(--seisbench-root "$SEISBENCH_ROOT")
    [[ -n $SEISBENCH_QUERY_ROOT ]] && command+=(--seisbench-query-root "$SEISBENCH_QUERY_ROOT")
    $NO_SIMD && command+=(--no-simd)
    if [[ $INDEX_TYPE == trie ]]; then
        [[ -n $TRIE_MBR_DIMS ]] && command+=(--trie-mbr-dims "$TRIE_MBR_DIMS")
        command+=(--n-segments "$TRIE_RECORD_LB_DIMS")
        [[ -n $TRIE_SPLIT_DIMS ]] && command+=(--trie-split-dims "$TRIE_SPLIT_DIMS")
        $TRIE_RECORD_MBR_SUFFIX_BOUND && command+=(--trie-record-mbr-suffix-bound)
        [[ $TRIE_RECORD_MBR_SUFFIX_BOUND == false ]] && command+=(--no-trie-record-mbr-suffix-bound)
        $TRIE_STREAMING_LEAF_SCAN && command+=(--trie-streaming-leaf-scan)
        [[ $TRIE_STREAMING_LEAF_SCAN == false ]] && command+=(--no-trie-streaming-leaf-scan)
        [[ $TRIE_LEAF_IVF != 0 ]] && command+=(--trie-leaf-ivf "$TRIE_LEAF_IVF")
        [[ $TRIE_LEAF_IVF == 0 ]] && command+=(--no-trie-leaf-ivf)
        [[ $TRIE_LEAF_IVF_RAW_BALL_BOUND == false ]] && command+=(--no-trie-leaf-ivf-raw-ball-bound)
        [[ $TRIE_LEAF_IVF_RADIAL_BOUND == true ]] && command+=(--trie-leaf-ivf-radial-bound)
        if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
            command+=(--trie-dynamic-alphabet --trie-min-fanout "$TRIE_MIN_FANOUT"
                      --trie-max-fanout "$TRIE_MAX_FANOUT"
                      --trie-alphabet-budget-bits "$TRIE_ALPHABET_BUDGET_BITS")
        else
            command+=(--trie-fanout "$TRIE_FANOUT")
        fi
    fi
    [[ -n $LEAF_SIZE ]] && command+=(--leaf-size "$LEAF_SIZE")
    [[ -n $MIN_LEAF_SIZE ]] && command+=(--min-leaf-size "$MIN_LEAF_SIZE")
    $TRIE_QUERY_PARALLEL && command+=(--trie-query-parallel)
    $TRIE_QUERY_BATCH && command+=(--trie-query-batch)
    [[ -n $QUERY_REPORT_INTERVAL ]] && command+=(--query-report-interval "$QUERY_REPORT_INTERVAL")
    $PROFILE_QUERY_PHASES && command+=(--profile-query-phases)
    [[ $DYNAMIC_ROOT_SPLIT_VARIANCE == true ]] && command+=(--dynamic-root-split-variance)
    [[ $DYNAMIC_ROOT_SPLIT_VARIANCE_SPECIFIED == true && $DYNAMIC_ROOT_SPLIT_VARIANCE == false ]] && command+=(--no-dynamic-root-split-variance)
    if [[ $INDEX_TYPE == isax ]]; then
        $ENABLE_SOFA_V2 && command+=(--enable-sofa-v2)
        $ISAX_NODE_MBR && command+=(--isax-node-mbr)
        $ISAX_RECORD_MBR_SUFFIX_BOUND && command+=(--isax-record-mbr-suffix-bound)
        [[ -n $ISAX_MBR_DIMS ]] && command+=(--isax-mbr-dimensions "$ISAX_MBR_DIMS")
        $ISAX_RECORD_LB_TABLE && command+=(--isax-record-lb-table)
        $TIGHT_BOUND && command+=(--tight-bound)
        [[ $TIGHT_BOUND == false ]] && command+=(--no-tight-bound)
    fi
    [[ -n $METHODS_OVERRIDE ]] && command+=(--methods "$METHODS_OVERRIDE")
    command+=("$@")
    $DRY_RUN && command+=(--dry-run)
    "${command[@]}"
    if [[ $DRY_RUN == false && $profile != high-frequency ]]; then
        "$SCRIPT_DIR/archive_results.sh" "$(dataset_label "$dataset")" "$result_value"
    fi
}

run_regular_suite() {
    local profile=$1 threads dataset value
    for threads in "${THREADS[@]}"; do
        case "$profile" in
            standard|high-frequency)
                for dataset in "${DATASETS[@]}"; do run_one "$dataset" "$profile" "$threads" "$threads"; done
                ;;
            knn)
                for value in "${K_VALUES[@]}"; do
                    for dataset in "${DATASETS[@]}"; do run_one "$dataset" knn "$threads" "$value" --k "$value"; done
                done
                ;;
            sampling)
                for value in "${SAMPLE_FACTORS[@]}"; do
                    for dataset in "${DATASETS[@]}"; do run_one "$dataset" sampling "$threads" "$value" --sample-factor "$value"; done
                done
                ;;
        esac
    done
}

run_query_suite() {
    local threads dataset query label entry
    local -a entries=()
    case "$SUITE" in
        generated-queries)
            entries=(
                'spacev1b|generated/spacev1B_noise_025.bin|SPACEV1B_ne_025'
                'spacev1b|generated/spacev1B_noise_05.bin|SPACEV1B_ne_05'
                'spacev1b|generated/spacev1B_noise_1.bin|SPACEV1B_ne_1'
                'text-to-image|generated/text-to-image_noise_005.bin|TEXTTOIMAGE_ne_005'
                'text-to-image|generated/text-to-image_noise_01.bin|TEXTTOIMAGE_ne_01'
                'text-to-image|generated/text-to-image_noise_025.bin|TEXTTOIMAGE_ne_025'
                'turinganns|generated/turingANNs_noise_01.bin|turingANNs_ne_01'
                'turinganns|generated/turingANNs_noise_025.bin|turingANNs_ne_025'
                'turinganns|generated/turingANNs_noise_05.bin|turingANNs_ne_05'
            )
            ;;
        hard-queries)
            entries=(
                'seismic|seismic_queries.bin|SEISMIC'
                'seismic|other_queries/queries_hard1p_seismic_len256_znorm.bin|queries_hard1p_seismic_len256_znorm'
                'seismic|other_queries/queries_hard2p_seismic_len256_znorm.bin|queries_hard2p_seismic_len256_znorm'
                'seismic|other_queries/queries_hard5p_seismic_len256_znorm.bin|queries_hard5p_seismic_len256_znorm'
                'seismic|other_queries/queries_hard10p_seismic_len256_znorm.bin|queries_hard10p_seismic_len256_znorm'
                'seismic|other_queries/queries_size100_seismic_len256_znorm.bin|queries_size100_seismic_len256_znorm'
                'deep1b|deep1b_queries.bin|DEEP1b'
                'deep1b|other_queries/queries_hard1p_deep1b_len96_znorm.bin|queries_hard1p_deep1b_len96_znorm'
                'deep1b|other_queries/queries_hard2p_deep1b_len96_znorm.bin|queries_hard2p_deep1b_len96_znorm'
                'deep1b|other_queries/queries_hard5p_deep1b_len96_znorm.bin|queries_hard5p_deep1b_len96_znorm'
                'deep1b|other_queries/queries_hard10p_deep1b_len96_znorm.bin|queries_hard10p_deep1b_len96_znorm'
                'deep1b|other_queries/queries_orig100_deep1b_len96_znorm.bin|queries_orig100_deep1b_len96_znorm'
                'sald|SALD_queries.bin|SALD'
                'sald|other_queries/queries_hard1p_sald_len128_znorm.bin|queries_hard1p_sald_len128_znorm'
                'sald|other_queries/queries_hard2p_sald_len128_znorm.bin|queries_hard2p_sald_len128_znorm'
                'sald|other_queries/queries_hard5p_sald_len128_znorm.bin|queries_hard5p_sald_len128_znorm'
                'sald|other_queries/queries_hard10p_sald_len128_znorm.bin|queries_hard10p_sald_len128_znorm'
                'sald|other_queries/queries_size100_sald_len128_znorm.bin|queries_size100_sald_len128_znorm'
            )
            ;;
        noise-workloads)
            entries=(
                'seismic|noiseSeismic001.bin|SEISMIC001' 'seismic|noiseSeismic002.bin|SEISMIC002'
                'seismic|noiseSeismic005.bin|SEISMIC005' 'seismic|noiseSeismic01.bin|SEISMIC01'
                'sald|noiseSALD001.bin|SALD001' 'sald|noiseSALD002.bin|SALD002'
                'sald|noiseSALD005.bin|SALD005' 'sald|noiseSALD01.bin|SALD01'
            )
            ;;
    esac

    for threads in "${THREADS[@]}"; do
        for entry in "${entries[@]}"; do
            IFS='|' read -r dataset query label <<< "$entry"
            if [[ $DRY_RUN == false && $RERUN_EXISTING == false ]] && archive_exists "$label"; then
                printf 'Skipping query workload=%s run=%s: archive directory already exists at %s/%s\n' \
                    "$label" "$threads" "$RESULTS_ROOT" "$label" >&2
                continue
            fi
            local methods=${METHODS_OVERRIDE:-sax,sfa-depth,sfa-width}
            [[ $INDEX_TYPE == trie && -z $METHODS_OVERRIDE ]] && methods=sfa-depth,sfa-width
            local queue_number=${QUEUE_NUMBER:-$threads}
            local -a command=("$SCRIPT_DIR/run_dataset.sh" "$dataset" standard --threads "$threads" \
                --queue-number "$queue_number" --numa "$NUMA_MODE" --index-type "$INDEX_TYPE")
            [[ -n $DATASET_FILE ]] && command+=(--dataset-file "$DATASET_FILE")
            [[ -n $DATASET_SIZE ]] && command+=(--dataset-size "$DATASET_SIZE")
            [[ -n $QUERY_SIZE ]] && command+=(--query-size "$QUERY_SIZE")
            [[ -n $SAMPLE_SIZE ]] && command+=(--sample-size "$SAMPLE_SIZE")
            command+=(--sample-type "$SAMPLE_TYPE")
            [[ -n $SAMPLING_SEED ]] && command+=(--sampling-seed "$SAMPLING_SEED")
            [[ -n $MESSI_EXECUTABLE ]] && command+=(--binary "$MESSI_EXECUTABLE")
            [[ -n $DATA_ROOT ]] && command+=(--data-root "$DATA_ROOT")
            [[ -n $QUERY_ROOT ]] && command+=(--query-root "$QUERY_ROOT")
            [[ -n $SEISBENCH_ROOT ]] && command+=(--seisbench-root "$SEISBENCH_ROOT")
            [[ -n $SEISBENCH_QUERY_ROOT ]] && command+=(--seisbench-query-root "$SEISBENCH_QUERY_ROOT")
            $NO_SIMD && command+=(--no-simd)
            if [[ $INDEX_TYPE == trie ]]; then
                [[ -n $TRIE_MBR_DIMS ]] && command+=(--trie-mbr-dims "$TRIE_MBR_DIMS")
                command+=(--n-segments "$TRIE_RECORD_LB_DIMS")
                [[ -n $TRIE_SPLIT_DIMS ]] && command+=(--trie-split-dims "$TRIE_SPLIT_DIMS")
                $TRIE_RECORD_MBR_SUFFIX_BOUND && command+=(--trie-record-mbr-suffix-bound)
                [[ $TRIE_RECORD_MBR_SUFFIX_BOUND == false ]] && command+=(--no-trie-record-mbr-suffix-bound)
                $TRIE_STREAMING_LEAF_SCAN && command+=(--trie-streaming-leaf-scan)
                [[ $TRIE_STREAMING_LEAF_SCAN == false ]] && command+=(--no-trie-streaming-leaf-scan)
                [[ $TRIE_LEAF_IVF != 0 ]] && command+=(--trie-leaf-ivf "$TRIE_LEAF_IVF")
                [[ $TRIE_LEAF_IVF == 0 ]] && command+=(--no-trie-leaf-ivf)
                [[ $TRIE_LEAF_IVF_RAW_BALL_BOUND == false ]] && command+=(--no-trie-leaf-ivf-raw-ball-bound)
                [[ $TRIE_LEAF_IVF_RADIAL_BOUND == true ]] && command+=(--trie-leaf-ivf-radial-bound)
                if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
                    command+=(--trie-dynamic-alphabet --trie-min-fanout "$TRIE_MIN_FANOUT"
                              --trie-max-fanout "$TRIE_MAX_FANOUT"
                              --trie-alphabet-budget-bits "$TRIE_ALPHABET_BUDGET_BITS")
                else
                    command+=(--trie-fanout "$TRIE_FANOUT")
                fi
            fi
            [[ -n $LEAF_SIZE ]] && command+=(--leaf-size "$LEAF_SIZE")
            [[ -n $MIN_LEAF_SIZE ]] && command+=(--min-leaf-size "$MIN_LEAF_SIZE")
            $TRIE_QUERY_PARALLEL && command+=(--trie-query-parallel)
            $TRIE_QUERY_BATCH && command+=(--trie-query-batch)
            [[ -n $QUERY_REPORT_INTERVAL ]] && command+=(--query-report-interval "$QUERY_REPORT_INTERVAL")
            $PROFILE_QUERY_PHASES && command+=(--profile-query-phases)
            [[ $DYNAMIC_ROOT_SPLIT_VARIANCE == true ]] && command+=(--dynamic-root-split-variance)
            [[ $DYNAMIC_ROOT_SPLIT_VARIANCE_SPECIFIED == true && $DYNAMIC_ROOT_SPLIT_VARIANCE == false ]] && command+=(--no-dynamic-root-split-variance)
            command+=(--query-file "$query" --methods "$methods")
            if [[ $INDEX_TYPE == isax ]]; then
                $ENABLE_SOFA_V2 && command+=(--enable-sofa-v2)
                $ISAX_NODE_MBR && command+=(--isax-node-mbr)
                $ISAX_RECORD_MBR_SUFFIX_BOUND && command+=(--isax-record-mbr-suffix-bound)
                [[ -n $ISAX_MBR_DIMS ]] && command+=(--isax-mbr-dimensions "$ISAX_MBR_DIMS")
                $ISAX_RECORD_LB_TABLE && command+=(--isax-record-lb-table)
                $TIGHT_BOUND && command+=(--tight-bound)
                [[ $TIGHT_BOUND == false ]] && command+=(--no-tight-bound)
            fi
            $DRY_RUN && command+=(--dry-run)
            "${command[@]}"
            if [[ $DRY_RUN == false ]]; then "$SCRIPT_DIR/archive_results.sh" "$label" "$threads"; fi
        done
    done
}

case "$SUITE" in
    standard|high-frequency|knn|sampling) run_regular_suite "$SUITE" ;;
    *) run_query_suite ;;
esac
