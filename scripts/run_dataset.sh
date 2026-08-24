#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=lib/datasets.sh
source "$SCRIPT_DIR/lib/datasets.sh"

usage() {
    cat <<'USAGE'
Usage: run_dataset.sh DATASET PROFILE --threads N|auto [--queue-number N] [OPTIONS]

Profiles: standard, high-frequency, knn, sampling

Options:
  --dataset-file PATH       Override the dataset filename/path
  --query-file PATH         Override the query filename/path
  --dataset-size N          Override dataset records; accepts 100m, 1mio, 20k
  --query-size N            Queries; accepts count suffixes (default: 100; high-frequency: 1)
  --leaf-size N             Maximum records per leaf; accepts count suffixes (default: 20000)
  --min-leaf-size N         iSAX query-leaf threshold; trie does not enforce it
  --k N                     Required by the knn profile
  --sample-factor F         Sampling fraction for sampling (default: 0.01)
  --sample-size N           Override the sample size; accepts count suffixes
  --sampling-seed N         Seed for random direct-CLI sampling (benchmark runners use uniform sampling)
  --no-simd                 Disable SIMD even when AVX2 is available
  --methods LIST            Comma-separated method names
  --index-type TYPE         Index layout: isax (default) or trie
  --trie-mbr-dims N         Trie MBR dimensions (default: 128; capped by series length)
  --n-segments N            Trie record-prefix lower-bound dimensions (default: 16; range: 16--64;
                            alias: --trie-record-lb-dims)
  --trie-split-dims N       Trie split-candidate dimensions (default: min(32, MBR dimensions))
  --trie-record-mbr-suffix-bound
                            Add leaf-MBR contributions outside the 16 record-bound dimensions
  --trie-fanout 2|4|8       Trie symbolic split fanout (default: 8)
  --trie-dynamic-alphabet  Use one global variance-weighted alphabet allocation
  --trie-min-fanout N      Minimum dynamic trie fanout (default: 2)
  --trie-max-fanout N      Maximum dynamic trie fanout (default: 16)
  --trie-alphabet-budget-bits N
                            Average dynamic alphabet budget in bits (default: 3)
  --dynamic-root-split-variance
                            Use variance-assigned root bits for iSAX SFA/PISA/SPARTAN
                            (the iSAX runner default)
  --no-dynamic-root-split-variance
                            Disable variance-assigned root bits for learned iSAX methods
  --trie-query-parallel     Parallelize each trie query (default; retained for compatibility)
  --trie-query-batch        Batch independent trie queries instead
  --query-report-interval N Print first, every Nth completed, and final query row (0=none; default: 10)
  --profile-query-phases    Measure traversal, lower-bound, and exact-distance work
  --no-tight-bound          Disable standard profile's tight-bound option
  --binary PATH             MESSI executable
  MESSI_SHELL_LOG_DIR       Directory for per-method shell output logs
  --data-root PATH          Main dataset root
  --query-root PATH         Main query root (default: data root)
  --seisbench-root PATH     SeisBench dataset root
  --seisbench-query-root P  SeisBench query root (default: dataset root)
  --numa MODE               auto (default), none, or number of NUMA nodes to use
  --dry-run                 Print commands without executing them
  -h, --help                Show this help

Methods: sax, sfa-depth, sfa-width, pisa-depth, pisa-width,
         spartan-depth, spartan-width

Note: trie benchmarks exclude SAX.  The current trie/SAX record bound uses
only a prefix of its 64-segment representation and is not comparable to the
whole-series 16-segment iSAX/SAX bound.
USAGE
}

die() { printf 'Error: %s\n' "$*" >&2; exit 2; }
is_positive_integer() { [[ $1 =~ ^[1-9][0-9]*$ ]]; }
is_nonnegative_integer() { [[ $1 =~ ^[0-9]+$ ]]; }
normalize_count() {
    local value=${1,,} number suffix multiplier=1
    if [[ $value =~ ^([1-9][0-9]*)(k|m|mio|g)?$ ]]; then
        number=${BASH_REMATCH[1]}
        suffix=${BASH_REMATCH[2]}
    else
        return 1
    fi
    case "$suffix" in
        k) multiplier=1000 ;;
        m|mio) multiplier=1000000 ;;
        g) multiplier=1000000000 ;;
    esac
    printf '%s\n' "$((number * multiplier))"
}
resolve_path() {
    local root=$1 path=$2
    if [[ $path == /* ]]; then printf '%s\n' "$path"; else printf '%s/%s\n' "${root%/}" "$path"; fi
}
print_command() {
    printf '%q ' "$@"
    printf '\n'
}
format_count() {
    local value=$1
    if (( value >= 1000000000 )); then
        awk -v value="$value" 'BEGIN { printf "%.3g G", value / 1000000000 }'
    elif (( value >= 1000000 )); then
        awk -v value="$value" 'BEGIN { printf "%.3g M", value / 1000000 }'
    elif (( value >= 1000 )); then
        awk -v value="$value" 'BEGIN { printf "%.3g K", value / 1000 }'
    else
        printf '%s' "$value"
    fi
}

# MESSI's legacy SIGINT handler is interactive and would otherwise prompt once
# for every command in a suite.  Run MESSI with SIGINT ignored and let this
# script own cancellation: one Ctrl-C terminates the active child and prevents
# subsequent methods from starting.
ACTIVE_MESSI_PID=
RUN_OUTPUT_FILE=
RUN_LOG_FILE=
stop_active_messi() {
    trap - INT TERM
    if [[ -n ${ACTIVE_MESSI_PID:-} ]]; then
        kill -TERM "$ACTIVE_MESSI_PID" 2>/dev/null || true
        wait "$ACTIVE_MESSI_PID" 2>/dev/null || true
        ACTIVE_MESSI_PID=
    fi
    printf '\nInterrupted; no further experiments will be started.\n' >&2
    exit 130
}
run_messi() {
    local capture_dir fifo tee_pid status=0
    capture_dir=$(mktemp -d "${TMPDIR:-/tmp}/messi-run.XXXXXX")
    fifo=$capture_dir/output
    mkfifo "$fifo"
    tee -a "$RUN_LOG_FILE" "$RUN_OUTPUT_FILE" < "$fifo" &
    tee_pid=$!
    # The wrapper execs MESSI, retaining SIGINT=ignore in the final process.
    bash -c 'trap "" INT; exec "$@"' messi "$@" > "$fifo" 2>&1 &
    ACTIVE_MESSI_PID=$!
    wait "$ACTIVE_MESSI_PID" || status=$?
    ACTIVE_MESSI_PID=
    wait "$tee_pid" || true
    rm -f -- "$fifo"
    rmdir "$capture_dir"
    return "$status"
}
trap stop_active_messi INT TERM

[[ $# -ge 2 ]] || { usage >&2; exit 2; }
DATASET_ARG=$1
PROFILE=$2
shift 2

case "$PROFILE" in
    standard|high-frequency|knn|sampling) ;;
    *) die "unknown profile '$PROFILE'" ;;
esac

THREADS=
NUMA_MODE=auto
QUEUE_NUMBER=
DATASET_OVERRIDE=
QUERY_OVERRIDE=
DATASET_SIZE_OVERRIDE=
QUERY_SIZE_OVERRIDE=
LEAF_SIZE=20000
MIN_LEAF_SIZE=
K_SIZE=
SAMPLE_FACTOR=0.01
SAMPLE_SIZE_OVERRIDE=
SAMPLING_SEED=1
NO_SIMD=false
METHODS_OVERRIDE=
INDEX_TYPE=isax
TRIE_QUERY_PARALLEL=false
TRIE_QUERY_BATCH=false
TRIE_MBR_DIMS=
TRIE_RECORD_LB_DIMS=16
TRIE_SPLIT_DIMS=
TRIE_RECORD_MBR_SUFFIX_BOUND=false
TRIE_FANOUT=8
TRIE_DYNAMIC_ALPHABET=false
TRIE_MIN_FANOUT=2
TRIE_MAX_FANOUT=16
TRIE_ALPHABET_BUDGET_BITS=3
PROFILE_QUERY_PHASES=false
QUERY_REPORT_INTERVAL=
# Resolved after parsing because the default depends on --index-type.
DYNAMIC_ROOT_SPLIT_VARIANCE=
TIGHT_BOUND=true
DRY_RUN=${MESSI_DRY_RUN:-false}
MESSI_EXECUTABLE=${MESSI_BINARY:-"$SCRIPT_DIR/../bin/MESSI"}
# Keep the shell transcript with MESSI's CSV logs by default.  An explicit
# MESSI_SHELL_LOG_DIR remains useful for CI or a separate transcript archive.
MESSI_SHELL_LOG_DIR=${MESSI_SHELL_LOG_DIR:-${MESSI_LOG_ROOT:-"$HOME/MESSI_logs"}}
DATA_ROOT=${MESSI_DATA_ROOT:-/vol/tmp/schaefpa/messi_datasets}
SEISBENCH_ROOT=${MESSI_SEISBENCH_ROOT:-/vol/tmp/schaefpa/seismic}
QUERY_ROOT=${MESSI_QUERY_ROOT:-}
SEISBENCH_QUERY_ROOT=${MESSI_SEISBENCH_QUERY_ROOT:-}

[[ $DRY_RUN == true || $DRY_RUN == false ]] || die 'MESSI_DRY_RUN must be true or false'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) [[ $# -ge 2 ]] || die "$1 requires a value"; THREADS=$2; shift 2 ;;
        --numa) [[ $# -ge 2 ]] || die "$1 requires a value"; NUMA_MODE=$2; shift 2 ;;
        --queue-number) [[ $# -ge 2 ]] || die "$1 requires a value"; QUEUE_NUMBER=$2; shift 2 ;;
        --dataset-file) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_OVERRIDE=$2; shift 2 ;;
        --query-file) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_OVERRIDE=$2; shift 2 ;;
        --dataset-size) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_SIZE_OVERRIDE=$2; shift 2 ;;
        --query-size) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_SIZE_OVERRIDE=$2; shift 2 ;;
        --leaf-size) [[ $# -ge 2 ]] || die "$1 requires a value"; LEAF_SIZE=$2; shift 2 ;;
        --min-leaf-size) [[ $# -ge 2 ]] || die "$1 requires a value"; MIN_LEAF_SIZE=$2; shift 2 ;;
        --k) [[ $# -ge 2 ]] || die "$1 requires a value"; K_SIZE=$2; shift 2 ;;
        --sample-factor) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_FACTOR=$2; shift 2 ;;
        --sample-size) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_SIZE_OVERRIDE=$2; shift 2 ;;
        --sampling-seed) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLING_SEED=$2; shift 2 ;;
        --no-simd) NO_SIMD=true; shift ;;
        --methods) [[ $# -ge 2 ]] || die "$1 requires a value"; METHODS_OVERRIDE=$2; shift 2 ;;
        --index-type) [[ $# -ge 2 ]] || die "$1 requires a value"; INDEX_TYPE=$2; shift 2 ;;
        --trie-query-parallel) TRIE_QUERY_PARALLEL=true; shift ;;
        --trie-query-batch) TRIE_QUERY_BATCH=true; shift ;;
        --trie-mbr-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MBR_DIMS=$2; shift 2 ;;
        --n-segments|--trie-record-lb-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_RECORD_LB_DIMS=$2; shift 2 ;;
        --trie-split-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_SPLIT_DIMS=$2; shift 2 ;;
        --trie-record-mbr-suffix-bound) TRIE_RECORD_MBR_SUFFIX_BOUND=true; shift ;;
        --trie-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_FANOUT=$2; shift 2 ;;
        --trie-dynamic-alphabet) TRIE_DYNAMIC_ALPHABET=true; shift ;;
        --trie-min-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MIN_FANOUT=$2; shift 2 ;;
        --trie-max-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MAX_FANOUT=$2; shift 2 ;;
        --trie-alphabet-budget-bits) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_ALPHABET_BUDGET_BITS=$2; shift 2 ;;
        --query-report-interval) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_REPORT_INTERVAL=$2; shift 2 ;;
        --profile-query-phases) PROFILE_QUERY_PHASES=true; shift ;;
        --dynamic-root-split-variance) DYNAMIC_ROOT_SPLIT_VARIANCE=true; shift ;;
        --no-dynamic-root-split-variance) DYNAMIC_ROOT_SPLIT_VARIANCE=false; shift ;;
        --no-tight-bound) TIGHT_BOUND=false; shift ;;
        --binary) [[ $# -ge 2 ]] || die "$1 requires a value"; MESSI_EXECUTABLE=$2; shift 2 ;;
        --data-root) [[ $# -ge 2 ]] || die "$1 requires a value"; DATA_ROOT=$2; shift 2 ;;
        --query-root) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_ROOT=$2; shift 2 ;;
        --seisbench-root) [[ $# -ge 2 ]] || die "$1 requires a value"; SEISBENCH_ROOT=$2; shift 2 ;;
        --seisbench-query-root) [[ $# -ge 2 ]] || die "$1 requires a value"; SEISBENCH_QUERY_ROOT=$2; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option '$1'" ;;
    esac
done

load_dataset "$DATASET_ARG" "$PROFILE"

[[ -n $DATASET_OVERRIDE ]] && DATASET_FILE=$DATASET_OVERRIDE
[[ -n $QUERY_OVERRIDE ]] && QUERY_FILE=$QUERY_OVERRIDE
[[ -n $DATASET_SIZE_OVERRIDE ]] && DATASET_SIZE=$DATASET_SIZE_OVERRIDE

[[ -n $THREADS ]] || die '--threads is required'
[[ $THREADS == auto ]] || is_positive_integer "$THREADS" || die '--threads must be a positive integer or auto'
[[ $NUMA_MODE == auto || $NUMA_MODE == none ]] || is_positive_integer "$NUMA_MODE" || die '--numa must be auto, none, or a positive integer'
[[ -z $QUEUE_NUMBER ]] || is_positive_integer "$QUEUE_NUMBER" || die '--queue-number must be a positive integer'
is_nonnegative_integer "$SAMPLING_SEED" || die '--sampling-seed must be a nonnegative integer'
[[ -z $QUERY_REPORT_INTERVAL ]] || is_nonnegative_integer "$QUERY_REPORT_INTERVAL" || die '--query-report-interval must be zero or a positive integer'
[[ -n $DATASET_FILE ]] || die '--dataset-file is required for this dataset'
[[ -n $QUERY_FILE ]] || die '--query-file is required for this dataset'
DATASET_SIZE=$(normalize_count "$DATASET_SIZE") || die '--dataset-size must be a positive integer or use k/m/mio/g'
[[ $INDEX_TYPE == isax || $INDEX_TYPE == trie ]] || die '--index-type must be isax or trie'
if [[ -z $DYNAMIC_ROOT_SPLIT_VARIANCE ]]; then
    [[ $INDEX_TYPE == isax ]] && DYNAMIC_ROOT_SPLIT_VARIANCE=true || DYNAMIC_ROOT_SPLIT_VARIANCE=false
fi
[[ $TRIE_QUERY_PARALLEL == false || $INDEX_TYPE == trie ]] || die '--trie-query-parallel requires --index-type trie'
[[ $TRIE_QUERY_BATCH == false || $INDEX_TYPE == trie ]] || die '--trie-query-batch requires --index-type trie'
[[ -z $TRIE_MBR_DIMS || $INDEX_TYPE == trie ]] || die '--trie-mbr-dims requires --index-type trie'
[[ $TRIE_RECORD_LB_DIMS == 16 || $INDEX_TYPE == trie ]] || die '--n-segments requires --index-type trie'
[[ $TRIE_RECORD_MBR_SUFFIX_BOUND == false || $INDEX_TYPE == trie ]] || die '--trie-record-mbr-suffix-bound requires --index-type trie'
[[ $TRIE_FANOUT == 8 || $INDEX_TYPE == trie ]] || die '--trie-fanout requires --index-type trie'
[[ $TRIE_DYNAMIC_ALPHABET == false || $INDEX_TYPE == trie ]] || die '--trie-dynamic-alphabet requires --index-type trie'
[[ $TRIE_QUERY_PARALLEL == false || $TRIE_QUERY_BATCH == false ]] || die 'choose at most one of --trie-query-parallel and --trie-query-batch'
[[ $DYNAMIC_ROOT_SPLIT_VARIANCE == false || $INDEX_TYPE == isax ]] || die '--dynamic-root-split-variance requires --index-type isax'

QUERY_SIZE=${QUERY_SIZE_OVERRIDE:-100}
[[ $PROFILE == high-frequency ]] && QUERY_SIZE=${QUERY_SIZE_OVERRIDE:-1}
QUERY_SIZE=$(normalize_count "$QUERY_SIZE") || die '--query-size must be a positive integer or use k/m/mio/g'
LEAF_SIZE=$(normalize_count "$LEAF_SIZE") || die '--leaf-size must be a positive integer or use k/m/mio/g'
MIN_LEAF_SIZE=${MIN_LEAF_SIZE:-$LEAF_SIZE}
MIN_LEAF_SIZE=$(normalize_count "$MIN_LEAF_SIZE") || die '--min-leaf-size must be a positive integer or use k/m/mio/g'
(( MIN_LEAF_SIZE <= LEAF_SIZE )) || die '--min-leaf-size cannot exceed --leaf-size'
if [[ $INDEX_TYPE == trie ]]; then
    is_positive_integer "$TRIE_RECORD_LB_DIMS" || die '--n-segments must be a positive integer'
    (( TRIE_RECORD_LB_DIMS >= 16 && TRIE_RECORD_LB_DIMS <= 64 )) || \
        die '--n-segments must be between 16 and 64 for trie runs'
    TRIE_MBR_DIMS=${TRIE_MBR_DIMS:-128}
    is_positive_integer "$TRIE_MBR_DIMS" || die '--trie-mbr-dims must be a positive integer'
    (( TRIE_MBR_DIMS <= 128 )) || TRIE_MBR_DIMS=128
    (( TRIE_MBR_DIMS <= TS_SIZE )) || TRIE_MBR_DIMS=$TS_SIZE
    (( TRIE_MBR_DIMS >= 16 )) || die '--trie-mbr-dims must be at least 16'
    (( TRIE_RECORD_LB_DIMS <= TRIE_MBR_DIMS )) || \
        die '--n-segments cannot exceed --trie-mbr-dims'
    TRIE_SPLIT_DIMS=${TRIE_SPLIT_DIMS:-$(( TRIE_RECORD_LB_DIMS > 32 ? TRIE_RECORD_LB_DIMS : (TRIE_MBR_DIMS < 32 ? TRIE_MBR_DIMS : 32) ))}
    is_positive_integer "$TRIE_SPLIT_DIMS" || die '--trie-split-dims must be a positive integer'
    (( TRIE_SPLIT_DIMS >= TRIE_RECORD_LB_DIMS && TRIE_SPLIT_DIMS <= TRIE_MBR_DIMS )) || \
        die '--trie-split-dims must be between --n-segments and --trie-mbr-dims'
    # SFA needs a training pool at least as wide as the trie word.  Keep the
    # historical iSAX coefficient defaults untouched; widen only a trie that
    # explicitly requests more MBR dimensions.  PISA uses the full real FFT.
    (( TRIE_MBR_DIMS <= COEFF_NUMBER )) || COEFF_NUMBER=$TRIE_MBR_DIMS
if [[ $TRIE_DYNAMIC_ALPHABET == false ]]; then
    [[ $TRIE_FANOUT == 2 || $TRIE_FANOUT == 4 || $TRIE_FANOUT == 8 ]] || die '--trie-fanout must be 2, 4, or 8'
else
    [[ $TRIE_FANOUT == 8 ]] || die '--trie-fanout cannot be combined with --trie-dynamic-alphabet'
    [[ $TRIE_MIN_FANOUT =~ ^(2|4|8|16|32|64|128|256)$ ]] || die '--trie-min-fanout must be a power of two between 2 and 256'
    [[ $TRIE_MAX_FANOUT =~ ^(2|4|8|16|32|64|128|256)$ ]] || die '--trie-max-fanout must be a power of two between 2 and 256'
    [[ $TRIE_ALPHABET_BUDGET_BITS =~ ^[1-8]$ ]] || die '--trie-alphabet-budget-bits must be between 1 and 8'
fi
fi

if [[ $PROFILE == knn ]]; then
    [[ -n $K_SIZE ]] || die '--k is required by the knn profile'
    is_positive_integer "$K_SIZE" || die '--k must be a positive integer'
fi

SAMPLE_SIZE=${SAMPLE_SIZE_OVERRIDE:-1000000}
if [[ $PROFILE == sampling && -z $SAMPLE_SIZE_OVERRIDE ]]; then
    [[ $SAMPLE_FACTOR =~ ^(0|[0-9]+)(\.[0-9]+)?$ ]] || die '--sample-factor must be a number in (0, 1]'
    awk -v factor="$SAMPLE_FACTOR" 'BEGIN { exit !(factor > 0 && factor <= 1) }' || die '--sample-factor must be in (0, 1]'
    SAMPLE_SIZE=$(awk -v size="$DATASET_SIZE" -v factor="$SAMPLE_FACTOR" 'BEGIN { printf "%.0f", size * factor }')
fi
SAMPLE_SIZE=$(normalize_count "$SAMPLE_SIZE") || die '--sample-size must be a positive integer or use k/m/mio/g'
if (( SAMPLE_SIZE > DATASET_SIZE )); then
    printf 'Warning: binning sample size %s exceeds dataset size %s; using dataset size.\n' \
        "$(format_count "$SAMPLE_SIZE")" "$(format_count "$DATASET_SIZE")" >&2
    SAMPLE_SIZE=$DATASET_SIZE
fi

ROOT=$DATA_ROOT
QUERY_BASE=${QUERY_ROOT:-$DATA_ROOT}
if [[ $DATASET_ROOT_KIND == seisbench ]]; then
    ROOT=$SEISBENCH_ROOT
    QUERY_BASE=${SEISBENCH_QUERY_ROOT:-$SEISBENCH_ROOT}
fi
DATASET_PATH=$(resolve_path "$ROOT" "$DATASET_FILE")
QUERY_PATH=$(resolve_path "$QUERY_BASE" "$QUERY_FILE")

case "$PROFILE" in
    standard) DEFAULT_METHODS=sax,sfa-depth,sfa-width,pisa-depth,pisa-width,spartan-depth,spartan-width ;;
    high-frequency) DEFAULT_METHODS=sfa-width ;;
    knn) DEFAULT_METHODS=sax,sfa-depth,sfa-width ;;
    sampling) DEFAULT_METHODS=sfa-depth,sfa-width ;;
esac
if [[ $INDEX_TYPE == trie ]]; then
    case "$PROFILE" in
        standard) DEFAULT_METHODS=sfa-depth,sfa-width,pisa-depth,pisa-width,spartan-depth,spartan-width ;;
        knn) DEFAULT_METHODS=sfa-depth,sfa-width ;;
    esac
fi
METHODS=${METHODS_OVERRIDE:-$DEFAULT_METHODS}

IFS=',' read -r -a METHOD_LIST <<< "$METHODS"
[[ ${#METHOD_LIST[@]} -gt 0 ]] || die 'at least one method is required'

COMMON_ARGS=(
    --dataset "$DATASET_PATH"
)
$APPLY_Z_NORM && COMMON_ARGS+=(--apply-z-norm)
$FILETYPE_INT && COMMON_ARGS+=(--filetype-int)
COMMON_ARGS+=(
    --in-memory
    --index-type "$INDEX_TYPE"
    --timeseries-size "$TS_SIZE"
    --dataset-size "$DATASET_SIZE"
    --flush-limit 300000
    --read-block 20000
    --sax-cardinality 8
    --queries "$QUERY_PATH"
    --queries-size "$QUERY_SIZE"
    --threads "$THREADS"
    --numa "$NUMA_MODE"
    --leaf-size "$LEAF_SIZE"
    --min-leaf-size "$MIN_LEAF_SIZE"
    --initial-lbl-size 20000
    --sampling-seed "$SAMPLING_SEED"
)
[[ $NO_SIMD == true ]] && COMMON_ARGS+=(--no-simd)
[[ -n $QUEUE_NUMBER ]] && COMMON_ARGS+=(--queue-number "$QUEUE_NUMBER")
[[ $INDEX_TYPE == trie ]] && COMMON_ARGS+=(--trie-mbr-dimensions "$TRIE_MBR_DIMS")
[[ $INDEX_TYPE == trie ]] && COMMON_ARGS+=(--n-segments "$TRIE_RECORD_LB_DIMS")
[[ $INDEX_TYPE == trie ]] && COMMON_ARGS+=(--trie-split-dimensions "$TRIE_SPLIT_DIMS")
[[ $TRIE_RECORD_MBR_SUFFIX_BOUND == true ]] && COMMON_ARGS+=(--trie-record-mbr-suffix-bound)
if [[ $INDEX_TYPE == trie ]]; then
    if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
        COMMON_ARGS+=(--trie-dynamic-alphabet
                      --trie-min-fanout "$TRIE_MIN_FANOUT"
                      --trie-max-fanout "$TRIE_MAX_FANOUT"
                      --trie-alphabet-budget-bits "$TRIE_ALPHABET_BUDGET_BITS")
    else
        COMMON_ARGS+=(--trie-fanout "$TRIE_FANOUT")
    fi
fi
if [[ $TRIE_QUERY_PARALLEL == true ]]; then
    COMMON_ARGS+=(--trie-query-parallel)
fi
if [[ $TRIE_QUERY_BATCH == true ]]; then
    COMMON_ARGS+=(--trie-query-batch)
fi
if [[ $PROFILE_QUERY_PHASES == true ]]; then
    COMMON_ARGS+=(--profile-query-phases)
fi
if [[ -n $QUERY_REPORT_INTERVAL ]]; then
    COMMON_ARGS+=(--query-report-interval "$QUERY_REPORT_INTERVAL")
fi

SUMMARY_ROWS=
collect_run_summary() {
    local method=$1 transcript=$2 fields method_name binning layout leaf_cap fanout wall lower lower_pct exact exact_pct
    fields=$(awk '
        /^  wall time[[:space:]]*:/ { wall = $4 " " $5 }
        /^  lower bounds[[:space:]]*:/ {
            line = $0; sub(/^.*:[[:space:]]*/, "", line)
            marker = index(line, " (")
            if (marker != 0) {
                lower = substr(line, 1, marker - 1)
                lower_pct = substr(line, marker + 2)
                sub(/[[:space:]].*$/, "", lower_pct)
            }
        }
        /^  exact distances[[:space:]]*:/ {
            line = $0; sub(/^.*:[[:space:]]*/, "", line)
            marker = index(line, " (")
            if (marker != 0) {
                exact = substr(line, 1, marker - 1)
                exact_pct = substr(line, marker + 2)
                sub(/[[:space:]].*$/, "", exact_pct)
            }
        }
        END {
            if (wall != "" && lower != "" && exact != "")
                print wall "\t" lower "\t" lower_pct "\t" exact "\t" exact_pct
        }
    ' "$transcript")
    if [[ -z $fields ]]; then
        printf 'warning: could not parse query summary for method=%s; omitting it from suite summary\n' "$method" |
            tee -a "$RUN_LOG_FILE" >&2
        return 0
    fi
    IFS=$'\t' read -r wall lower lower_pct exact exact_pct <<< "$fields"

    case "$method" in
        sax) method_name=SAX; binning=depth ;;
        sfa-depth) method_name=SFA; binning=depth ;;
        sfa-width) method_name=SFA; binning=width ;;
        pisa-depth) method_name=PISA; binning=depth ;;
        pisa-width) method_name=PISA; binning=width ;;
        spartan-depth) method_name=SPARTAN; binning=depth ;;
        spartan-width) method_name=SPARTAN; binning=width ;;
    esac
    [[ $INDEX_TYPE == trie ]] && layout=Trie || layout=iSAX
    leaf_cap=$(format_count "$LEAF_SIZE")
    if [[ $INDEX_TYPE == trie ]]; then
        if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
            fanout="dynamic ${TRIE_MIN_FANOUT}-${TRIE_MAX_FANOUT}"
        else
            fanout="${TRIE_FANOUT}-way"
        fi
    else
        fanout='-'
    fi
    SUMMARY_ROWS+="$layout|$method_name|$binning|$leaf_cap|$fanout|$lower|$lower_pct|$exact|$exact_pct|$wall"$'\n'
}

print_suite_summary() {
    [[ -n $SUMMARY_ROWS ]] || return 0
    printf '\n=== Benchmark summary: dataset=%s, profile=%s ===\n' "$DATASET_ID" "$PROFILE" |
        tee -a "$RUN_LOG_FILE" >&2
    printf '%-6s  %-8s  %-5s  %8s  %6s  %18s  %8s  %18s  %8s  %11s\n' \
        Layout Method Binning 'Leaf cap.' Fanout 'Lower bounds/query' 'LB %' 'Exact comps/query' 'Exact %' Walltime |
        tee -a "$RUN_LOG_FILE" >&2
    printf '%s\n' '------------------------------------------------------------------------------------------------------------------------' |
        tee -a "$RUN_LOG_FILE" >&2
    while IFS='|' read -r layout method binning leaf_cap fanout lower lower_pct exact exact_pct wall; do
        [[ -n $layout ]] || continue
        printf '%-6s  %-8s  %-5s  %8s  %6s  %18s  %8s  %18s  %8s  %11s\n' \
            "$layout" "$method" "$binning" "$leaf_cap" "$fanout" "$lower" "$lower_pct" "$exact" "$exact_pct" "$wall" |
            tee -a "$RUN_LOG_FILE" >&2
    done <<< "$SUMMARY_ROWS"
}

run_method() {
    local method=$1 function_type histogram_type=
    local -a args=("${COMMON_ARGS[@]}")

    case "$method" in
        sax) function_type=3 ;;
        sfa-depth) function_type=4; histogram_type=1 ;;
        sfa-width) function_type=4; histogram_type=2 ;;
        pisa-depth) function_type=6; histogram_type=1 ;;
        pisa-width) function_type=6; histogram_type=2 ;;
        spartan-depth) function_type=5; histogram_type=1 ;;
        spartan-width) function_type=5; histogram_type=2 ;;
        *) die "unknown method '$method'" ;;
    esac

    if [[ $INDEX_TYPE == trie && $method == sax ]]; then
        die 'trie/SAX is excluded from benchmark runs because its current record bound is prefix-only'
    fi

    args+=(--function-type "$function_type")
    # SAX has no learned transform variance.  Keep its legacy uniform root
    # split while applying the learned root-bit allocation to SFA/PISA/SPARTAN.
    if [[ $DYNAMIC_ROOT_SPLIT_VARIANCE == true && $function_type != 3 ]]; then
        args+=(--dynamic-root-split-variance)
    fi
    if [[ -n $histogram_type ]]; then
        # Uniformly spaced samples keep binning I/O monotonic.  Random sampling
        # previously issued one scattered seek/read for each sampled record.
        args+=(--sample-size "$SAMPLE_SIZE" --sample-type 2 --is-norm --histogram-type "$histogram_type")
        [[ $function_type == 4 ]] && args+=(--sfa-n-coefficients "$COEFF_NUMBER")
        if [[ $PROFILE == standard && $TIGHT_BOUND == true ]]; then args+=(--tight-bound); fi
        if [[ $DATASET_ID == sift1b && $histogram_type == 1 && ( $PROFILE == knn || $PROFILE == sampling ) ]]; then
            args+=(--tight-bound)
        fi
    fi
    if [[ $PROFILE == knn ]]; then args+=(--topk --k-size "$K_SIZE"); fi

    if [[ $DRY_RUN == false && -z $RUN_LOG_FILE ]]; then
        mkdir -p -- "$MESSI_SHELL_LOG_DIR"
        RUN_LOG_FILE="$MESSI_SHELL_LOG_DIR/MESSI_${DATASET_ID}_${PROFILE}_$(date +%Y%m%d_%H%M%S)_$$.log"
        : > "$RUN_LOG_FILE"
    fi
    if [[ $DRY_RUN == true ]]; then
        printf 'Running dataset=%s profile=%s method=%s\n' "$DATASET_ID" "$PROFILE" "$method" >&2
    else
        printf 'Running dataset=%s profile=%s method=%s\n' "$DATASET_ID" "$PROFILE" "$method" |
            tee -a "$RUN_LOG_FILE" >&2
    fi
    if [[ $DRY_RUN == true ]]; then
        print_command "$MESSI_EXECUTABLE" "${args[@]}"
    else
        [[ -x $MESSI_EXECUTABLE ]] || die "MESSI executable is not executable: $MESSI_EXECUTABLE"
        [[ -f $DATASET_PATH ]] || die "dataset file does not exist: $DATASET_PATH"
        [[ -f $QUERY_PATH ]] || die "query file does not exist: $QUERY_PATH"
        mkdir -p -- "$MESSI_SHELL_LOG_DIR"
        RUN_OUTPUT_FILE=$(mktemp "${TMPDIR:-/tmp}/messi-summary.XXXXXX")
        run_messi "$MESSI_EXECUTABLE" "${args[@]}"
        collect_run_summary "$method" "$RUN_OUTPUT_FILE"
        rm -f -- "$RUN_OUTPUT_FILE"
        RUN_OUTPUT_FILE=
    fi
}

for method in "${METHOD_LIST[@]}"; do
    [[ -n $method ]] || die 'method list contains an empty value'
    run_method "$method"
done
if [[ $DRY_RUN == false ]]; then
    print_suite_summary
fi
