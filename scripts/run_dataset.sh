#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=lib/datasets.sh
source "$SCRIPT_DIR/lib/datasets.sh"

usage() {
    cat <<'USAGE'
Usage: run_dataset.sh DATASET PROFILE --cpu-type N --queue-number N [OPTIONS]

Profiles: standard, high-frequency, knn, sampling

Options:
  --dataset-file PATH       Override the dataset filename/path
  --query-file PATH         Override the query filename/path
  --dataset-size N          Override the number of dataset records
  --query-size N            Number of queries (default: 100; high-frequency: 1)
  --k N                     Required by the knn profile
  --sample-factor F         Sampling fraction for sampling (default: 0.01)
  --sample-size N           Override the sample size directly
  --methods LIST            Comma-separated method names
  --no-tight-bound          Disable standard profile's tight-bound option
  --binary PATH             MESSI executable
  --data-root PATH          Main dataset root
  --query-root PATH         Main query root (default: data root)
  --seisbench-root PATH     SeisBench dataset root
  --seisbench-query-root P  SeisBench query root (default: dataset root)
  --dry-run                 Print commands without executing them
  -h, --help                Show this help

Methods: sax, sfa-depth, sfa-width, pisa-depth, pisa-width,
         spartan-depth, spartan-width
USAGE
}

die() { printf 'Error: %s\n' "$*" >&2; exit 2; }
is_positive_integer() { [[ $1 =~ ^[1-9][0-9]*$ ]]; }
resolve_path() {
    local root=$1 path=$2
    if [[ $path == /* ]]; then printf '%s\n' "$path"; else printf '%s/%s\n' "${root%/}" "$path"; fi
}
print_command() {
    printf '%q ' "$@"
    printf '\n'
}

[[ $# -ge 2 ]] || { usage >&2; exit 2; }
DATASET_ARG=$1
PROFILE=$2
shift 2

case "$PROFILE" in
    standard|high-frequency|knn|sampling) ;;
    *) die "unknown profile '$PROFILE'" ;;
esac

CPU_TYPE=
QUEUE_NUMBER=
DATASET_OVERRIDE=
QUERY_OVERRIDE=
DATASET_SIZE_OVERRIDE=
QUERY_SIZE_OVERRIDE=
K_SIZE=
SAMPLE_FACTOR=0.01
SAMPLE_SIZE_OVERRIDE=
METHODS_OVERRIDE=
TIGHT_BOUND=true
DRY_RUN=${MESSI_DRY_RUN:-false}
MESSI_EXECUTABLE=${MESSI_BINARY:-"$SCRIPT_DIR/../bin/MESSI"}
DATA_ROOT=${MESSI_DATA_ROOT:-/vol/tmp/schaefpa/messi_datasets}
SEISBENCH_ROOT=${MESSI_SEISBENCH_ROOT:-/vol/tmp/schaefpa/seismic}
QUERY_ROOT=${MESSI_QUERY_ROOT:-}
SEISBENCH_QUERY_ROOT=${MESSI_SEISBENCH_QUERY_ROOT:-}

[[ $DRY_RUN == true || $DRY_RUN == false ]] || die 'MESSI_DRY_RUN must be true or false'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu-type) [[ $# -ge 2 ]] || die "$1 requires a value"; CPU_TYPE=$2; shift 2 ;;
        --queue-number) [[ $# -ge 2 ]] || die "$1 requires a value"; QUEUE_NUMBER=$2; shift 2 ;;
        --dataset-file) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_OVERRIDE=$2; shift 2 ;;
        --query-file) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_OVERRIDE=$2; shift 2 ;;
        --dataset-size) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASET_SIZE_OVERRIDE=$2; shift 2 ;;
        --query-size) [[ $# -ge 2 ]] || die "$1 requires a value"; QUERY_SIZE_OVERRIDE=$2; shift 2 ;;
        --k) [[ $# -ge 2 ]] || die "$1 requires a value"; K_SIZE=$2; shift 2 ;;
        --sample-factor) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_FACTOR=$2; shift 2 ;;
        --sample-size) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_SIZE_OVERRIDE=$2; shift 2 ;;
        --methods) [[ $# -ge 2 ]] || die "$1 requires a value"; METHODS_OVERRIDE=$2; shift 2 ;;
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

[[ -n $CPU_TYPE ]] || die '--cpu-type is required'
[[ -n $QUEUE_NUMBER ]] || die '--queue-number is required'
is_positive_integer "$CPU_TYPE" || die '--cpu-type must be a positive integer'
is_positive_integer "$QUEUE_NUMBER" || die '--queue-number must be a positive integer'
[[ -n $DATASET_FILE ]] || die '--dataset-file is required for this dataset'
[[ -n $QUERY_FILE ]] || die '--query-file is required for this dataset'
is_positive_integer "$DATASET_SIZE" || die '--dataset-size must be a positive integer'

QUERY_SIZE=${QUERY_SIZE_OVERRIDE:-100}
[[ $PROFILE == high-frequency ]] && QUERY_SIZE=${QUERY_SIZE_OVERRIDE:-1}
is_positive_integer "$QUERY_SIZE" || die '--query-size must be a positive integer'

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
is_positive_integer "$SAMPLE_SIZE" || die '--sample-size must be a positive integer'

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
    --timeseries-size "$TS_SIZE"
    --dataset-size "$DATASET_SIZE"
    --flush-limit 300000
    --read-block 20000
    --sax-cardinality 8
    --queries "$QUERY_PATH"
    --queries-size "$QUERY_SIZE"
    --queue-number "$QUEUE_NUMBER"
    --cpu-type "$CPU_TYPE"
    --leaf-size 20000
    --min-leaf-size 20000
    --initial-lbl-size 20000
    --SIMD
)

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

    args+=(--function-type "$function_type")
    if [[ -n $histogram_type ]]; then
        args+=(--sample-size "$SAMPLE_SIZE" --sample-type 3 --is-norm --histogram-type "$histogram_type")
        [[ $function_type == 4 || $function_type == 6 ]] && args+=(--sfa-n-coefficients "$COEFF_NUMBER")
        if [[ $PROFILE == standard && $TIGHT_BOUND == true ]]; then args+=(--tight-bound); fi
        if [[ $DATASET_ID == sift1b && $histogram_type == 1 && ( $PROFILE == knn || $PROFILE == sampling ) ]]; then
            args+=(--tight-bound)
        fi
    fi
    if [[ $PROFILE == knn ]]; then args+=(--topk --k-size "$K_SIZE"); fi

    printf 'Running dataset=%s profile=%s method=%s\n' "$DATASET_ID" "$PROFILE" "$method" >&2
    if [[ $DRY_RUN == true ]]; then
        print_command "$MESSI_EXECUTABLE" "${args[@]}"
    else
        [[ -x $MESSI_EXECUTABLE ]] || die "MESSI executable is not executable: $MESSI_EXECUTABLE"
        [[ -f $DATASET_PATH ]] || die "dataset file does not exist: $DATASET_PATH"
        [[ -f $QUERY_PATH ]] || die "query file does not exist: $QUERY_PATH"
        "$MESSI_EXECUTABLE" "${args[@]}"
    fi
}

for method in "${METHOD_LIST[@]}"; do
    [[ -n $method ]] || die 'method list contains an empty value'
    run_method "$method"
done
