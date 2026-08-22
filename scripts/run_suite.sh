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
  --threads LIST          Comma-separated CPU/queue counts (default: 36)
  --k-values LIST         K values for knn (default: 20,50)
  --sample-factors LIST   Factors for sampling (default: 0.15,...,0.5)
  --datasets LIST         Limit regular suites to dataset IDs
  --methods LIST          Comma-separated methods to run
  --index-type TYPE       Index layout: isax (default) or trie
  --trie-mbr-dims N       Trie MBR/split dimensions (default: dataset coefficient count)
  --trie-fanout 2|4|8      Trie symbolic split fanout (default: 8)
  --trie-dynamic-alphabet Use one global variance-weighted alphabet allocation
  --trie-min-fanout N     Minimum dynamic trie fanout (default: 2)
  --trie-max-fanout N     Maximum dynamic trie fanout (default: 16)
  --trie-alphabet-budget-bits N
                          Average dynamic alphabet budget in bits (default: 3)
  --no-dynamic-root-split-variance
                          Disable variance-assigned root bits for learned iSAX methods
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

THREADS_CSV=36
THREADS_SET=false
K_VALUES_CSV=20,50
SAMPLE_FACTORS_CSV=0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5
DATASETS_CSV=
METHODS_OVERRIDE=
INDEX_TYPE=isax
TRIE_FANOUT=8
TRIE_MBR_DIMS=
TRIE_DYNAMIC_ALPHABET=false
TRIE_MIN_FANOUT=2
TRIE_MAX_FANOUT=16
TRIE_ALPHABET_BUDGET_BITS=3
NO_DYNAMIC_ROOT_SPLIT_VARIANCE=false
RERUN_EXISTING=false
DRY_RUN=false
RESULTS_ROOT=${MESSI_RESULTS_ROOT:-"$HOME/MESSI_SFA_logs"}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) [[ $# -ge 2 ]] || die "$1 requires a value"; THREADS_CSV=$2; THREADS_SET=true; shift 2 ;;
        --k-values) [[ $# -ge 2 ]] || die "$1 requires a value"; K_VALUES_CSV=$2; shift 2 ;;
        --sample-factors) [[ $# -ge 2 ]] || die "$1 requires a value"; SAMPLE_FACTORS_CSV=$2; shift 2 ;;
        --datasets) [[ $# -ge 2 ]] || die "$1 requires a value"; DATASETS_CSV=$2; shift 2 ;;
        --methods) [[ $# -ge 2 ]] || die "$1 requires a value"; METHODS_OVERRIDE=$2; shift 2 ;;
        --index-type) [[ $# -ge 2 ]] || die "$1 requires a value"; INDEX_TYPE=$2; shift 2 ;;
        --trie-mbr-dims) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MBR_DIMS=$2; shift 2 ;;
        --trie-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_FANOUT=$2; shift 2 ;;
        --trie-dynamic-alphabet) TRIE_DYNAMIC_ALPHABET=true; shift ;;
        --trie-min-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MIN_FANOUT=$2; shift 2 ;;
        --trie-max-fanout) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_MAX_FANOUT=$2; shift 2 ;;
        --trie-alphabet-budget-bits) [[ $# -ge 2 ]] || die "$1 requires a value"; TRIE_ALPHABET_BUDGET_BITS=$2; shift 2 ;;
        --no-dynamic-root-split-variance) NO_DYNAMIC_ROOT_SPLIT_VARIANCE=true; shift ;;
        --rerun-existing) RERUN_EXISTING=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option '$1'" ;;
    esac
done

if [[ $SUITE == noise-workloads && $THREADS_SET == false ]]; then THREADS_CSV=9,18,36; fi

case "$SUITE" in
    standard|high-frequency|knn|sampling|generated-queries|hard-queries|noise-workloads) ;;
    *) die "unknown suite '$SUITE'" ;;
esac
[[ $INDEX_TYPE == isax || $INDEX_TYPE == trie ]] || die '--index-type must be isax or trie'
[[ $TRIE_FANOUT == 2 || $TRIE_FANOUT == 4 || $TRIE_FANOUT == 8 ]] || \
    die '--trie-fanout must be 2, 4, or 8'
[[ $INDEX_TYPE == trie || $TRIE_FANOUT == 8 ]] || \
    die '--trie-fanout requires --index-type trie'
[[ -z $TRIE_MBR_DIMS || $INDEX_TYPE == trie ]] || \
    die '--trie-mbr-dims requires --index-type trie'
[[ $TRIE_DYNAMIC_ALPHABET == false || $INDEX_TYPE == trie ]] || \
    die '--trie-dynamic-alphabet requires --index-type trie'
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
    local -a command=("$SCRIPT_DIR/run_dataset.sh" "$dataset" "$profile" --threads "$threads" --queue-number "$threads" --index-type "$INDEX_TYPE")
    if [[ $INDEX_TYPE == trie ]]; then
        [[ -n $TRIE_MBR_DIMS ]] && command+=(--trie-mbr-dims "$TRIE_MBR_DIMS")
        if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
            command+=(--trie-dynamic-alphabet --trie-min-fanout "$TRIE_MIN_FANOUT"
                      --trie-max-fanout "$TRIE_MAX_FANOUT"
                      --trie-alphabet-budget-bits "$TRIE_ALPHABET_BUDGET_BITS")
        else
            command+=(--trie-fanout "$TRIE_FANOUT")
        fi
    fi
    $NO_DYNAMIC_ROOT_SPLIT_VARIANCE && command+=(--no-dynamic-root-split-variance)
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
            local -a command=("$SCRIPT_DIR/run_dataset.sh" "$dataset" standard --threads "$threads" --queue-number "$threads" --index-type "$INDEX_TYPE")
            if [[ $INDEX_TYPE == trie ]]; then
                [[ -n $TRIE_MBR_DIMS ]] && command+=(--trie-mbr-dims "$TRIE_MBR_DIMS")
                if [[ $TRIE_DYNAMIC_ALPHABET == true ]]; then
                    command+=(--trie-dynamic-alphabet --trie-min-fanout "$TRIE_MIN_FANOUT"
                              --trie-max-fanout "$TRIE_MAX_FANOUT"
                              --trie-alphabet-budget-bits "$TRIE_ALPHABET_BUDGET_BITS")
                else
                    command+=(--trie-fanout "$TRIE_FANOUT")
                fi
            fi
            $NO_DYNAMIC_ROOT_SPLIT_VARIANCE && command+=(--no-dynamic-root-split-variance)
            command+=(--query-file "$query" --methods "$methods" --no-tight-bound)
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
