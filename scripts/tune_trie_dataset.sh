#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
# shellcheck source=lib/datasets.sh
source "$SCRIPT_DIR/lib/datasets.sh"

usage() {
    cat <<'USAGE'
Usage: tune_trie_dataset.sh DATASET [OPTIONS]

Tune a SPARTAN trie independently for spartan-depth and spartan-width.
Completed runs are reused. If a run directory exists without a completion
marker, the script stops rather than overwriting possibly useful output.

Options:
  --threads N                 Worker threads (default: physical cores)
  --repeats N                 Final query repetitions per built finalist (default: 5)
  --output-root PATH          Output root (default: ./trie-tuning)
  --binary PATH               MESSI executable passed to run_suite.sh
  --data-root PATH            Main dataset root
  --query-root PATH           Main query root
  --seisbench-root PATH       SeisBench dataset root
  --seisbench-query-root PATH SeisBench query root
  --dataset-file PATH         Override the base-vector file
  --query-file PATH           Override the query file
  --dataset-size N            Override indexed record count
  --query-size N              Override query count
  --numa MODE                 NUMA mode passed to run_suite.sh (default: auto)
  --sample-size N             Override SPARTAN binning sample size
  --sampling-seed N           Binning sampling seed (default: runner default)
  -h, --help                  Show this help

Outputs under OUTPUT_ROOT/DATASET:
  all-runs.tsv/csv            One row per completed run
  METHOD/final-ranking.tsv    Repeated finalist statistics
  METHOD/best-config.env      Machine-readable winning configuration
  METHOD/best-command.sh      Reproducible run_suite.sh command
  best-configs.txt            Human-readable winners for both methods
USAGE
}

die() { printf 'Error: %s\n' "$*" >&2; exit 2; }
is_positive_integer() { [[ $1 =~ ^[1-9][0-9]*$ ]]; }
absolute_path() {
    local path=$1
    if [[ $path == /* ]]; then printf '%s\n' "$path"; return; fi
    printf '%s/%s\n' "$ORIGINAL_CWD" "$path"
}

[[ $# -ge 1 ]] || { usage >&2; exit 2; }
case "$1" in -h|--help) usage; exit 0 ;; esac

ORIGINAL_CWD=$PWD
DATASET_INPUT=$1
shift
THREADS=$(physical_core_count) || die 'unable to detect physical cores; pass --threads N'
REPEATS=5
OUTPUT_ROOT=$ORIGINAL_CWD/trie-tuning
NUMA_MODE=auto
RESOLVED_PATH=
declare -a RUNNER_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) [[ $# -ge 2 ]] || die "$1 requires a value"; THREADS=$2; shift 2 ;;
        --repeats) [[ $# -ge 2 ]] || die "$1 requires a value"; REPEATS=$2; shift 2 ;;
        --output-root) [[ $# -ge 2 ]] || die "$1 requires a value"; OUTPUT_ROOT=$(absolute_path "$2"); shift 2 ;;
        --binary)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            RESOLVED_PATH=$(absolute_path "$2")
            [[ -x $RESOLVED_PATH ]] || die "MESSI executable is not executable: $RESOLVED_PATH"
            RUNNER_ARGS+=(--binary "$RESOLVED_PATH")
            shift 2
            ;;
        --data-root|--query-root|--seisbench-root|--seisbench-query-root)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            RESOLVED_PATH=$(absolute_path "$2")
            [[ -d $RESOLVED_PATH ]] || die "directory does not exist: $RESOLVED_PATH"
            RUNNER_ARGS+=("$1" "$RESOLVED_PATH")
            shift 2
            ;;
        --dataset-file|--query-file)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            RESOLVED_PATH=$(absolute_path "$2")
            [[ -f $RESOLVED_PATH ]] || die "file does not exist: $RESOLVED_PATH"
            RUNNER_ARGS+=("$1" "$RESOLVED_PATH")
            shift 2
            ;;
        --dataset-size|--query-size|--sample-size|--sampling-seed)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            RUNNER_ARGS+=("$1" "$2")
            shift 2
            ;;
        --numa) [[ $# -ge 2 ]] || die "$1 requires a value"; NUMA_MODE=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option '$1'" ;;
    esac
done

is_positive_integer "$THREADS" || die '--threads must be a positive integer'
is_positive_integer "$REPEATS" || die '--repeats must be a positive integer'
load_dataset "$DATASET_INPUT" standard || exit 2
DATASET=$DATASET_ID
MAX_MBR_DIMS=$TS_SIZE
(( MAX_MBR_DIMS > 128 )) && MAX_MBR_DIMS=128
DATASET_ROOT=$OUTPUT_ROOT/$DATASET
DATASET_ROOT_EXISTED=false
[[ ! -e $DATASET_ROOT ]] || DATASET_ROOT_EXISTED=true
mkdir -p -- "$DATASET_ROOT"
cd -- "$REPO_ROOT"

MANIFEST=$DATASET_ROOT/tuning-invocation.txt
MANIFEST_TMP=$DATASET_ROOT/.tuning-invocation.tmp
{
    printf 'format=2\ndataset=%q\nthreads=%q\nrepeats=%q\nnuma=%q\n' \
        "$DATASET" "$THREADS" "$REPEATS" "$NUMA_MODE"
    printf 'runner_args='
    printf ' %q' "${RUNNER_ARGS[@]}"
    printf '\n'
} > "$MANIFEST_TMP"
if [[ -f $MANIFEST ]]; then
    if ! cmp -s "$MANIFEST_TMP" "$MANIFEST"; then
        rm -f -- "$MANIFEST_TMP"
        die "existing tuning directory used different options: $DATASET_ROOT (choose another --output-root)"
    fi
elif [[ $DATASET_ROOT_EXISTED == true ]]; then
    rm -f -- "$MANIFEST_TMP"
    die "existing directory is not a trie-tuner output directory: $DATASET_ROOT"
else
    mv -- "$MANIFEST_TMP" "$MANIFEST"
fi
rm -f -- "$MANIFEST_TMP"

SUMMARY_HEADER=$'dataset\tmethod\tphase\tconfig_id\trepeat\tleaf_size\tn_segments\tsplit_dims\tmbr_dims\tfanout\tivf_groups\tivf_min_size\tresidual\tresidual_order\tquery_s\tbuild_s\texact_per_query\tlower_bounds_per_query\teligible_leaves\tclusters\trun_dir\tradial_mode'

rebuild_summary() {
    local summary_tmp=$DATASET_ROOT/.all-runs.tsv.tmp metrics
    printf '%s\n' "$SUMMARY_HEADER" > "$summary_tmp"
    while IFS= read -r metrics; do
        sed -n '2,$p' "$metrics" >> "$summary_tmp"
    done < <(find "$DATASET_ROOT" -type f -name metrics.tsv -print | LC_ALL=C sort)
    mv -- "$summary_tmp" "$DATASET_ROOT/all-runs.tsv"
    awk 'BEGIN { FS="\t"; OFS="," } { for (i=1; i<=NF; ++i) { gsub(/"/, "\"\"", $i); $i="\"" $i "\"" } print }' \
        "$DATASET_ROOT/all-runs.tsv" > "$DATASET_ROOT/all-runs.csv"
}

human_count() {
    awk -v text="$1" 'BEGIN {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", text)
        n = split(text, p, /[[:space:]]+/)
        value = p[1]; suffix = p[2]
        sub(/\/query$/, "", value); sub(/\/query$/, "", suffix)
        multiplier = 1
        if (suffix == "K") multiplier = 1000
        else if (suffix == "M") multiplier = 1000000
        else if (suffix == "G") multiplier = 1000000000
        printf "%.0f", value * multiplier
    }'
}

parse_log() {
    local log=$1 value
    PARSED_QUERY_S=$(awk '/^>>> query wall time:/ { value=$5 } END { if (value != "") print value }' "$log")
    PARSED_REPEAT_TIMES=$(awk '/^>>> query repeat [0-9]+\/[0-9]+ wall time:/ { print $7 }' "$log")
    PARSED_BUILD_S=$(awk '
        /^>>> trie build timing/ { in_build=1; next }
        /^>>>/ && in_build { in_build=0 }
        in_build && /^[[:space:]]+total[[:space:]]*:/ { value=$3 }
        END { if (value != "") print value }
    ' "$log")
    value=$(awk '/^  exact distances[[:space:]]*:/ {
        line=$0; sub(/^.*:[[:space:]]*/, "", line); sub(/[[:space:]]*\(.*/, "", line); value=line
    } END { print value }' "$log")
    PARSED_EXACT=$(human_count "$value")
    value=$(awk '/^  symbolic record bounds[[:space:]]*:/ {
        line=$0; sub(/^.*:[[:space:]]*/, "", line); sub(/[[:space:]]*\(.*/, "", line); value=line
    } END { print value }' "$log")
    PARSED_LOWER=$(human_count "$value")
    PARSED_ELIGIBLE=$(awk '/^[[:space:]]+eligible leaves[[:space:]]*:/ { value=$4 } END { print value+0 }' "$log")
    PARSED_CLUSTERS=$(awk '/^[[:space:]]+clusters[[:space:]]*:/ { value=$3 } END { print value+0 }' "$log")
    [[ -n $PARSED_QUERY_S && -n $PARSED_BUILD_S ]] || return 1
}

print_config() {
    printf 'LEAF_SIZE=%q\n' "$LEAF_SIZE"
    printf 'N_SEGMENTS=%q\n' "$N_SEGMENTS"
    printf 'SPLIT_DIMS=%q\n' "$SPLIT_DIMS"
    printf 'MBR_DIMS=%q\n' "$MBR_DIMS"
    printf 'FANOUT=%q\n' "$FANOUT"
    printf 'IVF_GROUPS=%q\n' "$IVF_GROUPS"
    printf 'IVF_MIN_SIZE=%q\n' "$IVF_MIN_SIZE"
    printf 'RESIDUAL=%q\n' "$RESIDUAL"
    printf 'RESIDUAL_ORDER=%q\n' "$RESIDUAL_ORDER"
    printf 'RADIAL_MODE=%q\n' "$RADIAL_MODE"
}

write_config() {
    print_config > "$1"
}

run_config() {
    local method=$1 phase=$2 config_id=$3 query_repeats=${4:-1}
    local run_dir=$DATASET_ROOT/$method/$phase/$config_id
    local complete=$run_dir/.complete log=$run_dir/run.log
    local -a command=("$SCRIPT_DIR/run_suite.sh" standard
        --threads "$THREADS" --numa "$NUMA_MODE" --datasets "$DATASET"
        --index-type trie --methods "$method" --leaf-size "$LEAF_SIZE"
        --n-segments "$N_SEGMENTS" --trie-split-dims "$SPLIT_DIMS"
        --trie-mbr-dims "$MBR_DIMS" --trie-fanout "$FANOUT")

    if [[ -f $complete && -f $run_dir/metrics.tsv && -f $run_dir/config.env ]]; then
        if [[ $(print_config) != "$(<"$run_dir/config.env")" ]]; then
            die "completed run has a conflicting configuration: $run_dir"
        fi
        printf 'Reusing completed run: %s\n' "$run_dir"
        return 0
    fi
    if [[ -e $run_dir ]]; then
        die "incomplete run directory exists: $run_dir (move or remove it after inspection, then rerun)"
    fi
    mkdir -p -- "$(dirname -- "$run_dir")"
    if ! mkdir -- "$run_dir"; then
        die "run directory was created concurrently: $run_dir"
    fi
    write_config "$run_dir/config.env"

    if (( IVF_GROUPS == 0 )); then
        command+=(--no-trie-leaf-ivf)
    else
        command+=(--trie-leaf-ivf "$IVF_GROUPS" --trie-leaf-ivf-min-size "$IVF_MIN_SIZE")
        case "$RADIAL_MODE" in
            on) command+=(--trie-leaf-ivf-radial-bound) ;;
            auto) command+=(--trie-leaf-ivf-radial-bound-auto) ;;
            off) command+=(--no-trie-leaf-ivf-radial-bound) ;;
            *) die "invalid internal radial mode '$RADIAL_MODE'" ;;
        esac
    fi
    if [[ $RESIDUAL == true ]]; then
        command+=(--trie-residual-record-only --trie-residual-order "$RESIDUAL_ORDER")
    fi
    command+=(--query-repeats "$query_repeats")
    command+=("${RUNNER_ARGS[@]}")

    printf '\n[%s/%s] %s\n' "$method" "$phase" "$config_id"
    printf 'Command:'
    printf ' %q' "${command[@]}"
    printf '\n'
    if ! MESSI_RESULTS_ROOT="$run_dir/results" \
         MESSI_LOG_ROOT="$run_dir/live-logs" \
         MESSI_SHELL_LOG_DIR="$run_dir/live-logs" \
         "${command[@]}" 2>&1 | tee "$log"; then
        die "run failed; preserved incomplete output at $run_dir"
    fi
    if ! parse_log "$log"; then
        die "run completed but its timings could not be parsed; preserved output at $run_dir"
    fi
    local parsed_repeat_count
    parsed_repeat_count=$(awk '/^>>> query repeat [0-9]+\/[0-9]+ wall time:/ { count++ } END { print count+0 }' "$log")
    if (( query_repeats == 1 && parsed_repeat_count == 0 )); then
        PARSED_REPEAT_TIMES=$PARSED_QUERY_S
        parsed_repeat_count=1
    fi
    (( parsed_repeat_count == query_repeats )) ||
        die "expected $query_repeats query timings but parsed $parsed_repeat_count; preserved output at $run_dir"
    {
        printf '%s\n' "$SUMMARY_HEADER"
        local timing timing_repeat=0
        while IFS= read -r timing; do
            (( ++timing_repeat ))
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$DATASET" "$method" "$phase" "$config_id" "$timing_repeat" \
                "$LEAF_SIZE" "$N_SEGMENTS" "$SPLIT_DIMS" "$MBR_DIMS" "$FANOUT" \
                "$IVF_GROUPS" "$IVF_MIN_SIZE" "$RESIDUAL" "$RESIDUAL_ORDER" \
                "$timing" "$PARSED_BUILD_S" "$PARSED_EXACT" "$PARSED_LOWER" \
                "$PARSED_ELIGIBLE" "$PARSED_CLUSTERS" "$run_dir" "$RADIAL_MODE"
        done <<< "$PARSED_REPEAT_TIMES"
    } > "$run_dir/metrics.tsv"
    : > "$complete"
    rebuild_summary
}

best_run_dir() {
    local method=$1 phase=$2 rank=${3:-1}
    find "$DATASET_ROOT/$method/$phase" -type f -name metrics.tsv -print 2>/dev/null |
        while IFS= read -r file; do sed -n '2p' "$file"; done |
        LC_ALL=C sort -t $'\t' -k15,15g |
        sed -n "${rank}p" | awk -F '\t' '{ print $21 }'
}

load_run_config() {
    local run_dir=$1
    [[ -n $run_dir && -f $run_dir/config.env ]] || die 'could not locate the selected configuration'
    # This file is generated by write_config and contains only quoted scalar assignments.
    # shellcheck disable=SC1090
    source "$run_dir/config.env"
}

unique_mbr_values() {
    local prefix=$1 candidate seen=' '
    for candidate in 64 96 128; do
        (( candidate < prefix )) && candidate=$prefix
        (( candidate > MAX_MBR_DIMS )) && candidate=$MAX_MBR_DIMS
        [[ $seen == *" $candidate "* ]] && continue
        printf '%s\n' "$candidate"
        seen+="$candidate "
    done
}

rank_finalists() {
    local method=$1 phase_dir=$DATASET_ROOT/$method/06-final-query-only
    local output=$DATASET_ROOT/$method/final-ranking.tsv
    {
        printf 'candidate\truns\tmedian_query_s\tmean_query_s\tmin_query_s\tmax_query_s\tconfig_dir\n'
        find "$phase_dir" -type f -name metrics.tsv -print |
            while IFS= read -r file; do sed -n '2,$p' "$file"; done |
            awk -F '\t' '
                function sort_values(a, n, i, j, x) {
                    for (i=2; i<=n; ++i) { x=a[i]; j=i-1; while (j>=1 && a[j]>x) { a[j+1]=a[j]; --j } a[j+1]=x }
                }
                {
                    key=$4; n[key]++; values[key,n[key]]=$15+0; sum[key]+=$15
                    if (!(key in min) || $15<min[key]) min[key]=$15
                    if (!(key in max) || $15>max[key]) max[key]=$15
                    dir[key]=$21
                }
                END {
                    for (key in n) {
                        split("", a)
                        for (i=1; i<=n[key]; ++i) a[i]=values[key,i]
                        sort_values(a,n[key])
                        if (n[key]%2) median=a[(n[key]+1)/2]
                        else median=(a[n[key]/2]+a[n[key]/2+1])/2
                        printf "%s\t%d\t%.9f\t%.9f\t%.9f\t%.9f\t%s\n", key,n[key],median,sum[key]/n[key],min[key],max[key],dir[key]
                    }
                }' | LC_ALL=C sort -t $'\t' -k3,3g
    } > "$output"
}

write_recommendation() {
    local method=$1 ranking=$DATASET_ROOT/$method/final-ranking.tsv
    local winner_dir median
    winner_dir=$(sed -n '2p' "$ranking" | awk -F '\t' '{ print $7 }')
    median=$(sed -n '2p' "$ranking" | awk -F '\t' '{ print $3 }')
    load_run_config "$winner_dir"
    {
        print_config
        printf 'DATASET=%q\nMETHOD=%q\nMEDIAN_QUERY_S=%q\nREPEATS=%q\n' \
            "$DATASET" "$method" "$median" "$REPEATS"
    } > "$DATASET_ROOT/$method/best-config.env"

    local -a command=("$SCRIPT_DIR/run_suite.sh" standard
        --threads "$THREADS" --numa "$NUMA_MODE" --datasets "$DATASET"
        --index-type trie --methods "$method" --leaf-size "$LEAF_SIZE"
        --n-segments "$N_SEGMENTS" --trie-split-dims "$SPLIT_DIMS"
        --trie-mbr-dims "$MBR_DIMS" --trie-fanout "$FANOUT")
    if (( IVF_GROUPS == 0 )); then command+=(--no-trie-leaf-ivf)
    else
        command+=(--trie-leaf-ivf "$IVF_GROUPS" --trie-leaf-ivf-min-size "$IVF_MIN_SIZE")
        case "$RADIAL_MODE" in
            on) command+=(--trie-leaf-ivf-radial-bound) ;;
            auto) command+=(--trie-leaf-ivf-radial-bound-auto) ;;
            off) command+=(--no-trie-leaf-ivf-radial-bound) ;;
        esac
    fi
    [[ $RESIDUAL == false ]] || command+=(--trie-residual-record-only --trie-residual-order "$RESIDUAL_ORDER")
    command+=("${RUNNER_ARGS[@]}")
    {
        printf '#!/usr/bin/env bash\nset -Eeuo pipefail\n'
        printf 'cd %q\n' "$REPO_ROOT"
        printf 'exec'
        printf ' %q' "${command[@]}"
        printf ' "$@"\n'
    } > "$DATASET_ROOT/$method/best-command.sh"
    chmod +x "$DATASET_ROOT/$method/best-command.sh"
}

tune_method() {
    local method=$1 leaf prefix fanout mbr split groups min_size order radial candidate rank selected
    local -a radial_values=() finalist_dirs=()
    printf '\n=== Tuning %s on %s ===\n' "$method" "$DATASET"

    # Stage 1: structural parameters, with the current IVF baseline and no residual.
    for leaf in 10000 20000 40000; do
        for prefix in 32 48 64; do
            for fanout in 4 8; do
                LEAF_SIZE=$leaf; N_SEGMENTS=$prefix; SPLIT_DIMS=$prefix
                MBR_DIMS=$MAX_MBR_DIMS; FANOUT=$fanout
                IVF_GROUPS=16; IVF_MIN_SIZE=4096
                RESIDUAL=false; RESIDUAL_ORDER=symbolic-first
                RADIAL_MODE=on
                run_config "$method" 01-structure "leaf-${leaf}-prefix-${prefix}-fanout-${fanout}"
            done
        done
    done

    selected=$(best_run_dir "$method" 01-structure)
    load_run_config "$selected"
    local base_leaf=$LEAF_SIZE base_prefix=$N_SEGMENTS base_split=$SPLIT_DIMS base_fanout=$FANOUT

    # Stage 2: MBR width and the number of dimensions eligible for splitting.
    while IFS= read -r mbr; do
        for split in "$base_prefix" "$mbr"; do
            LEAF_SIZE=$base_leaf; N_SEGMENTS=$base_prefix; SPLIT_DIMS=$split
            MBR_DIMS=$mbr; FANOUT=$base_fanout
            IVF_GROUPS=16; IVF_MIN_SIZE=4096
            RESIDUAL=false; RESIDUAL_ORDER=symbolic-first; RADIAL_MODE=on
            run_config "$method" 02-mbr "mbr-${mbr}-split-${split}"
            [[ $mbr != "$base_prefix" ]] || break
        done
    done < <(unique_mbr_values "$base_prefix")

    selected=$(best_run_dir "$method" 02-mbr)
    load_run_config "$selected"
    local base_mbr=$MBR_DIMS
    base_split=$SPLIT_DIMS

    # Stage 3: choose IVF group count at the default 4096-record threshold.
    for groups in 0 8 16 32; do
        LEAF_SIZE=$base_leaf; N_SEGMENTS=$base_prefix; SPLIT_DIMS=$base_split
        MBR_DIMS=$base_mbr; FANOUT=$base_fanout
        IVF_GROUPS=$groups; IVF_MIN_SIZE=4096
        RESIDUAL=false; RESIDUAL_ORDER=symbolic-first; RADIAL_MODE=on
        run_config "$method" 03-ivf-groups "groups-$groups"
    done

    selected=$(best_run_dir "$method" 03-ivf-groups)
    load_run_config "$selected"
    local base_groups=$IVF_GROUPS

    # Stage 4: choose the IVF eligibility threshold. For IVF-off, retain one control.
    if (( base_groups == 0 )); then
        LEAF_SIZE=$base_leaf; N_SEGMENTS=$base_prefix; SPLIT_DIMS=$base_split
        MBR_DIMS=$base_mbr; FANOUT=$base_fanout
        IVF_GROUPS=0; IVF_MIN_SIZE=4096
        RESIDUAL=false; RESIDUAL_ORDER=symbolic-first; RADIAL_MODE=off
        run_config "$method" 04-ivf-min-size ivf-off
    else
        for min_size in 2048 4096 8192; do
            LEAF_SIZE=$base_leaf; N_SEGMENTS=$base_prefix; SPLIT_DIMS=$base_split
            MBR_DIMS=$base_mbr; FANOUT=$base_fanout
            IVF_GROUPS=$base_groups; IVF_MIN_SIZE=$min_size
            RESIDUAL=false; RESIDUAL_ORDER=symbolic-first; RADIAL_MODE=on
            run_config "$method" 04-ivf-min-size "min-$min_size"
        done
    fi

    selected=$(best_run_dir "$method" 04-ivf-min-size)
    load_run_config "$selected"
    local base_min=$IVF_MIN_SIZE

    # Stage 5: residual order and, when IVF is active, radial-filter policy.
    if (( base_groups == 0 )); then radial_values=(off); else radial_values=(on auto off); fi
    for radial in "${radial_values[@]}"; do
        for order in off symbolic-first residual-first; do
            LEAF_SIZE=$base_leaf; N_SEGMENTS=$base_prefix; SPLIT_DIMS=$base_split
            MBR_DIMS=$base_mbr; FANOUT=$base_fanout
            IVF_GROUPS=$base_groups; IVF_MIN_SIZE=$base_min; RADIAL_MODE=$radial
            if [[ $order == off ]]; then RESIDUAL=false; RESIDUAL_ORDER=symbolic-first
            else RESIDUAL=true; RESIDUAL_ORDER=$order
            fi
            run_config "$method" 05-residual "radial-${radial}-residual-${order}"
        done
    done

    # Stage 6: build each of the two fastest residual-stage configurations once,
    # then repeat only its query workload against the same in-memory trie.
    for rank in 1 2; do
        finalist_dirs[$rank]=$(best_run_dir "$method" 05-residual "$rank")
    done
    for rank in 1 2; do
        load_run_config "${finalist_dirs[$rank]}"
        candidate=$(printf 'candidate-%s' "$rank")
        run_config "$method" 06-final-query-only "$candidate" "$REPEATS"
    done
    rank_finalists "$method"
    write_recommendation "$method"
}

if [[ -f $DATASET_ROOT/.complete ]]; then
    rebuild_summary
    printf 'Tuning is already complete for dataset=%s; existing results were not overwritten.\n' "$DATASET"
else
    for method in spartan-depth spartan-width; do tune_method "$method"; done
    rebuild_summary
    : > "$DATASET_ROOT/.complete"
fi

{
    printf 'Dataset: %s\n' "$DATASET"
    for method in spartan-depth spartan-width; do
        printf '\n[%s]\n' "$method"
        if [[ -f $DATASET_ROOT/$method/best-config.env ]]; then
            sed -n '/^LEAF_SIZE=/,$p' "$DATASET_ROOT/$method/best-config.env"
            printf 'command=%s\n' "$DATASET_ROOT/$method/best-command.sh"
        else
            printf 'No completed recommendation.\n'
        fi
    done
} | tee "$DATASET_ROOT/best-configs.txt"

printf '\nAll-run summary: %s\n' "$DATASET_ROOT/all-runs.tsv"
printf 'Recommendations: %s\n' "$DATASET_ROOT/best-configs.txt"
