#!/usr/bin/env bash
# Paired CLI builds/searches on a bounded SALD fixture. Override paths/counts
# through the variables below for a server run; each run gets separate logs.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
BINARY=${MESSI_BINARY:-"$REPO_ROOT/build/bin/MESSI"}
DATASET=${RESIDUAL_DATASET:-"$REPO_ROOT/data_head/SALD_head.bin"}
QUERIES=${RESIDUAL_QUERIES:-"$REPO_ROOT/data_queries/SALD_queries.bin"}
RECORDS=${RESIDUAL_RECORDS:-20000}
LENGTH=${RESIDUAL_LENGTH:-128}
WORKERS=${RESIDUAL_WORKERS:-4}
REPEATS=${RESIDUAL_REPEATS:-3}
QUERY_MODE=${RESIDUAL_QUERY_MODE:-single}
HISTOGRAM=${RESIDUAL_HISTOGRAM:-1}
[[ $QUERY_MODE == single || $QUERY_MODE == batch ]] || { printf 'RESIDUAL_QUERY_MODE must be single or batch\n' >&2; exit 1; }
OUTPUT_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/messi-residual-bench.XXXXXX")
printf 'Logs: %s\n' "$OUTPUT_ROOT" >&2
printf 'prefix,repeat,residual,build_seconds,query_seconds,exact_per_query,record_checks,record_wins,record_prunes\n'
for prefix in 16 32 64; do
    for ((repeat=1; repeat<=REPEATS; ++repeat)); do
        for mode in 0 1; do
            run_dir="$OUTPUT_ROOT/k${prefix}-r${repeat}-m${mode}"
            mkdir -p -- "$run_dir"
            extra=()
            [[ $QUERY_MODE != batch ]] || extra+=(--trie-query-batch)
            ((mode == 0)) || extra+=(--trie-residual-record-only)
            MESSI_LOG_ROOT="$run_dir" "$BINARY" --in-memory --index-path "$run_dir/index-root" \
                --dataset "$DATASET" --dataset-size "$RECORDS" --timeseries-size "$LENGTH" \
                --queries "$QUERIES" --queries-size 100 --threads "$WORKERS" --numa none \
                --index-type trie --function-type 5 --sample-size 1000 --sample-type 2 --sampling-seed 1 \
                --is-norm --histogram-type "$HISTOGRAM" --sax-cardinality 8 --n-segments "$prefix" \
                --trie-mbr-dimensions "$((LENGTH < 128 ? LENGTH : 128))" \
                --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 \
                --trie-record-mbr-suffix-bound --trie-streaming-leaf-scan --profile-query-phases \
                --trie-leaf-ivf 16 --trie-leaf-ivf-radial-bound "${extra[@]}" \
                >"$run_dir/run.log" 2>&1
            awk -v k="$prefix" -v r="$repeat" -v m="$mode" '
                /    total      :/ { build=$3 }
                />>> query wall time:/ { query=$5 }
                END { printf "%d,%d,%d,%s,%s,",k,r,m,build,query }
            ' "$run_dir/run.log"
            # The final query CSV row contains means, including exact calls.
            awk -F, 'FNR>1 { exact=$13 } END { printf "%.3f,",exact }' "$run_dir"/query/*.csv
            awk '/>>> ResSPARTAN record:/ {
                split($4,a,"="); checks+=a[2]; split($5,a,"="); wins+=a[2]; split($6,a,"="); prunes+=a[2]
            } END { printf "%d,%d,%d\n",checks,wins,prunes }' "$run_dir/run.log"
            if ((mode == 1)); then
                baseline="$OUTPUT_ROOT/k${prefix}-r${repeat}-m0"
                awk -F, 'FNR==NR { if(FNR>1) reference[FNR]=$16; next }
                    FNR>1 && $16 != reference[FNR] { bad=1 }
                    END { if(bad) { print "Paired query distances differ" > "/dev/stderr"; exit 1 } }
                ' "$baseline"/query/*.csv "$run_dir"/query/*.csv
            fi
        done
    done
done
