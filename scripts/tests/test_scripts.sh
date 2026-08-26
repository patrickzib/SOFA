#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
TEST_COUNT=0

fail() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }
pass() { TEST_COUNT=$((TEST_COUNT + 1)); printf 'ok %d - %s\n' "$TEST_COUNT" "$1"; }
assert_contains() { [[ $1 == *"$2"* ]] || fail "expected output to contain: $2"; }
assert_not_contains() { [[ $1 != *"$2"* ]] || fail "expected output not to contain: $2"; }

while IFS= read -r script; do bash -n "$script"; done < <(find "$SCRIPT_DIR" -type f -name '*.sh' -not -path '*/old/*' -print)
pass 'all maintained scripts pass bash syntax validation'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --data-root '/tmp/data root' --binary /tmp/MESSI --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 7 ]] || fail 'standard profile should emit seven commands'
assert_contains "$OUTPUT" '/tmp/data\ root/astro.bin'
assert_contains "$OUTPUT" '--function-type 6'
assert_contains "$OUTPUT" '--function-type 5'
assert_contains "$OUTPUT" '--tight-bound'
assert_contains "$OUTPUT" '--sample-type 2'
assert_not_contains "$OUTPUT" '--sample-type 3'
[[ $(printf '%s\n' "$OUTPUT" | grep -c -- '--dynamic-root-split-variance') == 6 ]] || fail 'iSAX should enable variance root splits for learned methods by default'
pass 'standard profile emits the complete method matrix with uniform binning samples'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --data-root '/tmp/data root' --binary /tmp/MESSI --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--queue-number'
pass 'queue count is optional and defaults in MESSI'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type trie --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 6 ]] || fail 'trie standard profile should exclude SAX'
assert_not_contains "$OUTPUT" '--function-type 3'
assert_contains "$OUTPUT" '--function-type 4'
assert_contains "$OUTPUT" '--function-type 5'
assert_contains "$OUTPUT" '--function-type 6'
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_contains "$OUTPUT" '--trie-split-dimensions 32'
assert_contains "$OUTPUT" '--trie-fanout 8'
pass 'trie standard profile excludes SAX from its method matrix'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --n-segments 32 --trie-split-dims 32 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 100'
assert_contains "$OUTPUT" '--n-segments 32'
assert_contains "$OUTPUT" '--trie-split-dimensions 32'
if "$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 64 --trie-split-dims 65 --dry-run >/dev/null 2>&1; then
    fail 'trie benchmark accepted split dimensions wider than its MBR dimensions'
fi
pass 'trie split candidates are independent from MBR dimensions'

if "$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 32 --n-segments 48 --dry-run >/dev/null 2>&1; then
    fail 'trie benchmark accepted a record prefix wider than its MBR dimensions'
fi
pass 'trie record-prefix dimensions are forwarded and validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie --trie-fanout 4 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-fanout 4'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type trie --trie-fanout 3 --dry-run >/dev/null 2>&1; then
    fail 'trie benchmark accepted an invalid fanout'
fi
pass 'trie fanout is forwarded and validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type isax --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | grep -c -- '--dynamic-root-split-variance') == 6 ]] || fail 'variance root split should apply to each learned iSAX method'
pass 'iSAX variance root split is forwarded only to learned transforms'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type isax --no-dynamic-root-split-variance --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--dynamic-root-split-variance'
pass 'iSAX variance root split can be disabled explicitly'

if "$SCRIPT_DIR/run_dataset.sh" astro standard --threads 1 --queue-number 1 --index-type trie --methods sax --dry-run >/dev/null 2>&1; then
    fail 'trie benchmark accepted SAX explicitly'
fi
pass 'trie benchmark rejects SAX explicitly'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --queue-number 36 --index-type trie --trie-query-parallel --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-query-parallel'
pass 'trie query-parallel option is forwarded'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-record-mbr-suffix-bound --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-record-mbr-suffix-bound'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-record-mbr-suffix-bound --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie record/MBR suffix bound for iSAX'
fi
pass 'trie record/MBR suffix bound is forwarded and scoped to trie'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-leaf-kmeans 16 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-kmeans 16'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-leaf-kmeans 16 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie leaf k-means for iSAX'
fi
pass 'trie leaf k-means is forwarded and scoped to trie'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --query-report-interval 10 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--query-report-interval 10'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --query-report-interval -1 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted a negative query report interval'
fi
pass 'query report interval is forwarded and validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann high-frequency --threads 36 --queue-number 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--apply-z-norm'
assert_contains "$OUTPUT" '--filetype-int'
assert_contains "$OUTPUT" '--queries-size 1'
assert_contains "$OUTPUT" '--sfa-n-coefficients 50'
assert_contains "$OUTPUT" '--histogram-type 2'
assert_not_contains "$OUTPUT" '--histogram-type 1'
pass 'high-frequency profile preserves BigANN flags and coefficients'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --queue-number 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sfa-n-coefficients 64'
assert_contains "$OUTPUT" '--sampling-seed 1'
pass 'standard SFA uses the 64-coefficient training pool when permitted'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" deep1b standard --threads 1 --queue-number 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sfa-n-coefficients 48'
pass 'short series use the largest valid even coefficient pool'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 64 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 64'
assert_contains "$OUTPUT" '--sfa-n-coefficients 64'
pass 'trie MBR dimensions are capped by series length, not half length'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_contains "$OUTPUT" '--sfa-n-coefficients 128'
pass 'trie supports 128 MBR dimensions for 128-value series'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 100'
assert_contains "$OUTPUT" '--sfa-n-coefficients 100'
pass 'trie MBR dimensions remain capped by the 100-value series length'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --queue-number 1 --dataset-size 100k --dry-run 2>&1)
assert_contains "$OUTPUT" 'exceeds dataset size 100 K; using dataset size'
assert_contains "$OUTPUT" '--sample-size 100000'
pass 'binning sample size is capped for reduced datasets'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sift1b knn --threads 36 --queue-number 36 --k 20 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--topk --k-size 20'
assert_contains "$OUTPUT" '--histogram-type 1'
assert_contains "$OUTPUT" '--tight-bound'
pass 'KNN profile emits top-K and SIFT-specific options'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" seisbench sampling --threads 36 --queue-number 36 --dataset-file sample.bin --query-file queries.bin --dataset-size 1000 --sample-factor 0.25 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sample-size 250'
pass 'SeisBench sampling factor is wired into sample-size calculation'

OUTPUT=$(cd /tmp && MESSI_DRY_RUN=true "$SCRIPT_DIR/run_astro.sh" 9 18 2>/dev/null)
assert_contains "$OUTPUT" '--threads 9'
assert_contains "$OUTPUT" '--queue-number 18'
pass 'compatibility wrappers work outside the scripts directory'

if "$SCRIPT_DIR/run_dataset.sh" astro knn --threads 1 --queue-number 1 --dry-run >/dev/null 2>&1; then
    fail 'knn profile accepted a missing K value'
fi
if "$SCRIPT_DIR/run_dataset.sh" astro sampling --threads 1 --queue-number 1 --sample-factor 2 --dry-run >/dev/null 2>&1; then
    fail 'sampling profile accepted an out-of-range factor'
fi
pass 'invalid profile parameters are rejected'

TEMP_ROOT=$(mktemp -d)
trap 'rm -rf -- "$TEMP_ROOT"' EXIT
mkdir -p "$TEMP_ROOT/logs" "$TEMP_ROOT/results/ASTRO/36"
printf 'old\n' > "$TEMP_ROOT/results/ASTRO/36/old.log"
printf 'new\n' > "$TEMP_ROOT/logs/new.log"
MESSI_LOG_ROOT="$TEMP_ROOT/logs" MESSI_RESULTS_ROOT="$TEMP_ROOT/results" "$SCRIPT_DIR/archive_results.sh" ASTRO 36
[[ -f $TEMP_ROOT/results/ASTRO/36/new.log ]] || fail 'new result was not archived'
[[ ! -e $TEMP_ROOT/results/ASTRO/36/old.log ]] || fail 'old result was not replaced'
if MESSI_LOG_ROOT="$TEMP_ROOT/logs" MESSI_RESULTS_ROOT="$TEMP_ROOT/results" "$SCRIPT_DIR/archive_results.sh" ../escape 36 >/dev/null 2>&1; then
    fail 'archive helper accepted an escaping label'
fi
pass 'result replacement is bounded by the configured results root'

TEST_RUN_MESSI="$TEMP_ROOT/test_run_messi"
printf '%s\n' '#!/usr/bin/env bash' 'printf "%s\\n" "$*"' > "$TEST_RUN_MESSI"
chmod +x "$TEST_RUN_MESSI"
OUTPUT=$(MESSI_BIN="$TEST_RUN_MESSI" "$REPO_ROOT/test_run.sh" 4 isax)
assert_contains "$OUTPUT" '--threads 4'
assert_contains "$OUTPUT" '--index-type isax'
assert_contains "$OUTPUT" '--dynamic-root-split-variance'
pass 'test_run accepts index type as its second argument'

touch "$TEMP_ROOT/astro.bin" "$TEMP_ROOT/astro_queries.bin"
FAKE_MESSI="$TEMP_ROOT/fake_messi"
printf '%s\n' '#!/usr/bin/env bash' \
    'printf "%s\\n" "=== Query summary ===" >&2' \
    'printf "%s\\n" "  wall time        : 3.090 s (30.900 ms/query)" >&2' \
    'printf "%s\\n" "  lower bounds     : 24.35 M/query (24.35% of 100.00 M indexed series)" >&2' \
    'printf "%s\\n" "  exact distances  : 567.00 K/query (0.57% of 100.00 M indexed series)" >&2' \
    > "$FAKE_MESSI"
chmod +x "$FAKE_MESSI"
OUTPUT=$(MESSI_SHELL_LOG_DIR="$TEMP_ROOT/logs" "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --data-root "$TEMP_ROOT" \
    --binary "$FAKE_MESSI" --methods spartan-width 2>&1)
assert_contains "$OUTPUT" '=== Benchmark summary: dataset=astro, profile=high-frequency ==='
assert_contains "$OUTPUT" 'iSAX    SPARTAN   width'
assert_contains "$OUTPUT" '24.35 M/query'
[[ $OUTPUT =~ 24\.35\ M/query[[:space:]]+24\.35% ]] || fail 'suite summary columns are out of order'
pass 'runner prints a compact suite summary from query results'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" generated-queries --threads 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" 'spacev1B_noise_025.bin'
assert_contains "$OUTPUT" 'text-to-image_noise_01.bin'
assert_contains "$OUTPUT" 'turingANNs_noise_05.bin'
pass 'migrated generated-query suite expands expected workloads'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--index-type trie'
assert_not_contains "$OUTPUT" '--function-type 3'
pass 'suite forwards the selected index layout to dataset runs'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 64 --datasets astro --index-type trie \
    --methods spartan-depth,spartan-width --trie-mbr-dims 64 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 64'
pass 'suite forwards trie MBR dimensions to dataset runs'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 64 --datasets astro --index-type trie \
    --leaf-size 10k --min-leaf-size 5k --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--leaf-size 10000'
assert_contains "$OUTPUT" '--min-leaf-size 5000'
pass 'suite forwards leaf capacity and minimum occupancy'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-record-mbr-suffix-bound --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-record-mbr-suffix-bound'
pass 'suite forwards the trie record/MBR suffix bound'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-leaf-kmeans 16 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-kmeans 16'
pass 'suite forwards trie leaf k-means'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro \
    --index-type trie --trie-dynamic-alphabet --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-dynamic-alphabet'
assert_contains "$OUTPUT" '--trie-min-fanout 2'
assert_contains "$OUTPUT" '--trie-max-fanout 16'
assert_contains "$OUTPUT" '--trie-alphabet-budget-bits 3'
pass 'suite forwards the precomputed dynamic trie alphabet option'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type isax --no-dynamic-root-split-variance --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--dynamic-root-split-variance'
pass 'suite forwards the dynamic root split disable option'

mkdir -p "$TEMP_ROOT/results/ASTRO/36"
OUTPUT=$(MESSI_RESULTS_ROOT="$TEMP_ROOT/results" "$SCRIPT_DIR/run_suite.sh" standard \
    --threads 36 --datasets astro --index-type trie 2>&1)
assert_contains "$OUTPUT" 'Skipping dataset=astro profile=standard run=36'
pass 'suite skips workloads with existing result archives'

if command -v shellcheck >/dev/null 2>&1; then
    mapfile -t SHELL_SCRIPTS < <(find "$SCRIPT_DIR" -type f -name '*.sh' -not -path '*/old/*' -print)
    shellcheck -x "${SHELL_SCRIPTS[@]}"
    pass 'ShellCheck passes'
else
    printf '# SKIP ShellCheck is not installed\n'
fi

if [[ ${RUN_MESSI_INTEGRATION:-0} == 1 ]]; then
    mkdir -p "$TEMP_ROOT/home"
    HOME="$TEMP_ROOT/home" "$SCRIPT_DIR/run_dataset.sh" astro high-frequency \
        --threads 1 --queue-number 1 \
        --dataset-file "$REPO_ROOT/data_head/astro_head.bin" --dataset-size 1000 \
        --query-file "$REPO_ROOT/data_queries/astro_queries.bin" --query-size 1 \
        --sample-size 100 --binary "$REPO_ROOT/bin/MESSI"
    pass 'real MESSI fixture smoke test completes'
else
    printf '# SKIP set RUN_MESSI_INTEGRATION=1 for the real MESSI fixture test\n'
fi

printf '1..%d\n' "$TEST_COUNT"
