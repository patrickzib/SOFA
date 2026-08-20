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
pass 'standard profile emits the complete method matrix with quoted paths'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type trie --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 7 ]] || fail 'trie standard profile should emit all supported methods'
assert_contains "$OUTPUT" '--function-type 3'
assert_contains "$OUTPUT" '--function-type 4'
assert_contains "$OUTPUT" '--function-type 5'
assert_contains "$OUTPUT" '--function-type 6'
pass 'trie standard profile emits the complete method matrix'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann high-frequency --threads 36 --queue-number 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--apply-z-norm'
assert_contains "$OUTPUT" '--filetype-int'
assert_contains "$OUTPUT" '--queries-size 1'
assert_contains "$OUTPUT" '--sfa-n-coefficients 50'
assert_contains "$OUTPUT" '--histogram-type 2'
assert_not_contains "$OUTPUT" '--histogram-type 1'
pass 'high-frequency profile preserves BigANN flags and coefficients'

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

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" generated-queries --threads 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" 'spacev1B_noise_025.bin'
assert_contains "$OUTPUT" 'text-to-image_noise_01.bin'
assert_contains "$OUTPUT" 'turingANNs_noise_05.bin'
pass 'migrated generated-query suite expands expected workloads'

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
