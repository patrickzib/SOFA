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

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type isax --data-root '/tmp/data root' --binary /tmp/MESSI --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 5 ]] || fail 'standard profile should emit five default commands'
assert_contains "$OUTPUT" '/tmp/data\ root/astro.bin'
assert_not_contains "$OUTPUT" '--function-type 6'
assert_contains "$OUTPUT" '--function-type 5'
assert_contains "$OUTPUT" '--tight-bound'
assert_contains "$OUTPUT" '--sample-type 2'
assert_not_contains "$OUTPUT" '--sample-type 3'
assert_not_contains "$OUTPUT" '--dynamic-root-split-variance'
assert_contains "$OUTPUT" '--n-segments 16'
pass 'standard profile emits the complete method matrix with uniform binning samples'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 1 --index-type isax \
    --isax-n-segments 32 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--n-segments 32'
pass 'iSAX/SAX segment count defaults to 16 and can be overridden'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" simsearchnet standard --threads 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '/home/tmp/schaefpa/messi_datasets/SimSearchNet.bin'
assert_contains "$OUTPUT" '/home/tmp/schaefpa/messi_datasets/SimSearchNet_queries.bin'
assert_contains "$OUTPUT" '--timeseries-size 256'
assert_contains "$OUTPUT" '--dataset-size 100000000'
assert_contains "$OUTPUT" '--filetype-int'
assert_contains "$OUTPUT" '--dataset-header-bytes 8'
assert_contains "$OUTPUT" '--query-header-bytes 8'
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" seismic standard --threads 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '/home/tmp/schaefpa/messi_datasets/seismic.bin'
assert_contains "$OUTPUT" '--timeseries-size 256'
assert_contains "$OUTPUT" '--dataset-size 100000000'
assert_not_contains "$OUTPUT" '--dataset-header-bytes'
assert_not_contains "$OUTPUT" '--query-header-bytes'
pass 'seismic and headered SimSearchNet datasets have valid runner metadata'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--timeseries-size 128'
assert_contains "$OUTPUT" '--filetype-int'
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" spacev1b standard --threads 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--timeseries-size 100'
assert_contains "$OUTPUT" '--filetype-int8'
assert_not_contains "$OUTPUT" '--filetype-int '
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" text-to-image standard --threads 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--dataset-header-bytes 8'
assert_not_contains "$OUTPUT" '--query-header-bytes'
pass 'BigANN, SpaceV, and Text-to-Image encodings match their files'

OUTPUT=$(MESSI_PHYSICAL_CORES=7 "$SCRIPT_DIR/run_dataset.sh" astro standard --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 2 ]] || fail 'default trie standard profile should emit two SPARTAN commands'
assert_contains "$OUTPUT" '--threads 7'
assert_contains "$OUTPUT" '--index-type trie'
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_contains "$OUTPUT" '--n-segments 64'
assert_contains "$OUTPUT" '--trie-split-dimensions 64'
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
assert_contains "$OUTPUT" '--trie-streaming-leaf-scan'
pass 'runner defaults to the trie benchmark profile and physical-core thread count'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --datasets astro --threads 16,32,64 \
    --index-threads 64 --index-type trie --methods spartan-depth --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--threads 16'
assert_contains "$OUTPUT" '--threads 32'
assert_contains "$OUTPUT" '--threads 64'
assert_contains "$OUTPUT" '--index-threads 64'
pass 'suite separates fixed index workers from query-core scaling'

OUTPUT=$("$SCRIPT_DIR/run_segment_scaling_experiment.sh" --datasets astro --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 20 ]] || \
    fail 'segment scaling should emit five methods at each of four widths'
for segments in 16 32 48 64; do
    [[ $(printf '%s\n' "$OUTPUT" | grep -c -- "--n-segments $segments") == 5 ]] || \
        fail "segment scaling did not emit five commands at width $segments"
    assert_contains "$OUTPUT" "--function-type 5"
    printf '%s\n' "$OUTPUT" | grep -- "--function-type 5" | grep -- "--n-segments $segments" | \
        grep -q -- "--trie-split-dimensions $segments" || \
        fail "SPARTAN split dimensions do not track width $segments"
done
assert_contains "$OUTPUT" '--threads 64'
assert_contains "$OUTPUT" '--index-threads 64'
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_not_contains "$OUTPUT" '--enable-sofa-v2'
assert_not_contains "$OUTPUT" '--isax-node-mbr'
assert_not_contains "$OUTPUT" '--isax-record-mbr-suffix-bound'
assert_not_contains "$OUTPUT" '--isax-record-lb-table'
assert_contains "$(<"$SCRIPT_DIR/run_segment_scaling_experiment.sh")" 'segments-$segments'
pass 'segment-scaling runner fixes workers, isolates widths, and excludes SOFA-v2 bounds'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --sample-type 3 --binary /tmp/MESSI --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sample-type 3'
if "$SCRIPT_DIR/run_dataset.sh" astro standard --threads 1 --sample-type 4 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted an invalid binning sample type'
fi
pass 'binning sample type can be overridden and is validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --no-tight-bound --binary /tmp/MESSI --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--tight-bound'
pass 'iSAX tight-bound pruning can be disabled explicitly'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --data-root '/tmp/data root' --binary /tmp/MESSI --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--queue-number'
pass 'queue count is optional and defaults in MESSI'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type trie --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | wc -l | tr -d ' ') == 2 ]] || fail 'trie standard profile should run only SPARTAN by default'
assert_not_contains "$OUTPUT" '--function-type 3'
assert_contains "$OUTPUT" '--function-type 5'
assert_not_contains "$OUTPUT" '--function-type 6'
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_contains "$OUTPUT" '--trie-split-dimensions 64'
assert_contains "$OUTPUT" '--trie-fanout 8'
assert_contains "$OUTPUT" '--trie-record-mbr-suffix-bound'
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
assert_not_contains "$OUTPUT" '--tight-bound'
pass 'trie standard profile excludes SAX from its method matrix'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --n-segments 32 --trie-split-dims 32 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
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
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro standard --threads 36 --queue-number 36 --index-type isax --dynamic-root-split-variance --dry-run 2>/dev/null)
[[ $(printf '%s\n' "$OUTPUT" | grep -c -- '--dynamic-root-split-variance') == 4 ]] || fail 'variance root split should apply to each learned iSAX method'
pass 'iSAX variance root split is opt-in and forwarded only to learned transforms'

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
    --no-trie-record-mbr-suffix-bound --dry-run 2>/dev/null)
assert_not_contains "$OUTPUT" '--trie-record-mbr-suffix-bound'
assert_contains "$OUTPUT" '--no-trie-record-mbr-suffix-bound'
pass 'trie record-MBR suffix pruning is enabled by default and can be disabled'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-streaming-leaf-scan --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-streaming-leaf-scan'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-streaming-leaf-scan --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie streaming leaf scan for iSAX'
fi
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --no-trie-streaming-leaf-scan --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--no-trie-streaming-leaf-scan'
assert_not_contains "$OUTPUT" ' --trie-streaming-leaf-scan'
pass 'trie streaming leaf scan defaults on, supports heap opt-out, and is scoped to trie'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-leaf-ivf 16 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-leaf-ivf 8 --trie-leaf-ivf-min-size 2048 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf 8'
assert_contains "$OUTPUT" '--trie-leaf-ivf-min-size 2048'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type trie \
    --trie-leaf-ivf 8 --trie-leaf-ivf-min-size 0 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted a non-positive trie IVF minimum size'
fi
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-leaf-ivf 16 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie leaf IVF for iSAX'
fi
pass 'trie leaf IVF is forwarded and scoped to trie'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-leaf-ivf 16 --trie-leaf-ivf-radial-bound --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf-radial-bound'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-leaf-ivf-radial-bound --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie IVF radial bound for iSAX'
fi
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type trie \
    --no-trie-leaf-ivf --trie-leaf-ivf-radial-bound --dry-run >/dev/null 2>&1; then
    fail 'runner accepted trie IVF radial bound with IVF disabled'
fi
pass 'explicit trie IVF radial bound is forwarded and requires IVF'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 64 --queue-number 64 \
    --methods spartan-depth --trie-mbr-dims 128 --n-segments 64 --trie-split-dims 64 \
    --trie-leaf-ivf 16 --trie-leaf-ivf-min-size 4096 \
    --no-trie-leaf-ivf-raw-ball-bound --trie-leaf-ivf-radial-bound \
    --trie-record-mbr-suffix-bound --trie-streaming-leaf-scan \
    --trie-residual-record-only --trie-residual-order symbolic-first \
    --trie-pruning-curve --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-pruning-curve'
assert_contains "$OUTPUT" '--no-trie-leaf-ivf-raw-ball-bound'
assert_contains "$OUTPUT" '--trie-residual-order symbolic-first'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --methods spartan-depth \
    --trie-residual-record-only --trie-pruning-curve --dry-run >/dev/null 2>&1; then
    fail 'pruning trace accepted a configuration containing the non-paper raw-ball bound'
fi
pass 'paper pruning trace forwards and validates the five-bound cascade'

OUTPUT=$("$SCRIPT_DIR/run_paper_pruning_experiment.sh" --datasets astro \
    --threads 1 --queue-number 1 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--function-type 5'
assert_contains "$OUTPUT" '--apply-z-norm'
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_contains "$OUTPUT" '--n-segments 64'
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
assert_contains "$OUTPUT" '--no-trie-leaf-ivf-raw-ball-bound'
assert_contains "$OUTPUT" '--trie-pruning-curve'
pass 'paper pruning experiment fixes the LaTeX configuration and five bounds'

for runner in dataset suite; do
    if [[ $runner == dataset ]]; then
        RADIAL_COMMAND=("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type trie --dry-run)
    else
        RADIAL_COMMAND=("$SCRIPT_DIR/run_suite.sh" standard --datasets astro --threads 1 --index-type trie --dry-run)
    fi
    OUTPUT=$("${RADIAL_COMMAND[@]}" 2>/dev/null)
    assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
    assert_contains "$OUTPUT" ' --trie-leaf-ivf-radial-bound '
    OUTPUT=$("${RADIAL_COMMAND[@]}" --no-trie-leaf-ivf-radial-bound 2>/dev/null)
    assert_contains "$OUTPUT" '--no-trie-leaf-ivf-radial-bound'
    assert_not_contains "$OUTPUT" ' --trie-leaf-ivf-radial-bound '
    OUTPUT=$("${RADIAL_COMMAND[@]}" --no-trie-leaf-ivf 2>/dev/null)
    assert_not_contains "$OUTPUT" ' --trie-leaf-ivf-radial-bound '
    OUTPUT=$("${RADIAL_COMMAND[@]}" --trie-leaf-ivf-radial-bound-auto --no-trie-leaf-ivf-radial-bound 2>/dev/null)
    assert_not_contains "$OUTPUT" '--trie-leaf-ivf-radial-bound-auto'
    assert_contains "$OUTPUT" '--no-trie-leaf-ivf-radial-bound'
done
pass 'dataset and suite default to IVF16 plus radial pruning and preserve explicit opt-outs'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --trie-leaf-ivf 16 --trie-leaf-ivf-radial-bound-auto --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf-radial-bound-auto'
assert_not_contains "$OUTPUT" ' --trie-leaf-ivf-radial-bound '
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax \
    --trie-leaf-ivf-radial-bound-auto --dry-run >/dev/null 2>&1; then
    fail 'runner accepted automatic trie IVF radial bound for iSAX'
fi
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type trie \
    --no-trie-leaf-ivf --trie-leaf-ivf-radial-bound-auto --dry-run >/dev/null 2>&1; then
    fail 'runner accepted automatic trie IVF radial bound with IVF disabled'
fi
pass 'automatic trie IVF radial bound is forwarded and requires IVF'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --query-report-interval 10 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--query-report-interval 10'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --query-report-interval -1 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted a negative query report interval'
fi
pass 'query report interval is forwarded and validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 36 --index-type trie \
    --query-repeats 3 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--query-repeats 3'
if "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --query-repeats 0 --dry-run >/dev/null 2>&1; then
    fail 'runner accepted zero query repeats'
fi
pass 'query repeats are forwarded and validated'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann high-frequency --threads 36 --queue-number 36 --index-type isax --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--apply-z-norm'
assert_contains "$OUTPUT" '--filetype-int'
assert_contains "$OUTPUT" '--queries-size 1'
assert_contains "$OUTPUT" '--sfa-n-coefficients 64'
assert_contains "$OUTPUT" '--histogram-type 2'
assert_not_contains "$OUTPUT" '--histogram-type 1'
pass 'high-frequency profile preserves BigANN flags and coefficients'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --queue-number 1 --index-type isax --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sfa-n-coefficients 64'
assert_contains "$OUTPUT" '--sampling-seed 1'
pass 'standard SFA uses the 64-coefficient training pool when permitted'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" deep1b standard --threads 1 --queue-number 1 --index-type isax --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--sfa-n-coefficients 48'
pass 'short series use the largest valid even coefficient pool'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 64 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 64'
assert_not_contains "$OUTPUT" '--sfa-n-coefficients'
pass 'trie MBR dimensions are capped by series length, not half length and SPARTAN is the default'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_not_contains "$OUTPUT" '--sfa-n-coefficients'
pass 'trie supports 128 MBR dimensions for 128-value series'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" bigann standard --threads 1 --index-type trie \
    --trie-mbr-dims 128 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 128'
assert_not_contains "$OUTPUT" '--sfa-n-coefficients'
pass 'BigANN trie MBR dimensions use the corrected 128-value series length'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --threads 1 --queue-number 1 --dataset-size 100k --dry-run 2>&1)
assert_contains "$OUTPUT" 'exceeds dataset size 100 K; using dataset size'
assert_contains "$OUTPUT" '--sample-size 100000'
pass 'binning sample size is capped for reduced datasets'

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sift1b knn --threads 36 --queue-number 36 --index-type isax --k 20 --dry-run 2>/dev/null)
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

OUTPUT=$(cd /tmp && MESSI_DRY_RUN=true MESSI_PHYSICAL_CORES=7 "$SCRIPT_DIR/run_astro.sh" 2>/dev/null)
assert_contains "$OUTPUT" '--threads 7'
assert_contains "$OUTPUT" '--index-type trie'
pass 'compatibility wrappers inherit canonical defaults when counts are omitted'

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
    'printf "%s\\n" "  symbolic record bounds: 24.35 M/query (24.35% of 100.00 M indexed series)" >&2' \
    'printf "%s\\n" "  exact distances  : 567.00 K/query (0.57% of 100.00 M indexed series)" >&2' \
    > "$FAKE_MESSI"
chmod +x "$FAKE_MESSI"
OUTPUT=$(MESSI_SHELL_LOG_DIR="$TEMP_ROOT/logs" "$SCRIPT_DIR/run_dataset.sh" astro high-frequency --threads 1 --index-type isax --data-root "$TEMP_ROOT" \
    --binary "$FAKE_MESSI" --methods spartan-width 2>&1)
assert_contains "$OUTPUT" '=== Benchmark summary: dataset=astro, profile=high-frequency ==='
assert_contains "$OUTPUT" 'iSAX    SPARTAN   width'
assert_contains "$OUTPUT" '24.35 M/query'
[[ $OUTPUT =~ 24\.35\ M/query[[:space:]]+24\.35% ]] || fail 'suite summary columns are out of order'
pass 'runner parses the symbolic-record-bound label into its compact suite summary'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" generated-queries --threads 36 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" 'spacev1B_noise_025.bin'
assert_not_contains "$OUTPUT" 'text-to-image_noise_01.bin'
assert_not_contains "$OUTPUT" 'turingANNs_noise_05.bin'
assert_not_contains "$OUTPUT" '--query-header-bytes 8'
pass 'generated-query suite expands default workloads and excludes opt-in datasets'

OUTPUT=$(bash -c 'source "$1"; active_datasets' _ "$SCRIPT_DIR/lib/datasets.sh")
assert_not_contains "$OUTPUT" 'turinganns'
assert_not_contains "$OUTPUT" 'seismic'
assert_not_contains "$OUTPUT" 'simsearchnet'
assert_not_contains "$OUTPUT" 'text-to-image'
pass 'opt-in datasets are excluded from the default suite'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--index-type trie'
assert_not_contains "$OUTPUT" '--function-type 3'
pass 'suite forwards the selected index layout to dataset runs'

OUTPUT=$(MESSI_PHYSICAL_CORES=7 "$SCRIPT_DIR/run_suite.sh" high-frequency --datasets astro --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--threads 7'
assert_contains "$OUTPUT" '--queue-number 7'
assert_contains "$OUTPUT" '--index-type trie'
assert_contains "$OUTPUT" '--n-segments 64'
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
assert_contains "$OUTPUT" '--trie-streaming-leaf-scan'
pass 'suite defaults to trie and the available physical-core count'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 64 --datasets astro --index-type trie \
    --methods spartan-depth,spartan-width --trie-mbr-dims 64 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-mbr-dimensions 64'
pass 'suite forwards trie MBR dimensions to dataset runs'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 64 --datasets astro --index-type trie \
    --leaf-size 10k --min-leaf-size 5k --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--leaf-size 10000'
assert_contains "$OUTPUT" '--min-leaf-size 5000'
assert_contains "$OUTPUT" '--initial-lbl-size 10000'
OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 64 --datasets astro --index-type trie \
    --leaf-size 40k --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--initial-lbl-size 40000'
pass 'suite forwards leaf capacity and minimum occupancy'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-record-mbr-suffix-bound --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-record-mbr-suffix-bound'
pass 'suite forwards the trie record/MBR suffix bound'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-streaming-leaf-scan --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-streaming-leaf-scan'
OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --no-trie-streaming-leaf-scan --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--no-trie-streaming-leaf-scan'
assert_not_contains "$OUTPUT" ' --trie-streaming-leaf-scan'
pass 'suite defaults to streaming leaf refinement and forwards heap opt-out'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-leaf-ivf 16 --trie-leaf-ivf-min-size 8192 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf 16'
assert_contains "$OUTPUT" '--trie-leaf-ivf-min-size 8192'
pass 'suite forwards trie leaf IVF and its minimum eligible leaf size'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --query-repeats 3 --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--query-repeats 3'
if "$SCRIPT_DIR/run_suite.sh" standard --threads 1 --datasets astro --index-type isax \
    --query-repeats 3 --dry-run >/dev/null 2>&1; then
    fail 'suite accepted repeated queries for iSAX'
fi
pass 'suite forwards in-memory trie query repeats'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-leaf-ivf 16 --trie-leaf-ivf-radial-bound --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf-radial-bound'
if "$SCRIPT_DIR/run_suite.sh" standard --threads 1 --datasets astro --index-type trie \
    --no-trie-leaf-ivf --trie-leaf-ivf-radial-bound --dry-run >/dev/null 2>&1; then
    fail 'suite accepted trie IVF radial bound with IVF disabled'
fi
pass 'suite forwards and validates trie IVF radial bound'

OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --threads 36 --datasets astro --index-type trie \
    --trie-leaf-ivf 16 --trie-leaf-ivf-radial-bound-auto --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-leaf-ivf-radial-bound-auto'
assert_not_contains "$OUTPUT" ' --trie-leaf-ivf-radial-bound '
if "$SCRIPT_DIR/run_suite.sh" standard --threads 1 --datasets astro --index-type trie \
    --no-trie-leaf-ivf --trie-leaf-ivf-radial-bound-auto --dry-run >/dev/null 2>&1; then
    fail 'suite accepted automatic trie IVF radial bound with IVF disabled'
fi
pass 'suite forwards and validates automatic trie IVF radial bound'

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

DIDS_SCRIPT="$REPO_ROOT/Competitors/dids_scripts/run_dids.sh"
OUTPUT=$("$DIDS_SCRIPT" bigann --threads 1 --dry-run)
assert_contains "$OUTPUT" 'dids_benchmark_128'
assert_contains "$OUTPUT" '--dimension 128'
assert_contains "$OUTPUT" '--scalar-type uint8'
assert_contains "$OUTPUT" '/home/tmp/schaefpa/messi_datasets/bigANN.bin'

OUTPUT=$("$DIDS_SCRIPT" spacev1b --threads 1 --dry-run)
assert_contains "$OUTPUT" '--dimension 100'
assert_contains "$OUTPUT" '--scalar-type int8'

OUTPUT=$("$DIDS_SCRIPT" text-to-image --threads 1 --dry-run)
assert_contains "$OUTPUT" '--base-header u32-count-dim'
assert_contains "$OUTPUT" '--query-header raw'

OUTPUT=$("$DIDS_SCRIPT" simsearchnet --threads 1 --dry-run)
assert_contains "$OUTPUT" '--base-header u32-count-dim'
assert_contains "$OUTPUT" '--query-header u32-count-dim'

OUTPUT=$("$DIDS_SCRIPT" turinganns --threads 1 --dry-run)
assert_contains "$OUTPUT" '--scalar-type float32'
assert_contains "$OUTPUT" '--base-header u32-count-dim'
assert_contains "$OUTPUT" '--query-header u32-count-dim'

OUTPUT=$(bash -c 'source "$1"; active_datasets' _ "$REPO_ROOT/Competitors/dids_scripts/lib/datasets.sh")
assert_not_contains "$OUTPUT" 'seismic'
assert_not_contains "$OUTPUT" 'simsearchnet'
assert_not_contains "$OUTPUT" 'text-to-image'
assert_not_contains "$OUTPUT" 'turinganns'
pass 'DIDS dataset metadata and defaults match the evaluation loaders'

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

OUTPUT=$("$SCRIPT_DIR/run_dataset.sh" sald standard --methods spartan-depth --trie-residual-record-only --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-residual-record-only'
OUTPUT=$("$SCRIPT_DIR/run_suite.sh" standard --datasets SALD --methods spartan-depth --trie-residual-record-only --dry-run 2>/dev/null)
assert_contains "$OUTPUT" '--trie-residual-record-only'
if "$SCRIPT_DIR/run_dataset.sh" sald standard --methods pisa-depth --trie-residual-record-only --dry-run >/dev/null 2>&1; then
    fail 'residual bound must reject PISA'
fi
if "$SCRIPT_DIR/run_dataset.sh" sald standard --index-type isax --methods spartan-depth --trie-residual-record-only --dry-run >/dev/null 2>&1; then
    fail 'residual bound must reject iSAX'
fi
pass 'ResSPARTAN is forwarded and rejects unsupported method/layout combinations'

OUTPUT=$("$SCRIPT_DIR/tune_trie_dataset.sh" --help)
assert_contains "$OUTPUT" 'spartan-depth or spartan-width to win'
assert_contains "$OUTPUT" 'best-config.env'
assert_contains "$OUTPUT" 'Final query repetitions per built finalist'
assert_contains "$OUTPUT" 'DATASET[,DATASET...]'
pass 'dataset-specific trie tuner documents joint method selection and safe result outputs'

TUNER_FAKE_MESSI="$TEMP_ROOT/tuner_fake_messi"
printf '%s\n' '#!/usr/bin/env bash' \
    'histogram=1; repeats=1' \
    'while (( $# )); do case "$1" in --histogram-type) histogram=$2; shift 2 ;; --query-repeats) repeats=$2; shift 2 ;; *) shift ;; esac; done' \
    'query=1.000; (( histogram == 1 )) || query=0.500' \
    'printf "    eligible leaves  : 10\n    clusters         : 160\n"' \
    'printf ">>> trie build timing\n    total      : 2.000 s\n"' \
    'if (( repeats > 1 )); then for ((i=1; i<=repeats; ++i)); do printf ">>> query repeat %d/%d wall time: %s s\n" "$i" "$repeats" "$query"; done; fi' \
    'printf ">>> query wall time: %s s\n" "$query"' \
    'printf "=== Query summary ===\n  wall time        : %s s (1.000 ms/query)\n" "$query"' \
    'printf "  symbolic record bounds: 1.00 M/query (1.00%% of 100.00 M indexed series)\n"' \
    'printf "  exact distances  : 1.00 K/query (0.00%% of 100.00 M indexed series)\n"' \
    > "$TUNER_FAKE_MESSI"
chmod +x "$TUNER_FAKE_MESSI"
TUNER_ROOT="$TEMP_ROOT/tuner"
"$SCRIPT_DIR/tune_trie_dataset.sh" astro --threads 1 --repeats 3 \
    --output-root "$TUNER_ROOT" --binary "$TUNER_FAKE_MESSI" \
    --dataset-file "$TEMP_ROOT/astro.bin" --query-file "$TEMP_ROOT/astro_queries.bin" \
    --dataset-size 1 --query-size 1 >/dev/null
assert_contains "$(<"$TUNER_ROOT/astro/best-config.env")" 'METHOD=spartan-width'
[[ $(awk -F '\t' '$3 == "01-structure" && $2 == "spartan-depth" { found=1 } END { print found+0 }' "$TUNER_ROOT/astro/all-runs.tsv") == 1 ]] ||
    fail 'joint tuner did not screen spartan-depth'
[[ $(awk -F '\t' '$3 == "01-structure" && $2 == "spartan-width" { found=1 } END { print found+0 }' "$TUNER_ROOT/astro/all-runs.tsv") == 1 ]] ||
    fail 'joint tuner did not screen spartan-width'
for method in spartan-depth spartan-width; do
    for phase in 01-structure 02-mbr 03-ivf-groups 04-ivf-min-size 05-residual 06-final-query-only; do
        [[ -d $TUNER_ROOT/astro/$method/$phase ]] || fail "missing $method/$phase"
    done
done
[[ $(awk -F '\t' '$3 == "06-final-query-only" { count++ } END { print count+0 }' "$TUNER_ROOT/astro/all-runs.tsv") == 12 ]] ||
    fail 'joint tuner did not record three query repeats for both finalists'
[[ $(awk 'END {print NR}' "$TUNER_ROOT/astro/final-ranking.tsv") == 5 ]] ||
    fail 'final ranking must contain four distinct method/configuration finalists'
pass 'trie tuner completes both methods and compares all four finalists'

# Simulate an older completed output that stopped depth after stage 1.
LEGACY_SAVED="$TEMP_ROOT/legacy-saved"
mkdir -p "$LEGACY_SAVED"
mv "$TUNER_ROOT/astro/.complete-both-methods" "$LEGACY_SAVED/"
for phase in 02-mbr 03-ivf-groups 04-ivf-min-size 05-residual 06-final-query-only; do
    mv "$TUNER_ROOT/astro/spartan-depth/$phase" "$LEGACY_SAVED/"
done
cp "$TUNER_ROOT/astro/spartan-width/06-final-query-only/candidate-1/metrics.tsv" "$LEGACY_SAVED/width-metrics.tsv"
"$SCRIPT_DIR/tune_trie_dataset.sh" astro --threads 1 --repeats 3 \
    --output-root "$TUNER_ROOT" --binary "$TUNER_FAKE_MESSI" \
    --dataset-file "$TEMP_ROOT/astro.bin" --query-file "$TEMP_ROOT/astro_queries.bin" \
    --dataset-size 1 --query-size 1 >/dev/null
[[ -f $TUNER_ROOT/astro/.complete-both-methods ]] || fail 'legacy output was not upgraded'
[[ -d $TUNER_ROOT/astro/spartan-depth/06-final-query-only ]] || fail 'missing method was not resumed'
cmp -s "$LEGACY_SAVED/width-metrics.tsv" "$TUNER_ROOT/astro/spartan-width/06-final-query-only/candidate-1/metrics.tsv" ||
    fail 'completed finalist metrics changed during resume'
[[ $(awk 'END {print NR}' "$TUNER_ROOT/astro/final-ranking.tsv") == 5 ]] || fail 'resumed ranking missing finalists'
pass 'legacy single-method completion resumes missing stages and preserves completed runs'

printf '1..%d\n' "$TEST_COUNT"
