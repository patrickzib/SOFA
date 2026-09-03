#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="$REPO_ROOT/build/benchmarks"
CC_BIN=${CC:-gcc}

if [[ $(uname -m) != x86_64 ]]; then
    echo "error: this runner requires an x86-64 machine" >&2
    exit 1
fi
if ! grep -qw avx2 /proc/cpuinfo; then
    echo "error: CPU does not advertise AVX2" >&2
    exit 1
fi
if ! grep -qw avx512f /proc/cpuinfo; then
    echo "error: CPU does not advertise AVX-512F" >&2
    exit 1
fi

mkdir -p "$BUILD_DIR"

COMMON_FLAGS=(
    -O3
    -DNDEBUG
    -DADS_HAVE_AVX2=1
    -fcommon
    -fopenmp
    -march=native
    -I"$SCRIPT_DIR/benchmark_config"
    -I"$REPO_ROOT"
    -I"$REPO_ROOT/include"
)
SOURCES=(
    "$SCRIPT_DIR/bench_lb_vs_ed.c"
    "$REPO_ROOT/src/ads/sax/ts.c"
)

echo ">>> compiling AVX2-only benchmark"
"$CC_BIN" "${COMMON_FLAGS[@]}" -mavx2 -mfma -mno-avx512f \
    "${SOURCES[@]}" -lm -o "$BUILD_DIR/bench_lb_vs_ed_avx2"

echo ">>> compiling AVX-512 benchmark"
"$CC_BIN" "${COMMON_FLAGS[@]}" -mavx512f -mavx512bw -mavx512dq -mavx512vl \
    "${SOURCES[@]}" -lm -o "$BUILD_DIR/bench_lb_vs_ed_avx512"

export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}

echo ">>> running AVX2"
"$BUILD_DIR/bench_lb_vs_ed_avx2" "$@"
echo ">>> running AVX-512"
"$BUILD_DIR/bench_lb_vs_ed_avx512" "$@"
