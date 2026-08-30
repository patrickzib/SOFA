#!/usr/bin/env bash
set -euo pipefail

FFTW_PREFIX="/vol/home-vol3/wbi/schaefpa/fftw-3.3.10-install"
FFTW_LIBDIR="$FFTW_PREFIX/lib64"
[[ -d "$FFTW_LIBDIR" ]] || FFTW_LIBDIR="$FFTW_PREFIX/lib"

# Sonic's system OpenBLAS is built with GNU Fortran.  When this script is run
# from Conda, its linker does not automatically search the active compiler's
# libgfortran directory for OpenBLAS's transitive dependency.  Ask that
# compiler for the path instead of hard-coding a Conda version.  The override
# supports site-specific installations.
GFORTRAN_LIBDIR="${MESSI_GFORTRAN_LIBDIR:-}"
GFORTRAN_LINK_FLAG="-lgfortran"
if [[ -z "$GFORTRAN_LIBDIR" ]]; then
    GFORTRAN_LIBRARY=""
    for runtime_name in libgfortran.so libgfortran.so.5; do
        candidate="$("${CC:-gcc}" -print-file-name="$runtime_name")"
        if [[ "$candidate" != "$runtime_name" && -f "$candidate" ]]; then
            GFORTRAN_LIBRARY="$candidate"
            [[ "$runtime_name" == "libgfortran.so.5" ]] &&
                GFORTRAN_LINK_FLAG="-l:libgfortran.so.5"
            break
        fi
    done
    if [[ -z "$GFORTRAN_LIBRARY" ]]; then
        printf 'error: could not locate libgfortran for %s; set MESSI_GFORTRAN_LIBDIR explicitly\n' \
            "${CC:-gcc}" >&2
        exit 1
    fi
    GFORTRAN_LIBDIR="$(dirname "$GFORTRAN_LIBRARY")"
fi
[[ -d "$GFORTRAN_LIBDIR" ]] || {
    printf 'error: MESSI_GFORTRAN_LIBDIR is not a directory: %s\n' "$GFORTRAN_LIBDIR" >&2
    exit 1
}
if [[ ! -f "$GFORTRAN_LIBDIR/libgfortran.so" && ! -f "$GFORTRAN_LIBDIR/libgfortran.so.5" ]]; then
    printf 'error: no libgfortran runtime found in %s\n' "$GFORTRAN_LIBDIR" >&2
    exit 1
fi
[[ -f "$GFORTRAN_LIBDIR/libgfortran.so" ]] || GFORTRAN_LINK_FLAG="-l:libgfortran.so.5"

# Sanity check
test -f "$FFTW_PREFIX/include/fftw3.h"

# Standard autoconf vars
export CPPFLAGS="-I$FFTW_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$FFTW_LIBDIR -L$GFORTRAN_LIBDIR -Wl,-rpath,$GFORTRAN_LIBDIR -Wl,-rpath-link,$GFORTRAN_LIBDIR ${LDFLAGS:-}"
export LIBS="-lfftw3f ${LIBS:-}"

# OpenBLAS provides LAPACK and the optional CBLAS batch-projection path.
export LAPACK_LIBS="-L/usr/lib64 -l:libopenblas.so.0 $GFORTRAN_LINK_FLAG"
export CBLAS_LIBS="$LAPACK_LIBS"

# Optimization flags.  AVX-512 can downclock this workload on Sonic and is
# slower than AVX2 in practice, so AVX2 is the default.  Override with
# MESSI_SONIC_SIMD_FLAGS="-mavx512f" when explicitly benchmarking AVX-512.
SONIC_SIMD_FLAGS="${MESSI_SONIC_SIMD_FLAGS:--mavx2 -mno-avx512f}"
export CFLAGS="-O3 -DNDEBUG -march=native ${CFLAGS:-} ${SONIC_SIMD_FLAGS}"
export CXXFLAGS="-O3 -DNDEBUG -march=native ${CXXFLAGS:-} ${SONIC_SIMD_FLAGS}"

make clean || true
CPPFLAGS="$CPPFLAGS" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LAPACK_LIBS="$LAPACK_LIBS" CBLAS_LIBS="$CBLAS_LIBS" \
CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS" \
./configure --with-fftw="$FFTW_PREFIX"
make -j

# python3 -m pip install -e ./python --no-build-isolation
