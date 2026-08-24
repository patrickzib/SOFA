#!/usr/bin/env bash
set -euo pipefail

FFTW_PREFIX="/vol/home-vol3/wbi/schaefpa/fftw-3.3.10-install"
FFTW_LIBDIR="$FFTW_PREFIX/lib64"
[[ -d "$FFTW_LIBDIR" ]] || FFTW_LIBDIR="$FFTW_PREFIX/lib"

# Sanity check
test -f "$FFTW_PREFIX/include/fftw3.h"

# Standard autoconf vars
export CPPFLAGS="-I$FFTW_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$FFTW_LIBDIR ${LDFLAGS:-}"
export LIBS="-lfftw3f ${LIBS:-}"

# OpenBLAS provides LAPACK and the optional CBLAS batch-projection path.
export LAPACK_LIBS="-L/usr/lib64 -l:libopenblas.so.0"
export CBLAS_LIBS="$LAPACK_LIBS"

# Optimization flags
export CFLAGS="-O3 -DNDEBUG -march=native ${CFLAGS:-}"
export CXXFLAGS="-O3 -DNDEBUG -march=native ${CXXFLAGS:-}"

make clean || true
CPPFLAGS="$CPPFLAGS" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LAPACK_LIBS="$LAPACK_LIBS" CBLAS_LIBS="$CBLAS_LIBS" \
CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS" \
./configure --with-fftw="$FFTW_PREFIX"
make -j

# python3 -m pip install -e ./python --no-build-isolation
