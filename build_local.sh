#!/usr/bin/env bash
set -euo pipefail

FFTW_PREFIX="/opt/local"

# MacPorts usually uses /opt/local/include and /opt/local/lib
FFTW_LIBDIR="$FFTW_PREFIX/lib"

# Sanity check
test -f "$FFTW_PREFIX/include/fftw3.h"

# Standard autoconf vars
export CPPFLAGS="-I$FFTW_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$FFTW_LIBDIR ${LDFLAGS:-}"
export LIBS="-lfftw3f ${LIBS:-}"

# OpenBLAS provides both LAPACK and the optional CBLAS projection fast path.
export LAPACK_LIBS="-L$FFTW_PREFIX/lib -lopenblas"
export CBLAS_LIBS="$LAPACK_LIBS"

# Optimization flags
export CFLAGS="-O3 -DNDEBUG -march=native ${CFLAGS:-}"
export CXXFLAGS="-O3 -DNDEBUG -march=native ${CXXFLAGS:-}"

# Clean out-of-tree build
rm -rf build
mkdir -p build
cd build

# Configure
../configure --with-fftw="$FFTW_PREFIX"

# Build
mkdir -p bin
make -j

# Copy MESSI executable to top-level bin/
mkdir -p ../bin
cp bin/MESSI ../bin/MESSI

# Install Python package in editable mode
# python3 -m pip install -e ../python --no-build-isolation
