from pathlib import Path
from setuptools import setup, Extension
import numpy
import os

try:
    from Cython.Build import cythonize
except ImportError:  # pragma: no cover - fallback for runtime without Cython
    cythonize = None

ROOT = Path(__file__).resolve().parent.parent
FFTW_CFLAGS = os.environ.get("FFTW_CFLAGS", "")
FFTW_LIBS = os.environ.get("FFTW_LIBS", "")
SIMD_CFLAGS = os.environ.get("SIMD_CFLAGS", "")
COMMON_CFLAGS = [flag for flag in (FFTW_CFLAGS + " " + SIMD_CFLAGS).split() if flag]
COMMON_LDFLAGS = [flag for flag in FFTW_LIBS.split() if flag]

sources = [
    "src/ads/api.c",
    "src/ads/isax_file_loaders.c",
    "src/ads/isax_first_buffer_layer.c",
    "src/ads/isax_index.c",
    "src/ads/isax_node.c",
    "src/ads/isax_node_buffer.c",
    "src/ads/isax_node_record.c",
    "src/ads/isax_node_split.c",
    "src/ads/isax_query_engine.c",
    "src/ads/isax_visualize_index.c",
    "src/ads/pqueue.c",
    "src/ads/sax/sax.c",
    "src/ads/sax/ts.c",
    "src/ads/inmemory_query_engine.c",
    "src/ads/parallel_inmemory_query_engine.c",
    "src/ads/inmemory_index_engine.c",
    "src/ads/parallel_index_engine.c",
    "src/ads/parallel_query_engine.c",
    "src/ads/inmemory_topk_engine.c",
    "src/ads/sfa/dft.c",
    "src/ads/sfa/sfa.c",
    "src/ads/calc_utils.c",
]
sourcedirs = [str(ROOT / path) for path in sources]

compiler = os.environ.get("CC", "") or os.environ.get("CXX", "")
if compiler and "clang" in compiler:
    omp_compile = ["-Xpreprocessor", "-fopenmp"]
    omp_link = ["-lomp"]
else:
    omp_compile = ["-fopenmp"]
    omp_link = ["-fopenmp"]


def _env_flag(name: str, default: bool = False) -> bool:
    """Return True when the surrounding environment requested a feature."""
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


enable_asan = _env_flag("MESSI_ENABLE_ASAN")
asan_compile_flags = ["-fsanitize=address", "-fno-omit-frame-pointer"] if enable_asan else []
asan_link_flags = ["-fsanitize=address"] if enable_asan else []

source_file = "messi/_index.pyx" if cythonize is not None else "messi/_index.c"

extension = Extension(
    "messi._index",
    sources=[source_file] + sourcedirs,
    include_dirs=[str(ROOT), str(ROOT / "include"), numpy.get_include()],
    extra_compile_args=COMMON_CFLAGS + omp_compile + ["-fcommon", "-O0", "-g"] + asan_compile_flags,
    extra_link_args=COMMON_LDFLAGS + omp_link + asan_link_flags,
    language="c",
)

if cythonize is not None:
    extensions = cythonize([extension], compiler_directives={"language_level": 3})
else:
    extensions = [extension]

setup(
    name="messi",
    version="0.1.0",
    packages=["messi"],
    package_dir={"messi": "messi"},
    ext_modules=extensions,
)
