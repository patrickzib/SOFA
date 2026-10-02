# S³-Trie

S³-Trie provides fast and exact similarity search using SPARTAN transforms and
a symbolic trie. This repository includes the native engine, a Python API,
and benchmark scripts, as well as the MESSI (SAX+iSAX) and SOFA (SFA+iSAX)
configurations.

# Build S³-Trie

S³-Trie requires the single-precision FFTW library (`fftw3f`). OpenBLAS is
recommended: it provides the LAPACK routines required by SPARTAN and the
optional CBLAS acceleration used for bulk PCA projection. On systems where
FFTW is discoverable through `pkg-config`, build from the repository root:

```bash
./configure
make -j
```

If FFTW is installed outside the system paths, give its prefix to `configure`
and provide OpenBLAS for LAPACK/CBLAS. For a MacPorts installation this is:

```bash
export LAPACK_LIBS="-L/opt/local/lib -lopenblas"
export CBLAS_LIBS="$LAPACK_LIBS"
./configure --with-fftw=/opt/local
make -j
```

`build_local.sh` is the maintained MacPorts convenience build. It sets those
paths, enables native optimization, creates a fresh out-of-tree `build/`
directory, and copies the executable to `bin/MESSI` (the current executable
name):

```bash
./build_local.sh
```

`build_sonic.sh` is the equivalent site-specific helper for the Sonic cluster;
edit its FFTW prefix if that installation changes. It defaults to AVX2 because
AVX-512 can downclock this workload; set `MESSI_SONIC_SIMD_FLAGS="-mavx512f"`
to explicitly benchmark the AVX-512 path. Pass `--enable-simd=no` to
`configure` to disable SIMD detection entirely. ARM NEON uses its native
implementation automatically.

Run `autoreconf -fi` before `configure` only when modifying Autotools inputs
or when a clone does not include a usable generated `configure` script.

The detailed S³-Trie lower-bound profiler is compiled out by default. Build a
separate experiment binary with `./configure --enable-trie-pruning-trace`; its
`--trie-pruning-curve` runtime option is then available. Normal builds contain
neither the counters nor their timing calls.

# Build Python (Cython) API

The Python extension compiles the native engine directly. Use an environment
with NumPy and Cython installed, and provide the same FFTW/OpenBLAS paths when
they are not in the compiler defaults:

```bash
export FFTW_CFLAGS="-I/opt/local/include"
export FFTW_LIBS="-L/opt/local/lib -lfftw3f"
export LAPACK_LIBS="-L/opt/local/lib -lopenblas"

python3 -m pip install --no-build-isolation -e ./python
```

For a direct extension build, run the equivalent command from `python/`:
`python3 setup.py build_ext --inplace`.

The Python package is the maintained high-level binding; the legacy Node.js
binding and its vendored build artifacts are no longer included.

# Minimal Python API usage
The API accepts NumPy arrays as well as float32 binary datasets in the CLI
format.

```python
import numpy as np
from s3trie import Index

data = np.load("data.npy").astype(np.float32)
queries = np.load("queries.npy").astype(np.float32)

with Index(timeseries_size=data.shape[1], index="s3trie") as index:
    index.add(data)
    distances, ids = index.search(queries)
```

`Index.add(data)` accepts a two-dimensional NumPy array and creates an
owned temporary float32 raw-data snapshot for exact refinement.  The snapshot
is removed by `idx.close()` or by a context manager.  The native query engines
return exact 1-NN distances and zero-based sequential row IDs. `k` must be 1.
Use `Index.add_file(...)` to build directly from a CLI-format binary dataset.

`index="messi"` selects SAX+iSAX, `index="sofa"` selects SFA+iSAX, and
`index="s3trie"` selects the paper-oriented SPARTAN+trie configuration.
Advanced callers may still select `layout=` and `transform=` directly; these
must agree when combined with a preset.

The same names are accepted by the native CLI with `--index messi|sofa|s3trie`.
Explicit CLI options override the preset, so binning, segment counts, bounds,
and query settings can still be tuned per run.

The S³-Trie preset uses equi-width binning, a 64-dimensional record bound,
up to 128 MBR dimensions, 20K leaves, IVF-16 with raw-ball and radial pruning,
record-MBR suffix pruning, streaming refinement, and symbolic-first residual
record pruning. Each setting remains explicitly overrideable.

## Index defaults

The direct CLI and benchmark runners default to the trie layout, 64 record-LB
dimensions, node MBRs up to 128 dimensions, 16 leaf-IVF groups, and streaming
leaf refinement. Their automatic worker count uses available physical CPU cores
rather than SMT siblings. The Python `Index` defaults to the `s3trie` preset;
pass explicit `layout=`, `transform=`, or `function_type=` arguments for
lower-level configuration.
When iSAX is selected, the direct CLI, benchmark runners, and Python API all
enable tight-bound pruning by default.

| Setting | Direct CLI | Script runners | Python `Index` |
|---|---|---|---|
| Index layout | Trie by default; use `--index messi`, `--index sofa`, or `--index s3trie` | Trie by default; use `--index messi`, `--index sofa`, or `--index s3trie` | S³-Trie (`s3trie`) by default; use `index="messi"` or `index="sofa"` |
| Worker threads | Available physical cores; use `--threads N` | Available physical cores; use `--threads N` | 1; pass `max_query_threads=N` |
| iSAX tight-bound pruning | On; use `--no-tight-bound` to disable | On; use `--no-tight-bound` to disable | On; pass `tight_bound=False` to disable |
| iSAX variance root splitting | Off; use `--dynamic-root-split-variance` | Off; use `--dynamic-root-split-variance` | Off; pass `dynamic_root_split_variance=True` |
| Trie node-MBR width | Automatic `min(128, series length)`; use `--trie-mbr-dimensions N` | Automatic `min(128, series length)`; use `--trie-mbr-dims N` | Automatic `min(128, series length)`; pass `trie_mbr_dimensions=N` |
| Trie record-LB width | 64; use `--n-segments N` | 64; use `--n-segments N` | 64; pass `n_segments=N` or `trie_record_lb_dimensions=N` |
| Trie record-MBR suffix pruning | On; use `--no-trie-record-mbr-suffix-bound` | On; use `--no-trie-record-mbr-suffix-bound` | On; pass `trie_record_mbr_suffix_bound=False` |
| Trie leaf refinement | Streaming LB → ED; use `--no-trie-streaming-leaf-scan` for the heap | Streaming LB → ED; use `--no-trie-streaming-leaf-scan` for the heap | Streaming for trie; pass `trie_streaming_leaf_scan=False` for the heap |
| Trie leaf IVF groups | 16 for learned transforms; use `--no-trie-leaf-ivf` to disable | 16 for learned transforms; use `--no-trie-leaf-ivf` to disable | 16 with `index="s3trie"`; pass `trie_leaf_ivf=0` to disable |
| Trie IVF radial record bound | On with IVF; use `--no-trie-leaf-ivf-radial-bound`, or `--trie-leaf-ivf-radial-bound-auto` | On with IVF; use `--no-trie-leaf-ivf-radial-bound`, or `--trie-leaf-ivf-radial-bound-auto` | On with `index="s3trie"`; pass `trie_leaf_ivf_radial_bound=False` to disable |

Variance root splitting is valid only for learned iSAX transforms. Trie leaf
IVF is valid only for learned trie transforms, with `K` from 2 to 64.

### iSAX root dimensions

The fixed iSAX root table always has `2^16` entries. SAX uses 16 uniformly
spaced symbolic dimensions for its root key. Learned iSAX transforms select
the 16 highest-variance symbolic dimensions for the root key; all symbolic
dimensions remain available for deeper splits, node MBRs, record lower bounds,
and SIMD evaluation. The selected dimensions are persisted with index settings
and printed in the build and query configuration logs.

`--dynamic-root-split-variance` is a separate opt-in policy that allocates
root bits by variance. It does not change the fixed root-dimension map or
restrict the dimensions available below the root. Existing indexes without a
persisted map use the legacy uniform mapping.

# Scripts

See [the benchmark scripts](scripts/README.md) for the complete runner
documentation. Named presets keep the transform and layout together:

```bash
scripts/run_dataset.sh bigann standard --index messi
scripts/run_dataset.sh bigann standard --index sofa --methods depth
scripts/run_dataset.sh bigann standard --index sofa --methods width
scripts/run_dataset.sh bigann standard --index s3trie --methods depth
scripts/run_dataset.sh bigann standard --index s3trie --methods width
```

For iSAX, `--histogram-type 1` selects equi-depth binning and
`--histogram-type 2` selects equi-width binning. The legacy
`--dynamic-root-split-variance` option is separate from trie alphabet
allocation and is available for learned iSAX runs.

```bash
scripts/run_dataset.sh deep1b standard --index sofa --methods width
```

For help, please type:
```bash
./MESSI --help
```


# Datasets

Instruction for downloading the datasets is in the `datasets` folder. The size of the datasets is too large to provide a direct link.
Some datasets must be downloaded, others generated from seisbench.

## Table: Characteristics of benchmark datasets

| Dataset Name | Series       | Series Length |
|--------------|--------------|---------------|
| **Astro** [soldi2014long] | 100,000,000   | 256           |
| **BigANN** [simhadri2022results] | 100,000,000   | 100           |
| **Deep1b** [babenko2016efficient] | 100,000,000   | 96            |
| **ETHZ** [woollam2022seisbench] | 4,999,932     | 256           |
| **Iquique** [woollam2019convolutional] | 578,853       | 256           |
| **ISC_EHB_DepthPhases** [munchmeyer2024learning] | 100,000,000   | 256           |
| **LenDB** [magrini2020local] | 37,345,260    | 256           |
| **Meier2019JGR** [woollam2022seisbench] | 6,361,998     | 256           |
| **NEIC** [yeck2021leveraging] | 93,473,541    | 256           |
| **OBS** [bornstein2024pickblue] | 15,508,794    | 256           |
| **OBST2024** [niksejel2024obstransformer] | 4,160,286     | 256           |
| **PNW** [ni2023curated] | 31,982,766    | 256           |
| **SALD** [url:SALD] | 100,000,000   | 128           |
| **SCEDC** [center2013southern] | 100,000,000   | 256           |
| **Seismic** | 100,000,000 benchmark subset | 256 |
| **SIFT1b** [jegou2011searching] | 100,000,000   | 128           |
| **SimSearchNet++** [simhadri2022results] | 100,000,000 benchmark subset | 256 |
| **SpaceV1B** | 100,000,000 | 100 |
| **STEAD** [mousavi2019stanford] | 87,323,433    | 256           |
| **Text-to-image** | 100,000,000 | 200 |
| **TuringANNs** | 100,000,000 | 100 |
| **TXED** [chen2024txed] | 35,851,641    | 256           |

# Competitors

The competitors are stored within the `competitors` folder.
