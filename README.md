This is the supporting website for the paper "Fast and Exact Similarity Search in less than a Blink of an Eye".


# Build MESSI

MESSI requires the single-precision FFTW library (`fftw3f`). OpenBLAS is
recommended: it provides the LAPACK routines required by PISA/SPARTAN and the
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
directory, and copies the executable to `bin/MESSI`:

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

# Minimal Python API usage
The API consumes float32 binary datasets (same format as CLI). 

This mirrors `tests/cython_with_data.py`.

```python
import numpy as np
from messi import Index

ts_size = 256
idx = Index(timeseries_size=ts_size, transform="spartan", layout="trie",
            sample_size=1000, max_query_threads=8,
            trie_mbr_dimensions=128, trie_record_lb_dimensions=32,
            trie_split_dimensions=32)
idx.add_file("data_head/astro_head.bin", ts_num=1000)

queries = np.fromfile("data_queries/astro_queries.bin", dtype=np.float32, count=10 * ts_size)
queries = queries.reshape(10, ts_size)
distances, indices = idx.search(queries, k=1)
```

`Index.add(data)` accepts a two-dimensional NumPy array and creates an
owned temporary float32 raw-data snapshot for exact refinement.  The snapshot
is removed by `idx.close()` or by a context manager.  The native query engines
return exact 1-NN distances and zero-based sequential row IDs. `k` must be 1.

## Index defaults

The direct CLI and benchmark runners default to the trie layout, 64 record-LB
dimensions, node MBRs up to 128 dimensions, 16 leaf-IVF groups, and streaming
leaf refinement. Their
automatic worker count uses available physical CPU cores rather than SMT
siblings. The Python API retains its explicit, compatibility-oriented defaults.
When iSAX is selected, the direct CLI keeps tight-bound pruning off while the
benchmark runners and Python API enable it by default.

| Setting | Direct CLI | Script runners | Python `Index` |
|---|---|---|---|
| Index layout | Trie | Trie | iSAX; pass `layout="trie"` |
| Worker threads | Available physical cores | Available physical cores | 1; pass `max_query_threads=N` |
| iSAX tight-bound pruning | Off; use `--tight-bound` | On; use `--no-tight-bound` to disable | On; pass `tight_bound=False` to disable |
| iSAX variance root splitting | Off; use `--dynamic-root-split-variance` | Off; use `--dynamic-root-split-variance` | Off; pass `dynamic_root_split_variance=True` |
| Trie node-MBR width | Automatic `min(128, series length)` | Same | Same |
| Trie record-LB width | 64 | 64 | `n_segments` (16 by default) |
| Trie record-MBR suffix pruning | On; use `--no-trie-record-mbr-suffix-bound` | On; use `--no-trie-record-mbr-suffix-bound` | On; pass `trie_record_mbr_suffix_bound=False` |
| Trie leaf refinement | Streaming LB → ED; use `--no-trie-streaming-leaf-scan` for the heap | Same | Streaming for trie; pass `trie_streaming_leaf_scan=False` for the heap |
| Trie leaf IVF groups | 16 for learned transforms; use `--no-trie-leaf-ivf` to disable | 16 | Off; pass `trie_leaf_ivf=K` |
| Trie IVF radial record bound | Off; use `--trie-leaf-ivf-radial-bound-auto` for the adaptive 25% gate or `--trie-leaf-ivf-radial-bound` unconditionally | Same | Off; pass `trie_leaf_ivf_radial_bound_auto=True` or `trie_leaf_ivf_radial_bound=True` |

Variance root splitting is valid only for iSAX SFA, SPARTAN, and PISA. Trie
leaf IVF is valid only for learned trie transforms (SFA, SPARTAN, and
PISA), with `K` from 2 to 64.

# Scripts

See the provided scripts in the `scripts`-folder for examples to run SOFA with SFA summarization.

- SAX command is `--function-type 3`
- SFA/SOFA command is `--function-type 4`
- SPARTAN command is `--function-type 5`
- PISA command is `--function-type 6`

For trie indexes, `--trie-fanout 2|4|8` selects a fixed fanout. Learned
SFA/SPARTAN/PISA tries can instead use one global, precomputed dynamic
alphabet allocation with `--trie-dynamic-alphabet`; this uses a 3-bit average
budget by default and supports 1--4 bits per coefficient (fanouts 2--16).
The allocation is computed once from the training representation and reused
for all trie splits. These fixed and dynamic modes are mutually exclusive.

SAX tries remain on the fixed-fanout path. iSAX’s
`--dynamic-root-split-variance` is a separate legacy root-only allocation and
does not control trie alphabets.

```bash
FILE_PATH=/home/tmp/schaefpa/messi_datasets/deep1b.bin
QUERIES_PATH=/home/tmp/schaefpa/messi_datasets/$QUERY
TS_SIZE=96

COEFF_NUMBER=32
DATASET_SIZE=100000000
SAMPLE_SIZE=1000000
QUERY_SIZE=100

./MESSI 
  --dataset --dataset $FILE_PATH 
  --dataset-size $DATASET_SIZE 
  --queries $QUERIES_PATH 
  --queries-size $QUERY_SIZE 
  --timeseries-size $TS_SIZE  
  --function-type 4 
  --histogram-type 2 
  --sample-type 3 
  --sample-size $SAMPLE_SIZE 
  --sfa-n-coefficients $COEFF_NUMBER  
  --is-norm 
  --SIMD
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
| **SIFT1b** [jegou2011searching] | 100,000,000   | 128           |
| **SpaceV1B** | 100,000,000 | 100 |
| **STEAD** [mousavi2019stanford] | 87,323,433    | 256           |
| **Text-to-image** | 100,000,000 | 200 |
| **TuringANNs** | 100,000,000 | 100 |
| **TXED** [chen2024txed] | 35,851,641    | 256           |

# Competitors

The competitors are stored within the `competitors` folder.
