This is the supporting website for the paper "Fast and Exact Similarity Search in less than a Blink of an Eye".


# Build environment variables
Set these if FFTW/LAPACK aren’t in default paths
```bash
export FFTW_LIBS="-L/opt/local/lib -lfftw3f -lfftw3"
export LAPACK_LIBS="-L/opt/local/lib -llapack -lblas"
```

# To compile SOFA 
Using Autotools, call from repo root.
```bash
./configure
make
```


# Build Python (Cython) API 
Again, set environment-variables, if FFTW/LAPACK aren’t in default paths or to enable SIMD in the extension build.

Call from repo root.
```bash
export FFTW_CFLAGS="-I/opt/local/include"
export FFTW_LIBS="-L/opt/local/lib -lfftw3f -lfftw3"
export LAPACK_LIBS="-L/opt/local/lib -llapack -lblas"
export SIMD_CFLAGS="-mavx -mavx2 -msse3"

python3 -m pip install -e ./python
```

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

`Index.add_array(data)` accepts a two-dimensional NumPy array and creates an
owned temporary float32 raw-data snapshot for exact refinement.  The snapshot
is removed by `idx.close()` or by a context manager.  The native query engines
currently return exact 1-NN distances only: `k` must be 1 and `indices` is
`None` until stable raw-record offsets are propagated by the native backends.

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
FILE_PATH=/vol/tmp/schaefpa/messi_datasets/deep1b.bin
QUERIES_PATH=/vol/tmp/schaefpa/messi_datasets/$QUERY
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

## Table: Characteristics of the 17 datasets used

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
| **STEAD** [mousavi2019stanford] | 87,323,433    | 256           |
| **TXED** [chen2024txed] | 35,851,641    | 256           |

# Competitors

The competitors are stored within the `competitors` folder.
