# Benchmark scripts

`run_dataset.sh` is the canonical entrypoint for one dataset and profile. The
older `run_*.sh` names remain as compatibility wrappers. Both canonical
runners default to trie indexes with 128 MBR dimensions (capped by series
length), 64 record-bound dimensions, and 16 leaf-IVF groups. If `--threads`
is omitted, they use the physical cores available to the process; set
`MESSI_PHYSICAL_CORES=N` when platform or scheduler topology is unavailable.

```bash
scripts/run_dataset.sh astro standard --threads 36 --queue-number 36
scripts/run_dataset.sh bigann knn --threads 36 --queue-number 36 --k 20
scripts/run_dataset.sh sald sampling --threads 36 --queue-number 36 --sample-factor 0.2
scripts/run_dataset.sh astro high-frequency --threads 36 --queue-number 36 --dry-run
scripts/run_dataset.sh simsearchnet standard --threads 36
```

The profiles preserve the existing experiment matrices. With the default trie
layout, `standard` runs SFA and SPARTAN variants and `knn` runs SFA variants;
select `--index-type isax` for the wider legacy matrices below:

- `standard`: SAX, SFA, PISA, and SPARTAN variants
- `high-frequency`: one query with SFA equi-width
- `knn`: SAX and both SFA histogram variants with top-K search
- `sampling`: both SFA histogram variants at a configurable sample fraction

Use `run_suite.sh` for the complete benchmark matrices. It also exposes the
`generated-queries`, `hard-queries`, and `noise-workloads` suites migrated from
`scripts/old/`.

## Trie fanout and dynamic alphabets

Trie runs support either one fixed fanout or one globally precomputed dynamic
alphabet allocation; the two modes are mutually exclusive.

```bash
# Fixed fanout
scripts/run_suite.sh standard --threads 64 --index-type trie --trie-fanout 8

# Dynamic alphabet: average 3 bits, bounded to 1--4 bits per coefficient
scripts/run_suite.sh standard --threads 64 --index-type trie \
  --trie-dynamic-alphabet \
  --trie-min-fanout 2 --trie-max-fanout 16 \
  --trie-alphabet-budget-bits 3
```

Dynamic trie alphabets are computed once from the global training
representation for SFA, SPARTAN, and PISA. The resulting per-dimension bit
allocation is reused at every trie split; node-local alphabets are not
recomputed. A 1--4 bit range corresponds to fanouts 2, 4, 8, and 16.

SAX tries use the fixed-fanout path because SAX has no learned training
variance allocation. iSAX’s `--dynamic-root-split-variance` option is separate
and intentionally retains its legacy root-only bit allocation; it is not the
trie dynamic-alphabet mechanism.

## Paths

Command-line path options take precedence over these environment variables:

| Variable | Default | Purpose |
| --- | --- | --- |
| `MESSI_BINARY` | repository `bin/MESSI` | MESSI executable |
| `MESSI_DATA_ROOT` | `/home/tmp/schaefpa/messi_datasets` | main datasets and queries |
| `MESSI_QUERY_ROOT` | main data root | optional separate main query root |
| `MESSI_SEISBENCH_ROOT` | `/home/tmp/schaefpa/seismic` | SeisBench datasets and queries |
| `MESSI_SEISBENCH_QUERY_ROOT` | SeisBench root | optional separate SeisBench query root |
| `MESSI_LOG_ROOT` | `$HOME/MESSI_logs` | logs produced by MESSI |
| `MESSI_SHELL_LOG_DIR` | `MESSI_LOG_ROOT` | combined runner transcript for each dataset/profile run |
| `MESSI_RESULTS_ROOT` | `$HOME/MESSI_SFA_logs` | archived benchmark results |
| `MESSI_INDEX_ROOT` | `$HOME/index` | index directory cleared by `run_cleanup.sh` |
| `MESSI_DRY_RUN` | `false` | enable dry-run mode through compatibility wrappers |

Relative dataset and query overrides are resolved below the applicable data
root. Absolute overrides are used unchanged. `--dry-run` prints shell-escaped
commands and performs no benchmark or log archival.

The regular suite includes `seismic` and `simsearchnet`. Both use the first
100 million base records for comparability with the other suite datasets.
SimSearchNet is 256-dimensional uint8 data in the ANN benchmark binary format;
its dataset and query files have an 8-byte `(count, dimensions)` header, which
the runner skips automatically. BigANN is 128-dimensional uint8, while SpaceV
is 100-dimensional signed int8. Text-to-Image uses float32 with an 8-byte base
header and a raw local query file. `--dataset-header-bytes` and
`--query-header-bytes` override the offsets independently;
`--input-header-bytes` remains a shorthand for setting both.

TuringANNS is excluded from the default regular and generated-query suites
because the currently observed local `turingANNs.bin` is not the canonical
100-dimensional float32 collection. It remains selectable explicitly after
installing the correct headered base and query files.

For trie runs, `--trie-leaf-ivf 16` adds a flat, post-build 16-list IVF/MRB
directory inside terminal leaves with at least 4 K records. It is enabled by
default and can be disabled with `--no-trie-leaf-ivf` for A/B benchmarking.
Construction clusters eligible leaves independently in parallel, using the
existing `--threads` setting; the build log reports the active worker count.

`--trie-leaf-ivf-radial-bound` stores each record's distance from its raw-space
IVF centroid as a contiguous float32 array. Records remain in their existing
IVF-cluster order. During a leaf scan, AVX-512 classifies 16 radii per load and
AVX2 classifies 8; the survivor mask is processed in record order before the
symbolic record bound. Rejected records skip both the symbolic bound and exact
distance. The option adds one float per record in eligible leaves but performs
no radius sorting or two-front traversal.

Use `--trie-leaf-ivf-radial-bound-auto` for the conservative adaptive mode.
It measures at least 64 K radial candidates of each query (rounded to whole
concurrently scanned IVF ranges) and bypasses the bound for the remaining
ranges when fewer than 25% were rejected. The calibration is shared across
query workers and stops all adaptive accounting after its one query-level
decision. The existing
`--trie-leaf-ivf-radial-bound` option remains unconditional.

The query summary attributes pruning to the first successful bound: node MBR,
leaf-IVF symbolic MBR, leaf-IVF raw ball, per-record radial bound, then the
symbolic record bound. Each line reports an average count per query, its
stage-local pruning rate, and its share of all indexed records, so overlapping
bounds are not double-counted.

Trie leaf refinement streams by default: it computes each record's lower bound
and immediately runs exact distance when that bound passes, so an improved BSF
affects the very next record. Cluster and leaf traversal ordering is unchanged.
Use `--no-trie-streaming-leaf-scan` to restore the best-first record heap for
A/B benchmarks; `--trie-streaming-leaf-scan` remains as an explicit spelling
of the default.

Result archival intentionally preserves the historical behavior: an existing
`DATASET/RUN` directory is replaced. Labels are restricted to single path
components, and the destination is checked before removal.

## Validation

Run `scripts/tests/test_scripts.sh` for syntax, argument-generation, wrapper,
and result-archival tests. Set `RUN_MESSI_INTEGRATION=1` to additionally run a
small fixture test with the real `bin/MESSI`; no substitute executable is used.

Files below `scripts/old/` are historical and unsupported. Useful datasets and
query matrices from that directory are represented by the maintained runner.

## Lower-bound versus exact-distance microbenchmark

On Sonic, compile and run the record lower-bound and Euclidean-distance kernels
with explicit AVX2 and AVX-512 builds:

```bash
scripts/run_lb_vs_ed_sonic.sh --count 2000000 --threads 1 --trials 5
```

The default count is two million records (about 2.04 GiB of inputs), access is
randomized, and early abandonment is disabled so the output compares complete
kernel costs.  It reports LB16/LB32/LB48/LB64, ED100/ED256, and the minimum
pruning percentage at which each lower bound pays for itself.  For the
64-physical-core production case, interleave allocations across Sonic's two
NUMA nodes:

```bash
OMP_PLACES=cores OMP_PROC_BIND=close \
  numactl --interleave=all scripts/run_lb_vs_ed_sonic.sh \
  --count 2000000 --threads 64 --trials 5
```

Use `--sequential` to measure contiguous scans instead of the default randomized
record order.  Build products are written under `build/benchmarks/`.
