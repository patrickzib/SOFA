# Benchmark scripts

`run_dataset.sh` is the canonical entrypoint for one dataset and profile. The
older `run_*.sh` names remain as compatibility wrappers.

```bash
scripts/run_dataset.sh astro standard --threads 36 --queue-number 36
scripts/run_dataset.sh bigann knn --threads 36 --queue-number 36 --k 20
scripts/run_dataset.sh sald sampling --threads 36 --queue-number 36 --sample-factor 0.2
scripts/run_dataset.sh astro high-frequency --threads 36 --queue-number 36 --dry-run
```

The profiles preserve the existing experiment matrices:

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
| `MESSI_DATA_ROOT` | `/vol/tmp/schaefpa/messi_datasets` | main datasets and queries |
| `MESSI_QUERY_ROOT` | main data root | optional separate main query root |
| `MESSI_SEISBENCH_ROOT` | `/vol/tmp/schaefpa/seismic` | SeisBench datasets and queries |
| `MESSI_SEISBENCH_QUERY_ROOT` | SeisBench root | optional separate SeisBench query root |
| `MESSI_LOG_ROOT` | `$HOME/MESSI_logs` | logs produced by MESSI |
| `MESSI_SHELL_LOG_DIR` | `MESSI_LOG_ROOT` | combined runner transcript for each dataset/profile run |
| `MESSI_RESULTS_ROOT` | `$HOME/MESSI_SFA_logs` | archived benchmark results |
| `MESSI_INDEX_ROOT` | `$HOME/index` | index directory cleared by `run_cleanup.sh` |
| `MESSI_DRY_RUN` | `false` | enable dry-run mode through compatibility wrappers |

Relative dataset and query overrides are resolved below the applicable data
root. Absolute overrides are used unchanged. `--dry-run` prints shell-escaped
commands and performs no benchmark or log archival.

For trie runs, `--trie-leaf-ivf 16` adds a flat, post-build 16-list IVF/MRB
directory inside terminal leaves with at least 4 K records. It is disabled by
default and is intended for direct A/B benchmarking with the MBR suffix bound.
Construction clusters eligible leaves independently in parallel, using the
existing `--threads` setting; the build log reports the active worker count.

`--trie-streaming-leaf-scan` provides an A/B alternative to the default
best-first record heap. It computes each record's lower bound and immediately
runs exact distance when that bound passes, so an improved BSF affects the very
next record. Cluster and leaf traversal ordering is unchanged.

Result archival intentionally preserves the historical behavior: an existing
`DATASET/RUN` directory is replaced. Labels are restricted to single path
components, and the destination is checked before removal.

## Validation

Run `scripts/tests/test_scripts.sh` for syntax, argument-generation, wrapper,
and result-archival tests. Set `RUN_MESSI_INTEGRATION=1` to additionally run a
small fixture test with the real `bin/MESSI`; no substitute executable is used.

Files below `scripts/old/` are historical and unsupported. Useful datasets and
query matrices from that directory are represented by the maintained runner.
