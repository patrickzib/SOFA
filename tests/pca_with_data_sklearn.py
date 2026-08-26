#!/usr/bin/env python3
"""PCA test using scikit-learn (LAPACK SVD) on bundled astro data."""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np

try:
    from sklearn.decomposition import PCA
except ImportError as exc:
    raise SystemExit("scikit-learn is required: pip install scikit-learn") from exc

from messi import Index

TS_SIZE = 256
N_SEGMENTS = 16
SAMPLE_SIZE = 10_000
QUERY_COUNT = 100


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    data_path = repo_root / "data_head" / "astro_head.bin"
    if not data_path.exists():
        raise FileNotFoundError(f"Expected sample data at {data_path}")
    queries_path = repo_root / "data_queries" / "astro_queries.bin"
    if not queries_path.exists():
        raise FileNotFoundError(f"Expected query data at {queries_path}")

    idx = Index(
        timeseries_size=TS_SIZE,
        max_query_threads=8,
        function_type=5,
        sample_size=SAMPLE_SIZE,
        histogram_type=2,
        #sample_type=1,
        is_norm=1,
    )

    samples = np.fromfile(data_path, dtype=np.float32, count=SAMPLE_SIZE * TS_SIZE)
    if samples.size < SAMPLE_SIZE * TS_SIZE:
        raise RuntimeError(
            f"Expected at least {SAMPLE_SIZE * TS_SIZE} floats, got {samples.size}"
        )
    samples = samples.reshape(SAMPLE_SIZE, TS_SIZE)

    queries = np.fromfile(queries_path, dtype=np.float32, count=QUERY_COUNT * TS_SIZE)
    if queries.size < QUERY_COUNT * TS_SIZE:
        raise RuntimeError(
            f"Expected at least {QUERY_COUNT * TS_SIZE} floats, got {queries.size}"
        )
    queries = queries.reshape(QUERY_COUNT, TS_SIZE)

    start = time.perf_counter()
    idx.add(str(data_path), ts_num=SAMPLE_SIZE)
    index_seconds = time.perf_counter() - start

    start = time.perf_counter()
    _dists, _labels = idx.search(queries, k=1)
    query_seconds = time.perf_counter() - start

    c_pca_queries = idx.pca_transform(queries)

    # svd_solver="full" uses LAPACK via SciPy/NumPy for exact SVD.
    pca = PCA(n_components=N_SEGMENTS, svd_solver="full")

    start = time.perf_counter()
    pca.fit(samples)
    fit_seconds = time.perf_counter() - start

    start = time.perf_counter()
    query_pca = pca.transform(queries)
    transform_seconds = time.perf_counter() - start

    print(f"index_seconds: {index_seconds:.4f}")
    print(f"query_seconds: {query_seconds:.4f}")
    print(f"fit_seconds: {fit_seconds:.4f}")
    print(f"transform_seconds: {transform_seconds:.4f}")
    print("c_pca_queries:")
    for idx, row in enumerate(c_pca_queries):
        row_str = " ".join(f"{value:.6f}" for value in row)
        print(f"{idx}: {row_str}")
    print("sklearn_pca_queries:")
    for idx, row in enumerate(query_pca):
        row_str = " ".join(f"{value:.6f}" for value in row)
        print(f"{idx}: {row_str}")


if __name__ == "__main__":
    main()
