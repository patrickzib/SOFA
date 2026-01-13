#!/usr/bin/env python3
"""Smoke test that loads the bundled astro head dataset and runs a query."""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np

from messi import Index

TS_SIZE = 256
N_SEGMENTS = 64


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    data_path = repo_root / "data_head" / "astro_head.bin"
    if not data_path.exists():
        raise FileNotFoundError(f"Expected sample data at {data_path}")
    queries_path = repo_root / "data_queries" / "astro_queries.bin"
    if not queries_path.exists():
        raise FileNotFoundError(f"Expected query data at {queries_path}")

    sample_size = 10_000
    idx = Index(
        timeseries_size=TS_SIZE,
        max_query_threads=8,
        function_type=5,
        sample_size=sample_size,
        histogram_type=2,
        #sample_type=1,
        is_norm=1
        )

    samples = np.fromfile(data_path, dtype=np.float32, count=sample_size * TS_SIZE)
    if samples.size < sample_size * TS_SIZE:
        raise RuntimeError(f"Expected at least {sample_size * TS_SIZE} floats, got {samples.size}")    
    samples = samples.reshape(sample_size, TS_SIZE)
    
    """
    mean = samples.astype(np.float32).sum(axis=0, dtype=np.float32) / sample_size
    centered = samples.astype(np.float64) - mean.astype(np.float64)
    cov = (centered.T @ centered) / (sample_size - 1)
    evals, _ = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]

    print("numpy_pca_eigenvalues_sorted:")
    for idx_pos, eig in zip(order[:N_SEGMENTS], evals[order][:N_SEGMENTS]):
        print(f"{idx_pos}, ({eig:.4f})", end=" ")
    print()
    """

    start = time.perf_counter()
    idx.add(str(data_path), ts_num=sample_size)
    index_seconds = time.perf_counter() - start

    query_count = 100
    queries = np.fromfile(queries_path, dtype=np.float32, count=query_count*TS_SIZE)
    
    if queries.size < query_count*TS_SIZE:
        raise RuntimeError(f"Expected at least {query_count*TS_SIZE} floats for one query, got {queries.size}")

    queries = queries.reshape(query_count, TS_SIZE)

    start = time.perf_counter()
    dists, labels = idx.search(queries, k=1)
    query_seconds = time.perf_counter() - start
    # print("distance:", dists)
    print(f"index_seconds: {index_seconds:.4f}")
    print(f"query_seconds: {query_seconds:.4f}")
    #, "labels:", labels

if __name__ == "__main__":
    main()
