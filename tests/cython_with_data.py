#!/usr/bin/env python3
"""Smoke test that loads the bundled astro head dataset and runs a query."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from messi import Index

TS_SIZE = 256


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
        function_type=4,
        histogram_type=2,
        sample_size=1_000)

    idx.add(str(data_path), ts_num=1_000)

    query_count = 10
    queries = np.fromfile(queries_path, dtype=np.float32, count=query_count*TS_SIZE)
    
    if queries.size < query_count*TS_SIZE:
        raise RuntimeError(f"Expected at least {query_count*TS_SIZE} floats for one query, got {queries.size}")

    queries = queries.reshape(query_count, TS_SIZE)

    dists, labels = idx.search(queries, k=3)
    print("distance:", dists)
    #, "labels:", labels


if __name__ == "__main__":
    main()
