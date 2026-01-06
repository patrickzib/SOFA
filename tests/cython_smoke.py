#!/usr/bin/env python3
"""Minimal script to verify the Cython binding loads and runs."""

from __future__ import annotations

import numpy as np

from messi import Index


def main() -> None:
    idx = Index(
        timeseries_size=256,
        sax_bit_cardinality=4,
        max_leaf_size=16,
        min_leaf_size=16,
        initial_leaf_buffer_size=16,
        max_total_buffer_size=16,
        initial_fbl_buffer_size=16,
        root_directory="tmp/messi-smoke",
    )


    queries = np.zeros((1, idx._dim), dtype=np.float32)
    try:
        idx.search(queries, k=1)
    except RuntimeError as exc:
        print("Search failed (expected until data is indexed):", exc)
    else:
        print("Search completed – index contains data.")


if __name__ == "__main__":
    main()
