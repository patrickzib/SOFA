#!/usr/bin/env python3
"""Exercise the managed NumPy-array ingestion path when the extension is built."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np

from messi import Index


def main() -> None:
    data = np.arange(32 * 16, dtype=np.float32).reshape(32, 16)
    query = data[[7]].copy()
    with tempfile.TemporaryDirectory(prefix="messi-python-api-") as directory:
        with Index(timeseries_size=16, layout="trie", transform="sax",
                   n_segments=16, max_leaf_size=16,
                   root_directory=directory) as index:
            index.add_array(data[:, ::-1][:, ::-1], storage_dir=directory)
            raw_path = index.raw_data_path
            assert raw_path is not None and Path(raw_path).exists()
            distances, indices = index.search(query)
            assert distances.shape == (1, 1)
            assert distances[0, 0] == 0.0
            assert indices.shape == (1, 1)
            assert indices[0, 0] == 0
            distances, indices = index.search(query, k=2)
            assert distances.shape == (1, 2)
            assert indices.shape == (1, 2)
            assert np.array_equal(distances, [[0.0, 0.0]])
            assert np.array_equal(indices, [[0, 1]])
        assert not os.path.exists(raw_path)


if __name__ == "__main__":
    main()
