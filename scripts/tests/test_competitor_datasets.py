#!/usr/bin/env python3
"""Check that competitor dataset metadata stays aligned with MESSI."""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).parents[2]


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


FAISS = load_module(
    "competitor_faiss", "Competitors/faiss/run_benchmark_mini_batch.py"
)
SKLEARN = load_module(
    "competitor_sklearn",
    "Competitors/sklearn/run_benchmark_nearest_neighbors.py",
)


class CompetitorDatasetMetadataTests(unittest.TestCase):
    def test_new_datasets_match_messi_formats(self):
        self.assertEqual(
            FAISS.datasets["simsearchnet"],
            ["SimSearchNet.bin", "SimSearchNet_queries.bin", 256, 0, np.uint8, 8, 8],
        )
        self.assertEqual(
            SKLEARN.DATASETS["simsearchnet"],
            ("SimSearchNet.bin", "SimSearchNet_queries.bin", 256, False, np.uint8, 8, 8),
        )
        self.assertEqual(
            FAISS.datasets["seismic"],
            ["seismic.bin", "seismic_queries.bin", 256, 0, np.float32, 0, 0],
        )
        self.assertEqual(
            SKLEARN.DATASETS["seismic"],
            ("seismic.bin", "seismic_queries.bin", 256, False, np.float32, 0, 0),
        )

    def test_dataset_defaults_and_100m_prefixes(self):
        for name in ("seismic", "simsearchnet", "TEXTTOIMAGE"):
            self.assertNotIn(name, FAISS.DEFAULT_DATASETS)
            self.assertNotIn(name, SKLEARN.DEFAULT_DATASETS)
        for name in ("simsearchnet", "seismic"):
            self.assertEqual(FAISS.dataset_properties[name][0], 100_000_000)
            self.assertEqual(SKLEARN.DATASET_PROPERTIES[name][0], 100_000_000)


if __name__ == "__main__":
    unittest.main()
