#!/usr/bin/env python3
"""Regression tests for the public s3trie Python API."""

from __future__ import annotations

import unittest

import numpy as np

from s3trie import Index


class PresetTests(unittest.TestCase):
    def config_for(self, **kwargs):
        with Index(timeseries_size=128, **kwargs) as index:
            return index.config

    def test_named_presets(self):
        messi = self.config_for(index="messi")
        self.assertEqual((messi["transform"], messi["layout"]), ("sax", "isax"))

        sofa = self.config_for(index="sofa")
        self.assertEqual((sofa["transform"], sofa["layout"]), ("sfa", "isax"))

        s3trie = self.config_for(index="s3trie")
        self.assertEqual((s3trie["transform"], s3trie["layout"]), ("spartan", "trie"))
        self.assertEqual(s3trie["record_lb_dimensions"], 64)
        self.assertEqual(s3trie["transform_dimensions"], 128)
        self.assertEqual(s3trie["trie_split_dimensions"], 64)
        self.assertEqual(s3trie["max_leaf_size"], 20_000)
        self.assertEqual(s3trie["min_leaf_size"], 20_000)
        self.assertEqual(s3trie["trie_leaf_ivf"], 16)
        self.assertEqual(s3trie["trie_leaf_ivf_min_size"], 4096)
        self.assertTrue(s3trie["trie_leaf_ivf_raw_ball_bound"])
        self.assertTrue(s3trie["trie_leaf_ivf_radial_bound"])
        self.assertTrue(s3trie["trie_record_mbr_suffix_bound"])
        self.assertTrue(s3trie["trie_streaming_leaf_scan"])
        self.assertTrue(s3trie["trie_residual_record_only"])
        self.assertEqual(s3trie["trie_residual_order"], "symbolic-first")

    def test_s3trie_options_remain_overrideable(self):
        config = self.config_for(
            index="s3trie",
            trie_leaf_ivf=0,
            trie_residual_record_only=False,
            max_leaf_size=100,
            min_leaf_size=50,
        )
        self.assertEqual(config["trie_leaf_ivf"], 0)
        self.assertFalse(config["trie_leaf_ivf_raw_ball_bound"])
        self.assertFalse(config["trie_leaf_ivf_radial_bound"])
        self.assertFalse(config["trie_residual_record_only"])
        self.assertEqual(config["max_leaf_size"], 100)
        self.assertEqual(config["min_leaf_size"], 50)

    def test_invalid_or_conflicting_selection(self):
        with self.assertRaisesRegex(ValueError, "index must be"):
            Index(timeseries_size=128, index="unknown")
        with self.assertRaisesRegex(ValueError, "conflicts with layout"):
            Index(timeseries_size=128, index="sofa", layout="trie")
        with self.assertRaisesRegex(ValueError, "record bound width"):
            Index(
                timeseries_size=128,
                layout="trie",
                transform="spartan",
                trie_split_dimensions=32,
            )

    def test_advanced_selection_remains_available(self):
        config = self.config_for(
            layout="trie",
            transform="pisa",
            n_segments=32,
            trie_mbr_dimensions=64,
            trie_split_dimensions=32,
        )
        self.assertIsNone(config["index"])
        self.assertEqual((config["transform"], config["layout"]), ("pisa", "trie"))


class ArrayApiTests(unittest.TestCase):
    def test_s3trie_exact_self_query(self):
        data = np.random.default_rng(7).normal(size=(64, 64)).astype(np.float32)
        with Index(
            timeseries_size=64,
            index="s3trie",
            sample_size=64,
            max_leaf_size=64,
            min_leaf_size=64,
            initial_leaf_buffer_size=64,
            trie_leaf_ivf=0,
        ) as index:
            index.add(data)
            distances, ids = index.search(data[:4])

        self.assertEqual(distances.shape, (4, 1))
        self.assertEqual(ids.shape, (4, 1))
        np.testing.assert_allclose(distances, 0.0, atol=1e-5)
        np.testing.assert_array_equal(ids[:, 0], np.arange(4))


if __name__ == "__main__":
    unittest.main()
