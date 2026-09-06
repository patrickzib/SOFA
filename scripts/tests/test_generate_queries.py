import importlib.util
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).parents[2] / "notebooks" / "generate_queries.py"
SPEC = importlib.util.spec_from_file_location("generate_queries", MODULE_PATH)
generate_queries = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = generate_queries
SPEC.loader.exec_module(generate_queries)


class GenerateQueriesTest(unittest.TestCase):
    def test_text_to_image_is_not_generated_by_default(self):
        self.assertNotIn("text-to-image", generate_queries.DEFAULT_DATASETS)

    def test_spacev_preserves_signed_int8_and_writes_no_header(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            values = np.arange(-100, 100, dtype=np.int8)
            values.tofile(root / "spacev1B.bin")
            spec = generate_queries.DatasetSpec(100, np.dtype("i1"), 0, (0.0,))

            generate_queries.generate_dataset(root, "spacev1B", spec, seed=1)

            output = root / "generated" / "spacev1B_noise_0.bin"
            actual = np.fromfile(output, dtype=np.int8)
            np.testing.assert_array_equal(actual, values)
            self.assertEqual(output.stat().st_size, values.nbytes)

    def test_header_validation_allows_prefix_only_when_configured(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "vectors.bin"
            values = np.arange(200, dtype="<f4")
            with path.open("wb") as output:
                output.write(struct.pack("<II", 10, 100))
                values.tofile(output)

            prefix_spec = generate_queries.DatasetSpec(100, np.dtype("<f4"), 8, (0.1,))
            strict_spec = generate_queries.DatasetSpec(
                100, np.dtype("<f4"), 8, (0.1,), True
            )
            np.testing.assert_array_equal(
                generate_queries.read_source(path, prefix_spec), values
            )
            with self.assertRaisesRegex(ValueError, "declares 10 vectors"):
                generate_queries.read_source(path, strict_spec)


if __name__ == "__main__":
    unittest.main()
