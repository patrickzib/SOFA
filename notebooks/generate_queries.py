"""Generate reproducible noisy queries from the benchmark base collections.

The output files are dense, headerless arrays in the same dtype as their
source. This matches MESSI's generated-query suite, which passes
``--query-header-bytes 0`` explicitly.
"""

from __future__ import annotations

import argparse
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np


QUERY_VALUE_BUDGET = 1024 * 1024


@dataclass(frozen=True)
class DatasetSpec:
    dimensions: int
    dtype: np.dtype
    header_bytes: int
    noise_levels: tuple[float, ...]
    require_complete_declared_file: bool = False


DATASETS = {
    # Canonical Big-ANN encodings. The local TuringANNS file must be the
    # complete float32 xbin distribution; this deliberately rejects the known
    # 100 GB file whose payload is a copy of SpaceV.
    "turingANNs": DatasetSpec(
        100, np.dtype("<f4"), 8, (0.1, 0.25, 0.5), True
    ),
    # The local Text-to-Image base is a usable prefix of the official file, so
    # accept a shorter payload after validating that complete vectors exist.
    "text-to-image": DatasetSpec(
        200, np.dtype("<f4"), 8, (0.05, 0.1, 0.25)
    ),
    "spacev1B": DatasetSpec(
        100, np.dtype("i1"), 0, (0.25, 0.5, 1.0)
    ),
}
DEFAULT_DATASETS = [name for name in DATASETS if name != "turingANNs"]


def read_source(path: Path, spec: DatasetSpec) -> np.ndarray:
    declared_records = None
    with path.open("rb") as source:
        if spec.header_bytes:
            header = source.read(spec.header_bytes)
            if len(header) != spec.header_bytes or spec.header_bytes != 8:
                raise ValueError(f"{path}: expected a complete 8-byte xbin header")
            declared_records, dimensions = struct.unpack("<II", header)
            if dimensions != spec.dimensions:
                raise ValueError(
                    f"{path}: header declares dimension {dimensions}, "
                    f"expected {spec.dimensions}"
                )

        payload_bytes = path.stat().st_size - spec.header_bytes
        record_bytes = spec.dimensions * spec.dtype.itemsize
        complete_records, trailing_bytes = divmod(payload_bytes, record_bytes)
        if complete_records == 0:
            raise ValueError(f"{path}: contains no complete vectors")
        if trailing_bytes:
            print(
                f"warning: {path} ends with {trailing_bytes} trailing payload bytes; "
                "ignoring the incomplete final vector"
            )
        if declared_records is not None and complete_records != declared_records:
            message = (
                f"{path}: header declares {declared_records:,} vectors but the file "
                f"contains {complete_records:,} complete vectors"
            )
            if spec.require_complete_declared_file:
                raise ValueError(message)
            print(f"warning: {message}; generating from the available prefix")

        value_count = min(
            (QUERY_VALUE_BUDGET // spec.dimensions) * spec.dimensions,
            complete_records * spec.dimensions,
        )
        data = np.fromfile(source, dtype=spec.dtype, count=value_count)

    if data.size != value_count:
        raise ValueError(f"{path}: short read while loading query source vectors")
    return data


def add_gaussian_noise(
    data: np.ndarray, dtype: np.dtype, noise_level: float, rng: np.random.Generator
) -> np.ndarray:
    # Preserve the historical scale definition while making generation
    # deterministic. Float64 avoids precision loss during perturbation.
    values = data.astype(np.float64)
    noisy = values + rng.normal(0.0, np.std(values) * noise_level, size=values.size)
    if np.issubdtype(dtype, np.integer):
        limits = np.iinfo(dtype)
        noisy = np.clip(np.rint(noisy), limits.min, limits.max)
    return noisy.astype(dtype)


def generate_dataset(data_root: Path, name: str, spec: DatasetSpec, seed: int) -> None:
    input_path = data_root / f"{name}.bin"
    data = read_source(input_path, spec)
    output_dir = data_root / "generated"
    output_dir.mkdir(parents=True, exist_ok=True)

    for noise_index, noise_level in enumerate(spec.noise_levels):
        # Independent stable streams make a single noise level reproducible
        # even if levels are reordered or added later.
        rng = np.random.default_rng(np.random.SeedSequence([seed, noise_index]))
        noisy = add_gaussian_noise(data, spec.dtype, noise_level, rng)
        suffix = f"{noise_level:g}".replace(".", "")
        output_path = output_dir / f"{name}_noise_{suffix}.bin"
        with output_path.open("wb") as output:
            noisy.tofile(output)
        print(
            f"wrote {noisy.size // spec.dimensions:,} headerless "
            f"{spec.dtype.name} queries to {output_path}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root", type=Path, default=Path("."),
        help="directory containing base .bin files (default: current directory)",
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--datasets", nargs="+", choices=sorted(DATASETS), default=DEFAULT_DATASETS,
        help="collections to generate (default: all except TuringANNS)",
    )
    args = parser.parse_args()

    for name in args.datasets:
        spec = DATASETS[name]
        print(
            f"{name}: dimension={spec.dimensions} dtype={spec.dtype.name} "
            f"source_header={spec.header_bytes} bytes"
        )
        generate_dataset(args.data_root, name, spec, args.seed)


if __name__ == "__main__":
    main()
