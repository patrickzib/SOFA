#!/usr/bin/env python3
"""Run a scikit-learn NearestNeighbors squared-L2 benchmark.

Install dependencies with:
    python -m pip install numpy scikit-learn

``algorithm='auto'`` lets scikit-learn choose its backend.  For the project's
high-dimensional workloads that will normally be exact brute-force search.
``n_jobs`` parallelizes query work; scikit-learn does not provide an IVF index.
The raw data handling, z-normalization, and CSV result layout match the FAISS
benchmark in ../faiss/run_benchmark_mini_batch.py.
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np


NORMAL_PATH = "/vol/tmp/schaefpa/messi_datasets"
SEISBENCH_PATH = "/vol/tmp/schaefpa/seismic"

DATASETS = {
    "BIGANN": ("bigANN.bin", "bigANN_queries.bin", 100, False, np.uint8),
    "DEEP1b": ("deep1b.bin", "deep1b_queries.bin", 96, False, np.float32),
    "SIFT1b": ("sift1b.bin", "sift1b_queries.bin", 128, False, np.float32),
    "spacev1b": ("spacev1B.bin", "spacev1B_queries.bin", 100, False, np.uint8),
    "turinganns": ("turingANNs.bin", "turingANNs_queries.bin", 100, False, np.uint8),
    "ASTRO": ("astro.bin", "astro_queries.bin", 256, False, np.float32),
    "SALD": ("SALD.bin", "SALD_queries.bin", 128, False, np.float32),
    "SCEDC": ("SCEDC.bin", "SCEDC_queries.bin", 256, False, np.float32),
    "ETHC": ("ETHZ.bin", "ETHZ_queries.bin", 256, True, np.float32),
    "ISC_EHB_DepthPhases": ("ISC_EHB_DepthPhases.bin", "ISC_EHB_DepthPhases_queries.bin", 256, True, np.float32),
    "LenDB": ("LenDB.bin", "LenDB_queries.bin", 256, True, np.float32),
    "Iquique": ("Iquique.bin", "Iquique_queries.bin", 256, True, np.float32),
    "NEIC": ("NEIC.bin", "NEIC_queries.bin", 256, True, np.float32),
    "OBS": ("OBS.bin", "OBS_queries.bin", 256, True, np.float32),
    "OBST2024": ("OBST2024.bin", "OBST2024_queries.bin", 256, True, np.float32),
    "PNW": ("PNW.bin", "PNW_queries.bin", 256, True, np.float32),
    "Meier2019JGR": ("Meier2019JGR.bin", "Meier2019JGR_queries.bin", 256, True, np.float32),
    "STEAD": ("STEAD.bin", "STEAD_queries.bin", 256, True, np.float32),
    "TXED": ("TXED.bin", "TXED_queries.bin", 256, True, np.float32),
}

DATASET_PROPERTIES = {
    "BIGANN": (100_000_000, True), "DEEP1b": (100_000_000, False),
    "SIFT1b": (100_000_000, False), "spacev1b": (100_000_000, True),
    "turinganns": (100_000_000, True), "ASTRO": (100_000_000, False),
    "SALD": (100_000_000, False), "SCEDC": (100_000_000, False),
    "ETHC": (4_999_932, False), "ISC_EHB_DepthPhases": (100_000_000, False),
    "LenDB": (37_345_260, False), "Iquique": (578_853, False),
    "NEIC": (93_473_541, False), "OBS": (15_508_794, False),
    "OBST2024": (4_160_286, False), "PNW": (31_982_766, False),
    "Meier2019JGR": (6_361_998, False), "STEAD": (87_323_433, False),
    "TXED": (35_851_641, False),
}


def read_vectors(path, dimensions, dtype, count):
    values = np.fromfile(path, dtype=dtype, count=dimensions * count)
    if values.size == 0 or values.size % dimensions:
        raise ValueError(f"{path}: expected complete {dimensions}-value rows")
    return values.reshape(-1, dimensions).astype(np.float32, copy=False)


def z_normalize_rows(values, rows_per_batch=8192):
    for start in range(0, len(values), rows_per_batch):
        block = values[start:start + rows_per_batch]
        block -= block.mean(axis=1, dtype=np.float32)[:, None]
        stddev = np.sqrt(np.mean(block * block, axis=1, dtype=np.float32))
        stddev[stddev < 1e-8] = 1.0
        block /= stddev[:, None]


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark scikit-learn NearestNeighbors.")
    parser.add_argument("--data-root", default=os.environ.get("MESSI_DATA_ROOT", NORMAL_PATH))
    parser.add_argument("--seisbench-root", default=os.environ.get("MESSI_SEISBENCH_ROOT", SEISBENCH_PATH))
    parser.add_argument("--datasets", default=",".join(DATASETS), help="comma-separated dataset IDs")
    parser.add_argument("--threads", type=int, nargs="+", default=[64], metavar="N")
    parser.add_argument("--output-dir", default="logs_sklearn_nearest_neighbors")
    parser.add_argument("--record-count", type=int, help="override the configured indexed prefix")
    parser.add_argument("--query-count", type=int, default=100)
    parser.add_argument("--k", type=int, default=1)
    args = parser.parse_args()
    if any(value <= 0 for value in args.threads) or args.query_count <= 0 or args.k <= 0:
        parser.error("--threads, --query-count, and --k must be positive")
    args.datasets = [name.strip() for name in args.datasets.split(",") if name.strip()]
    unknown = sorted(set(args.datasets).difference(DATASETS))
    if unknown:
        parser.error("unknown dataset IDs: " + ", ".join(unknown))
    return args


def run_dataset(dataset, args):
    from sklearn.neighbors import NearestNeighbors

    data_file, query_file, dimensions, is_seisbench, dtype = DATASETS[dataset]
    root = Path(args.seisbench_root if is_seisbench else args.data_root)
    data_path, query_path = root / data_file, root / query_file
    missing = [str(path) for path in (data_path, query_path) if not path.is_file()]
    if missing:
        print(f"Skipping {dataset}: missing " + ", ".join(missing))
        return
    configured_count, apply_znorm = DATASET_PROPERTIES[dataset]
    data = read_vectors(data_path, dimensions, dtype, args.record_count or configured_count)
    queries = read_vectors(query_path, dimensions, dtype, args.query_count)
    if apply_znorm:
        z_normalize_rows(data)
        z_normalize_rows(queries)
    if args.k > len(data):
        raise ValueError(f"{dataset}: --k ({args.k}) exceeds {len(data)} loaded records")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Running {dataset}: data={data.shape}, queries={queries.shape}")
    for threads in args.threads:
        build_started = time.perf_counter()
        index = NearestNeighbors(
            n_neighbors=args.k, algorithm="auto", metric="euclidean", n_jobs=threads
        ).fit(data)
        build_ms = (time.perf_counter() - build_started) * 1000
        # NearestNeighbors retains a reference to this float32 data matrix.
        storage_mib = data.nbytes / 1024 ** 2

        query_started = time.perf_counter()
        distances, indices = index.kneighbors(queries, return_distance=True)
        query_ms = (time.perf_counter() - query_started) * 1000
        rows = [[query_id, query_ms / len(queries), distances[query_id, 0] ** 2, query_ms]
                for query_id in range(len(queries))]
        suffix = f"{dataset}_{threads}.csv"
        np.savetxt(output_dir / ("index_" + suffix), [[build_ms, storage_mib]], delimiter=",",
                   header="index creation time in ms, index storage in MiB", fmt="%s")
        np.savetxt(output_dir / ("queries_" + suffix), rows, delimiter=",",
                   header="query,time in ms,distance, total", fmt="%s")
        print(f"  threads={threads}: {query_ms:.2f} ms total; "
              f"{len(queries) / (query_ms / 1000):.2f} queries/s; "
              f"backend={getattr(index, '_fit_method', 'unknown')}")


def main():
    args = parse_args()
    for dataset in args.datasets:
        run_dataset(dataset, args)


if __name__ == "__main__":
    main()
