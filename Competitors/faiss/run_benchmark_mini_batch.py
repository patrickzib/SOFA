import argparse
import os
import time

import numpy as np

def read(fp, dim, data_type=np.float32, count=100, header_bytes=0):
    with open(fp, "rb") as source:
        source.seek(header_bytes)
        a = np.fromfile(source, dtype=data_type, count=dim * count)
    if data_type != np.float32:
        return a.reshape(-1, dim).copy().astype(np.float32, copy=False)
    return a.reshape(-1, dim)


def z_normalize_rows(values, rows_per_batch=8192):
    """Match MESSI's per-series zero-mean, unit-standard-deviation transform."""
    for start in range(0, len(values), rows_per_batch):
        block = values[start:start + rows_per_batch]
        mean = block.mean(axis=1, dtype=np.float32)
        block -= mean[:, None]
        stddev = np.sqrt(np.mean(block * block, axis=1, dtype=np.float32))
        stddev[stddev < 1e-8] = 1.0
        block /= stddev[:, None]


NORMAL_PATH = "/home/tmp/schaefpa/messi_datasets/"
SEISBENCH_PATH = "/home/tmp/schaefpa/seismic/"

datasets = {
    # Vector datasets.  Keep these definitions aligned with
    # scripts/lib/datasets.sh: the benchmark reads the configured 100M-record
    # prefix from each raw collection.
    "BIGANN": ["bigANN.bin", "bigANN_queries.bin", 128, 0, np.uint8, 0, 0],
    "DEEP1b": ["deep1b.bin", "deep1b_queries.bin", 96, 0, np.float32, 0, 0],
    "SIFT1b": ["sift1b.bin", "sift1b_queries.bin", 128, 0, np.float32, 0, 0],
    "simsearchnet": ["SimSearchNet.bin", "SimSearchNet_queries.bin", 256, 0, np.uint8, 8, 8],
    "spacev1b": ["spacev1B.bin", "spacev1B_queries.bin", 100, 0, np.int8, 0, 0],
    "TEXTTOIMAGE": ["text-to-image.bin", "text-to-image_queries.bin", 200, 0, np.float32, 8, 0],
    "turinganns": ["turingANNs.bin", "turingANNs_queries.bin", 100, 0, np.float32, 8, 8],

    # Time-series and SeisBench datasets intentionally remain opt-in here.
    "ASTRO": ["astro.bin", "astro_queries.bin", 256, 0, np.float32, 0, 0],
    "SALD": ["SALD.bin", "SALD_queries.bin", 128, 0, np.float32, 0, 0],
    "SCEDC": ["SCEDC.bin", "SCEDC_queries.bin", 256, 0, np.float32, 0, 0],
    "seismic": ["seismic.bin", "seismic_queries.bin", 256, 0, np.float32, 0, 0],

    # SeisBench
    "ETHC": ["ETHZ.bin", "ETHZ_queries.bin", 256, 1, np.float32, 0, 0],
    "ISC_EHB_DepthPhases": ["ISC_EHB_DepthPhases.bin", "ISC_EHB_DepthPhases_queries.bin", 256, 1, np.float32, 0, 0],
    "LenDB": ["LenDB.bin", "LenDB_queries.bin", 256, 1, np.float32, 0, 0],
    "Iquique": ["Iquique.bin", "Iquique_queries.bin", 256, 1, np.float32, 0, 0],
    "NEIC": ["NEIC.bin", "NEIC_queries.bin", 256, 1, np.float32, 0, 0],
    "OBS": ["OBS.bin", "OBS_queries.bin", 256, 1, np.float32, 0, 0],
    "OBST2024": ["OBST2024.bin", "OBST2024_queries.bin", 256, 1, np.float32, 0, 0],
    "PNW": ["PNW.bin", "PNW_queries.bin", 256, 1, np.float32, 0, 0],
    "Meier2019JGR": ["Meier2019JGR.bin", "Meier2019JGR_queries.bin", 256, 1, np.float32, 0, 0],
    "STEAD": ["STEAD.bin", "STEAD_queries.bin", 256, 1, np.float32, 0, 0],
    "TXED": ["TXED.bin", "TXED_queries.bin", 256, 1, np.float32, 0, 0],
}

DEFAULT_DATASETS = [
    name for name in datasets
    if name not in {"turinganns", "simsearchnet", "TEXTTOIMAGE"}
]

# Mirrors DATASET_SIZE and APPLY_Z_NORM in scripts/lib/datasets.sh.  The
# separate metadata keeps this competitor benchmark comparable to the runner
# even for SeisBench datasets that contain fewer than 100M records.
dataset_properties = {
    "BIGANN": (100_000_000, True), "DEEP1b": (100_000_000, False),
    "SIFT1b": (100_000_000, False), "simsearchnet": (100_000_000, True),
    "spacev1b": (100_000_000, True),
    "TEXTTOIMAGE": (100_000_000, True), "turinganns": (100_000_000, True),
    "ASTRO": (100_000_000, False), "SALD": (100_000_000, False),
    "SCEDC": (100_000_000, False), "seismic": (100_000_000, False),
    "ETHC": (4_999_932, False),
    "ISC_EHB_DepthPhases": (100_000_000, False), "LenDB": (37_345_260, False),
    "Iquique": (578_853, False), "NEIC": (93_473_541, False),
    "OBS": (15_508_794, False), "OBST2024": (4_160_286, False),
    "PNW": (31_982_766, False), "Meier2019JGR": (6_361_998, False),
    "STEAD": (87_323_433, False), "TXED": (35_851_641, False),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Run FAISS FlatL2 mini-batch benchmarks.")
    parser.add_argument(
        "--threads", "--cores", dest="threads", type=int, nargs="+", default=[64],
        metavar="N", help="FAISS/OpenMP worker count(s); default: 64",
    )
    parser.add_argument("--output-dir", default="logs_mini_batch2",
                        help="directory for CSV results (default: %(default)s)")
    parser.add_argument("--data-root", default=os.environ.get("MESSI_DATA_ROOT", NORMAL_PATH))
    parser.add_argument("--seisbench-root", default=os.environ.get("MESSI_SEISBENCH_ROOT", SEISBENCH_PATH))
    parser.add_argument("--datasets", default=",".join(DEFAULT_DATASETS),
                        help="comma-separated dataset IDs; TuringANNS, SimSearchNet, and Text-to-Image are opt-in")
    args = parser.parse_args()
    if any(threads <= 0 for threads in args.threads):
        parser.error("--threads values must be positive")
    args.datasets = [name.strip() for name in args.datasets.split(",") if name.strip()]
    unknown = sorted(set(args.datasets).difference(datasets))
    if unknown:
        parser.error("unknown dataset IDs: " + ", ".join(unknown))
    return args


def run_benchmark(args):
    all_threads = args.threads
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    faiss = None
    k = 1
    for dataset in args.datasets:
        file_data, file_queries, d, path_switch, data_type, data_header, query_header = datasets[dataset]
        path = args.data_root if path_switch == 0 else args.seisbench_root
        data_path = os.path.join(path, file_data)
        query_path = os.path.join(path, file_queries)
        missing = [candidate for candidate in (data_path, query_path) if not os.path.isfile(candidate)]
        if missing:
            print(f"Skipping {dataset}: missing " + ", ".join(missing))
            continue

        if faiss is None:
            try:
                import faiss
            except ImportError as exc:
                raise RuntimeError("FAISS is required; install the Python package before running this benchmark") from exc
        print("Running: ", dataset)
        record_count, apply_znorm = dataset_properties[dataset]
        data = read(data_path, dim=d, data_type=data_type, count=record_count,
                    header_bytes=data_header)
        queries = read(query_path, data_type=data_type, dim=d, count=100,
                       header_bytes=query_header)
        if apply_znorm:
            z_normalize_rows(data)
            z_normalize_rows(queries)

        print("\tData Shape", data.shape)
        print("\tQuery Shape", queries.shape)

        for threads in all_threads:
            faiss.omp_set_num_threads(threads)
            print(f"\t{dataset},{threads}\tBuilding indices...")
            tic = time.time_ns()

            index = faiss.IndexFlatL2(d)
            index.add(data)

            index_time = (time.time_ns() - tic) / 1e6
            # IndexFlatL2 stores an owned float32 copy of every vector.  This
            # is the index's actual payload size, unlike RSS deltas, which can
            # be zero or negative when a previous index is freed or allocator
            # pages are reused between benchmark configurations.
            memory_used = (
                index.ntotal * index.d * np.dtype(np.float32).itemsize
            ) / (1024 ** 2)
            print(f"\tRuntime: {index_time} ms")
            print(f"\tIndex storage: {memory_used:.2f} MiB (float32 payload)")
            index_creation = [[index_time, memory_used]]

            results = []
            start_tic = time.time_ns()
            print("Query, Runtime in ms, True Distance, Total Time in ms")
            for batch_start in range(0, len(queries), threads):
                batch = queries[batch_start:batch_start + threads]
                tic = time.time_ns()
                _, indices = index.search(batch, k)
                time_query = (time.time_ns() - tic) / 1e6
                total_time = (time.time_ns() - start_tic) / 1e6

                for offset, record_id in enumerate(indices[:, 0]):
                    query_id = batch_start + offset
                    distance = np.linalg.norm(data[record_id] - queries[query_id]) ** 2
                    results.append([query_id, time_query / len(batch), distance, total_time])
                    print(f"{query_id}\t {time_query / len(batch):.2f}\t {distance:0.2f}\t {total_time:.2f}")

            suffix = f"{dataset}_{threads}.csv"
            np.savetxt(os.path.join(output_dir, "index_" + suffix), index_creation,
                       delimiter=",", header="index creation time in ms, index storage in MiB", fmt="%s")
            np.savetxt(os.path.join(output_dir, "queries_" + suffix), results,
                       delimiter=",", header="query,time in ms,distance, total", fmt="%s")


def main():
    args = parse_args()
    run_benchmark(args)


if __name__ == "__main__":
    main()
