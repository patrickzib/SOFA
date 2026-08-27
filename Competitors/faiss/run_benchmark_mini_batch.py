import argparse
import os
import time

import psutil

import numpy as np

def read(fp, dim, data_type=np.float32, count=100):
    if data_type != np.float32:
        a = np.fromfile(fp, dtype=data_type, count=dim*count)
        return a.reshape(-1, dim).copy().astype(np.float32, copy=False)
    else:
        return np.fromfile(fp, dtype=np.float32, count=dim*count).reshape(-1, dim)


def z_normalize_rows(values, rows_per_batch=8192):
    """Match MESSI's per-series zero-mean, unit-standard-deviation transform."""
    for start in range(0, len(values), rows_per_batch):
        block = values[start:start + rows_per_batch]
        mean = block.mean(axis=1, dtype=np.float32)
        block -= mean[:, None]
        stddev = np.sqrt(np.mean(block * block, axis=1, dtype=np.float32))
        stddev[stddev < 1e-8] = 1.0
        block /= stddev[:, None]


NORMAL_PATH = "/vol/tmp/schaefpa/messi_datasets/"
SEISBENCH_PATH = "/vol/tmp/schaefpa/seismic/"

datasets = {
    # Vector datasets.  Keep these definitions aligned with
    # scripts/lib/datasets.sh: the benchmark reads the configured 100M-record
    # prefix from each raw collection.
    "BIGANN": ["bigANN.bin", "bigANN_queries.bin", 100, 0, np.uint8],
    "DEEP1B": ["deep1b.bin", "deep1b_queries.bin", 96, 0, np.float32],
    "SIFT1B": ["sift1b.bin", "sift1b_queries.bin", 128, 0, np.float32],
    "SPACEV1B": ["spacev1B.bin", "spacev1B_queries.bin", 100, 0, np.uint8],
    "TEXTTOIMAGE": ["text-to-image.bin", "text-to-image_queries.bin", 200, 0, np.float32],
    "TURINGANNS": ["turingANNs.bin", "turingANNs_queries.bin", 100, 0, np.uint8],

    # Time-series and SeisBench datasets intentionally remain opt-in here.
    "ASTRO": ["astro.bin", "astro_queries.bin", 256, 0, np.float32],
    "SALD": ["SALD.bin", "SALD_queries.bin", 128, 0, np.float32],
    "SCEDC": ["SCEDC.bin", "SCEDC_queries.bin", 256, 0, np.float32],

    # SeisBench
    "ETHZ": ["ETHZ.bin", "ETHZ_queries.bin", 256, 1, np.float32],
    "ISC_EHB_DepthPhases": ["ISC_EHB_DepthPhases.bin", "ISC_EHB_DepthPhases_queries.bin", 256, 1, np.float32],
    "LenDB": ["LenDB.bin", "LenDB_queries.bin", 256, 1, np.float32],
    "Iquique": ["Iquique.bin", "Iquique_queries.bin", 256, 1, np.float32],
    "NEIC": ["NEIC.bin", "NEIC_queries.bin", 256, 1, np.float32],
    "OBS": ["OBS.bin", "OBS_queries.bin", 256, 1, np.float32],
    "OBST2024": ["OBST2024.bin", "OBST2024_queries.bin", 256, 1, np.float32],
    "PNW": ["PNW.bin", "PNW_queries.bin", 256, 1, np.float32],
    "Meier2019JGR": ["Meier2019JGR.bin", "Meier2019JGR_queries.bin", 256, 1, np.float32],
    "STEAD": ["STEAD.bin", "STEAD_queries.bin", 256, 1, np.float32],
    "TXED": ["TXED.bin", "TXED_queries.bin", 256, 1, np.float32],
}

# Mirrors DATASET_SIZE and APPLY_Z_NORM in scripts/lib/datasets.sh.  The
# separate metadata keeps this competitor benchmark comparable to the runner
# even for SeisBench datasets that contain fewer than 100M records.
dataset_properties = {
    "BIGANN": (100_000_000, True), "DEEP1B": (100_000_000, False),
    "SIFT1B": (100_000_000, False), "SPACEV1B": (100_000_000, True),
    "TEXTTOIMAGE": (100_000_000, True), "TURINGANNS": (100_000_000, True),
    "ASTRO": (100_000_000, False), "SALD": (100_000_000, False),
    "SCEDC": (100_000_000, False), "ETHZ": (4_999_932, False),
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
    args = parser.parse_args()
    if any(threads <= 0 for threads in args.threads):
        parser.error("--threads values must be positive")
    return args


def run_benchmark(all_threads, output_dir):
    try:
        import faiss
    except ImportError as exc:
        raise RuntimeError("FAISS is required; install the Python package before running this benchmark") from exc
    os.makedirs(output_dir, exist_ok=True)
    k = 1
    for dataset, (file_data, file_queries, d, path_switch, data_type) in datasets.items():
        print("Running: ", dataset)
        path = NORMAL_PATH if path_switch == 0 else SEISBENCH_PATH
        record_count, apply_znorm = dataset_properties[dataset]
        data = read(path + file_data, dim=d, data_type=data_type, count=record_count)
        queries = read(path + file_queries, data_type=data_type, dim=d, count=100)
        if apply_znorm:
            z_normalize_rows(data)
            z_normalize_rows(queries)

        print("\tData Shape", data.shape)
        print("\tQuery Shape", queries.shape)

        for threads in all_threads:
            faiss.omp_set_num_threads(threads)
            print(f"\t{dataset},{threads}\tBuilding indices...")
            tic = time.time_ns()
            process = psutil.Process(os.getpid())
            mem_before = process.memory_info().rss

            index = faiss.IndexFlatL2(d)
            index.parallel_mode = 0
            index.add(data)

            mem_after = process.memory_info().rss
            index_time = (time.time_ns() - tic) / 1e6
            memory_used = (mem_after - mem_before) / 1048576
            print(f"\tRuntime: {index_time} ms")
            print(f"\tMemory: {memory_used} MB (rough estimate)")
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
                       delimiter=",", header="index creation time in ms, memory in MB", fmt="%s")
            np.savetxt(os.path.join(output_dir, "queries_" + suffix), results,
                       delimiter=",", header="query,time in ms,distance, total", fmt="%s")


def main():
    args = parse_args()
    run_benchmark(args.threads, args.output_dir)


if __name__ == "__main__":
    main()
