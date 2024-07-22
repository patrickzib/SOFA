import os
import time

import psutil

import faiss
import numpy as np

def read(fp, dim, data_type=np.float32, count=100):
    if data_type != np.float32:
        a = np.fromfile(fp, dtype=data_type, count=dim*count)
        return a.reshape(-1, dim).copy().astype(np.float32, copy=False)
    else:
        return np.fromfile(fp, dtype=np.float32, count=dim*count).reshape(-1, dim)


NORMAL_PATH = "/vol/tmp/schaefpa/messi_datasets/"
SEISBENCH_PATH = "/vol/tmp/schaefpa/seismic/"

datasets = {
    # Other DS
    #"ASTRO": ["astro.bin", "astro_queries.bin", 256, 0, np.float32],
    "BIGANN": ["bigANN.bin", "bigANN_queries.bin", 100, 0, np.int8],
    #"SALD": ["SALD.bin", "SALD_queries.bin", 128, 0, np.float32],
    #"SIFT1B": ["sift1b.bin", "sift1b_queries.bin", 128, 0, np.float32],
    #"DEPP1B": ["deep1b.bin", "deep1b_queries.bin", 96, 0, np.float32],
    #"SCEDC": ["SCEDC.bin", "SCEDC_queries.bin", 256, 0, np.float32],

    # Seisbench
    #"ETHZ": ["ETHZ.bin", "ETHZ_queries.bin", 256, 1, np.float32],
    #"ISC_EHB_DepthPhases": ["ISC_EHB_DepthPhases.bin", "ISC_EHB_DepthPhases_queries.bin", 256, 1, np.float32],
    #"LenDB": ["LenDB.bin", "LenDB_queries.bin", 256, 1, np.float32],
    #"Iquique": ["Iquique.bin", "Iquique_queries.bin", 256, 1, np.float32],
    #"NEIC": ["NEIC.bin", "NEIC_queries.bin", 256, 1, np.float32],
    #"OBS": ["OBS.bin", "OBS_queries.bin", 256, 1, np.float32],
    #"OBST2024": ["OBST2024.bin", "OBST2024_queries.bin", 256, 1, np.float32],
    #"PNW": ["PNW.bin", "PNW_queries.bin", 256, 1, np.float32],
    #"Meier2019JGR": ["Meier2019JGR.bin", "Meier2019JGR_queries.bin", 256, 1, np.float32],
    #"STEAD": ["STEAD.bin", "STEAD_queries.bin", 256, 1, np.float32],
    #"TXED": ["TXED.bin", "TXED_queries.bin", 256, 1, np.float32],
}
all_threads = [9, 18, 36]

k = 1

for dataset in datasets:
    print("Running: ", dataset)
    file_data, file_queries, d, path_switch, data_type = datasets[dataset]

    path = NORMAL_PATH if path_switch == 0 else SEISBENCH_PATH
    data = read(path + file_data, dim=d, data_type=data_type, count=100_000_000)  # data
    queries = read(path + file_queries, data_type=data_type, dim=d, count=100)  # queries

    print("\tData Shape", data.shape)
    print("\tQuery Shape", queries.shape)

    for threads in all_threads:
        faiss.omp_set_num_threads(threads)

        print(f"\t{dataset},{threads}\tBuilding indices...")
        tic = time.time_ns()

        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss

        # index = faiss.IndexFlatIP(d)
        index = faiss.IndexFlatL2(d)
        index.parallel_mode = 0
        index.add(data)

        mem_after = process.memory_info().rss

        index_time = (time.time_ns() - tic) / 1e6
        memory_used = (mem_after - mem_before) / (1048576)
        print(f"\tRuntime: {index_time} ms")
        print(f"\tMemory: {memory_used} MB (rough estimate)")

        index_creation = [[index_time, memory_used]]

        results = []
        start_tic = time.time_ns()

        print(f"Query, Runtime in ms, True Distance, Total Time in ms")

        # Batch mode
        """
        tic = time.time_ns()
        D, I = index.search(queries, k)
        end = time.time_ns()
        total_time = (end - start_tic) / 1e6

        for i in range(len(queries)):
            distance = np.linalg.norm(data[I[i]] - queries[i])**2
            results.append([i, total_time/len(queries), distance, total_time])
            print(f"{i}\t {0:.2f}\t {distance:0.2f}\t {total_time:.2f}")
        """

        # Small Batch query mode
        for i in range(np.int32(np.ceil(len(queries) / threads))):
            i_start = i*threads
            i_end = min((i+1)*threads, len(queries))
            batch = queries[i_start:i_end]

            tic = time.time_ns()
            D, I = index.search(batch, k)
            end = time.time_ns()
            time_query = (end - tic) / 1e6

            total_time = (end - start_tic) / 1e6

            for j in range(len(batch)):
                distance = np.linalg.norm(data[I[j]] - queries[i_start+j]) ** 2
                results.append([i_start+j, time_query / len(batch), distance, total_time])
                print(f"{i_start+j}\t {time_query / len(batch):.2f}\t {distance:0.2f}\t {total_time:.2f}")


        np.savetxt("logs_mini_batch2/index_" + dataset + "_" + str(threads) + ".csv", index_creation,
                   delimiter=",", header='index creation time in ms, memory in MB',
                   fmt="%s")

        np.savetxt("logs_mini_batch2/queries_" + dataset + "_" + str(threads) + ".csv", results,
                   delimiter=",", header='query,time in ms,distance, total',
                   fmt="%s")
