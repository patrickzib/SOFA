import os
import numpy as np
import pandas as pd
import fnmatch

import sys

sys.path.append("../")

default_config_names = [
    "MESSI+\nSAX+\nSIMD",
    "MESSI+\nSFA+EW",
    "MESSI+\nSFA+EW+\nSIMD",
    "MESSI+\nSFA+ED",
    "MESSI+\nSFA+ED+\nSIMD"]

default_path = "logs/MESSI_SFA_logs_DP"
logtypes = ["query", "index", "settings", "tree"]

def read_logs(
    log_type: str = "query",
    path: str = default_path,
    config_names=default_config_names
):
    if log_type not in logtypes:
        raise ValueError("Logtype must be one of ", logtypes)

    all_files = {}

    for d in np.sort(os.listdir(path)):
        dataset = path + "/" + d

        if not os.path.isdir(dataset):
            continue

        for c in np.sort(os.listdir(dataset)):
            config = dataset + "/" + c

            if not os.path.isdir(config):
                continue

            print("Config", config)

            queries = config + "/" + log_type + "/"

            if not os.path.isdir(queries):
                print("Missing query directory:", queries)
                continue

            key = d + " - gruenau1 - CPUs " + c

            csv_files = np.sort(
                fnmatch.filter(
                    os.listdir(queries),
                    "*.csv"
                )
            )

            # --------------------------------------------
            # Remove empty files FIRST
            # --------------------------------------------

            valid_files = []

            for q in csv_files:
                file_path = queries + "/" + q

                if os.path.getsize(file_path) == 0:
                    print("Skipping empty CSV:", file_path)
                    continue

                valid_files.append(q)

            # --------------------------------------------
            # Warn if number of configs does not match
            # --------------------------------------------

            if len(valid_files) != len(config_names):
                print(
                    f"WARNING: {queries} has "
                    f"{len(valid_files)} valid CSVs, "
                    f"but {len(config_names)} config names"
                )

                print("Files:", valid_files)
                print("Configs:", config_names)

            all_files[key] = {}

            # --------------------------------------------
            # Match only as many as we can safely assign
            # --------------------------------------------

            for name, q in zip(config_names, valid_files):
                all_files[key][name] = queries + "/" + q

    return all_files


def read_UCR_logs():
    path = "logs/UCR_SUITE_logs"
    all_files = {}
    for i, key in enumerate(np.sort(fnmatch.filter(os.listdir(path), "*.log"))):
        # print("Queries", i, key)
        config = key.split("_")[-1].split(".")[0]
        all_files[key] = {}
        all_files[key][config] = path + "/" + key
    return all_files

def read_UCR_logs_vdb():
    path = "logs/UCR_SUITE_logs_vdb"
    all_files = {}
    for i, key in enumerate(np.sort(fnmatch.filter(os.listdir(path), "*.log"))):
        # print("Queries", i, key)
        config = key.split("_")[-1].split(".")[0]
        all_files[key] = {}
        all_files[key][config] = path + "/" + key
    return all_files

def read_faiss_logs(log_type:str ="queries"):
    # path = "logs/FAISS_mini_batch_logs"
    path = "logs/FAISS_L2_logs_mini_batch"
    all_files = {}
    for i, key in enumerate(np.sort(fnmatch.filter(os.listdir(path), log_type+"*.csv"))):
        # print("Queries", i, key)
        config = key.split("_")[-1].split(".")[0]
        all_files[key] = {}
        all_files[key][config] = path + "/" + key
    return all_files


def read_faiss_logs_knn(log_type:str ="queries"):
    # path = "logs/FAISS_mini_batch_logs"
    path = "logs/FAISS_L2_logs_mini_batch_knn"
    all_files = {}
    for i, key in enumerate(np.sort(fnmatch.filter(os.listdir(path), log_type+"*.csv"))):
        # print("Queries", i, key)
        knns = key.split("_")[-2].split(".")[0]
        all_files[key] = {}
        all_files[key][knns] = path + "/" + key
    return all_files    


def read_faiss_logs_vdb(log_type:str ="queries"):
    # path = "logs/FAISS_mini_batch_logs"
    path = "logs/FAISS_L2_logs_mini_batch_vdb"
    all_files = {}
    for i, key in enumerate(np.sort(fnmatch.filter(os.listdir(path), log_type+"*.csv"))):
        # print("Queries", i, key)
        config = key.split("_")[-1].split(".")[0]
        all_files[key] = {}
        all_files[key][config] = path + "/" + key
    return all_files        