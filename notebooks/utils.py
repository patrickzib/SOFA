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

def read_logs(log_type: str = "query",
              path: str = default_path,
              config_names = default_config_names):
    if log_type not in logtypes:
        raise ValueError("Logtype must be one of ", logtypes)

    all_files = {}
    for d in np.sort(os.listdir(path)):
        dataset = path + "/" + d
        if (os.path.isdir(dataset)):
            # print("Dataset", dataset)
            for c in np.sort(os.listdir(dataset)):
                config = dataset + "/" + c
                if os.path.isdir(config):
                    print("Config", config)
                    queries = config + "/" + log_type + "/"
                    key = d + " - gruenau1 - CPUs " + c
                    all_files[key] = {}
                    for i, q in enumerate(
                            np.sort(fnmatch.filter(os.listdir(queries), "*.csv"))):
                        # print("Queries", i, q, key, q)
                        all_files[key][config_names[i]] = queries + "/" + q

            print("-----------------")

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