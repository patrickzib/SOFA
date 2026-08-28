import os
import numpy as np
import pandas as pd
import fnmatch
import atexit
import shutil
import tempfile
import zipfile

import sys
from pathlib import Path

sys.path.append("../")

default_config_names = [
    "MESSI+\nSAX+\nSIMD",
    "MESSI+\nSFA+EW",
    "MESSI+\nSFA+EW+\nSIMD",
    "MESSI+\nSFA+ED",
    "MESSI+\nSFA+ED+\nSIMD"]

default_path = "logs/MESSI_SFA_logs_DP"
logtypes = ["query", "index", "settings", "tree"]


# Archive members are extracted only when a caller requests them.  Returning
# regular temporary paths keeps existing notebooks compatible with
# ``pd.read_csv(file)`` and ``Path(file)``; neither needs ZIP-specific code.
_zip_log_caches = {}


def _zip_log_root_offset(members):
    """Find the optional archive prefix before dataset/config/log_type."""
    for offset in range(max((len(member.split("/")) for member in members), default=0)):
        if any(
            len(parts := member.split("/")) >= offset + 4
            and parts[offset + 2] in logtypes
            for member in members
        ):
            return offset
    return 0


def _read_logs_from_zip(log_type, archive_path, config_names):
    """Read the normal log hierarchy from a ZIP archive.

    The archive may contain either ``dataset/config/...`` directly or one
    enclosing directory (as produced by ``zip -r logs.zip logs``).
    """
    archive_path = os.path.abspath(archive_path)

    with zipfile.ZipFile(archive_path) as archive:
        members = sorted(
            info.filename for info in archive.infolist() if not info.is_dir()
        )
        root_offset = _zip_log_root_offset(members)
        entries = {}

        for member in members:
            parts = member.split("/")
            relative = parts[root_offset:]
            if (
                len(relative) != 4
                or relative[2] != log_type
                or not fnmatch.fnmatch(relative[3], "*.csv")
            ):
                continue

            dataset, config, _, filename = relative
            entries.setdefault(dataset, {}).setdefault(config, []).append(
                (filename, member, archive.getinfo(member).file_size, relative)
            )

        cache = _zip_log_caches.get(archive_path)
        if cache is None:
            cache = tempfile.TemporaryDirectory(prefix="messi-notebook-logs-")
            _zip_log_caches[archive_path] = cache

        all_files = {}
        for dataset_name in sorted(entries):
            for config_name in sorted(entries[dataset_name]):
                csv_entries = sorted(entries[dataset_name][config_name])
                display_path = f"{archive_path}!/{dataset_name}/{config_name}/{log_type}"
                print("Config", display_path)

                valid_files = []
                for filename, member, file_size, relative in csv_entries:
                    if file_size == 0:
                        print("Skipping empty CSV:", f"{archive_path}!/{member}")
                        continue

                    destination = os.path.join(cache.name, *relative)
                    if not os.path.exists(destination):
                        os.makedirs(os.path.dirname(destination), exist_ok=True)
                        with archive.open(member) as source, open(destination, "wb") as target:
                            shutil.copyfileobj(source, target)
                    valid_files.append((filename, destination))

                if len(valid_files) != len(config_names):
                    print(
                        f"WARNING: {display_path} has {len(valid_files)} valid CSVs, "
                        f"but {len(config_names)} config names"
                    )
                    print("Files:", [filename for filename, _ in valid_files])
                    print("Configs:", config_names)

                key = dataset_name + " - gruenau1 - CPUs " + config_name
                all_files[key] = {}
                for name, (_, filename) in zip(config_names, valid_files):
                    all_files[key][name] = filename

    return all_files


@atexit.register
def _cleanup_zip_log_caches():
    for cache in _zip_log_caches.values():
        cache.cleanup()

def read_logs(
    log_type: str = "query",
    path: str = default_path,
    config_names=default_config_names
):
    if log_type not in logtypes:
        raise ValueError("Logtype must be one of ", logtypes)

    # Allow an archived replacement of an existing directory without changing
    # notebook calls: read_logs("query", "logs/MESSI_SFA_logs") also accepts
    # ``logs/MESSI_SFA_logs.zip`` when the directory was removed after zipping.
    if not os.path.exists(path) and zipfile.is_zipfile(path + ".zip"):
        path += ".zip"

    if zipfile.is_zipfile(path):
        return _read_logs_from_zip(log_type, path, config_names)

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

def read_faiss_logs(log_type="queries"):
    path = Path("logs/FAISS_L2_logs_mini_batch/10")
    all_files = {}

    for file in sorted(path.glob(f"{log_type}_*.csv")):
        # queries_ASTRO_10.csv -> ASTRO, 10
        parts = file.stem.split("_")
        config = parts[-1]
        dataset = "_".join(parts[1:-1]).upper()

        all_files.setdefault(dataset, {})[config] = str(file)

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
