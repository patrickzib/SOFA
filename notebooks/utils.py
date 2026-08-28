import os
import sys
import atexit
import shutil
import tempfile
import zipfile
from pathlib import Path

import fnmatch
import numpy as np
import pandas as pd

sys.path.append("../")


# ============================================================
# DEFAULTS
# ============================================================

default_config_names = [
    "MESSI+\nSAX+\nSIMD",
    "MESSI+\nSFA+EW",
    "MESSI+\nSFA+EW+\nSIMD",
    "MESSI+\nSFA+ED",
    "MESSI+\nSFA+ED+\nSIMD",
]

default_path = "logs/MESSI_SFA_logs_DP"
logtypes = ["query", "index", "settings", "tree"]

_zip_log_caches = {}


# ============================================================
# GENERIC HELPERS
# ============================================================

def read_flat_logs(path, pattern, key_parser):
    path = Path(path)
    all_files = {}

    for file in sorted(path.glob(pattern)):
        config = key_parser(file)
        all_files.setdefault(file.name, {})[config] = str(file)

    return all_files


def align_columns(df, reference_columns):
    df = df.copy()

    for col in reference_columns:
        if col not in df.columns:
            df[col] = pd.NA

    return df[reference_columns]


# ============================================================
# ZIP LOG SUPPORT
# ============================================================

def _zip_log_root_offset(members):
    max_depth = max((len(member.split("/")) for member in members), default=0)

    for offset in range(max_depth):
        if any(
            len(parts := member.split("/")) >= offset + 4
            and parts[offset + 2] in logtypes
            for member in members
        ):
            return offset

    return 0


def _read_logs_from_zip(log_type, archive_path, config_names):
    archive_path = os.path.abspath(archive_path)

    with zipfile.ZipFile(archive_path) as archive:
        members = sorted(
            info.filename
            for info in archive.infolist()
            if not info.is_dir()
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
                (
                    filename,
                    member,
                    archive.getinfo(member).file_size,
                    relative,
                )
            )

        cache = _zip_log_caches.get(archive_path)

        if cache is None:
            cache = tempfile.TemporaryDirectory(prefix="messi-notebook-logs-")
            _zip_log_caches[archive_path] = cache

        all_files = {}

        for dataset_name in sorted(entries):
            for config_name in sorted(entries[dataset_name]):
                csv_entries = sorted(entries[dataset_name][config_name])

                display_path = (
                    f"{archive_path}!/"
                    f"{dataset_name}/"
                    f"{config_name}/"
                    f"{log_type}"
                )

                valid_files = []

                for filename, member, file_size, relative in csv_entries:
                    if file_size == 0:
                        print("Skipping empty CSV:", f"{archive_path}!/{member}")
                        continue

                    destination = os.path.join(cache.name, *relative)

                    if not os.path.exists(destination):
                        os.makedirs(
                            os.path.dirname(destination),
                            exist_ok=True,
                        )

                        with archive.open(member) as source, open(
                            destination, "wb"
                        ) as target:
                            shutil.copyfileobj(source, target)

                    valid_files.append((filename, destination))

                if len(valid_files) != len(config_names):
                    print(
                        f"WARNING: {display_path} has "
                        f"{len(valid_files)} valid CSVs, "
                        f"but {len(config_names)} config names"
                    )
                    print(
                        "Files:",
                        [filename for filename, _ in valid_files],
                    )
                    # print("Configs:", config_names)

                key = (
                    dataset_name
                    + " - gruenau1 - CPUs "
                    + config_name
                )

                all_files[key] = {}

                for name, (_, filename) in zip(config_names, valid_files):
                    all_files[key][name] = filename

    return all_files


@atexit.register
def _cleanup_zip_log_caches():
    for cache in _zip_log_caches.values():
        cache.cleanup()


# ============================================================
# MESSI / TRIE / iSAX LOG READERS
# ============================================================

def read_logs(
    log_type="query",
    path=default_path,
    config_names=default_config_names,
):
    if log_type not in logtypes:
        raise ValueError(f"log_type must be one of {logtypes}")

    if not os.path.exists(path) and zipfile.is_zipfile(path + ".zip"):
        path += ".zip"

    if zipfile.is_zipfile(path):
        return _read_logs_from_zip(
            log_type=log_type,
            archive_path=path,
            config_names=config_names,
        )

    all_files = {}

    for dataset_dir in sorted(Path(path).iterdir()):
        if not dataset_dir.is_dir():
            continue

        for config_dir in sorted(dataset_dir.iterdir()):
            if not config_dir.is_dir():
                continue

            # print("Config", config_dir)

            log_dir = config_dir / log_type

            if not log_dir.is_dir():
                print("Missing log directory:", log_dir)
                continue

            valid_files = []

            for file in sorted(log_dir.glob("*.csv")):
                if file.stat().st_size == 0:
                    print("Skipping empty CSV:", file)
                    continue

                valid_files.append(file)

            if len(valid_files) != len(config_names):
                print(
                    f"WARNING: {log_dir} has "
                    f"{len(valid_files)} valid CSVs, "
                    f"but {len(config_names)} config names"
                )
                print("Files:", [file.name for file in valid_files])
                print("Configs:", config_names)

            key = (
                dataset_dir.name
                + " - gruenau1 - CPUs "
                + config_dir.name
            )

            all_files[key] = {}

            for name, file in zip(config_names, valid_files):
                all_files[key][name] = str(file)

    return all_files


# ============================================================
# UCR READERS
# ============================================================

def read_UCR_logs():
    return read_flat_logs(
        path="logs/UCR_SUITE_logs",
        pattern="*.log",
        key_parser=lambda f: f.stem.split("_")[-1],
    )


def read_UCR_logs_vdb():
    return read_flat_logs(
        path="logs/UCR_SUITE_logs_vdb",
        pattern="*.log",
        key_parser=lambda f: f.stem.split("_")[-1],
    )


# ============================================================
# FAISS READERS
# ============================================================

def read_faiss_logs(log_type="queries"):
    path = Path("logs/FAISS_L2_logs_mini_batch/10")
    all_files = {}

    for file in sorted(path.glob(f"{log_type}_*.csv")):
        parts = file.stem.split("_")

        config = parts[-1]
        dataset = "_".join(parts[1:-1])

        all_files.setdefault(dataset, {})[config] = str(file)

    return all_files


def read_faiss_logs_knn(log_type="queries"):
    return read_flat_logs(
        path="logs/FAISS_L2_logs_mini_batch_knn",
        pattern=f"{log_type}*.csv",
        key_parser=lambda f: f.stem.split("_")[-2],
    )


def read_faiss_logs_vdb(log_type="queries"):
    return read_flat_logs(
        path="logs/FAISS_L2_logs_mini_batch_vdb",
        pattern=f"{log_type}*.csv",
        key_parser=lambda f: f.stem.split("_")[-1],
    )


# ============================================================
# QUERY LOADERS
# ============================================================

def load_logs(path, config_names, layout):
    all_files = read_logs(
        log_type="query",
        path=path,
        config_names=config_names,
    )

    dfs = []

    for ds_name, files in all_files.items():
        dataset = ds_name.split(" - ")[0]

        for method_name, file in files.items():
            df = pd.read_csv(file).iloc[:-1].copy()

            df["querying time"] /= 1_000_000
            df["method"] = method_name
            df["dataset"] = dataset
            df["layout"] = layout

            dfs.append(
                df[
                    [
                        "method",
                        "querying time",
                        "dataset",
                        "layout",
                    ]
                ]
            )

    columns = [
        "method",
        "querying time",
        "dataset",
        "layout",
    ]

    return (
        pd.concat(dfs, ignore_index=True)
        if dfs
        else pd.DataFrame(columns=columns)
    )


def load_dids_logs(path, layout="DIDS"):
    dfs = []

    for file in sorted(Path(path).rglob("queries.csv")):
        df = pd.read_csv(file).copy()

        df["querying time"] = df["query_time_ms"] / 1_000
        df["method"] = "DIDS"
        df["dids_mode"] = df["mode"]

        # DIDS_logs / dataset / threads-64 / queries.csv
        df["dataset"] = file.parent.parent.name
        df["threads"] = file.parent.name
        df["layout"] = layout

        dfs.append(
            df[
                [
                    "method",
                    "querying time",
                    "dataset",
                    "layout",
                    "dids_mode",
                    "threads",
                ]
            ]
        )

    columns = [
        "method",
        "querying time",
        "dataset",
        "layout",
        "dids_mode",
        "threads",
    ]

    return (
        pd.concat(dfs, ignore_index=True)
        if dfs
        else pd.DataFrame(columns=columns)
    )


def load_faiss_query_logs():
    dfs = []

    for dataset, files in read_faiss_logs().items():
        for cores, file in files.items():
            df = pd.read_csv(file).copy()

            df["querying time"] = df["time in ms"] / 1_000

            df["method"] = "FAISS"
            df["base_method"] = "FAISS"
            df["layout"] = "FAISS"
            df["Method"] = "FAISS"

            df["split"] = str(cores)
            df["Cores"] = str(cores)
            df["dataset"] = dataset

            dfs.append(
                df[
                    [
                        "method",
                        "base_method",
                        "layout",
                        "split",
                        "Cores",
                        "dataset",
                        "Method",
                        "querying time",
                    ]
                ]
            )

    columns = [
        "method",
        "base_method",
        "layout",
        "split",
        "Cores",
        "dataset",
        "Method",
        "querying time",
    ]

    return (
        pd.concat(dfs, ignore_index=True)
        if dfs
        else pd.DataFrame(columns=columns)
    )


# ============================================================
# METHOD METADATA
# ============================================================

def add_method_metadata(df, split):
    df = df.copy()

    df["split"] = split
    df["base_method"] = (
        df["method"]
        .str.replace(
            r"\+(ED|EW)$",
            "",
            regex=True,
        )
    )

    return df


# ============================================================
# CONFIG SELECTION
# ============================================================

def select_best_query_configs(query_times):
    group_cols = [
        "dataset",
        "layout",
        "split",
        "base_method",
        "method",
    ]

    mean_times = (
        query_times
        .groupby(
            group_cols,
            as_index=False,
            dropna=False,
        )["querying time"]
        .mean()
    )

    best_configs = mean_times.loc[
        mean_times.groupby(
            [
                "dataset",
                "layout",
                "split",
                "base_method",
            ],
            dropna=False,
        )["querying time"]
        .idxmin()
    ].reset_index(drop=True)

    best_overall = best_configs.loc[
        best_configs.groupby(
            [
                "dataset",
                "layout",
                "base_method",
            ],
            dropna=False,
        )["querying time"]
        .idxmin()
    ].reset_index(drop=True)

    result = query_times.merge(
        best_overall[group_cols],
        on=group_cols,
        how="inner",
    )

    result["selected_config"] = result["method"]
    result["Method"] = result["layout"] + " " + result["base_method"]
    result.loc[result["layout"].eq("DIDS"), "Method"] = "DIDS"

    return result


def select_best_faiss_config(df):
    mean_times = (
        df.groupby(
            ["dataset", "Cores"],
            as_index=False,
        )["querying time"]
        .mean()
    )

    best = mean_times.loc[
        mean_times.groupby("dataset")["querying time"]
        .idxmin()
    ].reset_index(drop=True)

    result = df.merge(
        best[["dataset", "Cores"]],
        on=["dataset", "Cores"],
        how="inner",
    )

    result["selected_config"] = (
        "FAISS "
        + result["Cores"].astype(str)
        + " cores"
    )

    return result