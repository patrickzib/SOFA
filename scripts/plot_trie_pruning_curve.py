from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Callable, Sequence

import matplotlib.pyplot as plt


STAGES = (
    "node_pruned",
    "cluster_pruned",
    "local_vq_pruned",
    "record_pruned",
)

BREAKDOWN_STAGES = (*STAGES, "exact_evaluated")
STAGE_LABELS = {
    "node_pruned": "node MBR",
    "cluster_pruned": "cluster bound",
    "local_vq_pruned": "local VQ",
    "record_pruned": "symbolic record",
    "exact_evaluated": "exact evaluation",
}
STAGE_COLORS = {
    "node_pruned": "tab:blue",
    "cluster_pruned": "tab:orange",
    "local_vq_pruned": "tab:green",
    "record_pruned": "tab:red",
    "exact_evaluated": "tab:purple",
}
STAGE_TIMING_COLUMNS = {
    "node_pruned": "node_mbr_us",
    "cluster_pruned": "cluster_bound_us",
    "local_vq_pruned": "local_vq_us",
    "record_pruned": "record_bound_us",
    "exact_evaluated": "exact_distance_us",
}
TIMING_COLUMNS = (*STAGE_TIMING_COLUMNS.values(), "other_us")

PAPER_STAGES = (
    "node_mbr",
    "group_mbr",
    "radial",
    "symbolic",
    "residual",
    "exact",
)
PAPER_LABELS = {
    "node_mbr": "1. Node MBR",
    "group_mbr": "2. Leaf MBR",
    "radial": "3. Per-Series Radial",
    "symbolic": "4. Symbolic Prefix + MBR suffix",
    "residual": "5. Symbolic Prefix + Residual",
    "exact": "6. Exact Distance",
}
PAPER_COLORS = {
    "node_mbr": "tab:blue",
    "group_mbr": "tab:orange",
    "radial": "tab:green",
    "symbolic": "tab:red",
    "residual": "tab:purple",
    "exact": "tab:brown",
}
PAPER_REQUIRED_COLUMNS = {
    "query", "total_records", "query_wall_us", "worker_search_us",
    "node_mbr_checks", "node_mbr_pruned", "node_mbr_us",
    "group_mbr_checks", "group_mbr_pruned", "group_mbr_us",
    "radial_checks", "radial_pruned", "radial_us",
    "symbolic_checks", "symbolic_pruned", "symbolic_us",
    "residual_checks", "residual_pruned", "residual_us",
    "exact_evaluated", "exact_us", "other_worker_us",
}

REQUIRED_COLUMNS = {
    "query",
    "elapsed_ms",
    "total_records",
    *BREAKDOWN_STAGES,
}


def read_paper_profile(path: Path | str) -> list[dict[str, float]]:
    """Read the final-per-query trace for the paper's five-bound cascade."""
    path = Path(path)
    rows: list[dict[str, float]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())
        missing = PAPER_REQUIRED_COLUMNS - fields
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            try:
                sample = {name: float(row[name]) for name in PAPER_REQUIRED_COLUMNS}
                sample["query"] = int(row["query"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"{path}:{line_number}: invalid profile row") from error
            if sample["total_records"] <= 0 or any(value < 0 for value in sample.values()):
                raise ValueError(f"{path}:{line_number}: counts and times must be nonnegative")
            resolved = sample["exact_evaluated"] + sum(
                sample[f"{stage}_pruned"] for stage in PAPER_STAGES[:-1]
            )
            if abs(resolved - sample["total_records"]) > 0.5:
                raise ValueError(
                    f"{path}:{line_number}: paper stages resolve {resolved:g} of "
                    f"{sample['total_records']:g} records"
                )
            rows.append(sample)
    if not rows:
        raise ValueError(f"{path}: no samples")
    return rows




def plot_paper_bound_profile(
    path: Path | str,
    *,
    aggregate: str = "mean",
    title: str = "S3-Trie paper-bound profile",
    output: Path | str | None = None,
    dpi: int = 180,
):
    """Plot the paper cascade as a clean time-weighted pruning staircase."""

    if aggregate not in {"mean", "median"}:
        raise ValueError("aggregate must be 'mean' or 'median'")

    rows = read_paper_profile(path)

    reducer = (
        statistics.mean
        if aggregate == "mean"
        else statistics.median
    )

    contributions = []
    worker_ms = []

    for stage in PAPER_STAGES:
        count_column = (
            "exact_evaluated"
            if stage == "exact"
            else f"{stage}_pruned"
        )

        contributions.append(
            reducer(
                100.0 * row[count_column] / row["total_records"]
                for row in rows
            )
        )

        worker_ms.append(
            reducer(
                row[f"{stage}_us"] / 1000.0
                for row in rows
            )
        )

    stage_worker_ms = sum(worker_ms)

    if stage_worker_ms <= 0:
        raise ValueError(
            f"{path}: total measured stage time must be positive"
        )

    # Smaller, more compact figure.
    fig, ax = plt.subplots(
        figsize=(7.0, 3.5)
    )

    previous_time = 0.0
    cumulative = 0.0

    for step, (stage, contribution, duration) in enumerate(
        zip(PAPER_STAGES, contributions, worker_ms),
        start=1,
    ):
        step_label = str(step) if stage != "exact" else "E"

        timestamp = (
            previous_time
            + 100.0 * duration / stage_worker_ms
        )

        next_cumulative = cumulative + contribution
        color = PAPER_COLORS[stage]

        ax.plot(
            [previous_time, timestamp, timestamp],
            [cumulative, cumulative, next_cumulative],
            color=color,
            linewidth=3.2,
            solid_capstyle="round",
            label=PAPER_LABELS[stage],
            zorder=2,
        )

        ax.scatter(
            timestamp,
            next_cumulative,
            color=color,
            s=70,
            edgecolor="white",
            linewidth=1.2,
            zorder=4,
        )

        ax.annotate(
            step_label,
            xy=(timestamp, next_cumulative),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11.5,
            fontweight="bold",
            color=color,
            bbox={
                "boxstyle": "circle,pad=0.15",
                "facecolor": "white",
                "edgecolor": color,
                "linewidth": 1.2,
                "alpha": 0.95,
            },
            zorder=5,
        )

        if contribution >= 0.01:
            if contribution < 0.1:
                contribution_text = f"+{contribution:.3f}%"
            else:
                contribution_text = f"+{contribution:.1f}%"

            midpoint_y = (cumulative + next_cumulative) / 2.0

            if timestamp > 90:
                xytext = (-10, -2)
                ha = "right"
            elif next_cumulative > 96:
                xytext = (8, 6)
                ha = "left"
            else:
                xytext = (8, 0)
                ha = "left"

            ax.annotate(
                contribution_text,
                xy=(timestamp, midpoint_y),
                xytext=xytext,
                textcoords="offset points",
                ha=ha,
                va="center",
                fontsize=11.5,
                fontweight="semibold",
                color=color,
                zorder=5,
            )

        previous_time = timestamp
        cumulative = next_cumulative

    # Title and labels.
    ax.set_title(
        title,
        fontsize=18,
        fontweight="semibold",
        pad=14,
    )

    ax.set_xlabel(
        f"{aggregate.capitalize()} time fraction",
        fontsize=15,
        labelpad=10,
    )

    ax.set_ylabel(
        f"{aggregate.capitalize()} series pruned",
        fontsize=15,
        labelpad=10,
    )

    ax.set_xlim(0, 103)
    ax.set_ylim(0, 105)

    ax.set_xticks(
        [0, 20, 40, 60, 80, 100],
        ["0%", "20%", "40%", "60%", "80%", "100%"],
    )

    ax.set_yticks(
        [0, 20, 40, 60, 80, 100],
        ["0%", "20%", "40%", "60%", "80%", "100%"],
    )

    ax.tick_params(
        axis="both",
        labelsize=13,
        length=5,
        width=1.0,
        pad=5,
    )

    ax.set_axisbelow(True)

    ax.grid(
        axis="y",
        alpha=0.16,
        linewidth=0.9,
    )

    ax.grid(
        axis="x",
        alpha=0.05,
        linewidth=0.7,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ------------------------------------------------------------
    # Compact legend below the axes
    # ------------------------------------------------------------

    handles, labels = ax.get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=3,
        fontsize=10,
        frameon=False,
        handlelength=2.4,
        handletextpad=0.6,
        columnspacing=1.6,
        labelspacing=0.6,
    )

    # Reserve just enough space for the legend.
    fig.subplots_adjust(
        left=0.12,
        right=0.98,
        top=0.88,
        bottom=0.23,
    )

    if output is not None:
        output = Path(output)
        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig.savefig(
            output,
            dpi=dpi,
            bbox_inches="tight",
        )

    return fig, ax




def newest_dataset_curve(dataset_directory: Path) -> Path | None:
    curves = list(
        dataset_directory.glob(
            "**/trie_pruning_curve/MESSI_TRIE_PRUNING_CURVE_*.csv"
        )
    )
    return max(curves, key=lambda path: path.stat().st_mtime) if curves else None
