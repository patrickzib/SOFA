from __future__ import annotations

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

REQUIRED_COLUMNS = {
    "query",
    "elapsed_ms",
    "total_records",
    *BREAKDOWN_STAGES,
}


def read_curve(
    path: Path | str,
) -> dict[int, list[dict[str, float]]]:
    """Read and validate a pruning curve.

    Samples must already be ordered within each query.  Equal timestamps are
    allowed because several cumulative updates can happen between clock
    ticks; only the last sample at an equal timestamp is retained.
    """
    path = Path(path)

    by_query: dict[int, list[dict[str, float]]] = defaultdict(list)

    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())

        missing = REQUIRED_COLUMNS - fields

        if missing:
            raise ValueError(
                f"{path}: missing columns: {', '.join(sorted(missing))}"
            )

        present_timing = set(TIMING_COLUMNS) & fields
        if present_timing and present_timing != set(TIMING_COLUMNS):
            missing_timing = set(TIMING_COLUMNS) - fields
            raise ValueError(
                f"{path}: incomplete stage-timing columns: "
                f"{', '.join(sorted(missing_timing))}"
            )

        for line_number, row in enumerate(reader, start=2):
            try:
                sample = {
                    name: float(row[name])
                    for name in (
                        "elapsed_ms",
                        "total_records",
                        *BREAKDOWN_STAGES,
                        *(TIMING_COLUMNS if present_timing else ()),
                    )
                }

                query = int(row["query"])

            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid curve row"
                ) from error

            if sample["total_records"] <= 0:
                raise ValueError(
                    f"{path}:{line_number}: total_records must be positive"
                )

            by_query[query].append(sample)

    if not by_query:
        raise ValueError(f"{path}: no samples")

    cleaned: dict[int, list[dict[str, float]]] = {}

    for query, samples in by_query.items():
        deduplicated: list[dict[str, float]] = []
        total_records = samples[0]["total_records"]
        previous_elapsed = -1.0
        previous_counts = {
            stage: 0.0
            for stage in BREAKDOWN_STAGES
        }

        for sample in samples:
            elapsed = sample["elapsed_ms"]

            if elapsed < previous_elapsed:
                raise ValueError(
                    f"{path}: query {query}: elapsed_ms decreases "
                    f"from {previous_elapsed:g} to {elapsed:g}"
                )

            if sample["total_records"] != total_records:
                raise ValueError(
                    f"{path}: query {query}: total_records changes "
                    "within the query"
                )

            for stage in BREAKDOWN_STAGES:
                count = sample[stage]
                if count < previous_counts[stage]:
                    raise ValueError(
                        f"{path}: query {query}: {stage} decreases "
                        f"from {previous_counts[stage]:g} to {count:g}"
                    )
                previous_counts[stage] = count

            resolved = sum(sample[stage] for stage in BREAKDOWN_STAGES)
            if resolved > total_records + 0.5:
                raise ValueError(
                    f"{path}: query {query}: cumulative stage counters "
                    f"exceed total_records ({resolved:g} > {total_records:g})"
                )

            if (
                deduplicated
                and elapsed == deduplicated[-1]["elapsed_ms"]
            ):
                deduplicated[-1] = sample
            else:
                deduplicated.append(sample)

            previous_elapsed = elapsed

        final_resolved = sum(
            deduplicated[-1][stage]
            for stage in BREAKDOWN_STAGES
        )
        if abs(final_resolved - total_records) > 0.5:
            raise ValueError(
                f"{path}: query {query}: final stage counters account for "
                f"{final_resolved:g} of {total_records:g} records; the trace "
                "cannot produce a certified 100% staircase"
            )

        cleaned[query] = deduplicated

    return cleaned


def cumulative_percent(
    sample: dict[str, float],
    stage: str | None = None,
) -> float:
    """
    Calculate cumulative pruning percentage.

    If stage is None:
        sum all pruning stages.

    If stage is given:
        sum all stages up to and including that stage.
    """
    if stage is None:
        value = sum(
            sample[name]
            for name in STAGES
        )

    else:
        stage_index = BREAKDOWN_STAGES.index(stage)

        value = sum(
            sample[name]
            for name in BREAKDOWN_STAGES[: stage_index + 1]
        )

    return 100.0 * value / sample["total_records"]


def interpolate(
    samples: list[dict[str, float]],
    elapsed_ms: float,
    value: Callable[[dict[str, float]], float],
) -> float:
    """Step interpolation: pruning happens at discrete search events."""
    current = value(samples[0])

    for sample in samples:
        if sample["elapsed_ms"] > elapsed_ms:
            break

        current = value(sample)

    return current


def aggregate_curve(
    curves: dict[int, list[dict[str, float]]],
    aggregate: str,
    value: Callable[[dict[str, float]], float],
) -> tuple[list[float], list[float]]:
    """Aggregate query curves on a shared elapsed-time grid."""

    max_elapsed = max(
        samples[-1]["elapsed_ms"]
        for samples in curves.values()
    )

    points = max(
        2,
        min(
            500,
            int(max_elapsed) + 1,
        ),
    )

    x_values = [
        max_elapsed * index / (points - 1)
        for index in range(points)
    ]

    reducer = (
        statistics.median
        if aggregate == "median"
        else statistics.mean
    )

    observations: list[list[float]] = [[] for _ in x_values]
    for samples in curves.values():
        sample_index = 0
        current = value(samples[0])
        for grid_index, elapsed_ms in enumerate(x_values):
            while (
                sample_index + 1 < len(samples)
                and samples[sample_index + 1]["elapsed_ms"] <= elapsed_ms
            ):
                sample_index += 1
                current = value(samples[sample_index])
            observations[grid_index].append(current)

    y_values = [reducer(values) for values in observations]

    return x_values, y_values


def selected_curve(
    curves: dict[int, list[dict[str, float]]],
    query: int,
    value: Callable[[dict[str, float]], float],
) -> tuple[list[float], list[float]]:
    """Return the curve for one particular query."""

    if query not in curves:
        available = ", ".join(
            str(number)
            for number in sorted(curves)
        )

        raise ValueError(
            f"query {query} is absent "
            f"(available: {available})"
        )

    samples = curves[query]

    x_values = [
        sample["elapsed_ms"]
        for sample in samples
    ]

    y_values = [
        value(sample)
        for sample in samples
    ]

    return x_values, y_values


def annotate_numeric_column(
    ax,
    samples: list[dict[str, float]],
    value: Callable[[dict[str, float]], float],
    column: str,
    *,
    fontsize: int = 8,
    offset: tuple[int, int] = (5, 5),
) -> None:
    """
    Annotate plot points using values from an existing numeric CSV column.

    Example columns:
        node_pruned
        cluster_pruned
        local_vq_pruned
        record_pruned
        total_records
        elapsed_ms
    """

    if column not in samples[0]:
        available = ", ".join(
            sorted(samples[0].keys())
        )

        raise ValueError(
            f"Unknown annotation column {column!r}. "
            f"Available columns: {available}"
        )

    for sample in samples:
        x = sample["elapsed_ms"]
        y = value(sample)

        label_value = sample[column]

        ax.annotate(
            f"{label_value:g}",
            xy=(x, y),
            xytext=offset,
            textcoords="offset points",
            fontsize=fontsize,
        )


def stage_timing_breakdown(
    curves: dict[int, list[dict[str, float]]],
    aggregate: str,
    query: int | None,
) -> tuple[list[tuple[str, float, float]], float, float]:
    """Return fixed-order stage durations and final contributions."""
    if query is None:
        selected = curves
    else:
        if query not in curves:
            available = ", ".join(str(number) for number in sorted(curves))
            raise ValueError(
                f"query {query} is absent (available: {available})"
            )
        selected = {query: curves[query]}

    reducer = statistics.median if aggregate == "median" else statistics.mean
    stage_summaries: list[tuple[str, float, float]] = []

    missing_timing = [
        query_number
        for query_number, samples in selected.items()
        if any(column not in samples[-1] for column in TIMING_COLUMNS)
    ]
    if missing_timing:
        raise ValueError(
            "This breakdown requires a new per-query timing log. Rerun with "
            "--trie-pruning-curve after rebuilding the index executable."
        )

    for stage in BREAKDOWN_STAGES:
        contributions: list[float] = []
        durations: list[float] = []

        for samples in selected.values():
            total_records = samples[0]["total_records"]
            maximum = max(sample[stage] for sample in samples)
            contributions.append(100.0 * maximum / total_records)
            durations.append(samples[-1][STAGE_TIMING_COLUMNS[stage]] / 1000.0)

        contribution = reducer(contributions)
        duration = reducer(durations)
        if contribution > 1e-12 or duration > 1e-12:
            stage_summaries.append(
                (stage, duration, contribution)
            )

    other_time = reducer(
        [samples[-1]["other_us"] / 1000.0 for samples in selected.values()]
    )
    full_query_time = reducer(
        [samples[-1]["elapsed_ms"] for samples in selected.values()]
    )
    return stage_summaries, other_time, full_query_time


def plot_stage_staircase(
    ax,
    stage_summaries: Sequence[tuple[str, float, float]],
    other_time: float,
    full_query_time: float,
    *,
    markers: bool,
    annotated_stages: Sequence[str],
) -> None:
    """Plot a fixed-order staircase from measured stage durations."""
    if not stage_summaries:
        raise ValueError("The curve contains no stage contributions.")

    previous_time = 0.0
    cumulative = 0.0
    for stage, duration, contribution in stage_summaries:
        timestamp = previous_time + duration
        next_cumulative = cumulative + contribution
        ax.plot(
            [previous_time, timestamp, timestamp],
            [cumulative, cumulative, next_cumulative],
            color=STAGE_COLORS[stage],
            linewidth=2.5,
            solid_capstyle="butt",
            label=STAGE_LABELS[stage],
        )
        if stage in annotated_stages:
            contribution_text = (
                f"{contribution:.4f}%"
                if contribution < 0.01
                else f"{contribution:.2f}%"
            )
            near_endpoint = timestamp >= 0.82 * full_query_time
            offset_x = -8 if near_endpoint else 8
            offset_y = -12 if next_cumulative > 96.0 else 0
            ax.annotate(
                f"{STAGE_LABELS[stage]}\n+{contribution_text}",
                xy=(timestamp, 0.5 * (cumulative + next_cumulative)),
                xytext=(offset_x, offset_y),
                textcoords="offset points",
                ha="right" if near_endpoint else "left",
                va="top" if offset_y < 0 else "center",
                color=STAGE_COLORS[stage],
                fontsize=8,
            )
        if markers:
            ax.scatter(
                [timestamp],
                [next_cumulative],
                color=STAGE_COLORS[stage],
                s=25,
                zorder=3,
            )
        previous_time = timestamp
        cumulative = next_cumulative

    endpoint = max(previous_time + other_time, full_query_time)
    ax.plot(
        [previous_time, endpoint],
        [cumulative, cumulative],
        color=STAGE_COLORS[stage_summaries[-1][0]],
        linewidth=2.5,
    )
    ax.annotate(
        f"full query: {cumulative:.2f}%",
        xy=(endpoint, cumulative),
        xytext=(8, -10),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
    )
    ax.set_xlim(left=0, right=max(1.0, endpoint * 1.30))


def plot_trie_pruning_curve(
    curves: Sequence[Path | str],
    *,
    labels: Sequence[str] | None = None,
    query: int | None = None,
    aggregate: str = "median",
    breakdown: bool = False,
    title: str = "Trie pruning curve",
    output: Path | str | None = None,
    figsize: tuple[float, float] = (8.0, 4.8),
    dpi: int = 180,
    markers: bool = True,
    annotation_column: str | None = None,
    annotation_fontsize: int = 8,
):
    """
    Plot one or more trie-pruning curves.

    Parameters
    ----------
    curves:
        CSV files written by --trie-pruning-curve.

    labels:
        Legend labels, one per CSV.

    query:
        Plot one particular query.

        If None, all queries are aggregated.

    aggregate:
        "median" or "mean".

    breakdown:
        Show one fixed-order cumulative staircase from the final per-query
        stage counters and timings. Each stage's horizontal segment is its
        measured duration; its vertical jump is the final contribution.

        Requires exactly one CSV.

    title:
        Plot title.

    output:
        Optional path for PNG/PDF/SVG output.

    figsize:
        Matplotlib figure size.

    dpi:
        DPI when saving.

    markers:
        Show circles at individual plotted points.

    annotation_column:
        In breakdown mode, optionally annotate only this stage instead of all
        stages. In a normal comparison plot, annotate this numeric
        column at every sampled point.

        Examples:
            "node_pruned"
            "cluster_pruned"
            "local_vq_pruned"
            "record_pruned"

        Normal-plot annotations are only available when query is set.

    annotation_fontsize:
        Font size for point annotations.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """

    curves = [
        Path(path)
        for path in curves
    ]

    if not curves:
        raise ValueError(
            "At least one curve must be supplied."
        )

    if aggregate not in {"median", "mean"}:
        raise ValueError(
            f"aggregate must be 'median' or 'mean', "
            f"got {aggregate!r}"
        )

    if labels is not None and len(labels) != len(curves):
        raise ValueError(
            f"Got {len(labels)} labels for "
            f"{len(curves)} CSV files. "
            "There must be exactly one label per CSV."
        )

    if breakdown and len(curves) != 1:
        raise ValueError(
            "breakdown=True requires exactly one CSV file."
        )

    if annotation_column is not None and query is None and not breakdown:
        raise ValueError(
            "annotation_column requires query=<query number>. "
            "Annotations cannot be mapped uniquely onto an "
            "aggregated curve."
        )

    labels = (
        list(labels)
        if labels is not None
        else [
            path.stem
            for path in curves
        ]
    )

    loaded = [
        read_curve(path)
        for path in curves
    ]

    fig, ax = plt.subplots(
        figsize=figsize,
        constrained_layout=True,
    )

    if breakdown:
        curve = loaded[0]
        stages_to_annotate = (
            (annotation_column,)
            if annotation_column is not None
            else BREAKDOWN_STAGES
        )
        unknown = set(stages_to_annotate) - set(BREAKDOWN_STAGES)
        if unknown:
            raise ValueError(
                "breakdown annotation_column must be a stage column: "
                + ", ".join(sorted(BREAKDOWN_STAGES))
            )
        stage_summaries, other_time, full_query_time = stage_timing_breakdown(
            curve,
            aggregate,
            query,
        )
        plot_stage_staircase(
            ax,
            stage_summaries,
            other_time,
            full_query_time,
            markers=markers,
            annotated_stages=stages_to_annotate,
        )

    else:
        for curve, label in zip(
            loaded,
            labels,
        ):

            if query is None:
                x_values, y_values = aggregate_curve(
                    curve,
                    aggregate,
                    cumulative_percent,
                )

            else:
                x_values, y_values = selected_curve(
                    curve,
                    query,
                    cumulative_percent,
                )

            ax.plot(
                x_values,
                y_values,
                linewidth=2,
                marker="o" if markers else None,
                label=label,
            )

            if (
                annotation_column is not None
                and query is not None
            ):
                annotate_numeric_column(
                    ax,
                    curve[query],
                    cumulative_percent,
                    annotation_column,
                    fontsize=annotation_fontsize,
                )

    ax.set_title(title)

    ax.set_xlabel(
        "Elapsed query time (ms)"
    )

    ax.set_ylabel(
        (
            f"{aggregate.capitalize()} cumulative records resolved (%)"
            if query is None
            else "Cumulative records resolved (%)"
        )
        if breakdown
        else "Cumulative database records pruned (%)"
    )

    ax.set_ylim(
        0,
        100,
    )

    ax.set_xlim(
        left=0,
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        loc="upper left" if breakdown else "lower right",
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
        )

    return fig, ax


def newest_curve(directory: Path) -> Path:
    curves = list(directory.glob("MESSI_TRIE_PRUNING_CURVE_*.csv"))
    if not curves:
        raise FileNotFoundError(f"No pruning-curve CSV found in {directory}")
    return max(curves, key=lambda path: path.stat().st_mtime)


def newest_dataset_curve(dataset_directory: Path) -> Path | None:
    curves = list(
        dataset_directory.glob(
            "**/trie_pruning_curve/MESSI_TRIE_PRUNING_CURVE_*.csv"
        )
    )
    return max(curves, key=lambda path: path.stat().st_mtime) if curves else None


def main() -> None:
    repository = Path(__file__).resolve().parents[1]
    local_curve_root = Path("trie_pruning_curves")
    curve_root = (
        local_curve_root
        if local_curve_root.is_dir()
        else repository / "notebooks" / "trie_pruning_curves"
    )
    if not curve_root.is_dir():
        raise FileNotFoundError(f"No trie-pruning-curves directory found at {curve_root}")

    output_directory = curve_root / "plots"
    for dataset_directory in sorted(curve_root.iterdir()):
        if not dataset_directory.is_dir() or dataset_directory.name == "plots":
            continue
        curve_file = newest_dataset_curve(dataset_directory)
        if curve_file is None:
            continue
        dataset = dataset_directory.name
        output = output_directory / f"{dataset}_pruning_staircase.png"
        fig, _ = plot_trie_pruning_curve(
            [curve_file],
            breakdown=True,
            aggregate="mean",
            markers=True,
            title=f"{dataset}: mean pruning progression",
            output=output,
        )
        plt.close(fig)
        print(output)


if __name__ == "__main__":
    main()
