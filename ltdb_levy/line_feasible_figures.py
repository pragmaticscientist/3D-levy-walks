"""Figures for the condition-specific line-feasible target sweep."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .nk_experiment import FINITE_PLOT_ORDER, UNBOUNDED_CLASSES
from .publication_figures import (
    CLASS_STYLES,
    COMPACT_CONDITION_LABELS,
    DETAILED_SEARCH_HEIGHT_INCHES,
    FINITE_ROW_LABELS,
    SELECTED_CONDITIONS,
    SHAPE_ORDER,
    SHAPE_STYLES,
    _panel_label,
    _save_figure,
)


def condition_c_label(
    biological_label: str,
    crossover_um: float,
) -> str:
    """Return the only numeric condition annotation used in the figures."""
    return "{}\n{}".format(
        biological_label,
        r"$c={:.2f}\,\mu\mathrm{{m}}$".format(float(crossover_um)),
    )


def normalized_area_ticks(maximum_area: float) -> Tuple[float, ...]:
    """Sparse powers-of-two tick locations for the normalized log area axis."""
    candidates = (4.0, 8.0, 16.0, 32.0, 64.0, 128.0)
    return tuple(value for value in candidates if value <= maximum_area + 1e-12)


def _set_normalized_area_axis(axis: object, areas: Sequence[float]) -> None:
    values = np.asarray(sorted(set(float(value) for value in areas)), dtype=float)
    if len(values) == 0:
        raise ValueError("Cannot configure an area axis without data")
    axis.set_xscale("log", base=2)
    lower = float(values[0]) / 1.12
    upper = float(values[-1]) * 1.12
    axis.set_xlim(lower, upper)
    ticks = normalized_area_ticks(float(values[-1]))
    axis.set_xticks(ticks)
    axis.set_xticklabels(["{:g}".format(value) for value in ticks])


def _condition_maps(fits: pd.DataFrame) -> Tuple[Dict[str, float], Dict[str, str]]:
    crossovers = fits.set_index("condition")["crossover"].astype(float).to_dict()
    labels = fits.set_index("condition")["condition_label"].astype(str).to_dict()
    return crossovers, labels


def _finite_upper_limit(frame: pd.DataFrame) -> float:
    values = 100.0 * frame["detection_probability_ci_high"].to_numpy(dtype=float)
    return max(1e-3, 1.13 * float(np.nanmax(values)))


def _positive_log_errors(
    means: np.ndarray,
    standard_errors: np.ndarray,
) -> np.ndarray:
    """Return positive asymmetric normal intervals suitable for a log axis."""
    lower_endpoint = np.maximum(
        means - 1.96 * standard_errors,
        np.maximum(np.finfo(float).tiny, means * 1e-3),
    )
    upper_endpoint = means + 1.96 * standard_errors
    return np.vstack((means - lower_endpoint, upper_endpoint - means))


def _plot_all_condition_unbounded(
    fits: pd.DataFrame,
    unbounded: pd.DataFrame,
    figure_directory: Path,
) -> List[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    crossovers, labels = _condition_maps(fits)
    row_order = (
        "Natural Killer Cells",
        "Neutrophils",
        "T Cells",
        "B Cells",
    )
    row_cells = [value for value in row_order if value in set(unbounded["cell_type"])]
    conditions = sorted(
        unbounded["condition"].astype(str).unique(),
        key=lambda value: labels[value],
    )
    palette = plt.get_cmap("tab20")
    condition_colors = {
        condition: palette(index / max(1, len(conditions) - 1))
        for index, condition in enumerate(conditions)
    }
    class_styles: Mapping[str, Tuple[str, str, str]] = {
        "vector_uniform_rotated": (
            "-",
            "o",
            "Empirical run pool",
        ),
        "fitted_mu": (
            "--",
            "s",
            r"Fitted $\mu$ model",
        ),
    }

    fig, axes = plt.subplots(
        len(row_cells),
        len(SHAPE_ORDER),
        figsize=(12.0, 2.65 * len(row_cells)),
        sharex=True,
        sharey="row",
        squeeze=False,
    )
    for row_index, cell_type in enumerate(row_cells):
        for column_index, shape in enumerate(SHAPE_ORDER):
            axis = axes[row_index, column_index]
            selected = unbounded[
                (unbounded["cell_type"] == cell_type)
                & (unbounded["shape"] == shape)
            ]
            for condition in sorted(
                selected["condition"].astype(str).unique(),
                key=lambda value: labels[value],
            ):
                for class_name in UNBOUNDED_CLASSES:
                    group = selected[
                        (selected["condition"] == condition)
                        & (selected["trajectory_class"] == class_name)
                    ].sort_values("projected_area")
                    if group.empty:
                        continue
                    linestyle, marker, _ = class_styles[class_name]
                    axis.plot(
                        group["projected_area"],
                        group["mean_detection_distance_model_units"],
                        color=condition_colors[condition],
                        linestyle=linestyle,
                        marker=marker,
                        markersize=3.0,
                        linewidth=0.95,
                    )
            _set_normalized_area_axis(
                axis,
                unbounded["projected_area"].astype(float).unique(),
            )
            axis.set_yscale("log")
            axis.grid(alpha=0.20)
            if row_index == 0:
                axis.set_title(shape, fontweight="bold")
            if column_index == 0:
                axis.set_ylabel(
                    "{}\nMean distance $D/c$".format(
                        {
                            "Natural Killer Cells": "NK",
                            "Neutrophils": "Neutrophils",
                            "T Cells": "T cells",
                            "B Cells": "B cells",
                        }[cell_type]
                    )
                )
            if row_index == len(row_cells) - 1:
                axis.set_xlabel(r"Projected area $A/c^2$")

    condition_handles = [
        Line2D(
            [0],
            [0],
            color=condition_colors[condition],
            linewidth=2.0,
            label=condition_c_label(labels[condition], crossovers[condition]).replace(
                "\n", " · "
            ),
        )
        for condition in conditions
    ]
    class_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=style[0],
            marker=style[1],
            markersize=3.5,
            linewidth=1.1,
            label=style[2],
        )
        for style in class_styles.values()
    ]
    fig.legend(
        handles=condition_handles,
        loc="upper center",
        ncol=3,
        fontsize=6.5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.972),
        columnspacing=1.1,
    )
    fig.legend(
        handles=class_handles,
        loc="lower center",
        ncol=2,
        fontsize=7.0,
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )
    fig.suptitle(
        "Unbounded first-passage distance on each condition's line-feasible grid",
        y=0.997,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.82))
    paths = _save_figure(
        fig,
        figure_directory,
        "unbounded_detection_distance_by_condition",
    )
    plt.close(fig)
    return paths


def _plot_detailed_search_performance(
    fits: pd.DataFrame,
    finite: pd.DataFrame,
    unbounded: pd.DataFrame,
    ratios: pd.DataFrame,
    figure_directory: Path,
) -> List[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    fit_index = fits.set_index("condition")
    missing = [condition for condition in SELECTED_CONDITIONS if condition not in fit_index.index]
    if missing:
        raise ValueError("Selected fitted conditions are missing")
    finite = finite[finite["condition"].isin(SELECTED_CONDITIONS)].copy()
    unbounded = unbounded[unbounded["condition"].isin(SELECTED_CONDITIONS)].copy()
    ratios = ratios[ratios["condition"].isin(SELECTED_CONDITIONS)].copy()

    areas_by_condition = {
        condition: tuple(
            sorted(
                finite.loc[
                    finite["condition"] == condition,
                    "projected_area",
                ].astype(float).unique()
            )
        )
        for condition in SELECTED_CONDITIONS
    }

    def title(condition: str) -> str:
        return condition_c_label(
            COMPACT_CONDITION_LABELS[condition].replace("\n", ", "),
            float(fit_index.loc[condition, "crossover"]),
        )

    fig = plt.figure(figsize=(7.01, DETAILED_SEARCH_HEIGHT_INCHES))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=(3.35, 1.30, 1.15),
        left=0.12,
        right=0.99,
        bottom=0.085,
        top=0.865,
        hspace=0.58,
    )
    top = outer[0].subgridspec(3, 3, wspace=0.25, hspace=0.30)
    raw = outer[1].subgridspec(1, 3, wspace=0.31)
    bottom = outer[2].subgridspec(1, 3, wspace=0.27)

    panel_index = 0
    for row, condition in enumerate(SELECTED_CONDITIONS):
        condition_data = finite[finite["condition"] == condition]
        row_upper = _finite_upper_limit(condition_data)
        for column, shape in enumerate(SHAPE_ORDER):
            axis = fig.add_subplot(top[row, column])
            cell = condition_data[condition_data["shape"] == shape]
            for class_name in FINITE_PLOT_ORDER:
                style = CLASS_STYLES[class_name]
                group = cell[cell["trajectory_class"] == class_name].sort_values(
                    "projected_area"
                )
                x = group["projected_area"].to_numpy(dtype=float)
                y = 100.0 * group["detection_probability"].to_numpy(dtype=float)
                low = 100.0 * group["detection_probability_ci_low"].to_numpy(
                    dtype=float
                )
                high = 100.0 * group["detection_probability_ci_high"].to_numpy(
                    dtype=float
                )
                axis.errorbar(
                    x,
                    y,
                    yerr=np.vstack((y - low, high - y)),
                    color=style["color"],
                    linestyle=style["linestyle"],
                    marker=style["marker"],
                    markersize=2.6,
                    markeredgewidth=0.5,
                    linewidth=0.85,
                    elinewidth=0.4,
                    capsize=1.0,
                    alpha=0.96,
                )
            _set_normalized_area_axis(axis, areas_by_condition[condition])
            axis.set_ylim(0.0, row_upper)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
            axis.grid(axis="y", color="0.90", linewidth=0.45)
            if row < 2:
                axis.tick_params(labelbottom=False)
            else:
                axis.set_xlabel(r"Projected area $A/c^2$")
            if column == 0:
                axis.set_ylabel(
                    "{}\nDetection (%)".format(FINITE_ROW_LABELS[condition])
                )
            else:
                axis.tick_params(labelleft=False)
            if row == 0:
                axis.set_title(shape, pad=2.0, fontweight="bold")
            _panel_label(axis, chr(ord("A") + panel_index))
            panel_index += 1

    class_handles = [
        Line2D(
            [0],
            [0],
            color=CLASS_STYLES[class_name]["color"],
            linestyle=CLASS_STYLES[class_name]["linestyle"],
            marker=CLASS_STYLES[class_name]["marker"],
            markersize=3.5,
            linewidth=1.0,
            label=str(CLASS_STYLES[class_name]["label"]),
        )
        for class_name in FINITE_PLOT_ORDER
    ]
    fig.legend(
        handles=class_handles,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.982),
        ncol=3,
        frameon=False,
        columnspacing=1.25,
        handlelength=2.0,
    )
    fig.text(
        0.55,
        0.897,
        "Finite search: means and 95% whole-track bootstrap intervals",
        ha="center",
        va="center",
        fontsize=6.5,
    )

    raw_process_styles = {
        "vector_uniform_rotated": ("-", "Empirical run pool"),
        "fitted_mu": ("--", r"Fitted $\mu$ model"),
    }
    raw_axes: List[object] = []
    for column, condition in enumerate(SELECTED_CONDITIONS):
        axis = fig.add_subplot(raw[0, column])
        raw_axes.append(axis)
        condition_data = unbounded[unbounded["condition"] == condition]
        for class_name in UNBOUNDED_CLASSES:
            process_linestyle, _ = raw_process_styles[class_name]
            for shape in SHAPE_ORDER:
                shape_style = SHAPE_STYLES[shape]
                group = condition_data[
                    (condition_data["trajectory_class"] == class_name)
                    & (condition_data["shape"] == shape)
                ].sort_values("projected_area")
                x = group["projected_area"].to_numpy(dtype=float)
                y = group["mean_detection_distance_model_units"].to_numpy(dtype=float)
                se = group["detection_distance_mc_se_model_units"].to_numpy(
                    dtype=float
                )
                axis.errorbar(
                    x,
                    y,
                    yerr=_positive_log_errors(y, se),
                    color=shape_style["color"],
                    linestyle=process_linestyle,
                    marker=shape_style["marker"],
                    markersize=2.7,
                    markeredgewidth=0.45,
                    linewidth=0.85,
                    elinewidth=0.4,
                    capsize=1.0,
                    alpha=0.96,
                )
        axis.set_yscale("log")
        _set_normalized_area_axis(axis, areas_by_condition[condition])
        axis.tick_params(labelbottom=False)
        axis.grid(which="major", axis="y", color="0.90", linewidth=0.45)
        axis.set_title(title(condition), pad=2.0, fontsize=8.0)
        if column == 0:
            axis.set_ylabel(r"Mean distance $D/c$")
        _panel_label(axis, chr(ord("J") + column))

    raw_means = unbounded["mean_detection_distance_model_units"].to_numpy(
        dtype=float
    )
    raw_errors = unbounded[
        "detection_distance_mc_se_model_units"
    ].to_numpy(dtype=float)
    raw_lower = np.maximum(
        raw_means - 1.96 * raw_errors,
        np.maximum(np.finfo(float).tiny, raw_means * 1e-3),
    )
    raw_upper = raw_means + 1.96 * raw_errors
    common_raw_low = float(np.min(raw_lower))
    common_raw_high = float(np.max(raw_upper))
    log_padding = float((common_raw_high / common_raw_low) ** 0.04)
    for axis in raw_axes:
        axis.set_ylim(
            common_raw_low / log_padding,
            common_raw_high * log_padding,
        )

    ratio_low = float(ratios["ratio_ci_low"].min())
    ratio_high = float(ratios["ratio_ci_high"].max())
    padding = max(0.025, 0.08 * (ratio_high - ratio_low))
    common_low = max(0.0, min(1.0, ratio_low) - padding)
    common_high = max(1.0, ratio_high) + padding
    bottom_axes: List[object] = []
    for column, condition in enumerate(SELECTED_CONDITIONS):
        axis = fig.add_subplot(bottom[0, column])
        bottom_axes.append(axis)
        condition_data = ratios[ratios["condition"] == condition]
        for shape in SHAPE_ORDER:
            style = SHAPE_STYLES[shape]
            group = condition_data[condition_data["shape"] == shape].sort_values(
                "projected_area"
            )
            x = group["projected_area"].to_numpy(dtype=float)
            y = group["fitted_to_empirical_ratio"].to_numpy(dtype=float)
            low = group["ratio_ci_low"].to_numpy(dtype=float)
            high = group["ratio_ci_high"].to_numpy(dtype=float)
            axis.errorbar(
                x,
                y,
                yerr=np.vstack((y - low, high - y)),
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                markersize=2.8,
                linewidth=0.85,
                elinewidth=0.5,
                capsize=1.2,
            )
        axis.axhline(1.0, color="0.50", linestyle=":", linewidth=0.75, zorder=0)
        _set_normalized_area_axis(axis, areas_by_condition[condition])
        axis.set_ylim(common_low, common_high)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.grid(axis="y", color="0.90", linewidth=0.45)
        axis.set_xlabel(r"Projected area $A/c^2$")
        if column == 0:
            axis.set_ylabel(r"Fitted / empirical $D$")
        else:
            axis.tick_params(labelleft=False)
        _panel_label(axis, chr(ord("M") + column))

    shape_handles = [
        Line2D(
            [0],
            [0],
            color=SHAPE_STYLES[shape]["color"],
            linestyle=SHAPE_STYLES[shape]["linestyle"],
            marker=SHAPE_STYLES[shape]["marker"],
            markersize=3.5,
            linewidth=1.0,
            label=shape,
        )
        for shape in SHAPE_ORDER
    ]
    process_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=style[0],
            linewidth=1.0,
            label=style[1],
        )
        for style in raw_process_styles.values()
    ]
    fig.legend(
        handles=process_handles + shape_handles,
        loc="lower center",
        bbox_to_anchor=(0.55, 0.006),
        ncol=5,
        frameon=False,
        columnspacing=1.0,
        handlelength=2.0,
    )
    ratio_heading_y = max(axis.get_position().y1 for axis in bottom_axes) + 0.020
    fig.text(
        0.55,
        ratio_heading_y,
        "Normalized comparison: fitted / empirical mean distance",
        ha="center",
        va="center",
        fontsize=6.5,
    )
    paths = _save_figure(fig, figure_directory, "figure2_search_performance")
    plt.close(fig)
    return paths


def plot_line_feasible_figures(
    output_directory: Path,
    fits: pd.DataFrame,
    finite: pd.DataFrame,
    unbounded: pd.DataFrame,
    ratios: pd.DataFrame,
) -> List[Path]:
    """Write both requested PDF/PNG pairs from validated normalized tables."""
    import matplotlib

    matplotlib.use("Agg")
    figure_directory = Path(output_directory) / "figures"
    paths: List[Path] = []
    paths.extend(_plot_all_condition_unbounded(fits, unbounded, figure_directory))
    paths.extend(
        _plot_detailed_search_performance(
            fits,
            finite,
            unbounded,
            ratios,
            figure_directory,
        )
    )
    return paths


__all__ = [
    "condition_c_label",
    "normalized_area_ticks",
    "plot_line_feasible_figures",
]
