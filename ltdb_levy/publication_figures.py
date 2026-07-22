"""Publication-ready figures for three prespecified LTDB conditions.

This is an additive presentation stage: it reads the signed all-condition
replay and the independently rotated-vector augmentation, validates the full
selected grid, and writes compact PNAS-width figures plus auditable data
subsets.  It never modifies either source stage.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactStore
from .dataio.artifacts import sha256_file
from .dist.mixture import cdf as mixture_cdf
from .logging_utils import get_logger
from .nk_experiment import (
    FINITE_CONTRASTS,
    FINITE_PLOT_ORDER,
    UNBOUNDED_CLASSES,
)

logger = get_logger(__name__)

BASELINE_STAGE = "all_conditions_replay"
AUGMENTED_STAGE = "ordered_length_rotation_replay"
OUTPUT_STAGE = "publication_selected_conditions"
SELECTION_VERSION = 2
ONE_COLUMN_WIDTH_INCHES = 8.7 / 2.54
FINITE_COMPACT_HEIGHT_INCHES = 2.72
UNBOUNDED_COMPACT_HEIGHT_INCHES = 2.12
DETAILED_SEARCH_HEIGHT_INCHES = 9.25
FINITE_COMPACT_BOOTSTRAP_REPLICATES = 10_000
FINITE_COMPACT_BOOTSTRAP_SEED_OFFSET = 71_903

SELECTED_CONDITIONS: Tuple[str, ...] = (
    (
        "cell_type=Natural Killer Cells | organ=popliteal lymph node | "
        "stimulus=Influenza Vaccine | dt_seconds=30.0"
    ),
    (
        "cell_type=Neutrophils | organ=popliteal lymph node | "
        "stimulus=Influenza Vaccine | dt_seconds=15.0"
    ),
    (
        "cell_type=T Cells | organ=popliteal lymph node | "
        "stimulus=HIV-infected humanized T cell | dt_seconds=15.0"
    ),
)

SELECTED_LABELS: Tuple[str, ...] = (
    "NK · PLN · influenza · 30 s",
    "Neutrophils · PLN · influenza · 15 s",
    "T · PLN · HIV model · 15 s",
)

SHAPE_ORDER: Tuple[str, ...] = ("Ball", "Disk", "Line")

CLASS_STYLES: Mapping[str, Mapping[str, object]] = {
    "exact_orientation": {
        "label": "Observed orientation",
        "color": "#000000",
        "linestyle": "-",
        "marker": "o",
    },
    "exact_global_rotation": {
        "label": "Whole-track rotation",
        "color": "#777777",
        "linestyle": "--",
        "marker": "s",
    },
    "ordered_length_rotated": {
        "label": "Per-vector rotation",
        "color": "#0072B2",
        "linestyle": "-.",
        "marker": "D",
    },
    "vector_uniform_rotated": {
        "label": "Empirical run pool",
        "color": "#009E73",
        "linestyle": ":",
        "marker": "^",
    },
    "fitted_mu": {
        "label": r"Fitted $\mu$ model",
        "color": "#D55E00",
        "linestyle": (0, (5, 1, 1, 1)),
        "marker": "X",
    },
}

COMPACT_CLASS_STYLES: Mapping[str, Mapping[str, object]] = {
    "exact_orientation": {
        "label": "Observed",
        "color": "#FFFFFF",
        "hatch": "////",
    },
    "exact_global_rotation": {
        "label": "Whole-track",
        "color": "#8C8C8C",
        "hatch": "",
    },
    "ordered_length_rotated": {
        "label": "Per-vector",
        "color": "#56B4E9",
        "hatch": "...",
    },
    "vector_uniform_rotated": {
        "label": "Run pool",
        "color": "#009E73",
        "hatch": "\\\\\\",
    },
    "fitted_mu": {
        "label": r"Fitted $\mu$",
        "color": "#D55E00",
        "hatch": "xx",
    },
}

COMPACT_CONDITION_LABELS: Mapping[str, str] = {
    SELECTED_CONDITIONS[0]: "NK\n30 s",
    SELECTED_CONDITIONS[1]: "Neutrophils\n15 s",
    SELECTED_CONDITIONS[2]: "T–HIV\n15 s",
}

SHAPE_STYLES: Mapping[str, Mapping[str, object]] = {
    "Ball": {
        "color": "#000000",
        "marker": "o",
        "linestyle": "-",
    },
    "Disk": {
        "color": "#0072B2",
        "marker": "s",
        "linestyle": "--",
    },
    "Line": {
        "color": "#D55E00",
        "marker": "^",
        "linestyle": "-.",
    },
}

SHORT_CONDITION_LABELS: Mapping[str, str] = {
    SELECTED_CONDITIONS[0]: "NK, influenza, 30 s",
    SELECTED_CONDITIONS[1]: "Neutrophils, influenza, 15 s",
    SELECTED_CONDITIONS[2]: "T cells, HIV model, 15 s",
}

FINITE_ROW_LABELS: Mapping[str, str] = {
    SELECTED_CONDITIONS[0]: "NK (30 s)",
    SELECTED_CONDITIONS[1]: "Neutrophils (15 s)",
    SELECTED_CONDITIONS[2]: "T–HIV (15 s)",
}


def _condition_order(frame: pd.DataFrame) -> pd.DataFrame:
    order = {condition: index for index, condition in enumerate(SELECTED_CONDITIONS)}
    result = frame.copy()
    result["_condition_order"] = result["condition"].map(order)
    sort_columns = ["_condition_order"]
    for name in (
        "trajectory_class",
        "shape",
        "projected_area",
        "run_length_model_units",
    ):
        if name in result.columns:
            sort_columns.append(name)
    return (
        result.sort_values(sort_columns, kind="mergesort")
        .drop(columns="_condition_order")
        .reset_index(drop=True)
    )


def _select_conditions(frame: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """Select all and only the prespecified conditions, failing if one is absent."""
    if "condition" not in frame.columns:
        raise ValueError("{} has no condition column".format(table_name))
    selected = frame[frame["condition"].isin(SELECTED_CONDITIONS)].copy()
    present = set(selected["condition"].astype(str))
    missing = [
        condition for condition in SELECTED_CONDITIONS if condition not in present
    ]
    if missing:
        raise ValueError(
            "{} is missing selected condition(s): {}".format(
                table_name, "; ".join(missing)
            )
        )
    return _condition_order(selected)


def _validate_unique_grid(
    frame: pd.DataFrame,
    classes: Sequence[str],
    shapes: Sequence[str],
    areas: Sequence[float],
    table_name: str,
) -> None:
    expected_conditions = set(SELECTED_CONDITIONS)
    if set(frame["condition"].astype(str)) != expected_conditions:
        raise ValueError("{} has an unexpected condition set".format(table_name))
    if set(frame["trajectory_class"].astype(str)) != set(classes):
        raise ValueError("{} has an unexpected trajectory-class set".format(table_name))
    if set(frame["shape"].astype(str)) != set(shapes):
        raise ValueError("{} has an unexpected target-shape set".format(table_name))
    actual_areas = sorted(frame["projected_area"].astype(float).unique())
    expected_areas = sorted(float(value) for value in areas)
    if len(actual_areas) != len(expected_areas) or not np.allclose(
        actual_areas, expected_areas, rtol=0.0, atol=1e-12
    ):
        raise ValueError("{} has an unexpected projected-area grid".format(table_name))
    keys = ["condition", "trajectory_class", "shape", "projected_area"]
    if frame.duplicated(keys).any():
        raise ValueError("{} contains duplicate grid cells".format(table_name))
    expected_rows = (
        len(SELECTED_CONDITIONS) * len(classes) * len(shapes) * len(areas)
    )
    if len(frame) != expected_rows:
        raise ValueError(
            "{} has {} rows; expected {}".format(
                table_name, len(frame), expected_rows
            )
        )


def _validate_selected_data(
    fits: pd.DataFrame,
    empirical_runs: pd.DataFrame,
    finite: pd.DataFrame,
    unbounded: pd.DataFrame,
) -> Tuple[float, ...]:
    """Validate selection, normalization, uncertainty columns, and plot grids."""
    if list(fits["condition"].astype(str)) != list(SELECTED_CONDITIONS):
        raise ValueError("fit_summary does not have the required condition order")
    if fits.duplicated("condition").any() or len(fits) != len(SELECTED_CONDITIONS):
        raise ValueError("fit_summary must contain one row per selected condition")

    fit_by_condition = fits.set_index("condition")
    if empirical_runs["run_uid"].duplicated().any():
        raise ValueError("selected empirical runs contain duplicate run_uid values")
    for condition in SELECTED_CONDITIONS:
        fit = fit_by_condition.loc[condition]
        condition_runs = empirical_runs[empirical_runs["condition"] == condition]
        if len(condition_runs) != int(fit["n_runs"]):
            raise ValueError(
                "{} empirical runs do not match fit n_runs".format(condition)
            )
        expected = (
            condition_runs["run_length_um"].to_numpy(dtype=float)
            / float(fit["crossover"])
        )
        actual = condition_runs["run_length_model_units"].to_numpy(dtype=float)
        if not np.allclose(actual, expected, rtol=1e-10, atol=1e-12):
            raise ValueError(
                "{} empirical run normalization is inconsistent".format(condition)
            )

    shapes = tuple(
        shape for shape in SHAPE_ORDER if shape in set(finite["shape"].astype(str))
    )
    if shapes != SHAPE_ORDER:
        raise ValueError("finite summary must contain Ball, Disk, and Line")
    areas = tuple(sorted(finite["projected_area"].astype(float).unique()))
    _validate_unique_grid(
        finite,
        FINITE_PLOT_ORDER,
        shapes,
        areas,
        "selected finite summary",
    )
    _validate_unique_grid(
        unbounded,
        UNBOUNDED_CLASSES,
        shapes,
        areas,
        "selected unbounded summary",
    )

    for condition in SELECTED_CONDITIONS:
        fit = fit_by_condition.loc[condition]
        finite_tracks = finite.loc[
            finite["condition"] == condition, "n_tracks"
        ].to_numpy(dtype=int)
        if not np.all(finite_tracks == int(fit["n_tracks"])):
            raise ValueError(
                "{} finite n_tracks does not match fit n_tracks".format(condition)
            )

    probabilities = finite["detection_probability"].to_numpy(dtype=float)
    lower = finite["detection_probability_ci_low"].to_numpy(dtype=float)
    upper = finite["detection_probability_ci_high"].to_numpy(dtype=float)
    if (
        np.any(~np.isfinite(probabilities))
        or np.any((probabilities < 0.0) | (probabilities > 1.0))
        or np.any(lower > probabilities)
        or np.any(upper < probabilities)
    ):
        raise ValueError("finite probability intervals are invalid")

    means = unbounded["mean_detection_distance_model_units"].to_numpy(dtype=float)
    standard_errors = unbounded[
        "detection_distance_mc_se_model_units"
    ].to_numpy(dtype=float)
    if (
        np.any(~np.isfinite(means))
        or np.any(means <= 0.0)
        or np.any(~np.isfinite(standard_errors))
        or np.any(standard_errors < 0.0)
    ):
        raise ValueError("unbounded distance means or standard errors are invalid")
    return areas


def _validate_selected_auxiliary_data(
    fits: pd.DataFrame,
    finite_contrasts: pd.DataFrame,
    unbounded_contrasts: pd.DataFrame,
    track_inputs: pd.DataFrame,
    areas: Sequence[float],
) -> None:
    """Validate the selected contrast and per-track tables used by the report."""
    expected_conditions = set(SELECTED_CONDITIONS)
    expected_areas = sorted(float(value) for value in areas)
    if set(track_inputs["condition"].astype(str)) != expected_conditions:
        raise ValueError("selected track_inputs has an unexpected condition set")
    if track_inputs.duplicated(["condition", "track_uid"]).any():
        raise ValueError("selected track_inputs contains duplicate tracks")
    expected_track_counts = fits.set_index("condition")["n_tracks"].astype(int)
    actual_track_counts = track_inputs.groupby("condition")["track_uid"].nunique()
    if any(
        int(actual_track_counts.loc[condition])
        != int(expected_track_counts.loc[condition])
        for condition in SELECTED_CONDITIONS
    ):
        raise ValueError("selected track_inputs does not match fitted track counts")

    finite_keys = [
        "condition",
        "left_class",
        "right_class",
        "shape",
        "projected_area",
    ]
    if finite_contrasts.duplicated(finite_keys).any():
        raise ValueError("selected finite_contrasts contains duplicate cells")
    actual_pairs = set(
        finite_contrasts[["left_class", "right_class"]].itertuples(
            index=False, name=None
        )
    )
    if actual_pairs != set(FINITE_CONTRASTS):
        raise ValueError("selected finite_contrasts has an unexpected contrast set")
    expected_finite_rows = (
        len(SELECTED_CONDITIONS)
        * len(FINITE_CONTRASTS)
        * len(SHAPE_ORDER)
        * len(expected_areas)
    )
    if (
        len(finite_contrasts) != expected_finite_rows
        or set(finite_contrasts["condition"].astype(str)) != expected_conditions
        or set(finite_contrasts["shape"].astype(str)) != set(SHAPE_ORDER)
        or not np.allclose(
            sorted(finite_contrasts["projected_area"].astype(float).unique()),
            expected_areas,
            rtol=0.0,
            atol=1e-12,
        )
    ):
        raise ValueError("selected finite_contrasts has an incomplete grid")

    unbounded_keys = ["condition", "shape", "projected_area"]
    expected_unbounded_rows = (
        len(SELECTED_CONDITIONS) * len(SHAPE_ORDER) * len(expected_areas)
    )
    if (
        unbounded_contrasts.duplicated(unbounded_keys).any()
        or len(unbounded_contrasts) != expected_unbounded_rows
        or set(unbounded_contrasts["condition"].astype(str))
        != expected_conditions
        or set(unbounded_contrasts["shape"].astype(str)) != set(SHAPE_ORDER)
        or not np.allclose(
            sorted(
                unbounded_contrasts["projected_area"].astype(float).unique()
            ),
            expected_areas,
            rtol=0.0,
            atol=1e-12,
        )
    ):
        raise ValueError("selected unbounded_contrasts has an incomplete grid")


def _empirical_ccdf(lengths: Sequence[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Return ordered lengths and the non-strict empirical Pr[L >= l]."""
    values = np.sort(np.asarray(lengths, dtype=float))
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("CCDF requires at least one length")
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("CCDF lengths must be finite and positive")
    survival = (len(values) - np.arange(len(values), dtype=float)) / len(values)
    return values, survival


def _build_ccdf_table(
    fits: pd.DataFrame, empirical_runs: pd.DataFrame
) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    fit_by_condition = fits.set_index("condition")
    for condition in SELECTED_CONDITIONS:
        subset = empirical_runs[empirical_runs["condition"] == condition]
        lengths, survival = _empirical_ccdf(
            subset["run_length_model_units"].to_numpy(dtype=float)
        )
        fit = fit_by_condition.loc[condition]
        rows.append(
            pd.DataFrame(
                {
                    "condition": condition,
                    "condition_id": fit["condition_id"],
                    "condition_label": fit["condition_label"],
                    "order_index": np.arange(len(lengths), dtype=int),
                    "run_length_over_c": lengths,
                    "empirical_pr_length_ge": survival,
                }
            )
        )
    return _condition_order(pd.concat(rows, ignore_index=True))


def _unbounded_ratio_table(unbounded: pd.DataFrame) -> pd.DataFrame:
    """Derive fitted/empirical mean-distance ratios and 95% Monte Carlo CIs."""
    index = ["condition", "condition_id", "condition_label", "shape", "projected_area"]
    means = unbounded.pivot(index=index, columns="trajectory_class", values="mean_detection_distance_model_units")
    errors = unbounded.pivot(
        index=index,
        columns="trajectory_class",
        values="detection_distance_mc_se_model_units",
    )
    fitted = means["fitted_mu"].to_numpy(dtype=float)
    empirical = means["vector_uniform_rotated"].to_numpy(dtype=float)
    fitted_se = errors["fitted_mu"].to_numpy(dtype=float)
    empirical_se = errors["vector_uniform_rotated"].to_numpy(dtype=float)
    ratio = fitted / empirical
    # Independent-sample delta method.  Trial streams were not designed as
    # paired draws, so no unverified covariance reduction is introduced.
    ratio_se = np.sqrt(
        np.square(fitted_se / empirical)
        + np.square(fitted * empirical_se / np.square(empirical))
    )
    log_ratio_se = np.sqrt(
        np.square(fitted_se / fitted)
        + np.square(empirical_se / empirical)
    )
    table = means.reset_index()[index].copy()
    table["fitted_mean_distance_model_units"] = fitted
    table["fitted_mc_se_model_units"] = fitted_se
    table["empirical_mean_distance_model_units"] = empirical
    table["empirical_mc_se_model_units"] = empirical_se
    table["fitted_to_empirical_ratio"] = ratio
    table["ratio_mc_se_delta"] = ratio_se
    table["log_ratio_mc_se_delta"] = log_ratio_se
    table["ratio_ci_low"] = np.exp(np.log(ratio) - 1.96 * log_ratio_se)
    table["ratio_ci_high"] = np.exp(np.log(ratio) + 1.96 * log_ratio_se)
    table["ci_method"] = (
        "independent-sample log-ratio delta method, normal 95% interval"
    )
    return _condition_order(table)


def _validate_selected_finite_tracks(
    fits: pd.DataFrame,
    finite_tracks: pd.DataFrame,
    areas: Sequence[float],
) -> None:
    """Validate the complete track-by-class-by-target grid used for bars."""
    required = {
        "condition",
        "condition_id",
        "condition_label",
        "track_uid",
        "trajectory_class",
        "shape",
        "projected_area",
        "n_trials",
        "n_detected",
        "detection_probability",
    }
    missing = sorted(required - set(finite_tracks.columns))
    if missing:
        raise ValueError(
            "selected finite track summaries are missing columns: {}".format(
                ", ".join(missing)
            )
        )
    if set(finite_tracks["condition"].astype(str)) != set(SELECTED_CONDITIONS):
        raise ValueError(
            "selected finite track summaries have an unexpected condition set"
        )
    if set(finite_tracks["trajectory_class"].astype(str)) != set(
        FINITE_PLOT_ORDER
    ):
        raise ValueError(
            "selected finite track summaries have an unexpected class set"
        )
    if set(finite_tracks["shape"].astype(str)) != set(SHAPE_ORDER):
        raise ValueError(
            "selected finite track summaries have an unexpected shape set"
        )
    actual_areas = sorted(
        finite_tracks["projected_area"].astype(float).unique()
    )
    expected_areas = sorted(float(value) for value in areas)
    if len(actual_areas) != len(expected_areas) or not np.allclose(
        actual_areas,
        expected_areas,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(
            "selected finite track summaries have an unexpected area grid"
        )

    keys = [
        "condition",
        "track_uid",
        "trajectory_class",
        "shape",
        "projected_area",
    ]
    if finite_tracks.duplicated(keys).any():
        raise ValueError("selected finite track summaries contain duplicate cells")
    expected_cells = len(SHAPE_ORDER) * len(expected_areas)
    cell_counts = finite_tracks.groupby(
        ["condition", "track_uid", "trajectory_class"], sort=False
    ).size()
    if not np.all(cell_counts.to_numpy(dtype=int) == expected_cells):
        raise ValueError(
            "selected finite track summaries contain incomplete target grids"
        )

    fitted_track_counts = fits.set_index("condition")["n_tracks"].astype(int)
    for condition in SELECTED_CONDITIONS:
        subset = finite_tracks[finite_tracks["condition"] == condition]
        expected_tracks = int(fitted_track_counts.loc[condition])
        class_track_counts = subset.groupby("trajectory_class")[
            "track_uid"
        ].nunique()
        if (
            set(class_track_counts.index.astype(str)) != set(FINITE_PLOT_ORDER)
            or np.any(class_track_counts.to_numpy(dtype=int) != expected_tracks)
        ):
            raise ValueError(
                "{} compact finite data do not contain every fitted track "
                "in every class".format(condition)
            )
        class_track_sets = [
            frozenset(
                subset.loc[
                    subset["trajectory_class"] == class_name, "track_uid"
                ].astype(str)
            )
            for class_name in FINITE_PLOT_ORDER
        ]
        if len(set(class_track_sets)) != 1:
            raise ValueError(
                "{} compact finite classes do not share the same tracks".format(
                    condition
                )
            )

    probabilities = finite_tracks["detection_probability"].to_numpy(dtype=float)
    trials = finite_tracks["n_trials"].to_numpy(dtype=float)
    detected = finite_tracks["n_detected"].to_numpy(dtype=float)
    if (
        np.any(~np.isfinite(probabilities))
        or np.any((probabilities < 0.0) | (probabilities > 1.0))
        or np.any(~np.isfinite(trials))
        or np.any(trials <= 0.0)
        or np.any(detected < 0.0)
        or np.any(detected > trials)
        or not np.allclose(probabilities, detected / trials, atol=1e-12)
    ):
        raise ValueError("selected finite track probabilities are invalid")


def _finite_compact_ratio_table(
    finite_tracks: pd.DataFrame,
    bootstrap_replicates: int = FINITE_COMPACT_BOOTSTRAP_REPLICATES,
    bootstrap_seed: int = 12345,
) -> pd.DataFrame:
    """Aggregate the target grid and bootstrap ratios by whole track.

    Each track first receives equal weight in each of the nine prespecified
    target cells.  Tracks are then averaged equally within a condition.  A
    paired cluster bootstrap resamples tracks and recomputes all five class
    means before division by the whole-track-rotation reference.
    """
    if bootstrap_replicates < 100:
        raise ValueError("bootstrap_replicates must be at least 100")
    reference_class = "exact_global_rotation"
    present = set(finite_tracks["condition"].astype(str))
    conditions = [
        condition for condition in SELECTED_CONDITIONS if condition in present
    ]
    rows: List[Dict[str, object]] = []
    seed_sequence = np.random.SeedSequence(int(bootstrap_seed))
    condition_seeds = seed_sequence.spawn(len(conditions))

    for condition_index, (condition, condition_seed) in enumerate(
        zip(conditions, condition_seeds)
    ):
        subset = finite_tracks[finite_tracks["condition"] == condition]
        metadata = subset.iloc[0]
        track_means = (
            subset.groupby(
                ["track_uid", "trajectory_class"], sort=True
            )["detection_probability"]
            .mean()
            .unstack("trajectory_class")
        )
        if set(track_means.columns.astype(str)) != set(FINITE_PLOT_ORDER):
            raise ValueError(
                "{} compact finite table has an incomplete class set".format(
                    condition
                )
            )
        track_means = track_means.loc[:, list(FINITE_PLOT_ORDER)]
        if track_means.isna().any().any():
            raise ValueError(
                "{} compact finite table has incomplete paired tracks".format(
                    condition
                )
            )
        values = track_means.to_numpy(dtype=float)
        class_means = values.mean(axis=0)
        reference_index = list(FINITE_PLOT_ORDER).index(reference_class)
        reference_mean = float(class_means[reference_index])
        if reference_mean <= 0.0:
            raise ValueError(
                "{} has zero grid-averaged whole-track detection".format(condition)
            )

        generator = np.random.default_rng(condition_seed)
        bootstrap_ratios = np.full(
            (bootstrap_replicates, len(FINITE_PLOT_ORDER)),
            np.nan,
            dtype=float,
        )
        batch_size = min(1000, bootstrap_replicates)
        for start in range(0, bootstrap_replicates, batch_size):
            stop = min(start + batch_size, bootstrap_replicates)
            indices = generator.integers(
                0,
                len(values),
                size=(stop - start, len(values)),
            )
            sampled_means = values[indices, :].mean(axis=1)
            sampled_reference = sampled_means[:, reference_index]
            valid = sampled_reference > 0.0
            batch_ratios = np.full(
                (stop - start, len(FINITE_PLOT_ORDER)),
                np.nan,
                dtype=float,
            )
            batch_ratios[valid, :] = (
                sampled_means[valid, :]
                / sampled_reference[valid, np.newaxis]
            )
            bootstrap_ratios[start:stop, :] = batch_ratios

        valid_draws = np.isfinite(bootstrap_ratios[:, reference_index])
        valid_count = int(valid_draws.sum())
        if valid_count < max(100, int(np.ceil(0.95 * bootstrap_replicates))):
            raise ValueError(
                "{} has too many zero-reference bootstrap draws".format(condition)
            )
        lower, upper = np.nanpercentile(
            bootstrap_ratios,
            [2.5, 97.5],
            axis=0,
        )
        ratios = class_means / reference_mean
        for class_index, class_name in enumerate(FINITE_PLOT_ORDER):
            rows.append(
                {
                    "condition": condition,
                    "condition_id": metadata["condition_id"],
                    "condition_label": metadata["condition_label"],
                    "trajectory_class": class_name,
                    "experiment_label": COMPACT_CLASS_STYLES[class_name][
                        "label"
                    ],
                    "reference_class": reference_class,
                    "n_tracks": int(len(track_means)),
                    "n_target_cells": int(
                        subset[
                            subset["trajectory_class"] == class_name
                        ][["shape", "projected_area"]]
                        .drop_duplicates()
                        .shape[0]
                    ),
                    "target_grid_mean_detection_probability": float(
                        class_means[class_index]
                    ),
                    "reference_target_grid_mean_detection_probability": (
                        reference_mean
                    ),
                    "normalized_detection_probability": float(
                        ratios[class_index]
                    ),
                    "ratio_ci_low": float(lower[class_index]),
                    "ratio_ci_high": float(upper[class_index]),
                    "bootstrap_replicates_requested": int(
                        bootstrap_replicates
                    ),
                    "bootstrap_replicates_valid": valid_count,
                    "bootstrap_seed": int(bootstrap_seed),
                    "bootstrap_condition_index": int(condition_index),
                    "grid_weighting": (
                        "equal target-cell weight within track; equal track "
                        "weight within condition"
                    ),
                    "interval_method": (
                        "paired whole-track percentile bootstrap, 95%"
                    ),
                }
            )
    return pd.DataFrame(rows)


def _unbounded_compact_ratio_table(unbounded: pd.DataFrame) -> pd.DataFrame:
    """Aggregate nine target cells and normalize to the empirical run pool."""
    reference_class = "vector_uniform_rotated"
    present = set(unbounded["condition"].astype(str))
    conditions = [
        condition for condition in SELECTED_CONDITIONS if condition in present
    ]
    rows: List[Dict[str, object]] = []
    for condition in conditions:
        subset = unbounded[unbounded["condition"] == condition]
        metadata = subset.iloc[0]
        aggregates: Dict[str, Tuple[float, float, int, int]] = {}
        for class_name in UNBOUNDED_CLASSES:
            class_data = subset[
                subset["trajectory_class"] == class_name
            ]
            means = class_data[
                "mean_detection_distance_model_units"
            ].to_numpy(dtype=float)
            errors = class_data[
                "detection_distance_mc_se_model_units"
            ].to_numpy(dtype=float)
            if len(means) == 0:
                raise ValueError(
                    "{} lacks unbounded class {}".format(condition, class_name)
                )
            # The shape-area trial streams are separate simulation tasks.
            # Their independent MC variances therefore add when the nine
            # equally weighted cell means are averaged.
            aggregate_mean = float(means.mean())
            aggregate_se = float(np.sqrt(np.square(errors).sum()) / len(means))
            aggregates[class_name] = (
                aggregate_mean,
                aggregate_se,
                int(len(means)),
                int(class_data["n_trials"].sum()),
            )

        reference_mean, reference_se, _, _ = aggregates[reference_class]
        if reference_mean <= 0.0:
            raise ValueError(
                "{} has non-positive empirical mean distance".format(condition)
            )
        for class_name in UNBOUNDED_CLASSES:
            mean, standard_error, n_cells, n_trials = aggregates[class_name]
            ratio = mean / reference_mean
            if class_name == reference_class:
                log_ratio_se = 0.0
                ratio_low = 1.0
                ratio_high = 1.0
            else:
                log_ratio_se = float(
                    np.sqrt(
                        np.square(standard_error / mean)
                        + np.square(reference_se / reference_mean)
                    )
                )
                ratio_low = float(
                    np.exp(np.log(ratio) - 1.96 * log_ratio_se)
                )
                ratio_high = float(
                    np.exp(np.log(ratio) + 1.96 * log_ratio_se)
                )
            rows.append(
                {
                    "condition": condition,
                    "condition_id": metadata["condition_id"],
                    "condition_label": metadata["condition_label"],
                    "trajectory_class": class_name,
                    "experiment_label": COMPACT_CLASS_STYLES[class_name][
                        "label"
                    ],
                    "reference_class": reference_class,
                    "n_target_cells": n_cells,
                    "n_trials_total": n_trials,
                    "target_grid_mean_detection_distance_model_units": mean,
                    "target_grid_mean_mc_se_model_units": standard_error,
                    "reference_target_grid_mean_detection_distance_model_units": (
                        reference_mean
                    ),
                    "reference_target_grid_mean_mc_se_model_units": reference_se,
                    "normalized_detection_distance": float(ratio),
                    "ratio_ci_low": ratio_low,
                    "ratio_ci_high": ratio_high,
                    "log_ratio_mc_se_delta": log_ratio_se,
                    "grid_weighting": "equal weight across nine target cells",
                    "interval_method": (
                        "independent-stream log-ratio delta method, 95%"
                        if class_name != reference_class
                        else "normalized reference fixed at one"
                    ),
                }
            )
    return pd.DataFrame(rows)


def _configure_matplotlib() -> None:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.0,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7.0,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "lines.linewidth": 0.9,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _panel_label(axis: object, label: str) -> None:
    axis.text(
        0.015,
        0.985,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
        zorder=20,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.7, "alpha": 0.82},
    )


def _draw_segmentation_schematic(axis: object) -> None:
    """Draw a vector schematic without relying on an external image asset."""
    from matplotlib.patches import Arc

    points = np.array(
        [
            [0.05, 0.20],
            [0.42, 0.23],
            [0.78, 0.30],
            [0.58, 0.67],
            [0.38, 0.95],
            [0.72, 0.79],
            [0.98, 0.59],
        ]
    )
    for start, stop in ((0, 2), (4, 6)):
        axis.plot(
            points[start : stop + 1, 0],
            points[start : stop + 1, 1],
            color="0.58",
            linestyle="--",
            marker="o",
            markersize=2.5,
            linewidth=1.1,
        )
    axis.plot(
        points[2:5, 0],
        points[2:5, 1],
        color="#0072B2",
        marker="o",
        markersize=2.8,
        linewidth=1.6,
    )
    for start, stop in zip(points[:-1], points[1:]):
        delta = stop - start
        axis.annotate(
            "",
            xy=stop,
            xytext=start,
            arrowprops={
                "arrowstyle": "-|>",
                "color": "#0072B2" if np.array_equal(start, points[2]) or np.array_equal(start, points[3]) else "0.58",
                "linewidth": 0.8,
                "mutation_scale": 5,
                "shrinkA": 2,
                "shrinkB": 2,
            },
        )
    axis.add_patch(
        Arc(
            points[2],
            0.27,
            0.27,
            angle=0.0,
            theta1=118.0,
            theta2=185.0,
            color="#D55E00",
            linewidth=1.0,
        )
    )
    axis.text(0.69, 0.43, r"$\theta>90^\circ$", color="#D55E00", fontsize=6.7)
    axis.text(
        0.06,
        0.08,
        "boundary run\nexcluded",
        color="0.38",
        fontsize=6.5,
        ha="left",
    )
    axis.text(
        0.35,
        1.00,
        "complete run retained",
        color="#0072B2",
        fontsize=6.5,
        ha="center",
    )
    axis.text(
        0.98,
        0.73,
        "boundary\nexcluded",
        color="0.38",
        fontsize=6.5,
        ha="right",
    )
    axis.text(
        0.50,
        -0.04,
        r"turn $>90^\circ$ starts a new directed run",
        ha="center",
        va="top",
        fontsize=6.5,
    )
    axis.set_xlim(0.0, 1.02)
    axis.set_ylim(0.0, 1.12)
    axis.set_title("Directed-run segmentation", pad=2.0)
    axis.set_axis_off()


def _save_figure(fig: object, figure_directory: Path, stem: str) -> List[Path]:
    figure_directory.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for extension, kwargs in (
        ("pdf", {}),
        ("png", {"dpi": 300}),
    ):
        path = figure_directory / "{}.{}".format(stem, extension)
        temporary = figure_directory / ".{}.{}.tmp".format(stem, extension)
        fig.savefig(
            temporary,
            format=extension,
            facecolor="white",
            edgecolor="none",
            **kwargs
        )
        os.replace(str(temporary), str(path))
        paths.append(path)
    return paths


def _plot_figure1(
    fits: pd.DataFrame,
    empirical_runs: pd.DataFrame,
    figure_directory: Path,
) -> List[Path]:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.01, 2.8))
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=(1.12, 1.0, 1.0, 1.0),
        left=0.035,
        right=0.995,
        bottom=0.19,
        top=0.90,
        wspace=0.35,
    )
    schematic = fig.add_subplot(grid[0, 0])
    _draw_segmentation_schematic(schematic)
    _panel_label(schematic, "A")

    fit_by_condition = fits.set_index("condition")
    for index, condition in enumerate(SELECTED_CONDITIONS):
        axis = fig.add_subplot(grid[0, index + 1])
        fit = fit_by_condition.loc[condition]
        lengths = empirical_runs.loc[
            empirical_runs["condition"] == condition,
            "run_length_model_units",
        ].to_numpy(dtype=float)
        lengths, empirical = _empirical_ccdf(lengths)
        lmax_over_c = float(fit["lmax"]) / float(fit["crossover"])
        grid_values = np.geomspace(
            max(float(lengths.min()), 1e-5),
            lmax_over_c * (1.0 - 1e-9),
            500,
        )
        fitted = 1.0 - mixture_cdf(
            grid_values,
            float(fit["mu_hat"]),
            lmax_over_c,
            crossover=1.0,
        )
        axis.step(
            lengths,
            empirical,
            where="post",
            color="black",
            linewidth=0.85,
            label="Empirical",
        )
        axis.plot(
            grid_values,
            np.maximum(fitted, 0.25 / len(lengths)),
            color="#D55E00",
            linewidth=1.25,
            label="Fitted mixture",
        )
        axis.axvline(
            1.0,
            color="0.48",
            linestyle=":",
            linewidth=0.75,
            zorder=0,
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_ylim(1e-4, 1.15)
        axis.set_xlim(
            max(0.75 * float(lengths.min()), 1e-3),
            1.06 * max(float(lengths.max()), lmax_over_c),
        )
        axis.grid(which="major", color="0.90", linewidth=0.45)
        axis.set_title(SHORT_CONDITION_LABELS[condition], pad=2.0)
        if index == 0:
            axis.set_ylabel(r"$\Pr(L\geq \ell)$")
            axis.legend(
                loc="lower left",
                frameon=False,
                handlelength=1.5,
                borderaxespad=0.2,
            )
        else:
            axis.tick_params(labelleft=False)
        annotation = (
            "{:d} tracks; {:d} runs\n"
            r"$\mu={:.2f}$ [{:.2f}, {:.2f}]"
            "\n"
            r"$c={:.2f}\,\mu\mathrm{{m}}$; "
            r"$L_{{max}}={:.1f}\,\mu\mathrm{{m}}$"
            "\n"
            r"KS $p={:.3f}$"
        ).format(
            int(fit["n_tracks"]),
            int(fit["n_runs"]),
            float(fit["mu_hat"]),
            float(fit["mu_ci_low"]),
            float(fit["mu_ci_high"]),
            float(fit["crossover"]),
            float(fit["lmax"]),
            float(fit["ks_bootstrap_p"]),
        )
        axis.text(
            0.98,
            0.98,
            annotation,
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=6.5,
            linespacing=1.15,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.85",
                "linewidth": 0.4,
                "pad": 1.5,
                "alpha": 0.92,
            },
        )
        _panel_label(axis, chr(ord("B") + index))

    fig.text(
        0.625,
        0.055,
        r"Normalized directed-run length $\ell/c$",
        ha="center",
        va="center",
        fontsize=7.0,
    )
    paths = _save_figure(fig, figure_directory, "figure1_data_and_fit")
    plt.close(fig)
    return paths


def _finite_row_upper_limit(frame: pd.DataFrame) -> float:
    values = 100.0 * frame["detection_probability_ci_high"].to_numpy(dtype=float)
    maximum = float(np.nanmax(values))
    return max(1e-3, 1.13 * maximum)


def _plot_figure2(
    fits: pd.DataFrame,
    finite: pd.DataFrame,
    unbounded: pd.DataFrame,
    ratios: pd.DataFrame,
    areas: Sequence[float],
    figure_directory: Path,
) -> List[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

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
    fitted_mu_by_condition = (
        fits.set_index("condition")["mu_hat"].astype(float).to_dict()
    )

    def condition_title(condition: str) -> str:
        return "{}\n{}".format(
            SHORT_CONDITION_LABELS[condition],
            r"$\mu={:.3f}$".format(fitted_mu_by_condition[condition]),
        )

    top_axes: List[object] = []
    panel_index = 0

    for row, condition in enumerate(SELECTED_CONDITIONS):
        condition_data = finite[finite["condition"] == condition]
        row_upper = _finite_row_upper_limit(condition_data)
        for column, shape in enumerate(SHAPE_ORDER):
            axis = fig.add_subplot(top[row, column])
            top_axes.append(axis)
            cell = condition_data[condition_data["shape"] == shape]
            for class_name in FINITE_PLOT_ORDER:
                style = CLASS_STYLES[class_name]
                group = cell[
                    cell["trajectory_class"] == class_name
                ].sort_values("projected_area")
                x = group["projected_area"].to_numpy(dtype=float)
                y = 100.0 * group["detection_probability"].to_numpy(dtype=float)
                low = 100.0 * group[
                    "detection_probability_ci_low"
                ].to_numpy(dtype=float)
                high = 100.0 * group[
                    "detection_probability_ci_high"
                ].to_numpy(dtype=float)
                axis.errorbar(
                    x,
                    y,
                    yerr=np.vstack((y - low, high - y)),
                    color=style["color"],
                    linestyle=style["linestyle"],
                    marker=style["marker"],
                    markersize=2.8,
                    markeredgewidth=0.55,
                    linewidth=0.9,
                    elinewidth=0.45,
                    capsize=1.1,
                    alpha=0.96,
                )
            axis.set_xlim(min(areas) - 0.3, max(areas) + 0.3)
            axis.set_xticks(list(areas))
            axis.set_ylim(0.0, row_upper)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
            axis.grid(axis="y", color="0.90", linewidth=0.45)
            if row < 2:
                axis.tick_params(labelbottom=False)
            else:
                axis.set_xlabel(r"Projected area $A/c^2$")
            if column == 0:
                axis.set_ylabel(
                    "{}\nDetection (%)".format(
                        FINITE_ROW_LABELS[condition]
                    )
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
        "Finite search: points are means; bars are 95% whole-track bootstrap intervals",
        ha="center",
        va="center",
        fontsize=6.5,
    )

    raw_axes: List[object] = []
    raw_process_styles = {
        "vector_uniform_rotated": ("-", "Empirical run pool (raw)"),
        "fitted_mu": ("--", r"Fitted $\mu$ (raw)"),
    }
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
                y = group[
                    "mean_detection_distance_model_units"
                ].to_numpy(dtype=float)
                standard_error = group[
                    "detection_distance_mc_se_model_units"
                ].to_numpy(dtype=float)
                axis.errorbar(
                    x,
                    y,
                    yerr=1.96 * standard_error,
                    color=shape_style["color"],
                    linestyle=process_linestyle,
                    marker=shape_style["marker"],
                    markersize=2.9,
                    markeredgewidth=0.5,
                    linewidth=0.9,
                    elinewidth=0.45,
                    capsize=1.1,
                    alpha=0.96,
                )
        axis.set_yscale("log")
        axis.set_xlim(min(areas) - 0.3, max(areas) + 0.3)
        axis.set_xticks(list(areas))
        axis.tick_params(labelbottom=False)
        axis.grid(which="major", axis="y", color="0.90", linewidth=0.45)
        axis.set_title(condition_title(condition), pad=2.0)
        if column == 0:
            axis.set_ylabel(r"Mean distance ($c$ units)")
        _panel_label(axis, chr(ord("J") + column))

    ratio_low = float(ratios["ratio_ci_low"].min())
    ratio_high = float(ratios["ratio_ci_high"].max())
    padding = max(0.025, 0.08 * (ratio_high - ratio_low))
    common_low = min(1.0, ratio_low) - padding
    common_high = max(1.0, ratio_high) + padding
    bottom_axes: List[object] = []
    for column, condition in enumerate(SELECTED_CONDITIONS):
        axis = fig.add_subplot(bottom[0, column])
        bottom_axes.append(axis)
        condition_data = ratios[ratios["condition"] == condition]
        for shape in SHAPE_ORDER:
            style = SHAPE_STYLES[shape]
            group = condition_data[
                condition_data["shape"] == shape
            ].sort_values("projected_area")
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
                markersize=3.0,
                linewidth=0.9,
                elinewidth=0.55,
                capsize=1.3,
            )
        axis.axhline(1.0, color="0.50", linestyle=":", linewidth=0.75, zorder=0)
        axis.set_xlim(min(areas) - 0.3, max(areas) + 0.3)
        axis.set_xticks(list(areas))
        axis.set_ylim(common_low, common_high)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.grid(axis="y", color="0.90", linewidth=0.45)
        axis.set_title(condition_title(condition), pad=2.0)
        axis.set_xlabel(r"Projected area $A/c^2$")
        if column == 0:
            axis.set_ylabel("Fitted / empirical mean distance")
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
    raw_heading_y = max(
        axis.get_position().y1 for axis in raw_axes
    ) + 0.050
    fig.text(
        0.55,
        raw_heading_y,
        "Unbounded search: absolute mean distance (95% Monte Carlo intervals)",
        ha="center",
        va="center",
        fontsize=6.5,
    )
    ratio_heading_y = max(
        axis.get_position().y1 for axis in bottom_axes
    ) + 0.034
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


def _plot_compact_finite_figure(
    compact_finite: pd.DataFrame,
    figure_directory: Path,
) -> List[Path]:
    """Plot five finite-search experiments in a one-column grouped bar chart."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    fig, axis = plt.subplots(
        figsize=(ONE_COLUMN_WIDTH_INCHES, FINITE_COMPACT_HEIGHT_INCHES)
    )
    fig.subplots_adjust(left=0.175, right=0.995, bottom=0.19, top=0.79)
    centers = np.arange(len(SELECTED_CONDITIONS), dtype=float)
    bar_width = 0.145
    offsets = (
        np.arange(len(FINITE_PLOT_ORDER), dtype=float)
        - (len(FINITE_PLOT_ORDER) - 1.0) / 2.0
    ) * bar_width

    maximum = 100.0
    for class_index, class_name in enumerate(FINITE_PLOT_ORDER):
        style = COMPACT_CLASS_STYLES[class_name]
        class_data = (
            compact_finite[
                compact_finite["trajectory_class"] == class_name
            ]
            .set_index("condition")
            .loc[list(SELECTED_CONDITIONS)]
        )
        values = (
            100.0
            * class_data["normalized_detection_probability"].to_numpy(
                dtype=float
            )
        )
        low = 100.0 * class_data["ratio_ci_low"].to_numpy(dtype=float)
        high = 100.0 * class_data["ratio_ci_high"].to_numpy(dtype=float)
        maximum = max(maximum, float(np.nanmax(high)))
        axis.bar(
            centers + offsets[class_index],
            values,
            width=bar_width * 0.88,
            color=style["color"],
            edgecolor="black",
            linewidth=0.45,
            hatch=style["hatch"],
            zorder=2,
        )
        if class_name != "exact_global_rotation":
            axis.errorbar(
                centers + offsets[class_index],
                values,
                yerr=np.vstack((values - low, high - values)),
                fmt="none",
                ecolor="black",
                elinewidth=0.55,
                capsize=1.4,
                capthick=0.55,
                zorder=4,
            )

    axis.axhline(
        100.0,
        color="0.35",
        linestyle=(0, (2.0, 1.5)),
        linewidth=0.7,
        zorder=1,
    )
    axis.set_xlim(-0.48, len(SELECTED_CONDITIONS) - 0.52)
    axis.set_ylim(0.0, 1.12 * maximum)
    axis.set_xticks(centers)
    axis.set_xticklabels(
        [COMPACT_CONDITION_LABELS[value] for value in SELECTED_CONDITIONS]
    )
    axis.set_ylabel("Relative finite detection\nprobability (%)")
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axis.grid(axis="y", color="0.90", linewidth=0.45, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    handles = [
        Patch(
            facecolor=COMPACT_CLASS_STYLES[class_name]["color"],
            edgecolor="black",
            linewidth=0.45,
            hatch=COMPACT_CLASS_STYLES[class_name]["hatch"],
            label=str(COMPACT_CLASS_STYLES[class_name]["label"]),
        )
        for class_name in FINITE_PLOT_ORDER
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.585, 0.985),
        ncol=3,
        frameon=False,
        handlelength=1.2,
        handleheight=0.8,
        columnspacing=0.85,
        borderaxespad=0.0,
        fontsize=6.1,
    )
    paths = _save_figure(
        fig,
        figure_directory,
        "figure1_finite_search_compact",
    )
    plt.close(fig)
    return paths


def _plot_compact_unbounded_figure(
    compact_unbounded: pd.DataFrame,
    figure_directory: Path,
) -> List[Path]:
    """Plot empirical and fitted unbounded searches at one-column width."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    classes = ("vector_uniform_rotated", "fitted_mu")
    fig, axis = plt.subplots(
        figsize=(ONE_COLUMN_WIDTH_INCHES, UNBOUNDED_COMPACT_HEIGHT_INCHES)
    )
    fig.subplots_adjust(left=0.175, right=0.995, bottom=0.245, top=0.82)
    centers = np.arange(len(SELECTED_CONDITIONS), dtype=float)
    bar_width = 0.25
    offsets = (-0.14, 0.14)
    maximum = 100.0

    for class_index, class_name in enumerate(classes):
        style = COMPACT_CLASS_STYLES[class_name]
        class_data = (
            compact_unbounded[
                compact_unbounded["trajectory_class"] == class_name
            ]
            .set_index("condition")
            .loc[list(SELECTED_CONDITIONS)]
        )
        values = (
            100.0
            * class_data["normalized_detection_distance"].to_numpy(dtype=float)
        )
        low = 100.0 * class_data["ratio_ci_low"].to_numpy(dtype=float)
        high = 100.0 * class_data["ratio_ci_high"].to_numpy(dtype=float)
        maximum = max(maximum, float(np.nanmax(high)))
        positions = centers + offsets[class_index]
        axis.bar(
            positions,
            values,
            width=bar_width,
            color=style["color"],
            edgecolor="black",
            linewidth=0.45,
            hatch=style["hatch"],
            zorder=2,
        )
        if class_name != "vector_uniform_rotated":
            axis.errorbar(
                positions,
                values,
                yerr=np.vstack((values - low, high - values)),
                fmt="none",
                ecolor="black",
                elinewidth=0.6,
                capsize=1.5,
                capthick=0.6,
                zorder=4,
            )
            for x_value, y_value, high_value in zip(
                positions, values, high
            ):
                axis.text(
                    x_value,
                    max(y_value, high_value) + 1.5,
                    "{:.1f}".format(y_value),
                    ha="center",
                    va="bottom",
                    fontsize=6.0,
                )

    axis.axhline(
        100.0,
        color="0.35",
        linestyle=(0, (2.0, 1.5)),
        linewidth=0.7,
        zorder=1,
    )
    axis.set_xlim(-0.48, len(SELECTED_CONDITIONS) - 0.52)
    axis.set_ylim(0.0, max(112.0, 1.09 * maximum))
    axis.set_xticks(centers)
    axis.set_xticklabels(
        [COMPACT_CONDITION_LABELS[value] for value in SELECTED_CONDITIONS]
    )
    axis.set_ylabel("Relative mean detection\ndistance (%)")
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axis.grid(axis="y", color="0.90", linewidth=0.45, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    handles = [
        Patch(
            facecolor=COMPACT_CLASS_STYLES[class_name]["color"],
            edgecolor="black",
            linewidth=0.45,
            hatch=COMPACT_CLASS_STYLES[class_name]["hatch"],
            label=str(COMPACT_CLASS_STYLES[class_name]["label"]),
        )
        for class_name in classes
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.585, 0.98),
        ncol=2,
        frameon=False,
        handlelength=1.2,
        handleheight=0.8,
        columnspacing=1.2,
        borderaxespad=0.0,
        fontsize=6.2,
    )
    paths = _save_figure(
        fig,
        figure_directory,
        "figure2_unbounded_search_compact",
    )
    plt.close(fig)
    return paths


def _source_entry(store: ArtifactStore, stage: str, name: str) -> Dict[str, object]:
    path = store.frame_path(stage, name)
    sidecar = path.with_suffix(".meta.json")
    with sidecar.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    return {
        "stage": stage,
        "artifact": name,
        "path": str(path.resolve()),
        "sha256": metadata["sha256"],
        "rows": metadata["rows"],
    }


def run_selected_publication_figures(
    config: Config,
    baseline_stage: str = BASELINE_STAGE,
    augmented_stage: str = AUGMENTED_STAGE,
    output_stage: str = OUTPUT_STAGE,
) -> Dict[str, object]:
    """Validate selected replay results and build the two publication figures."""
    if not baseline_stage.strip() or not augmented_stage.strip() or not output_stage.strip():
        raise ValueError("stage names must be non-empty")
    if output_stage in {baseline_stage, augmented_stage}:
        raise ValueError("output_stage must differ from both source stages")

    store = ArtifactStore(config.output.root)
    source_names = (
        (baseline_stage, "fit_summary"),
        (baseline_stage, "empirical_pool_inputs"),
        (baseline_stage, "track_inputs"),
        (baseline_stage, "finite_track_summaries"),
        (augmented_stage, "finite_track_summaries"),
        (augmented_stage, "augmented_finite_summary"),
        (augmented_stage, "augmented_finite_contrasts"),
        (baseline_stage, "unbounded_summary"),
        (baseline_stage, "unbounded_contrasts"),
    )
    source_entries = [
        _source_entry(store, stage, name) for stage, name in source_names
    ]
    fits = _select_conditions(
        store.read_frame(baseline_stage, "fit_summary"),
        "fit_summary",
    )
    empirical_runs = _select_conditions(
        store.read_frame(baseline_stage, "empirical_pool_inputs"),
        "empirical_pool_inputs",
    )
    finite = _select_conditions(
        store.read_frame(augmented_stage, "augmented_finite_summary"),
        "augmented_finite_summary",
    )
    finite_contrasts = _select_conditions(
        store.read_frame(augmented_stage, "augmented_finite_contrasts"),
        "augmented_finite_contrasts",
    )
    unbounded = _select_conditions(
        store.read_frame(baseline_stage, "unbounded_summary"),
        "unbounded_summary",
    )
    unbounded_contrasts = _select_conditions(
        store.read_frame(baseline_stage, "unbounded_contrasts"),
        "unbounded_contrasts",
    )
    track_inputs = _select_conditions(
        store.read_frame(baseline_stage, "track_inputs"),
        "track_inputs",
    )
    baseline_finite_tracks = _select_conditions(
        store.read_frame(baseline_stage, "finite_track_summaries"),
        "baseline finite_track_summaries",
    )
    augmented_finite_tracks = _select_conditions(
        store.read_frame(augmented_stage, "finite_track_summaries"),
        "augmented finite_track_summaries",
    )
    finite_tracks = _condition_order(
        pd.concat(
            [baseline_finite_tracks, augmented_finite_tracks],
            ignore_index=True,
        )
    )
    areas = _validate_selected_data(fits, empirical_runs, finite, unbounded)
    _validate_selected_auxiliary_data(
        fits,
        finite_contrasts,
        unbounded_contrasts,
        track_inputs,
        areas,
    )
    _validate_selected_finite_tracks(fits, finite_tracks, areas)
    ccdf = _build_ccdf_table(fits, empirical_runs)
    ratios = _unbounded_ratio_table(unbounded)
    compact_finite = _finite_compact_ratio_table(
        finite_tracks,
        bootstrap_replicates=FINITE_COMPACT_BOOTSTRAP_REPLICATES,
        bootstrap_seed=(
            int(config.replay.random_seed)
            + FINITE_COMPACT_BOOTSTRAP_SEED_OFFSET
        ),
    )
    compact_unbounded = _unbounded_compact_ratio_table(unbounded)
    ratio_check = ratios.merge(
        unbounded_contrasts[
            [
                "condition",
                "shape",
                "projected_area",
                "fitted_to_empirical_distance_ratio",
                "ratio_ci_low",
                "ratio_ci_high",
            ]
        ],
        on=["condition", "shape", "projected_area"],
        suffixes=("_derived", "_stored"),
        validate="one_to_one",
    )
    for derived, stored in (
        ("fitted_to_empirical_ratio", "fitted_to_empirical_distance_ratio"),
        ("ratio_ci_low_derived", "ratio_ci_low_stored"),
        ("ratio_ci_high_derived", "ratio_ci_high_stored"),
    ):
        if not np.allclose(
            ratio_check[derived].to_numpy(dtype=float),
            ratio_check[stored].to_numpy(dtype=float),
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                "derived unbounded ratios do not match signed contrasts"
            )

    common_metadata = {
        "config_hash": config.hash(),
        "selection_version": SELECTION_VERSION,
        "selected_conditions": list(SELECTED_CONDITIONS),
        "source_stages": {
            "baseline": baseline_stage,
            "augmented": augmented_stage,
        },
    }
    output_paths: List[Path] = []
    for name, frame, source in (
        ("selected_fit_summary", fits, "{}:fit_summary".format(baseline_stage)),
        (
            "selected_empirical_runs",
            empirical_runs,
            "{}:empirical_pool_inputs".format(baseline_stage),
        ),
        (
            "selected_empirical_ccdf",
            ccdf,
            "derived: non-strict empirical survival",
        ),
        (
            "selected_finite_summary",
            finite,
            "{}:augmented_finite_summary".format(augmented_stage),
        ),
        (
            "selected_finite_contrasts",
            finite_contrasts,
            "{}:augmented_finite_contrasts".format(augmented_stage),
        ),
        (
            "selected_unbounded_summary",
            unbounded,
            "{}:unbounded_summary".format(baseline_stage),
        ),
        (
            "selected_unbounded_contrasts",
            unbounded_contrasts,
            "{}:unbounded_contrasts".format(baseline_stage),
        ),
        (
            "selected_track_inputs",
            track_inputs,
            "{}:track_inputs".format(baseline_stage),
        ),
        (
            "selected_finite_track_summaries",
            finite_tracks,
            (
                "{}:finite_track_summaries + "
                "{}:finite_track_summaries"
            ).format(baseline_stage, augmented_stage),
        ),
        (
            "selected_unbounded_ratios",
            ratios,
            "derived: fitted_mu / vector_uniform_rotated",
        ),
        (
            "selected_compact_finite_ratios",
            compact_finite,
            (
                "derived: equal target-cell means, equal track means, "
                "normalized to exact_global_rotation"
            ),
        ),
        (
            "selected_compact_unbounded_ratios",
            compact_unbounded,
            (
                "derived: equal target-cell mean distances, normalized to "
                "vector_uniform_rotated"
            ),
        ),
    ):
        metadata = dict(common_metadata)
        metadata["source"] = source
        output_paths.append(
            store.write_frame(
                output_stage,
                name,
                frame,
                extra_metadata=metadata,
            )
        )

    _configure_matplotlib()
    figure_directory = store.stage_dir(output_stage) / "figures"
    figure_paths = _plot_compact_finite_figure(
        compact_finite,
        figure_directory,
    )
    figure_paths.extend(
        _plot_compact_unbounded_figure(
            compact_unbounded,
            figure_directory,
        )
    )
    # Retain the larger diagnostic panels as supplementary material.
    figure_paths.extend(_plot_figure1(fits, empirical_runs, figure_directory))
    figure_paths.extend(
        _plot_figure2(
            fits,
            finite,
            unbounded,
            ratios,
            areas,
            figure_directory,
        )
    )
    output_paths.extend(figure_paths)

    figure_dimensions = {
        "figure1_finite_search_compact": (
            ONE_COLUMN_WIDTH_INCHES,
            FINITE_COMPACT_HEIGHT_INCHES,
        ),
        "figure2_unbounded_search_compact": (
            ONE_COLUMN_WIDTH_INCHES,
            UNBOUNDED_COMPACT_HEIGHT_INCHES,
        ),
        "figure1_data_and_fit": (7.01, 2.8),
        "figure2_search_performance": (
            7.01,
            DETAILED_SEARCH_HEIGHT_INCHES,
        ),
    }
    figure_manifest = [
        {
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
            "format": path.suffix.lstrip("."),
            "stem": path.stem,
            "width_inches": figure_dimensions[path.stem][0],
            "height_inches": figure_dimensions[path.stem][1],
            "role": (
                "main_one_column"
                if path.stem
                in {
                    "figure1_finite_search_compact",
                    "figure2_unbounded_search_compact",
                }
                else "supplementary_diagnostic"
            ),
        }
        for path in figure_paths
    ]
    table_manifest = [
        {
            "artifact": path.stem,
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
            "rows": int(
                json.loads(path.with_suffix(".meta.json").read_text())["rows"]
            ),
        }
        for path in output_paths
        if path.suffix == ".csv"
    ]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_version": SELECTION_VERSION,
        "config_hash": config.hash(),
        "source_artifacts": source_entries,
        "selected_conditions": [
            {
                "condition": condition,
                "condition_label": label,
                "order": index + 1,
            }
            for index, (condition, label) in enumerate(
                zip(SELECTED_CONDITIONS, SELECTED_LABELS)
            )
        ],
        "selection_rationale": (
            "Primary NK example plus the two remaining eligible conditions "
            "with the largest complete-run counts."
        ),
        "normalization": "each condition divided by its fitted crossover c",
        "segmentation": (
            "a turning angle greater than 90 degrees starts a new directed "
            "run; first and last boundary runs are excluded from fitting"
        ),
        "empirical_ccdf_definition": (
            "Pr[L >= sorted_L[i]] = (N - i) / N for zero-based i"
        ),
        "finite_classes": list(FINITE_PLOT_ORDER),
        "target_shapes": list(SHAPE_ORDER),
        "projected_areas": list(areas),
        "finite_interval": (
            "95% whole-track bootstrap interval from augmented finite summary"
        ),
        "compact_finite_estimand": (
            "equal mean over target cells within each track, then equal mean "
            "over tracks; normalized to whole-track rotation"
        ),
        "compact_finite_interval": (
            "{}-replicate paired whole-track percentile bootstrap, 95%"
        ).format(FINITE_COMPACT_BOOTSTRAP_REPLICATES),
        "unbounded_ratio_interval": (
            "95% independent-sample log-ratio delta-method Monte Carlo interval"
        ),
        "compact_unbounded_estimand": (
            "equal mean over target-cell MC means; normalized to empirical "
            "run-pool search"
        ),
        "tables": table_manifest,
        "figures": figure_manifest,
    }
    manifest_path = store.write_json(output_stage, "manifest", manifest)
    output_paths.append(manifest_path)
    summary: Dict[str, object] = {
        "conditions": len(SELECTED_CONDITIONS),
        "tracks": int(fits["n_tracks"].sum()),
        "complete_runs": int(fits["n_runs"].sum()),
        "finite_classes": len(FINITE_PLOT_ORDER),
        "finite_cells": int(len(finite)),
        "finite_contrast_cells": int(len(finite_contrasts)),
        "unbounded_cells": int(len(unbounded)),
        "unbounded_contrast_cells": int(len(unbounded_contrasts)),
        "unbounded_ratio_cells": int(len(ratios)),
        "compact_finite_rows": int(len(compact_finite)),
        "compact_unbounded_rows": int(len(compact_unbounded)),
        "selected_table_artifacts": [entry["artifact"] for entry in table_manifest],
        "figures": [str(path.resolve()) for path in figure_paths],
        "output_stage": output_stage,
        "validated": True,
    }
    summary_path = store.write_json(output_stage, "summary", summary)
    output_paths.append(summary_path)
    logger.info(
        "Publication figure summary: %s",
        json.dumps(summary, sort_keys=True),
    )
    return summary


__all__ = [
    "AUGMENTED_STAGE",
    "BASELINE_STAGE",
    "OUTPUT_STAGE",
    "SELECTED_CONDITIONS",
    "SELECTED_LABELS",
    "_empirical_ccdf",
    "_finite_compact_ratio_table",
    "_unbounded_compact_ratio_table",
    "_unbounded_ratio_table",
    "_validate_selected_data",
    "_validate_selected_finite_tracks",
    "run_selected_publication_figures",
]
