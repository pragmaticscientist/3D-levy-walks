"""Track-clustered biological uncertainty for replay outcomes."""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

ADJACENT_CLASS_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("empirical", "ordered_length_rotated"),
    ("ordered_length_rotated", "pool_isotropic"),
    ("pool_isotropic", "fitted_mu"),
    ("fitted_mu", "mu2"),
)


def _interval(values: np.ndarray) -> Tuple[float, float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.nan, np.nan, np.nan
    low, median, high = np.quantile(finite, [0.025, 0.5, 0.975])
    return float(low), float(median), float(high)


def clustered_outcome_intervals(
    track_summaries: pd.DataFrame,
    n_bootstraps: int,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Bootstrap whole tracks with paired class/shape outcomes kept together.

    Returns:
        ``(outcome_intervals, adjacent_class_contrasts, shape_intervals)``.
    """
    if n_bootstraps < 1:
        raise ValueError("n_bootstraps must be positive")
    cell_fields = ["trajectory_class", "shape", "projected_area"]
    outcome_rows: List[Dict[str, object]] = []
    contrast_rows: List[Dict[str, object]] = []
    shape_rows: List[Dict[str, object]] = []

    for condition, condition_frame in track_summaries.groupby("condition", sort=True):
        tracks = sorted(condition_frame["track_uid"].unique())
        if not tracks:
            continue
        cells = (
            condition_frame[cell_fields]
            .drop_duplicates()
            .sort_values(cell_fields, kind="mergesort")
            .reset_index(drop=True)
        )
        cell_keys = [tuple(row) for row in cells.to_numpy()]
        cell_index = {key: index for index, key in enumerate(cell_keys)}
        track_index = {track: index for index, track in enumerate(tracks)}
        detection = np.full((len(tracks), len(cells)), np.nan, dtype=float)
        rmdd = np.full_like(detection, np.nan)
        for _, row in condition_frame.iterrows():
            key = (
                row["trajectory_class"],
                row["shape"],
                row["projected_area"],
            )
            i = track_index[row["track_uid"]]
            j = cell_index[key]
            detection[i, j] = float(row["detection_probability"])
            rmdd[i, j] = float(row["restricted_mean_detection_distance"])
        if np.any(~np.isfinite(detection)) or np.any(~np.isfinite(rmdd)):
            raise ValueError(
                "{} has incomplete track-by-cell replay summaries; paired bootstrap "
                "requires every selected track in every cell".format(condition)
            )

        sampled = rng.integers(0, len(tracks), size=(n_bootstraps, len(tracks)))
        detection_draws = np.mean(detection[sampled, :], axis=1)
        rmdd_draws = np.mean(rmdd[sampled, :], axis=1)
        for j, key in enumerate(cell_keys):
            det_low, det_median, det_high = _interval(detection_draws[:, j])
            rmdd_low, rmdd_median, rmdd_high = _interval(rmdd_draws[:, j])
            outcome_rows.append(
                {
                    "condition": condition,
                    "trajectory_class": key[0],
                    "shape": key[1],
                    "projected_area": key[2],
                    "n_tracks": len(tracks),
                    "detection_probability_ci_low": det_low,
                    "detection_probability_bootstrap_median": det_median,
                    "detection_probability_ci_high": det_high,
                    "rmdd_ci_low": rmdd_low,
                    "rmdd_bootstrap_median": rmdd_median,
                    "rmdd_ci_high": rmdd_high,
                    "n_bootstraps": n_bootstraps,
                }
            )

        classes = set(cells["trajectory_class"])
        shapes = sorted(cells["shape"].unique())
        areas = sorted(cells["projected_area"].unique())
        for left, right in ADJACENT_CLASS_PAIRS:
            if left not in classes or right not in classes:
                continue
            for shape in shapes:
                for area in areas:
                    left_index = cell_index.get((left, shape, area))
                    right_index = cell_index.get((right, shape, area))
                    if left_index is None or right_index is None:
                        continue
                    detection_difference = (
                        detection_draws[:, left_index] - detection_draws[:, right_index]
                    )
                    log_rmdd_difference = (
                        np.log(rmdd_draws[:, left_index])
                        - np.log(rmdd_draws[:, right_index])
                    )
                    det_low, det_median, det_high = _interval(detection_difference)
                    log_low, log_median, log_high = _interval(log_rmdd_difference)
                    contrast_rows.append(
                        {
                            "condition": condition,
                            "left_class": left,
                            "right_class": right,
                            "shape": shape,
                            "projected_area": area,
                            "detection_probability_difference_ci_low": det_low,
                            "detection_probability_difference_median": det_median,
                            "detection_probability_difference_ci_high": det_high,
                            "log_rmdd_difference_ci_low": log_low,
                            "log_rmdd_difference_median": log_median,
                            "log_rmdd_difference_ci_high": log_high,
                            "n_bootstraps": n_bootstraps,
                        }
                    )

        for trajectory_class in sorted(classes):
            for area in areas:
                indices = [
                    cell_index[(trajectory_class, shape, area)]
                    for shape in shapes
                    if (trajectory_class, shape, area) in cell_index
                ]
                if len(indices) < 2:
                    continue
                selected = rmdd_draws[:, indices]
                sensitivity = np.max(selected, axis=1) / np.min(selected, axis=1)
                low, median, high = _interval(sensitivity)
                shape_rows.append(
                    {
                        "condition": condition,
                        "trajectory_class": trajectory_class,
                        "projected_area": area,
                        "shape_sensitivity_ci_low": low,
                        "shape_sensitivity_bootstrap_median": median,
                        "shape_sensitivity_ci_high": high,
                        "n_bootstraps": n_bootstraps,
                    }
                )

    return (
        pd.DataFrame(outcome_rows),
        pd.DataFrame(contrast_rows),
        pd.DataFrame(shape_rows),
    )

