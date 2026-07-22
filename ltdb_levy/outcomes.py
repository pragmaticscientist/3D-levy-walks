"""Replay outcome summaries with censoring-aware primary estimands."""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
import pandas as pd


def summarise_trials(
    trials: pd.DataFrame,
    group_fields: Sequence[str],
) -> pd.DataFrame:
    """Summarise detection probability and restricted mean detection distance."""
    required = {"detected", "detection_distance", "budget"}
    missing = sorted(required.difference(trials.columns))
    if missing:
        raise ValueError("Trial table is missing columns: {}".format(", ".join(missing)))
    rows: List[Dict[str, object]] = []
    for keys, group in trials.groupby(list(group_fields), sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        detected = group["detected"].astype(bool).to_numpy()
        distance = group["detection_distance"].to_numpy(dtype=float)
        budget = group["budget"].to_numpy(dtype=float)
        restricted = np.where(detected, distance, budget)
        row: Dict[str, object] = dict(zip(group_fields, keys))
        row.update(
            {
                "n_trials": int(len(group)),
                "n_detected": int(detected.sum()),
                "detection_probability": float(detected.mean()),
                "detection_probability_mc_se": float(
                    np.sqrt(detected.mean() * (1.0 - detected.mean()) / len(group))
                ),
                "restricted_mean_detection_distance": float(np.mean(restricted)),
                "restricted_mean_detection_distance_mc_se": (
                    float(np.std(restricted, ddof=1) / np.sqrt(len(group)))
                    if len(group) > 1
                    else np.nan
                ),
                "conditional_mean_detection_distance": (
                    float(np.mean(distance[detected])) if detected.any() else np.nan
                ),
                "mean_budget": float(np.mean(budget)),
                "_restricted_sum": float(np.sum(restricted)),
                "_restricted_sum_squares": float(np.sum(np.square(restricted))),
                "_detected_distance_sum": float(np.nansum(distance[detected])),
                "_budget_sum": float(np.sum(budget)),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def combine_trial_summaries(
    summaries: pd.DataFrame,
    group_fields: Sequence[str],
) -> pd.DataFrame:
    """Combine disjoint trial summaries without reloading their source trials."""
    required = {
        "n_trials",
        "n_detected",
        "_restricted_sum",
        "_restricted_sum_squares",
        "_detected_distance_sum",
        "_budget_sum",
    }
    missing = sorted(required.difference(summaries.columns))
    if missing:
        raise ValueError(
            "Summaries cannot be combined; missing columns: {}".format(", ".join(missing))
        )
    rows: List[Dict[str, object]] = []
    for keys, group in summaries.groupby(list(group_fields), sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        n_trials = int(group["n_trials"].sum())
        n_detected = int(group["n_detected"].sum())
        restricted_sum = float(group["_restricted_sum"].sum())
        restricted_sum_squares = float(group["_restricted_sum_squares"].sum())
        rmdd = restricted_sum / n_trials
        variance = (
            max(0.0, (restricted_sum_squares - n_trials * rmdd * rmdd) / (n_trials - 1))
            if n_trials > 1
            else np.nan
        )
        probability = n_detected / n_trials
        row: Dict[str, object] = dict(zip(group_fields, keys))
        row.update(
            {
                "n_trials": n_trials,
                "n_detected": n_detected,
                "detection_probability": probability,
                "detection_probability_mc_se": float(
                    np.sqrt(probability * (1.0 - probability) / n_trials)
                ),
                "restricted_mean_detection_distance": rmdd,
                "restricted_mean_detection_distance_mc_se": (
                    float(np.sqrt(variance / n_trials)) if n_trials > 1 else np.nan
                ),
                "conditional_mean_detection_distance": (
                    float(group["_detected_distance_sum"].sum() / n_detected)
                    if n_detected
                    else np.nan
                ),
                "mean_budget": float(group["_budget_sum"].sum() / n_trials),
                "_restricted_sum": restricted_sum,
                "_restricted_sum_squares": restricted_sum_squares,
                "_detected_distance_sum": float(group["_detected_distance_sum"].sum()),
                "_budget_sum": float(group["_budget_sum"].sum()),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def shape_sensitivity(
    summaries: pd.DataFrame,
    group_fields: Sequence[str],
    shape_field: str = "shape",
) -> pd.DataFrame:
    """Compute ``G=max_shape(RMDD)/min_shape(RMDD)`` and its log form."""
    value = "restricted_mean_detection_distance"
    rows: List[Dict[str, object]] = []
    for keys, group in summaries.groupby(list(group_fields), sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = group[value].to_numpy(dtype=float)
        finite_positive = values[np.isfinite(values) & (values > 0.0)]
        row: Dict[str, object] = dict(zip(group_fields, keys))
        if len(finite_positive) and group[shape_field].nunique() >= 2:
            minimum = float(np.min(finite_positive))
            maximum = float(np.max(finite_positive))
            row["shape_sensitivity_G"] = maximum / minimum
            row["log_shape_sensitivity"] = float(np.log(maximum) - np.log(minimum))
            row["best_shape"] = str(group.loc[group[value].idxmin(), shape_field])
            row["worst_shape"] = str(group.loc[group[value].idxmax(), shape_field])
        else:
            row["shape_sensitivity_G"] = np.nan
            row["log_shape_sensitivity"] = np.nan
            row["best_shape"] = ""
            row["worst_shape"] = ""
        rows.append(row)
    return pd.DataFrame(rows)
