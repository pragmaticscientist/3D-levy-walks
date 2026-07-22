"""Empirical relocation pools with track-aware weights and audit diagnostics."""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .logging_utils import get_logger

logger = get_logger(__name__)

POOL_COLUMNS = (
    "condition",
    "condition_source",
    "file",
    "video_id",
    "track_uid",
    "fragment_uid",
    "run_uid",
    "run_length_um",
    "run_duration_s",
    "run_vector_x",
    "run_vector_y",
    "run_vector_z",
    "cell_balanced_weight",
    "step_weighted_weight",
)


def _eligible_run_mask(runs: pd.DataFrame, complete_only: bool) -> pd.Series:
    if "included" not in runs.columns:
        return pd.Series(True, index=runs.index)
    if complete_only:
        return runs["included"].astype(bool)
    reason = runs.get("exclusion_reason", pd.Series("", index=runs.index)).fillna("")
    return runs["included"].astype(bool) | (reason == "first_or_last_run_censored")


def attach_pool_weights(runs: pd.DataFrame) -> pd.DataFrame:
    """Attach primary cell-balanced and sensitivity step-weighted probabilities.

    Cell-balanced sampling first chooses a track uniformly and then a run
    uniformly within that track. Thus every track contributes total mass
    ``1 / n_tracks`` regardless of its duration. Step-weighted sampling chooses
    every run uniformly.
    """
    if runs.empty:
        out = runs.copy()
        out["cell_balanced_weight"] = pd.Series(dtype=float)
        out["step_weighted_weight"] = pd.Series(dtype=float)
        return out
    required = {"condition", "track_uid", "run_length_um"}
    missing = sorted(required.difference(runs.columns))
    if missing:
        raise ValueError("Cannot weight pools; missing columns: {}".format(", ".join(missing)))

    pieces: List[pd.DataFrame] = []
    for _, grp in runs.groupby("condition", sort=True):
        grp = grp.copy()
        n_tracks = int(grp["track_uid"].nunique())
        per_track = grp.groupby("track_uid")["run_uid"].transform("size").astype(float)
        grp["cell_balanced_weight"] = 1.0 / (float(n_tracks) * per_track)
        grp["step_weighted_weight"] = 1.0 / float(len(grp))
        pieces.append(grp)
    return pd.concat(pieces, ignore_index=True)


def build_pools(
    runs: pd.DataFrame,
    minimum_tracks: int,
    minimum_runs: int,
    complete_only: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build nested per-condition pools and a full eligibility audit table."""
    selected = runs.loc[_eligible_run_mask(runs, complete_only)].copy()
    if "condition" not in selected.columns:
        raise ValueError("Runs must be assigned a condition before pools are built")
    selected = selected[selected["condition"].fillna("").astype(str) != ""]
    selected = selected[
        np.isfinite(selected["run_length_um"].to_numpy(dtype=float))
        & (selected["run_length_um"].to_numpy(dtype=float) >= 0.0)
    ]
    weighted = attach_pool_weights(selected)

    audit_rows: List[Dict[str, object]] = []
    all_conditions = sorted(runs["condition"].dropna().astype(str).unique())
    for condition in all_conditions:
        grp = weighted[weighted["condition"].astype(str) == condition]
        n_runs = int(len(grp))
        n_tracks = int(grp["track_uid"].nunique()) if n_runs else 0
        reasons: List[str] = []
        if n_tracks < minimum_tracks:
            reasons.append("fewer than {} tracks".format(minimum_tracks))
        if n_runs < minimum_runs:
            reasons.append("fewer than {} runs".format(minimum_runs))
        lengths = grp["run_length_um"].to_numpy(dtype=float)
        weights = (
            grp["cell_balanced_weight"].to_numpy(dtype=float)
            if n_runs
            else np.array([], dtype=float)
        )
        track_mass = (
            grp.groupby("track_uid")["cell_balanced_weight"].sum().to_numpy(dtype=float)
            if n_runs
            else np.array([], dtype=float)
        )
        audit_rows.append(
            {
                "condition": condition,
                "condition_source": (
                    str(grp["condition_source"].iloc[0])
                    if n_runs and "condition_source" in grp.columns
                    else ""
                ),
                "n_tracks": n_tracks,
                "n_runs": n_runs,
                "n_videos": int(grp["video_id"].nunique()) if n_runs else 0,
                "length_min": float(np.min(lengths)) if n_runs else np.nan,
                "length_q25": float(np.quantile(lengths, 0.25)) if n_runs else np.nan,
                "length_median": float(np.median(lengths)) if n_runs else np.nan,
                "length_q75": float(np.quantile(lengths, 0.75)) if n_runs else np.nan,
                "length_q95": float(np.quantile(lengths, 0.95)) if n_runs else np.nan,
                "length_max": float(np.max(lengths)) if n_runs else np.nan,
                "effective_sample_size": (
                    float(1.0 / np.sum(np.square(weights))) if n_runs else 0.0
                ),
                "maximum_track_weight_share": (
                    float(np.max(track_mass)) if track_mass.size else np.nan
                ),
                "eligible": not reasons,
                "ineligible_reason": "; ".join(reasons),
            }
        )
    audit = pd.DataFrame(audit_rows)
    if not audit.empty:
        audit = audit.sort_values(
            ["eligible", "n_runs"], ascending=[False, False], kind="mergesort"
        ).reset_index(drop=True)
    eligible = set(audit.loc[audit["eligible"], "condition"]) if not audit.empty else set()
    weighted["pool_eligible"] = weighted["condition"].isin(eligible)
    logger.info(
        "Built %d relocation pools; %d meet the configured minimums",
        len(audit),
        int(audit["eligible"].sum()) if not audit.empty else 0,
    )
    return weighted, audit


def sample_pool(
    pool: pd.DataFrame,
    size: int,
    rng: np.random.Generator,
    weighting: str = "cell_balanced",
) -> pd.DataFrame:
    """Draw pool rows with replacement using the requested weighting."""
    if size < 0:
        raise ValueError("size must be non-negative")
    if pool.empty and size:
        raise ValueError("Cannot sample from an empty pool")
    column = "{}_weight".format(weighting)
    if column not in pool.columns:
        raise ValueError("Unknown pool weighting {!r}".format(weighting))
    probabilities = pool[column].to_numpy(dtype=float)
    probabilities = probabilities / probabilities.sum()
    indices = rng.choice(len(pool), size=size, replace=True, p=probabilities)
    return pool.iloc[indices].reset_index(drop=True)

