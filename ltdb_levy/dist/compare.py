"""Track-held-out grouped cross-validation for relocation-length models."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..errors import FitError
from .families import FAMILY_NAMES, fit_family, log_density as family_log_density
from .mixture import log_density as mixture_log_density
from .mle import fit_mixture

MODEL_NAMES = ("c_mixture",) + FAMILY_NAMES


def _comparison_weights(frame: pd.DataFrame, weighting: str) -> np.ndarray:
    if weighting == "cell_balanced":
        n_tracks = frame["track_uid"].nunique()
        per_track = (
            frame.groupby("track_uid")["run_uid"]
            .transform("size")
            .to_numpy(dtype=float)
        )
        return 1.0 / (float(n_tracks) * per_track)
    if weighting == "step_weighted":
        return np.full(len(frame), 1.0 / len(frame), dtype=float)
    raise ValueError("Unknown weighting {!r}".format(weighting))


def _equal_run_oof_summary(
    score_pieces: Sequence[np.ndarray],
    track_pieces: Sequence[np.ndarray],
) -> tuple[float, float]:
    """Pool out-of-fold run scores and return a track-clustered mean SE.

    The point estimand is the mean log density of a randomly selected run, so
    every held-out run receives mass ``1 / N`` even when grouped folds contain
    different numbers of runs.  The standard error treats tracks as the
    independent sampling clusters.  It is the finite-cluster-corrected
    sandwich standard error for the run-weighted ratio estimator:

        sqrt(J / (J - 1) * sum_j S_j**2) / N,

    where ``S_j = sum_i(y_ij - y_bar)`` within track ``j``.
    """
    if not score_pieces:
        return float("-inf"), np.nan

    scores = np.concatenate(score_pieces).astype(float, copy=False)
    tracks = np.concatenate(track_pieces)
    mean_score = float(np.mean(scores))
    unique_tracks = pd.unique(tracks)
    if len(unique_tracks) < 2:
        return mean_score, np.nan

    residuals = scores - mean_score
    cluster_totals = (
        pd.DataFrame({"track_uid": tracks, "residual": residuals})
        .groupby("track_uid", sort=False, dropna=False)["residual"]
        .sum()
        .to_numpy(dtype=float)
    )
    correction = float(len(cluster_totals) / (len(cluster_totals) - 1))
    standard_error = float(
        np.sqrt(correction * np.dot(cluster_totals, cluster_totals)) / len(scores)
    )
    return mean_score, standard_error


def grouped_model_comparison(
    runs: pd.DataFrame,
    lmax: float,
    crossover: float,
    mu_bounds: Sequence[float],
    n_folds: int,
    rng: np.random.Generator,
    weighting: str = "cell_balanced",
    estimate_crossover: bool = False,
    crossover_bounds: Optional[Sequence[float]] = None,
) -> pd.DataFrame:
    """Compare models with folds that hold out complete tracked cells."""
    if weighting not in ("cell_balanced", "step_weighted"):
        raise ValueError("Unknown weighting {!r}".format(weighting))
    tracks = runs["track_uid"].drop_duplicates().to_numpy()
    if len(tracks) < 2:
        raise FitError("Grouped cross-validation requires at least two tracks")
    folds = min(max(2, int(n_folds)), len(tracks))
    shuffled = tracks.copy()
    rng.shuffle(shuffled)
    assignments: Dict[object, int] = {
        track: index % folds for index, track in enumerate(shuffled)
    }
    fold_scores: Dict[str, List[float]] = {name: [] for name in MODEL_NAMES}
    oof_score_pieces: Dict[str, List[np.ndarray]] = {
        name: [] for name in MODEL_NAMES
    }
    oof_track_pieces: Dict[str, List[np.ndarray]] = {
        name: [] for name in MODEL_NAMES
    }

    for fold in range(folds):
        test_mask = runs["track_uid"].map(assignments).to_numpy() == fold
        train = runs.loc[~test_mask]
        test = runs.loc[test_mask]
        train_values = train["run_length_um"].to_numpy(dtype=float)
        test_values = test["run_length_um"].to_numpy(dtype=float)
        test_tracks = test["track_uid"].to_numpy()
        train_weights = _comparison_weights(train, weighting)
        test_weights = _comparison_weights(test, weighting)
        test_weights = test_weights / test_weights.sum()

        fold_crossover = float(crossover)
        try:
            mixture = fit_mixture(
                train_values,
                lmax=lmax,
                crossover=crossover,
                mu_bounds=mu_bounds,
                weights=train_weights,
                estimate_crossover=estimate_crossover,
                crossover_bounds=crossover_bounds,
            )
            fold_crossover = mixture.crossover
            logs = mixture_log_density(
                test_values, mixture.mu, lmax, fold_crossover
            )
            fold_score = float(np.dot(test_weights, logs))
            fold_scores["c_mixture"].append(fold_score)
            if weighting == "step_weighted" and np.isfinite(fold_score):
                oof_score_pieces["c_mixture"].append(np.asarray(logs, dtype=float))
                oof_track_pieces["c_mixture"].append(test_tracks)
        except FitError:
            fold_scores["c_mixture"].append(float("-inf"))

        for name in FAMILY_NAMES:
            try:
                fitted = fit_family(
                    name,
                    train_values,
                    lmax=lmax,
                    weights=train_weights,
                    crossover=fold_crossover,
                )
                logs = family_log_density(
                    name,
                    test_values,
                    fitted.parameters,
                    lmax,
                    fold_crossover,
                )
                fold_score = float(np.dot(test_weights, logs))
                fold_scores[name].append(fold_score)
                if weighting == "step_weighted" and np.isfinite(fold_score):
                    oof_score_pieces[name].append(np.asarray(logs, dtype=float))
                    oof_track_pieces[name].append(test_tracks)
            except FitError:
                fold_scores[name].append(float("-inf"))

    rows: List[Dict[str, object]] = []
    for name in MODEL_NAMES:
        scores = np.asarray(fold_scores[name], dtype=float)
        finite = scores[np.isfinite(scores)]
        if weighting == "step_weighted":
            mean_score, standard_error = _equal_run_oof_summary(
                oof_score_pieces[name],
                oof_track_pieces[name],
            )
            score_aggregation = "pooled_out_of_fold_runs"
            se_method = "track_clustered_sandwich"
        else:
            mean_score = float(np.mean(finite)) if finite.size else -np.inf
            standard_error = (
                float(np.std(finite, ddof=1) / np.sqrt(len(finite)))
                if len(finite) > 1
                else np.nan
            )
            score_aggregation = "equal_fold_means"
            se_method = "standard_error_across_folds"
        rows.append(
            {
                "model": name,
                "mean_test_log_density": mean_score,
                "se_test_log_density": standard_error,
                "score_aggregation": score_aggregation,
                "se_method": se_method,
                "folds_successful": int(finite.size),
                "folds_total": int(folds),
            }
        )
    result = pd.DataFrame(rows)
    result["rank"] = (
        result["mean_test_log_density"].rank(method="min", ascending=False).astype(int)
    )
    return result.sort_values("rank", kind="mergesort").reset_index(drop=True)
