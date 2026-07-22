from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad

import ltdb_levy.dist.compare as compare_module
from ltdb_levy.dist.compare import MODEL_NAMES, grouped_model_comparison
from ltdb_levy.dist.families import FAMILY_NAMES, fit_family, log_density
from ltdb_levy.dist.mixture import sample


@pytest.mark.parametrize("name", FAMILY_NAMES)
def test_alternative_family_fits_normalized_density(name):
    rng = np.random.default_rng(102)
    values = sample(2.0, 20.0, 500, rng)
    fitted = fit_family(name, values, lmax=20.0)
    integral, _ = quad(
        lambda value: float(
            np.exp(
                log_density(
                    name,
                    np.array([value]),
                    fitted.parameters,
                    fitted.lmax,
                    fitted.crossover,
                )[0]
            )
        ),
        0.0,
        20.0,
        points=[1.0],
        limit=100,
    )
    assert integral == pytest.approx(1.0, abs=2e-6)


def test_grouped_comparison_returns_every_model_and_holds_out_tracks():
    rng = np.random.default_rng(103)
    pieces = []
    for track_index in range(6):
        lengths = sample(2.1, 20.0, 40, rng)
        pieces.append(
            pd.DataFrame(
                {
                    "track_uid": "track{}".format(track_index),
                    "run_uid": [
                        "track{}:{}".format(track_index, i) for i in range(len(lengths))
                    ],
                    "run_length_um": lengths,
                }
            )
        )
    runs = pd.concat(pieces, ignore_index=True)
    comparison = grouped_model_comparison(
        runs,
        lmax=20.0,
        crossover=1.0,
        mu_bounds=(1.001, 6.0),
        n_folds=3,
        rng=rng,
    )
    assert set(comparison["model"]) == set(MODEL_NAMES)
    assert comparison["folds_total"].eq(3).all()
    assert comparison["folds_successful"].eq(3).all()
    assert np.isfinite(comparison["mean_test_log_density"]).all()


def test_grouped_comparison_estimates_crossover_inside_training_folds():
    rng = np.random.default_rng(104)
    pieces = []
    for track_index in range(8):
        lengths = sample(2.4, 30.0, 60, rng, crossover=2.5)
        pieces.append(
            pd.DataFrame(
                {
                    "track_uid": "track{}".format(track_index),
                    "run_uid": [
                        "track{}:{}".format(track_index, i)
                        for i in range(len(lengths))
                    ],
                    "run_length_um": lengths,
                }
            )
        )
    comparison = grouped_model_comparison(
        pd.concat(pieces, ignore_index=True),
        lmax=30.0,
        crossover=3.75,
        mu_bounds=(1.001, 6.0),
        n_folds=4,
        rng=rng,
        estimate_crossover=True,
    )
    assert comparison["folds_successful"].eq(4).all()
    assert np.isfinite(comparison["mean_test_log_density"]).all()


def test_grouped_comparison_accepts_equal_run_weighting():
    rng = np.random.default_rng(105)
    pieces = []
    for track_index, count in enumerate((20, 35, 50, 65)):
        lengths = sample(2.2, 20.0, count, rng, crossover=2.0)
        pieces.append(
            pd.DataFrame(
                {
                    "track_uid": "track{}".format(track_index),
                    "run_uid": [
                        "track{}:{}".format(track_index, i)
                        for i in range(len(lengths))
                    ],
                    "run_length_um": lengths,
                }
            )
        )
    comparison = grouped_model_comparison(
        pd.concat(pieces, ignore_index=True),
        lmax=20.0,
        crossover=2.0,
        mu_bounds=(1.001, 6.0),
        n_folds=2,
        rng=rng,
        weighting="step_weighted",
        estimate_crossover=True,
    )
    assert set(comparison["model"]) == set(MODEL_NAMES)
    assert comparison["folds_successful"].eq(2).all()


def _patch_log_scores_to_equal_run_lengths(monkeypatch):
    monkeypatch.setattr(
        compare_module,
        "fit_mixture",
        lambda values, lmax, crossover, **kwargs: SimpleNamespace(
            mu=2.0,
            crossover=crossover,
        ),
    )
    monkeypatch.setattr(
        compare_module,
        "mixture_log_density",
        lambda values, *args, **kwargs: np.asarray(values, dtype=float),
    )
    monkeypatch.setattr(
        compare_module,
        "fit_family",
        lambda *args, **kwargs: SimpleNamespace(parameters={}),
    )
    monkeypatch.setattr(
        compare_module,
        "family_log_density",
        lambda name, values, *args, **kwargs: np.asarray(values, dtype=float),
    )


def _unequal_track_runs():
    pieces = []
    for track_uid, count, score in (
        ("short", 1, 1.0),
        ("medium", 2, 2.0),
        ("long", 9, 10.0),
    ):
        pieces.append(
            pd.DataFrame(
                {
                    "track_uid": track_uid,
                    "run_uid": [
                        "{}:{}".format(track_uid, index) for index in range(count)
                    ],
                    "run_length_um": np.full(count, score),
                }
            )
        )
    return pd.concat(pieces, ignore_index=True)


def test_equal_run_comparison_pools_runs_and_uses_track_clustered_se(monkeypatch):
    _patch_log_scores_to_equal_run_lengths(monkeypatch)
    runs = _unequal_track_runs()

    comparison = grouped_model_comparison(
        runs,
        lmax=20.0,
        crossover=2.0,
        mu_bounds=(1.001, 6.0),
        n_folds=3,
        rng=np.random.default_rng(106),
        weighting="step_weighted",
    )

    expected_mean = float(runs["run_length_um"].mean())
    track_residual_sums = (
        (runs["run_length_um"] - expected_mean)
        .groupby(runs["track_uid"])
        .sum()
        .to_numpy()
    )
    n_tracks = runs["track_uid"].nunique()
    expected_se = float(
        np.sqrt(
            n_tracks
            / (n_tracks - 1)
            * np.dot(track_residual_sums, track_residual_sums)
        )
        / len(runs)
    )
    mixture = comparison.set_index("model").loc["c_mixture"]
    assert mixture["mean_test_log_density"] == pytest.approx(expected_mean)
    assert mixture["mean_test_log_density"] != pytest.approx(
        np.mean([1.0, 2.0, 10.0])
    )
    assert mixture["se_test_log_density"] == pytest.approx(expected_se)
    assert mixture["score_aggregation"] == "pooled_out_of_fold_runs"
    assert mixture["se_method"] == "track_clustered_sandwich"


def test_cell_balanced_comparison_preserves_equal_fold_summary(monkeypatch):
    _patch_log_scores_to_equal_run_lengths(monkeypatch)
    runs = _unequal_track_runs()
    seed = 107
    tracks = runs["track_uid"].drop_duplicates().to_numpy()
    shuffled = tracks.copy()
    np.random.default_rng(seed).shuffle(shuffled)
    assignments = {
        track_uid: index % 2 for index, track_uid in enumerate(shuffled)
    }
    expected_fold_scores = []
    for fold in range(2):
        test = runs[runs["track_uid"].map(assignments) == fold]
        per_track_means = test.groupby("track_uid")["run_length_um"].mean()
        expected_fold_scores.append(float(per_track_means.mean()))

    comparison = grouped_model_comparison(
        runs,
        lmax=20.0,
        crossover=2.0,
        mu_bounds=(1.001, 6.0),
        n_folds=2,
        rng=np.random.default_rng(seed),
        weighting="cell_balanced",
    )

    mixture = comparison.set_index("model").loc["c_mixture"]
    assert mixture["mean_test_log_density"] == pytest.approx(
        np.mean(expected_fold_scores)
    )
    assert mixture["se_test_log_density"] == pytest.approx(
        np.std(expected_fold_scores, ddof=1) / np.sqrt(2)
    )
    assert mixture["score_aggregation"] == "equal_fold_means"
    assert mixture["se_method"] == "standard_error_across_folds"
