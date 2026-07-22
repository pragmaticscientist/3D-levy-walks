from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ltdb_levy.classes import (
    block_resampled,
    fitted_isotropic,
    ordered_length_rotated,
    pool_isotropic,
)
from ltdb_levy.outcomes import (
    combine_trial_summaries,
    shape_sensitivity,
    summarise_trials,
)


def test_ordered_rotation_preserves_every_length_in_order():
    vectors = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    rotated = ordered_length_rotated(vectors, np.random.default_rng(77))
    assert np.allclose(np.linalg.norm(rotated, axis=1), [1.0, 2.0, 3.0])
    assert not np.allclose(rotated, vectors)


def test_pool_and_fitted_surrogates_cover_budget_and_respect_support():
    rng = np.random.default_rng(78)
    support = np.array([0.5, 1.5, 4.0])
    vectors = pool_isotropic(
        support, np.array([0.2, 0.3, 0.5]), budget=30.0, rng=rng
    )
    lengths = np.linalg.norm(vectors, axis=1)
    assert lengths.sum() >= 30.0
    assert set(np.round(lengths, 12)).issubset(set(support))

    fitted = fitted_isotropic(2.0, 12.0, budget=30.0, rng=rng)
    fitted_lengths = np.linalg.norm(fitted, axis=1)
    assert fitted_lengths.sum() >= 30.0
    assert np.all((fitted_lengths >= 0.0) & (fitted_lengths <= 12.0))


def test_block_resampling_keeps_within_block_vectors_under_one_rotation():
    sequence = np.array(
        [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]]
    )
    result = block_resampled(
        [sequence], block_length=2, budget=20.0, rng=np.random.default_rng(79)
    )
    for start in range(0, len(result), 2):
        block = result[start : start + 2]
        if len(block) == 2:
            # Every candidate pair in this fixture is collinear with a known
            # ratio except the middle [2x, y] block; rotations preserve dot angle.
            candidate_angles = []
            for original_start in range(3):
                original = sequence[original_start : original_start + 2]
                cosine = np.dot(original[0], original[1]) / np.prod(
                    np.linalg.norm(original, axis=1)
                )
                candidate_angles.append(cosine)
            observed = np.dot(block[0], block[1]) / np.prod(
                np.linalg.norm(block, axis=1)
            )
            assert any(observed == pytest.approx(value) for value in candidate_angles)


def test_outcomes_use_budget_for_censored_trials_and_compute_shape_G():
    trials = pd.DataFrame(
        {
            "class": ["empirical"] * 4 + ["empirical"] * 4,
            "shape": ["Ball"] * 4 + ["Line"] * 4,
            "detected": [True, False, True, False, True, True, False, False],
            "detection_distance": [2.0, np.nan, 6.0, np.nan, 1.0, 3.0, np.nan, np.nan],
            "budget": [10.0] * 8,
        }
    )
    summary = summarise_trials(trials, ["class", "shape"])
    ball = summary[summary["shape"] == "Ball"].iloc[0]
    line = summary[summary["shape"] == "Line"].iloc[0]
    assert ball["restricted_mean_detection_distance"] == pytest.approx(7.0)
    assert line["restricted_mean_detection_distance"] == pytest.approx(6.0)
    sensitivity = shape_sensitivity(summary, ["class"])
    assert sensitivity.iloc[0]["shape_sensitivity_G"] == pytest.approx(7.0 / 6.0)
    assert sensitivity.iloc[0]["best_shape"] == "Line"


def test_disjoint_trial_summaries_combine_exactly():
    trials = pd.DataFrame(
        {
            "cell": ["x"] * 6,
            "detected": [True, False, False, True, True, False],
            "detection_distance": [1.0, np.nan, np.nan, 3.0, 5.0, np.nan],
            "budget": [10.0] * 6,
        }
    )
    direct = summarise_trials(trials, ["cell"])
    chunks = pd.concat(
        [
            summarise_trials(trials.iloc[:2], ["cell"]),
            summarise_trials(trials.iloc[2:], ["cell"]),
        ],
        ignore_index=True,
    )
    combined = combine_trial_summaries(chunks, ["cell"])
    for column in (
        "detection_probability",
        "restricted_mean_detection_distance",
        "conditional_mean_detection_distance",
        "restricted_mean_detection_distance_mc_se",
    ):
        assert combined.iloc[0][column] == pytest.approx(direct.iloc[0][column])
