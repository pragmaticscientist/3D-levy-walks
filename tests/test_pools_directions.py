from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ltdb_levy.directions import (
    independently_rotate_vectors,
    isotropic_directions,
    random_quaternions,
    rotate_vectors,
)
from ltdb_levy.pools import attach_pool_weights, sample_pool


def _pool():
    return pd.DataFrame(
        {
            "condition": ["c"] * 6,
            "track_uid": ["a", "a", "b", "b", "b", "b"],
            "run_uid": ["a0", "a1", "b0", "b1", "b2", "b3"],
            "run_length_um": np.arange(1.0, 7.0),
        }
    )


def test_cell_balanced_weights_give_every_track_equal_mass():
    weighted = attach_pool_weights(_pool())
    masses = weighted.groupby("track_uid")["cell_balanced_weight"].sum()
    assert masses["a"] == pytest.approx(0.5)
    assert masses["b"] == pytest.approx(0.5)
    assert np.allclose(weighted["step_weighted_weight"], 1.0 / 6.0)

    sampled = sample_pool(
        weighted, 30_000, np.random.default_rng(7), weighting="cell_balanced"
    )
    fraction_a = (sampled["track_uid"] == "a").mean()
    assert fraction_a == pytest.approx(0.5, abs=0.015)


def test_isotropic_directions_have_unit_norm_and_uniform_cosine():
    directions = isotropic_directions(60_000, np.random.default_rng(8))
    assert np.max(np.abs(np.linalg.norm(directions, axis=1) - 1.0)) < 1e-12
    assert directions[:, 2].mean() == pytest.approx(0.0, abs=0.01)
    assert np.mean(directions[:, 2] ** 2) == pytest.approx(1.0 / 3.0, abs=0.01)


def test_random_rotation_preserves_lengths_and_one_rotation_preserves_block():
    rng = np.random.default_rng(9)
    vectors = np.array([[2.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
    quaternion = random_quaternions(1, rng)
    rotated_block = rotate_vectors(vectors, quaternion)
    assert np.allclose(np.linalg.norm(rotated_block, axis=1), [2.0, 4.0])
    assert np.allclose(rotated_block[1], 2.0 * rotated_block[0])

    independent = independently_rotate_vectors(vectors, rng)
    assert np.allclose(np.linalg.norm(independent, axis=1), [2.0, 4.0])
    assert not np.allclose(independent[1], 2.0 * independent[0])

