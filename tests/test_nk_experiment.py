from __future__ import annotations

import numpy as np
import pytest

from ltdb_levy.config import (
    Config,
    ConditionsConfig,
    EmpiricalPoolConfig,
    FittingConfig,
    InputConfig,
    OutputConfig,
    PreprocessingConfig,
    ReplayConfig,
    SegmentationConfig,
    SelectionConfig,
)
from ltdb_levy.nk_experiment import (
    _draw_renewable_batch,
    _finite_sequences,
    _make_targets,
    _mixture_mean,
    _rotate_vector_batch,
    _unbounded_contrasts,
)


def _target_test_config(tmp_path):
    return Config(
        schema_version=1,
        run_id="test",
        input=InputConfig(tracks_directory=tmp_path),
        selection=SelectionConfig(),
        preprocessing=PreprocessingConfig(),
        segmentation=SegmentationConfig(),
        conditions=ConditionsConfig(),
        empirical_pool=EmpiricalPoolConfig(),
        fitting=FittingConfig(),
        replay=ReplayConfig(shapes=["Ball", "Disk", "Line"]),
        output=OutputConfig(root=tmp_path),
    )


def test_batch_global_rotation_preserves_path_geometry():
    vectors = np.array(
        [
            [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        ]
    )
    quaternions = np.array(
        [
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, np.sqrt(0.5), np.sqrt(0.5)],
        ]
    )
    rotated = _rotate_vector_batch(vectors, quaternions)
    assert np.allclose(
        np.linalg.norm(rotated, axis=2), np.linalg.norm(vectors, axis=2)
    )
    original_dot = np.sum(vectors[:, 0] * vectors[:, 1], axis=1)
    rotated_dot = np.sum(rotated[:, 0] * rotated[:, 1], axis=1)
    assert np.allclose(rotated_dot, original_dot)


def test_vector_uniform_class_draws_only_empirical_lengths():
    pool = np.array(
        [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 4.0]]
    )
    drawn = _draw_renewable_batch(
        "vector_uniform_rotated",
        n_trials=20,
        n_steps=30,
        rng=np.random.default_rng(101),
        pool_vectors=pool,
        fitted_mu=2.0,
        fitted_crossover=1.0,
        fitted_lmax=10.0,
    )
    lengths = set(np.round(np.linalg.norm(drawn, axis=2).ravel(), 12))
    assert lengths.issubset({1.0, 2.0, 4.0})


def test_exact_classes_keep_source_vector_lengths_and_cover_budget():
    source = np.array(
        [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    pool = source.copy()
    budget = 6.0
    for class_name in (
        "exact_orientation",
        "exact_global_rotation",
        "ordered_length_rotated",
    ):
        sequences = _finite_sequences(
            class_name,
            source,
            pool,
            budget,
            n_trials=8,
            rng=np.random.default_rng(102),
            fitted_mu=2.2,
            fitted_crossover=1.0,
            fitted_lmax=10.0,
        )
        assert np.allclose(
            np.linalg.norm(sequences, axis=2),
            np.broadcast_to([1.0, 2.0, 3.0], (8, 3)),
        )


def test_ordered_length_rotation_destroys_internal_path_geometry():
    source = np.array(
        [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    sequences = _finite_sequences(
        "ordered_length_rotated",
        source,
        source,
        budget=6.0,
        n_trials=32,
        rng=np.random.default_rng(104),
        fitted_mu=2.2,
        fitted_crossover=1.0,
        fitted_lmax=10.0,
    )
    assert np.allclose(
        np.linalg.norm(sequences, axis=2),
        np.broadcast_to([1.0, 2.0, 3.0], (32, 3)),
    )
    adjacent_cosines = np.sum(
        sequences[:, :-1] * sequences[:, 1:], axis=2
    ) / (
        np.linalg.norm(sequences[:, :-1], axis=2)
        * np.linalg.norm(sequences[:, 1:], axis=2)
    )
    assert not np.allclose(adjacent_cosines, 1.0)


def test_mixture_mean_matches_monte_carlo():
    rng = np.random.default_rng(103)
    mu = 2.3
    crossover = 0.6
    lmax = 10.7
    uniforms = rng.random(500_000)
    sampled = __import__(
        "ltdb_levy.dist.mixture", fromlist=["ppf"]
    ).ppf(uniforms, mu, lmax, crossover)
    np.testing.assert_allclose(
        np.mean(sampled),
        _mixture_mean(mu, lmax, crossover),
        rtol=0.01,
    )


def test_unbounded_contrast_reports_ratio_and_direction():
    import pandas as pd

    summary = pd.DataFrame(
        {
            "trajectory_class": ["fitted_mu", "vector_uniform_rotated"],
            "shape": ["Ball", "Ball"],
            "projected_area": [8.0, 8.0],
            "mean_detection_distance_model_units": [80.0, 100.0],
            "detection_distance_mc_se_model_units": [2.0, 3.0],
        }
    )
    result = _unbounded_contrasts(summary).iloc[0]
    assert result["fitted_minus_empirical_distance"] == -20.0
    assert result["fitted_to_empirical_distance_ratio"] == 0.8
    assert result["relative_distance_change_percent"] == pytest.approx(-20.0)


def test_target_grid_requires_full_neighbourhood_to_fit_torus(tmp_path):
    config = _target_test_config(tmp_path)
    targets = _make_targets(config, [32.0], side_length=21.356)
    assert len(targets) == 3
    assert all(
        target.dimension + 2.0 * target.detection_radius <= target.side_length
        for _, target in targets
    )
    with pytest.raises(ValueError, match="periodic self-overlap"):
        _make_targets(config, [64.0], side_length=21.356)
