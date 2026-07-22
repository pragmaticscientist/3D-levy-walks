from __future__ import annotations

import numpy as np
import pytest

from ltdb_levy.replay import (
    replay_batch,
    replay_targets,
    replay_targets_batch,
    replay_trial,
)
from ltdb_levy.targets import (
    Target,
    dimension_from_neighbourhood_volume,
    dimension_from_projected_area,
    minimum_image,
    projected_area,
    wrap_positions,
)


def test_minimum_image_and_wrapping_adversarial_inputs():
    side = 10.0
    values = np.array([-20.0, -10.0, -1e-12, 0.0, 10.0, 20.0, 25.0])
    wrapped = wrap_positions(values, side)
    assert np.all((wrapped >= 0.0) & (wrapped < side))
    assert np.allclose(wrapped[[0, 1, 3, 4, 5]], 0.0)
    image = minimum_image(values, side)
    assert np.all(image >= -side / 2.0)
    assert np.all(image < side / 2.0)
    assert image[2] == pytest.approx(-1e-12)
    assert image[-1] == pytest.approx(-5.0)


@pytest.mark.parametrize("shape", ["Ball", "Disk", "Line"])
def test_projected_area_roundtrip_and_membership_boundaries(shape):
    area = 25.0
    dimension = dimension_from_projected_area(area, shape, detection_radius=1.0)
    assert projected_area(dimension, shape, 1.0) == pytest.approx(area)
    target = Target(shape, dimension, side_length=100.0)
    center = np.asarray(target.center)
    if shape == "Ball":
        boundary = center + [dimension / 2.0 + 1.0, 0.0, 0.0]
        outside = boundary + [1e-8, 0.0, 0.0]
    elif shape == "Disk":
        boundary = center + [dimension / 2.0 + 1.0, 0.0, 1.0]
        outside = boundary + [0.0, 0.0, 1e-8]
    else:
        boundary = center + [1.0, dimension / 2.0, 0.0]
        outside = boundary + [1e-8, 0.0, 0.0]
    assert target.contains(boundary)
    assert not target.contains(outside)


def test_endpoint_only_detection_and_overshoot_policies():
    target = Target("Ball", dimension=0.0, side_length=10.0)
    start = np.array([3.5, 5.0, 5.0])
    vectors = np.array([[2.0, 0.0, 0.0]])

    discarded = replay_trial(vectors, start, target, budget=0.5, overshoot_policy="discard")
    assert not discarded.detected
    assert discarded.endpoints_checked == 0
    assert discarded.restricted_detection_distance == 0.5

    truncated = replay_trial(vectors, start, target, budget=0.5, overshoot_policy="truncate")
    assert truncated.detected
    assert truncated.endpoints_checked == 1
    assert truncated.detection_distance == pytest.approx(0.5)


def test_initial_position_is_not_checked_by_default():
    target = Target("Ball", dimension=0.0, side_length=10.0)
    start = np.array([5.0, 5.0, 5.0])
    vectors = np.array([[3.0, 0.0, 0.0]])
    primary = replay_trial(vectors, start, target, budget=3.0)
    sensitivity = replay_trial(
        vectors, start, target, budget=3.0, check_initial_position=True
    )
    assert not primary.detected
    assert sensitivity.detected
    assert sensitivity.detection_distance == 0.0


def test_multi_target_endpoint_reuse_matches_independent_replay():
    vectors = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    start = np.array([2.0, 5.0, 5.0])
    targets = [
        Target("Ball", 0.0, 10.0),
        Target("Disk", 2.0, 10.0),
        Target("Line", 3.0, 10.0),
    ]
    reused = replay_targets(vectors, start, targets, budget=3.0)
    independent = [
        replay_trial(vectors, start, target, budget=3.0) for target in targets
    ]
    for reused_result, independent_result in zip(reused, independent):
        assert reused_result.detected == independent_result.detected
        assert reused_result.first_hit_step == independent_result.first_hit_step
        assert reused_result.endpoints_checked == independent_result.endpoints_checked
        assert reused_result.budget == independent_result.budget
        assert np.allclose(
            reused_result.detection_distance,
            independent_result.detection_distance,
            equal_nan=True,
        )
        assert reused_result.restricted_detection_distance == pytest.approx(
            independent_result.restricted_detection_distance
        )

    padded = np.zeros((2, 2, 3))
    padded[0] = vectors
    padded[1, 0] = vectors[0]
    valid = np.array([[True, True], [True, False]])
    batched = replay_targets_batch(
        padded,
        np.vstack((start, start)),
        targets,
        np.array([3.0, 1.0]),
        valid_steps=valid,
    )
    for target_index, target in enumerate(targets):
        expected = replay_batch(
            padded,
            np.vstack((start, start)),
            target,
            np.array([3.0, 1.0]),
            valid_steps=valid,
        )
        assert np.array_equal(batched[target_index].detected, expected.detected)
        assert np.allclose(
            batched[target_index].restricted_detection_distance,
            expected.restricted_detection_distance,
        )


@pytest.mark.parametrize("shape", ["Ball", "Disk", "Line"])
def test_analytic_single_step_detection_probability(shape):
    """Translation invariance gives a known one-step answer on the torus."""
    rng = np.random.default_rng(450 + len(shape))
    n_trials = 250_000
    side = 20.0
    target = Target(shape, dimension=2.0, side_length=side)
    starts = rng.uniform(0.0, side, size=(n_trials, 3))
    step = np.array([3.7, -2.1, 4.4])
    vectors = np.broadcast_to(step, (n_trials, 1, 3)).copy()
    budgets = np.full(n_trials, np.linalg.norm(step) + 1e-12)
    result = replay_batch(vectors, starts, target, budgets)
    estimate = float(result.detected.mean())
    expected = target.neighbourhood_volume() / side ** 3
    standard_error = np.sqrt(expected * (1.0 - expected) / n_trials)
    assert abs(estimate - expected) <= 4.0 * standard_error


@pytest.mark.parametrize("shape", ["Ball", "Disk", "Line"])
@pytest.mark.parametrize("volume", [7.0, 40.0, 925.0])
def test_neighbourhood_volume_dimension_round_trip(shape, volume):
    dimension = dimension_from_neighbourhood_volume(volume, shape, 1.0)
    target = Target(
        shape=shape,
        dimension=dimension,
        side_length=2048.0,
        detection_radius=1.0,
    )
    assert target.neighbourhood_volume() == pytest.approx(volume)


def test_common_volume_must_exceed_disk_point_target_volume():
    with pytest.raises(Exception, match="below"):
        dimension_from_neighbourhood_volume(1.0, "Disk", 1.0)
