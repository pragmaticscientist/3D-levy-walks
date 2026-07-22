from __future__ import annotations

import numpy as np
import pandas as pd
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
from ltdb_levy.ordered_rotation_experiment import (
    BASELINE_FINITE_CLASSES,
    NEW_CLASS,
    NEW_FINITE_CONTRASTS,
    OrderedRotationReplaySettings,
    _baseline_signature_payload,
    _build_condition_tasks,
    _paired_track_contrasts,
)
from ltdb_levy.targets import Target


def _config(tmp_path):
    return Config(
        schema_version=1,
        run_id="test",
        input=InputConfig(tracks_directory=tmp_path),
        selection=SelectionConfig(),
        preprocessing=PreprocessingConfig(),
        segmentation=SegmentationConfig(),
        conditions=ConditionsConfig(),
        empirical_pool=EmpiricalPoolConfig(primary_weighting="step_weighted"),
        fitting=FittingConfig(clustered_bootstraps=100),
        replay=ReplayConfig(
            shapes=["Ball", "Disk", "Line"],
            projected_areas=[4.0],
        ),
        output=OutputConfig(root=tmp_path),
    )


def test_baseline_signature_payload_is_pinned_to_old_class_ladder(tmp_path):
    config = _config(tmp_path)
    baseline_settings = {
        "projected_areas": [4.0, 6.0],
        "finite_replays_per_track": 10_000,
        "unbounded_trials_per_cell": 2_000,
        "finite_batch_size": 1_000,
        "unbounded_chunk_trials": 25,
        "unbounded_step_block": 256,
        "random_seed": 24_680,
    }
    payload = _baseline_signature_payload(
        config,
        baseline_settings,
        "condition-a",
        {"mu_hat": 2.3, "crossover": 4.0, "lmax": 40.0},
        {"runs": "a" * 64, "pools": "b" * 64, "mixture_fits": "c" * 64},
    )
    assert tuple(payload["finite_classes"]) == BASELINE_FINITE_CLASSES
    assert NEW_CLASS not in payload["finite_classes"]
    assert payload["fitted_crossover_model_units"] == 1.0
    assert payload["fitted_lmax_model_units"] == 10.0


def test_condition_tasks_preserve_exact_length_order_and_budget(tmp_path):
    config = _config(tmp_path)
    settings = OrderedRotationReplaySettings(
        finite_replays_per_track=7,
        workers=1,
        finite_batch_size=4,
        random_seed=123,
    )
    condition = "condition-a"
    runs = pd.DataFrame(
        {
            "condition": [condition] * 3,
            "track_uid": ["track-1"] * 3,
            "video_id": ["video-1"] * 3,
            "fragment_uid": ["f"] * 3,
            "run_id": [2, 0, 1],
            "frame_from": [2, 0, 1],
            "frame_to": [3, 1, 2],
            "run_length_um": [6.0, 2.0, 4.0],
            "run_vector_x": [0.0, 2.0, 0.0],
            "run_vector_y": [0.0, 0.0, 4.0],
            "run_vector_z": [6.0, 0.0, 0.0],
        }
    )
    pool = runs.iloc[[1, 2]].copy()
    baseline_tracks = pd.DataFrame(
        {
            "condition": [condition],
            "video_id": ["video-1"],
            "track_uid": ["track-1"],
            "n_runs": [3],
            "budget_model_units": [6.0],
        }
    )
    target_rows = []
    for shape in config.replay.shapes:
        target = Target(
            shape=shape,
            dimension=0.25,
            side_length=20.0,
            detection_radius=1.0,
        )
        target_rows.append(
            {
                "shape": shape,
                "projected_area": 4.0,
                "dimension": target.dimension,
                "side_length": target.side_length,
                "detection_radius": target.detection_radius,
            }
        )
    tasks = _build_condition_tasks(
        config,
        settings,
        condition,
        {"mu_hat": 2.2, "crossover": 2.0, "lmax": 20.0},
        runs,
        pool,
        baseline_tracks,
        pd.DataFrame(target_rows),
    )
    assert len(tasks) == 1
    task = tasks[0]
    np.testing.assert_allclose(
        np.linalg.norm(task["source_vectors"], axis=1),
        [1.0, 2.0, 3.0],
    )
    assert np.linalg.norm(task["source_vectors"], axis=1).sum() == pytest.approx(
        6.0
    )
    assert task["finite_classes"] == (NEW_CLASS,)
    assert task["n_replays"] == 7


def test_new_contrasts_are_paired_over_tracks_and_have_expected_direction():
    rows = []
    probabilities = {
        "t1": {
            "exact_global_rotation": 0.10,
            NEW_CLASS: 0.20,
            "vector_uniform_rotated": 0.25,
        },
        "t2": {
            "exact_global_rotation": 0.30,
            NEW_CLASS: 0.35,
            "vector_uniform_rotated": 0.50,
        },
    }
    for track_uid, values in probabilities.items():
        for class_name, probability in values.items():
            rows.append(
                {
                    "condition": "condition-a",
                    "track_uid": track_uid,
                    "trajectory_class": class_name,
                    "shape": "Ball",
                    "projected_area": 4.0,
                    "detection_probability": probability,
                }
            )
    result = _paired_track_contrasts(
        pd.DataFrame(rows), n_bootstraps=200, random_seed=99
    )
    pairs = {
        (row.left_class, row.right_class): row
        for row in result.itertuples(index=False)
    }
    assert set(pairs) == set(NEW_FINITE_CONTRASTS)
    assert pairs[(NEW_CLASS, "exact_global_rotation")].n_tracks == 2
    assert pairs[
        (NEW_CLASS, "exact_global_rotation")
    ].detection_probability_difference == pytest.approx(0.075)
    assert pairs[
        ("vector_uniform_rotated", NEW_CLASS)
    ].detection_probability_difference == pytest.approx(0.10)


def test_new_contrast_rejects_missing_track_pair():
    frame = pd.DataFrame(
        {
            "condition": ["condition-a"] * 5,
            "track_uid": ["t1", "t1", "t1", "t2", "t2"],
            "trajectory_class": [
                "exact_global_rotation",
                NEW_CLASS,
                "vector_uniform_rotated",
                "exact_global_rotation",
                NEW_CLASS,
            ],
            "shape": ["Ball"] * 5,
            "projected_area": [4.0] * 5,
            "detection_probability": [0.1, 0.2, 0.3, 0.1, 0.2],
        }
    )
    with pytest.raises(ValueError, match="every track"):
        _paired_track_contrasts(frame, n_bootstraps=10, random_seed=2)
