from __future__ import annotations

import pandas as pd
import pytest

from ltdb_levy.multicondition_experiment import (
    AllConditionsReplaySettings,
    REPLAY_SIGNATURE_VERSION,
    _canonical_signature,
    _expected_replay_row_counts,
    _ordered_shape_area_pairs,
    _stored_signature_matches,
    _validate_replay_frames,
)
from ltdb_levy.nk_experiment import (
    FINITE_CLASSES,
    FINITE_CONTRASTS,
    UNBOUNDED_CLASSES,
)


def test_replay_signature_is_canonical_and_rejects_legacy_or_changed_settings():
    first = {
        "condition": "test",
        "projected_areas": [4.0, 8.0],
        "finite_replays": 100,
    }
    reordered = {
        "finite_replays": 100,
        "projected_areas": [4.0, 8.0],
        "condition": "test",
    }
    signature = _canonical_signature(first)
    assert signature == _canonical_signature(reordered)
    assert signature != _canonical_signature(
        dict(first, finite_replays=101)
    )
    assert not _stored_signature_matches({}, signature)
    assert not _stored_signature_matches(
        {
            "replay_signature_version": REPLAY_SIGNATURE_VERSION - 1,
            "replay_signature": signature,
        },
        signature,
    )
    assert _stored_signature_matches(
        {
            "replay_signature_version": REPLAY_SIGNATURE_VERSION,
            "replay_signature": signature,
        },
        signature,
    )


def test_expected_replay_counts_are_exact():
    settings = AllConditionsReplaySettings(
        projected_areas=(4.0, 8.0),
        finite_replays_per_track=7,
        unbounded_trials_per_cell=11,
    )
    assert _expected_replay_row_counts(
        n_tracks=3,
        n_pool_vectors=17,
        n_shapes=3,
        settings=settings,
    ) == {
        "track_inputs": 3,
        "empirical_pool_inputs": 17,
        "finite_track_summaries": 90,
        "finite_summary": 30,
        "finite_contrasts": 30,
        "unbounded_trials": 132,
        "unbounded_summary": 12,
        "unbounded_contrasts": 6,
        "targets": 6,
    }


def test_shape_area_plot_order_is_canonical_and_numeric():
    frame = pd.DataFrame(
        {
            "shape": ["Line", "Ball", "Disk", "Ball", "Line", "Disk"],
            "projected_area": [4.0, 16.0, 24.0, 4.0, 24.0, 4.0],
        }
    )
    assert _ordered_shape_area_pairs(frame) == [
        ("Ball", 4.0),
        ("Ball", 16.0),
        ("Disk", 4.0),
        ("Disk", 24.0),
        ("Line", 4.0),
        ("Line", 24.0),
    ]


def _tiny_complete_replay():
    settings = AllConditionsReplaySettings(
        projected_areas=(4.0,),
        finite_replays_per_track=2,
        unbounded_trials_per_cell=3,
    )
    finite_tracks = pd.DataFrame(
        [
            {
                "track_uid": "t1",
                "trajectory_class": class_name,
                "shape": "Ball",
                "projected_area": 4.0,
                "n_trials": 2,
            }
            for class_name in FINITE_CLASSES
        ]
    )
    finite_summary = pd.DataFrame(
        [
            {
                "trajectory_class": class_name,
                "shape": "Ball",
                "projected_area": 4.0,
                "n_tracks": 1,
                "n_trials": 2,
            }
            for class_name in FINITE_CLASSES
        ]
    )
    finite_contrasts = pd.DataFrame(
        [
            {
                "left_class": left,
                "right_class": right,
                "shape": "Ball",
                "projected_area": 4.0,
            }
            for left, right in FINITE_CONTRASTS
        ]
    )
    unbounded_trials = pd.DataFrame(
        [
            {
                "trajectory_class": class_name,
                "shape": "Ball",
                "projected_area": 4.0,
                "trial_index": trial,
            }
            for class_name in UNBOUNDED_CLASSES
            for trial in range(3)
        ]
    )
    unbounded_summary = pd.DataFrame(
        [
            {
                "trajectory_class": class_name,
                "shape": "Ball",
                "projected_area": 4.0,
                "n_trials": 3,
            }
            for class_name in UNBOUNDED_CLASSES
        ]
    )
    frames = {
        "track_inputs": pd.DataFrame({"track_uid": ["t1"]}),
        "empirical_pool_inputs": pd.DataFrame(
            {
                "run_uid": ["r1", "r2"],
                "uniform_run_sampling_probability": [0.5, 0.5],
            }
        ),
        "finite_track_summaries": finite_tracks,
        "finite_summary": finite_summary,
        "finite_contrasts": finite_contrasts,
        "unbounded_trials": unbounded_trials,
        "unbounded_summary": unbounded_summary,
        "unbounded_contrasts": pd.DataFrame(
            {"shape": ["Ball"], "projected_area": [4.0]}
        ),
        "targets": pd.DataFrame(
            {
                "shape": ["Ball"],
                "projected_area": [4.0],
                "fits_without_periodic_self_overlap": [True],
            }
        ),
    }
    summary = {
        "tracks": 1,
        "pool_vectors": 2,
        "finite_trial_count": 10,
        "unbounded_trial_count": 6,
    }
    return settings, frames, summary


def test_replay_frame_validator_fails_closed_on_bad_counts():
    settings, frames, summary = _tiny_complete_replay()
    _validate_replay_frames(
        frames,
        summary,
        settings,
        n_tracks=1,
        n_pool_vectors=2,
        shapes=["Ball"],
    )

    damaged = dict(frames)
    damaged["unbounded_trials"] = frames["unbounded_trials"].iloc[:-1].copy()
    with pytest.raises(ValueError, match="unbounded_trials has 5 rows"):
        _validate_replay_frames(
            damaged,
            summary,
            settings,
            n_tracks=1,
            n_pool_vectors=2,
            shapes=["Ball"],
        )
