from __future__ import annotations

import numpy as np
import pandas as pd

from ltdb_levy.preprocess import DISPLACEMENT_COLUMNS
from ltdb_levy.segment import SegmentParams, segment_runs


def _displacements(vectors):
    vectors = np.asarray(vectors, dtype=float)
    lengths = np.linalg.norm(vectors, axis=1)
    angles = np.full(len(vectors), np.nan)
    if len(vectors) > 1:
        angles[1:] = np.arctan2(
            np.linalg.norm(np.cross(vectors[:-1], vectors[1:]), axis=1),
            np.einsum("ij,ij->i", vectors[:-1], vectors[1:]),
        )
    n = len(vectors)
    return pd.DataFrame(
        {
            "file": "TEST.csv",
            "cohort": "TEST",
            "video": "TEST001",
            "population": "a",
            "video_id": "TEST001_a",
            "track_id": 1,
            "track_uid": "TEST001_a#1",
            "fragment_id": 0,
            "fragment_uid": "TEST001_a#1/0",
            "step_index": np.arange(n),
            "frame_from": np.arange(n),
            "frame_to": np.arange(1, n + 1),
            "time_s": np.arange(n, dtype=float),
            "dt_s": np.ones(n),
            "dx_um": vectors[:, 0],
            "dy_um": vectors[:, 1],
            "dz_um": vectors[:, 2],
            "length_um": lengths,
            "speed_um_s": lengths,
            "turn_angle_rad": angles,
        },
        columns=DISPLACEMENT_COLUMNS,
    )


def test_turning_displacement_starts_new_run_and_threshold_is_strict():
    displacements = _displacements(
        [[1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0]]
    )
    split, _ = segment_runs(
        displacements,
        SegmentParams(turning_angle_degrees=89, exclude_first_and_last_run=False),
    )
    assert list(zip(split["step_start"], split["step_end"])) == [(0, 1), (2, 3)]

    tied, _ = segment_runs(
        displacements,
        SegmentParams(turning_angle_degrees=90, exclude_first_and_last_run=False),
    )
    assert len(tied) == 1
    assert tied.iloc[0]["n_steps"] == 4


def test_first_and_last_runs_are_flagged_not_dropped():
    displacements = _displacements([[1, 0, 0], [0, 1, 0], [1, 0, 0]])
    runs, audit = segment_runs(
        displacements,
        SegmentParams(turning_angle_degrees=45, exclude_first_and_last_run=True),
    )
    assert len(runs) == 3
    assert runs["included"].tolist() == [False, True, False]
    assert audit.excluded_first_last == 2


def test_speed_threshold_omits_slow_steps_and_is_inclusive():
    displacements = _displacements(
        [[1.0, 0, 0], [0.5, 0, 0], [0.75, 0, 0], [1.0, 0, 0]]
    )
    runs, audit = segment_runs(
        displacements,
        SegmentParams(
            turning_angle_degrees=90,
            speed_threshold=0.75,
            exclude_first_and_last_run=False,
        ),
    )

    assert list(zip(runs["step_start"], runs["step_end"])) == [(0, 0), (2, 3)]
    assert runs["n_steps"].tolist() == [1, 2]
    assert audit.non_moving_steps == 1
