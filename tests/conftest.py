"""Shared pytest helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def make_observations(frames, positions, dt=1.0, track_uid="TEST001_a#1"):
    """Create a minimal tidy observation table for preprocessing tests."""
    rows = []
    for frame, position in zip(frames, positions):
        rows.append(
            {
                "file": "TEST001_a_GT.csv",
                "cohort": "TEST",
                "video": "TEST001",
                "population": "a",
                "video_id": "TEST001_a",
                "track_id": 1,
                "track_uid": track_uid,
                "x_um": float(position[0]),
                "y_um": float(position[1]),
                "z_um": float(position[2]),
                "frame": int(frame),
                "time_s": float(frame) * dt,
            }
        )
    return pd.DataFrame(rows)

