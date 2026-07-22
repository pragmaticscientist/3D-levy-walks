from __future__ import annotations

import numpy as np
import pytest

from conftest import make_observations
from ltdb_levy.preprocess import (
    PreprocessAudit,
    preprocess,
    resolve_duplicates,
    turning_angles,
)


def test_missing_frame_splits_fragments_without_bridging():
    observations = make_observations(
        [0, 1, 3, 4],
        [(0, 0, 0), (1, 0, 0), (100, 0, 0), (101, 0, 0)],
    )
    cleaned, displacements, audit = preprocess(
        observations, minimum_observations_per_track=2
    )
    assert cleaned["fragment_uid"].nunique() == 2
    assert len(displacements) == 2
    assert np.allclose(displacements["length_um"], 1.0)
    assert audit.fragments_created_by_gaps == 1


def test_turning_angles_use_stable_full_range_geometry():
    vectors = np.array(
        [[1.0, 0.0, 0.0], [1.0, 1e-15, 0.0], [-1.0, 0.0, 0.0]]
    )
    angles = turning_angles(vectors)
    assert np.isnan(angles[0])
    assert angles[1] == pytest.approx(1e-15, abs=1e-18)
    assert angles[2] == pytest.approx(np.pi, abs=1e-14)


def test_duplicate_mean_policy_preserves_schema_and_averages_coordinates():
    observations = make_observations(
        [0, 0, 1], [(0, 0, 0), (2, 4, 6), (3, 4, 5)]
    )
    audit = PreprocessAudit()
    resolved = resolve_duplicates(observations, policy="mean", audit=audit)
    assert list(resolved.columns) == list(observations.columns)
    first = resolved.loc[resolved["frame"] == 0].iloc[0]
    assert (first["x_um"], first["y_um"], first["z_um"]) == (1.0, 2.0, 3.0)
    assert audit.duplicate_rows_dropped == 1

