from __future__ import annotations

import numpy as np
import pandas as pd

from ltdb_levy.stats import clustered_outcome_intervals


def test_clustered_outcome_bootstrap_preserves_paired_track_cells():
    rows = []
    for track, base in (("a", 2.0), ("b", 8.0), ("c", 5.0)):
        for trajectory_class, multiplier in (
            ("empirical", 1.0),
            ("ordered_length_rotated", 2.0),
        ):
            for shape, shape_multiplier in (("Ball", 1.0), ("Line", 1.5)):
                rows.append(
                    {
                        "condition": "condition",
                        "track_uid": track,
                        "trajectory_class": trajectory_class,
                        "shape": shape,
                        "projected_area": 100.0,
                        "detection_probability": base / 10.0,
                        "restricted_mean_detection_distance": (
                            base * multiplier * shape_multiplier
                        ),
                    }
                )
    frame = pd.DataFrame(rows)
    outcomes, contrasts, shapes = clustered_outcome_intervals(
        frame, 500, np.random.default_rng(91)
    )
    assert len(outcomes) == 4
    assert len(contrasts) == 2
    assert len(shapes) == 2
    # Pairing makes empirical/rotated log-RMDD contrast exactly log(1/2)
    # under every possible track resample.
    assert np.allclose(
        contrasts["log_rmdd_difference_ci_low"], -np.log(2.0)
    )
    assert np.allclose(
        contrasts["log_rmdd_difference_ci_high"], -np.log(2.0)
    )

