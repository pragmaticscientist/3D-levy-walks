from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ltdb_levy.nk_experiment import FINITE_PLOT_ORDER, UNBOUNDED_CLASSES
from ltdb_levy.publication_figures import (
    SELECTED_CONDITIONS,
    SHAPE_ORDER,
    _empirical_ccdf,
    _finite_compact_ratio_table,
    _select_conditions,
    _unbounded_compact_ratio_table,
    _unbounded_ratio_table,
    _validate_selected_data,
)


def test_empirical_ccdf_is_non_strict_and_keeps_largest_observation_mass():
    lengths, survival = _empirical_ccdf([3.0, 1.0, 2.0, 2.0])
    np.testing.assert_allclose(lengths, [1.0, 2.0, 2.0, 3.0])
    np.testing.assert_allclose(survival, [1.0, 0.75, 0.50, 0.25])
    assert survival[-1] == pytest.approx(1.0 / len(lengths))


def test_selected_condition_filter_fails_when_a_prespecified_condition_is_absent():
    incomplete = pd.DataFrame(
        {
            "condition": [SELECTED_CONDITIONS[0], SELECTED_CONDITIONS[1]],
            "value": [1, 2],
        }
    )
    with pytest.raises(ValueError, match="missing selected condition"):
        _select_conditions(incomplete, "test table")


def _complete_selected_frames():
    fit_rows = []
    run_rows = []
    finite_rows = []
    unbounded_rows = []
    for condition_index, condition in enumerate(SELECTED_CONDITIONS):
        fit_rows.append(
            {
                "condition": condition,
                "condition_id": "condition-{}".format(condition_index),
                "condition_label": "Condition {}".format(condition_index),
                "n_runs": 2,
                "n_tracks": 1,
                "crossover": 2.0,
            }
        )
        for run_index, length in enumerate((2.0, 4.0)):
            run_rows.append(
                {
                    "condition": condition,
                    "run_uid": "{}-run-{}".format(condition_index, run_index),
                    "run_length_um": length,
                    "run_length_model_units": length / 2.0,
                }
            )
        for class_name in FINITE_PLOT_ORDER:
            for shape in SHAPE_ORDER:
                finite_rows.append(
                    {
                        "condition": condition,
                        "trajectory_class": class_name,
                        "shape": shape,
                        "projected_area": 4.0,
                        "n_tracks": 1,
                        "detection_probability": 0.2,
                        "detection_probability_ci_low": 0.1,
                        "detection_probability_ci_high": 0.3,
                    }
                )
        for class_name in UNBOUNDED_CLASSES:
            for shape in SHAPE_ORDER:
                unbounded_rows.append(
                    {
                        "condition": condition,
                        "trajectory_class": class_name,
                        "shape": shape,
                        "projected_area": 4.0,
                        "mean_detection_distance_model_units": 10.0,
                        "detection_distance_mc_se_model_units": 0.5,
                    }
                )
    return (
        pd.DataFrame(fit_rows),
        pd.DataFrame(run_rows),
        pd.DataFrame(finite_rows),
        pd.DataFrame(unbounded_rows),
    )


def test_selected_grid_validator_requires_all_five_finite_classes():
    fits, runs, finite, unbounded = _complete_selected_frames()
    assert _validate_selected_data(fits, runs, finite, unbounded) == (4.0,)

    damaged = finite[
        finite["trajectory_class"] != "ordered_length_rotated"
    ].copy()
    with pytest.raises(ValueError, match="trajectory-class set"):
        _validate_selected_data(fits, runs, damaged, unbounded)


def test_unbounded_ratio_uses_independent_sample_delta_interval():
    condition = SELECTED_CONDITIONS[0]
    frame = pd.DataFrame(
        [
            {
                "condition": condition,
                "condition_id": "nk",
                "condition_label": "NK",
                "shape": "Ball",
                "projected_area": 4.0,
                "trajectory_class": "fitted_mu",
                "mean_detection_distance_model_units": 20.0,
                "detection_distance_mc_se_model_units": 2.0,
            },
            {
                "condition": condition,
                "condition_id": "nk",
                "condition_label": "NK",
                "shape": "Ball",
                "projected_area": 4.0,
                "trajectory_class": "vector_uniform_rotated",
                "mean_detection_distance_model_units": 10.0,
                "detection_distance_mc_se_model_units": 1.0,
            },
        ]
    )
    result = _unbounded_ratio_table(frame).iloc[0]
    expected_se = np.sqrt((2.0 / 10.0) ** 2 + (20.0 * 1.0 / 10.0**2) ** 2)
    expected_log_se = np.sqrt((2.0 / 20.0) ** 2 + (1.0 / 10.0) ** 2)
    assert result["fitted_to_empirical_ratio"] == pytest.approx(2.0)
    assert result["ratio_mc_se_delta"] == pytest.approx(expected_se)
    assert result["log_ratio_mc_se_delta"] == pytest.approx(expected_log_se)
    assert result["ratio_ci_low"] == pytest.approx(
        np.exp(np.log(2.0) - 1.96 * expected_log_se)
    )
    assert result["ratio_ci_high"] == pytest.approx(
        np.exp(np.log(2.0) + 1.96 * expected_log_se)
    )


def test_compact_finite_ratios_use_paired_whole_track_bootstrap():
    condition = SELECTED_CONDITIONS[0]
    multipliers = {
        "exact_orientation": 1.1,
        "exact_global_rotation": 1.0,
        "ordered_length_rotated": 1.4,
        "vector_uniform_rotated": 1.5,
        "fitted_mu": 1.6,
    }
    rows = []
    for track_uid, baseline in (("track-a", 0.01), ("track-b", 0.02)):
        for class_name in FINITE_PLOT_ORDER:
            for shape in SHAPE_ORDER:
                for area in (4.0, 6.0, 8.0):
                    rows.append(
                        {
                            "condition": condition,
                            "condition_id": "nk",
                            "condition_label": "NK",
                            "track_uid": track_uid,
                            "trajectory_class": class_name,
                            "shape": shape,
                            "projected_area": area,
                            "detection_probability": (
                                baseline * multipliers[class_name]
                            ),
                        }
                    )
    result = _finite_compact_ratio_table(
        pd.DataFrame(rows),
        bootstrap_replicates=200,
        bootstrap_seed=17,
    ).set_index("trajectory_class")
    assert len(result) == len(FINITE_PLOT_ORDER)
    assert set(result["n_target_cells"].astype(int)) == {9}
    assert set(result["n_tracks"].astype(int)) == {2}
    for class_name, multiplier in multipliers.items():
        assert result.loc[
            class_name, "normalized_detection_probability"
        ] == pytest.approx(multiplier)
        assert result.loc[class_name, "ratio_ci_low"] == pytest.approx(
            multiplier
        )
        assert result.loc[class_name, "ratio_ci_high"] == pytest.approx(
            multiplier
        )


def test_compact_unbounded_ratio_propagates_target_cell_mc_errors():
    condition = SELECTED_CONDITIONS[0]
    rows = []
    values = {
        "vector_uniform_rotated": ((10.0, 14.0), (1.0, 2.0)),
        "fitted_mu": ((20.0, 28.0), (2.0, 4.0)),
    }
    for class_name, (means, errors) in values.items():
        for area, mean, error in zip((4.0, 6.0), means, errors):
            rows.append(
                {
                    "condition": condition,
                    "condition_id": "nk",
                    "condition_label": "NK",
                    "trajectory_class": class_name,
                    "shape": "Ball",
                    "projected_area": area,
                    "n_trials": 2000,
                    "mean_detection_distance_model_units": mean,
                    "detection_distance_mc_se_model_units": error,
                }
            )
    result = _unbounded_compact_ratio_table(pd.DataFrame(rows)).set_index(
        "trajectory_class"
    )
    empirical = result.loc["vector_uniform_rotated"]
    fitted = result.loc["fitted_mu"]
    expected_empirical_se = np.sqrt(1.0**2 + 2.0**2) / 2.0
    expected_fitted_se = np.sqrt(2.0**2 + 4.0**2) / 2.0
    expected_log_se = np.sqrt(
        (expected_fitted_se / 24.0) ** 2
        + (expected_empirical_se / 12.0) ** 2
    )
    assert empirical[
        "target_grid_mean_detection_distance_model_units"
    ] == pytest.approx(12.0)
    assert empirical["normalized_detection_distance"] == pytest.approx(1.0)
    assert empirical["ratio_ci_low"] == pytest.approx(1.0)
    assert fitted[
        "target_grid_mean_detection_distance_model_units"
    ] == pytest.approx(24.0)
    assert fitted["normalized_detection_distance"] == pytest.approx(2.0)
    assert fitted["log_ratio_mc_se_delta"] == pytest.approx(expected_log_se)
    assert fitted["ratio_ci_low"] == pytest.approx(
        np.exp(np.log(2.0) - 1.96 * expected_log_se)
    )
