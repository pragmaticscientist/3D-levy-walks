from __future__ import annotations

import math

import pandas as pd
import pytest

from ltdb_levy.arrest_scale_experiment import (
    ArrestScaleExperimentSettings,
    _target_manifest_for_fit,
    line_neighbourhood_volume_limit,
    line_projected_area_limit,
)


@pytest.mark.parametrize(
    "side,area_limit,volume_limit",
    [
        (80.086, 159.3135926535898, 249.50316643653757),
        (533.786, 1066.7135926535898, 1674.8437790466987),
        (295.189, 589.5195926535898, 925.2691980267808),
    ],
)
def test_selected_line_geometry_limits(side, area_limit, volume_limit):
    assert line_projected_area_limit(side) == pytest.approx(area_limit)
    assert line_neighbourhood_volume_limit(side) == pytest.approx(volume_limit)


def test_nk_manifest_keeps_only_common_three_shape_values():
    settings = ArrestScaleExperimentSettings(
        goodness_of_fit_bootstraps=1,
        clustered_bootstraps=1,
    )
    fit = {
        "condition": "nk",
        "condition_id": "nk",
        "condition_label": "NK",
        "c_arrest_um": 1.0,
        "torus_side_model_units": 80.086,
    }
    manifest = _target_manifest_for_fit(fit, settings)
    areas = sorted(
        manifest.loc[
            manifest["matching_rule"] == "projected_area",
            "matched_value_model_units",
        ].unique()
    )
    volumes = sorted(
        manifest.loc[
            manifest["matching_rule"] == "neighbourhood_volume",
            "matched_value_model_units",
        ].unique()
    )
    assert areas == [4.0, 8.0, 16.0, 32.0, 64.0, 125.0]
    assert volumes == [7.0, 10.0, 20.0, 40.0, 80.0, 160.0]
    assert (manifest.groupby("matching_rule")["shape"].nunique() == 3).all()
    assert manifest["fits_without_periodic_self_overlap"].all()


def test_volume_targets_match_volume_not_projected_area():
    settings = ArrestScaleExperimentSettings(
        area_values=(4.0,),
        volume_values=(40.0,),
        goodness_of_fit_bootstraps=1,
        clustered_bootstraps=1,
    )
    fit = {
        "condition": "t",
        "condition_id": "t",
        "condition_label": "T",
        "c_arrest_um": 0.5,
        "torus_side_model_units": 295.189,
    }
    manifest = _target_manifest_for_fit(fit, settings)
    volume = manifest[manifest["matching_rule"] == "neighbourhood_volume"]
    assert set(volume["neighbourhood_volume_model_units3"].round(10)) == {40.0}
    assert volume["projected_area_model_units2"].nunique() == 3
    assert set(volume["matched_value_physical"].round(10)) == {5.0}
