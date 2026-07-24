from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from ltdb_levy.arrest_scale_experiment import (
    SELECTED_CONDITIONS,
    ArrestScaleExperimentSettings,
    _make_figures,
    _target_manifest_for_fit,
    line_neighbourhood_volume_limit,
    line_projected_area_limit,
    load_arrest_scale_settings,
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


def test_arrest_scale_settings_resolve_paths_from_config(tmp_path):
    source = tmp_path / "source.yaml"
    source.write_text("placeholder\n", encoding="utf-8")
    config = tmp_path / "arrest.yaml"
    config.write_text(
        "\n".join(
            [
                "schema_version: 1",
                "source_config: source.yaml",
                "output_root: generated",
                "area_values: [4.0]",
                "volume_values: [7.0]",
                "finite_replays_per_track: 2",
                "finite_batch_size: 1",
                "workers: 1",
                "clustered_bootstraps: 1",
                "goodness_of_fit_bootstraps: 1",
                "random_seed: 42",
                "resume: false",
            ]
        ),
        encoding="utf-8",
    )

    settings = load_arrest_scale_settings(config)

    assert settings.source_config == source
    assert settings.output_root == tmp_path / "generated"
    assert settings.area_values == (4.0,)
    assert settings.volume_values == (7.0,)
    assert not settings.resume


def test_arrest_scale_plotting_generates_all_figure_layouts(tmp_path):
    rows = []
    for condition_index, spec in enumerate(SELECTED_CONDITIONS):
        for shape_index, shape in enumerate(("Ball", "Disk", "Line")):
            for trajectory_class in ("ordered_length_rotated", "fitted_mu"):
                for matched_value in (4.0, 8.0):
                    probability = (
                        0.0001
                        * (condition_index + 1)
                        * (shape_index + 1)
                        * (1.0 if matched_value == 4.0 else 1.5)
                    )
                    rows.append(
                        {
                            "matching_rule": "projected_area",
                            "trajectory_class": trajectory_class,
                            "condition_id": spec.condition_id,
                            "shape": shape,
                            "matched_value_model_units": matched_value,
                            "detection_probability": probability,
                            "detection_probability_ci_low": probability * 0.9,
                            "detection_probability_ci_high": probability * 1.1,
                        }
                    )
    summary = pd.DataFrame(rows)
    fits = pd.DataFrame(
        {"condition_id": [spec.condition_id for spec in SELECTED_CONDITIONS]}
    )

    paths = _make_figures(Path(tmp_path), summary, fits)

    assert len(paths) == 5
    assert all(path.is_file() for path in paths)
    assert all(path.with_suffix(".png").is_file() for path in paths)
