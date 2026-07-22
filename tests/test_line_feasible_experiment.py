from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ltdb_levy.line_feasible_experiment import (
    DEFAULT_CANDIDATE_AREAS,
    _enhance_targets,
    line_feasible_areas,
    line_projected_area_limit,
)
from ltdb_levy.line_feasible_figures import (
    condition_c_label,
    normalized_area_ticks,
)


def test_line_feasibility_is_inclusive_at_periodic_fit_boundary():
    area = 12.0
    line_dimension = (area - math.pi) / 2.0
    exact_side = line_dimension + 2.0
    assert line_projected_area_limit(exact_side) == pytest.approx(area)
    assert line_feasible_areas((4.0, area), exact_side) == (4.0, area)


def test_line_grid_is_applied_as_the_common_shape_limit():
    nk_side_over_c = 34.339
    assert line_feasible_areas(DEFAULT_CANDIDATE_AREAS, nk_side_over_c) == (
        4.0,
        6.0,
        8.0,
        12.0,
        16.0,
        24.0,
        32.0,
        48.0,
        64.0,
    )
    # A=96 would fit the NK ball/disk geometry, but not its line capsule.  It
    # must therefore be excluded from every shape in the common comparison.
    ball_disk_limit = math.pi * nk_side_over_c**2 / 4.0
    assert ball_disk_limit > 96.0
    assert line_projected_area_limit(nk_side_over_c) < 96.0


def test_selected_condition_grid_maxima_match_the_approved_values():
    grids = {
        "NK": line_feasible_areas(DEFAULT_CANDIDATE_AREAS, 34.339),
        "Neutrophils": line_feasible_areas(DEFAULT_CANDIDATE_AREAS, 85.061),
        "T-HIV": line_feasible_areas(DEFAULT_CANDIDATE_AREAS, 35.249),
    }
    assert {name: values[-1] for name, values in grids.items()} == {
        "NK": 64.0,
        "Neutrophils": 128.0,
        "T-HIV": 64.0,
    }


def test_target_physical_columns_are_exact_c_rescalings():
    targets = pd.DataFrame(
        {
            "shape": ["Ball"],
            "projected_area": [8.0],
            "dimension": [1.25],
            "side_length": [34.0],
        }
    )
    result = _enhance_targets(targets, crossover_um=2.5).iloc[0]
    assert result["projected_area_model_units2"] == 8.0
    assert result["projected_area_um2"] == pytest.approx(50.0)
    assert result["dimension_um"] == pytest.approx(3.125)
    assert result["torus_side_um"] == pytest.approx(85.0)
    assert result["crossover_um"] == pytest.approx(2.5)


def test_figure_condition_annotation_contains_only_c_as_numeric_parameter():
    label = condition_c_label("NK, influenza, 30 s", 2.3322093)
    assert "NK, influenza, 30 s" in label
    assert "c=2.33" in label
    assert "\\mu\\mathrm{m}" in label
    assert "L_" not in label
    assert "S=" not in label
    assert "mu=" not in label


def test_normalized_area_ticks_are_sparse_powers_of_two():
    assert normalized_area_ticks(64.0) == (4.0, 8.0, 16.0, 32.0, 64.0)
    assert normalized_area_ticks(128.0) == (
        4.0,
        8.0,
        16.0,
        32.0,
        64.0,
        128.0,
    )
