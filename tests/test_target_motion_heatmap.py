from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from conftest import REPOSITORY_ROOT


NOTEBOOK_PATH = REPOSITORY_ROOT / "simulation-plots.ipynb"


def _target_motion_plotting_namespace():
    notebook = json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))
    definition_source = next(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
        and "def build_target_motion_comparison(" in "".join(cell["source"])
    )
    namespace = {"Path": Path, "np": np, "pd": pd, "plt": plt}
    exec(compile(definition_source, str(NOTEBOOK_PATH), "exec"), namespace)
    return namespace


def _experiment_rows(
    motion: str, target_step_length: float, detection_times: tuple[float, float]
) -> pd.DataFrame:
    rows = []
    for shape in ("Ball", "Disk", "Line"):
        for detection_time in detection_times:
            rows.append(
                {
                    "mu": 1.0,
                    "surface": 64.0,
                    "TargetShape": shape,
                    "detection_time": detection_time,
                    "target_motion": motion,
                    "target_step_length": target_step_length,
                }
            )
    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    ("moving_motion", "target_step_length"),
    (
        ("moving_continuous_unit_run", 1.0),
        ("uniform_random_relocation", 0.0),
    ),
)
def test_notebook_builds_target_motion_comparison_pdf(
    tmp_path: Path, moving_motion: str, target_step_length: float
):
    namespace = _target_motion_plotting_namespace()
    fixed_path = tmp_path / "fixed.csv"
    moving_path = tmp_path / "moving.csv"
    ratio_path = tmp_path / "ratios.csv"
    figure_path = tmp_path / "heatmap.pdf"
    _experiment_rows("fixed", 0.0, (2.0, 4.0)).to_csv(
        fixed_path, index=False
    )
    _experiment_rows(
        moving_motion, target_step_length, (1.0, 2.0)
    ).to_csv(moving_path, index=False)

    ratios = namespace["build_target_motion_comparison"](
        fixed_path=fixed_path,
        moving_path=moving_path,
        moving_motion=moving_motion,
        output_csv=ratio_path,
        output_pdf=figure_path,
    )

    assert len(ratios) == 3
    assert ratios["comparison_motion"].eq(moving_motion).all()
    assert ratios["ratio_definition"].eq("mean_moving/mean_fixed").all()
    assert ratios["detection_time_ratio"].eq(0.5).all()
    assert ratio_path.is_file()
    assert figure_path.read_bytes().startswith(b"%PDF-")


def test_notebook_target_motion_heatmap_requires_pdf(tmp_path: Path):
    namespace = _target_motion_plotting_namespace()

    with pytest.raises(ValueError, match="must be a PDF"):
        namespace["plot_target_motion_ratio_heatmap"](
            pd.DataFrame(),
            "moving_continuous_unit_run",
            tmp_path / "heatmap.png",
        )
