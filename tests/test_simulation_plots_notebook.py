from __future__ import annotations

import json

from conftest import REPOSITORY_ROOT


NOTEBOOK_PATH = REPOSITORY_ROOT / "simulation-plots.ipynb"


def _notebook():
    return json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))


def _code_source() -> str:
    return "\n\n".join(
        "".join(cell["source"])
        for cell in _notebook()["cells"]
        if cell["cell_type"] == "code"
    )


def test_simulation_plot_notebook_contains_code_only():
    notebook = _notebook()
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            assert cell["execution_count"] is None
            assert cell["outputs"] == []

    encoded = NOTEBOOK_PATH.read_text(encoding="utf-8")
    assert '"image/png"' not in encoded
    assert '"image/jpeg"' not in encoded
    assert '"text/html"' not in encoded


def test_simulation_plot_notebook_compiles_and_reads_merged_outputs():
    source = _code_source()
    compile(source, str(NOTEBOOK_PATH), "exec")

    expected_configs = (
        "detection_time_cauchy_projected_surface.conf",
        "detection_time_mu.conf",
        "detection_time_fixed_volume.conf",
        "detection_time_fixed_surface.conf",
        "detection_time_small_delta_large_mu.conf",
        "detection_time_small_delta_small_mu.conf",
        "detection_time_ratio_fixed_surface.conf",
    )
    for config_name in expected_configs:
        config_path = (
            REPOSITORY_ROOT / "experiments" / "configs" / config_name
        )
        values = {}
        for raw_line in config_path.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
        output_path = "{}{}".format(
            values["save_directory"].removeprefix("./"),
            values["file_name"],
        )
        assert output_path in source

    assert (
        "glob.glob('./results/detection_time_ratio_fixed_surface/*')"
        not in source
    )
