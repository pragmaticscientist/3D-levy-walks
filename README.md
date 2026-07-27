# LTDB Lévy-walk simulations

This repository contains code and reproducible configuration, but no source
trajectories, generated tables, figures, or simulation results. The tracked
`simulation-plots.ipynb` notebook contains plotting code only: its outputs and
execution counts are cleared.

## Local data setup

Place the LTDB ground-truth trajectory files in:

```text
datasets/GT_TRACKS/
```

The YAML pipeline configurations expect the biological annotation table at:

```text
configs/ltdb_metadata.csv
```

Both locations are ignored by Git. The trajectory format is the one described
by Pizzagalli et al. (2018), DOI `10.1038/sdata.2018.129`. If the annotation
table is absent, the inspect stage creates a template from the installed
trajectory files:

```bash
python -m ltdb_levy inspect \
  --config configs/equal_run_weight_run_eligible_conditions.yaml
```

Complete and verify that template before preprocessing. Generated artifacts are
written beneath `outputs/` or `results/`; both are ignored.

## Arrest-scale experiment

Install the package in a virtual environment, including test dependencies:

```bash
python -m pip install -e '.[test]'
```

Build the source preprocessing artifacts required by the arrest-scale replay:

```bash
python -m ltdb_levy preprocess \
  --config configs/equal_run_weight_run_eligible_conditions.yaml
```

Run the experiment from its tracked configuration:

```bash
python -m ltdb_levy.arrest_scale_cli run \
  --config configs/arrest_scale.yaml
```

The equivalent repository wrapper is:

```bash
python scripts/run_arrest_scale_experiment.py run \
  --config configs/arrest_scale.yaml
```

The experiment validates the source artifact hashes, writes all tables and
figures under the configured `output_root`, and resumes an already completed
matching run by default.

## Regenerating plots

Plot generation is independent of the raw trajectory files once the experiment
tables exist:

```bash
python -m ltdb_levy.arrest_scale_cli plot \
  --config configs/arrest_scale.yaml
```

Use `--output-root PATH` with either command to override the configured artifact
directory. Run `python -m ltdb_levy.arrest_scale_cli --help` for the complete
set of reproducibility overrides.

## C-simulation plots

The non-empirical figures are reproduced by the output-free
`simulation-plots.ipynb` notebook. First run the experiment scripts needed for
the desired figure from the repository root, for example:

```bash
bash experiments/detection_time_cauchy_projected_surface.sh
bash experiments/detection_time_mu.sh
bash experiments/detection_time_fixed_volume.sh
bash experiments/detection_time_fixed_surface.sh
bash experiments/detection_time_small_delta_large_mu.sh
bash experiments/detection_time_small_delta_small_mu.sh
bash experiments/detection_time_ratio_fixed_surface.sh
bash experiments/detection_time_fixed_target_projected_surface.sh
bash experiments/detection_time_moving_target_projected_surface.sh
bash experiments/detection_time_uniform_relocation_projected_surface.sh
```

These commands compile the simulator, create their ignored result directories,
and write a merged CSV at the path used by the notebook. Launch Jupyter from
the repository root and execute:

```bash
jupyter lab simulation-plots.ipynb
```

The generated PDFs are written beneath `plots/pdf_figures/`, which is also
ignored by Git.
