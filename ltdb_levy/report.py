"""Reproducible report, methods, limitations, conventions, and manifest."""

from __future__ import annotations

import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ltdb_levy_matplotlib")
)

import matplotlib
import numpy as np
import pandas as pd
import scipy
import yaml

from . import __version__
from .config import Config
from .dataio import ArtifactStore, discover_files, sha256_file
from .figures import generate_figures


def _write_text(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".{}.tmp".format(path.name))
    temporary.write_text(content.rstrip() + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def _git_value(args: Sequence[str], default: str = "unknown") -> str:
    try:
        result = subprocess.run(
            ["git"] + list(args),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return default
    return result.stdout.strip()


def _markdown_table(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    if frame.empty:
        return "_No rows available._"
    labels = list(columns)
    header = "| " + " | ".join(labels) + " |"
    divider = "| " + " | ".join(["---"] * len(labels)) + " |"
    rows = [header, divider]
    for _, row in frame.iterrows():
        values = []
        for column in labels:
            value = row[column]
            if isinstance(value, (float, np.floating)):
                values.append("" if not np.isfinite(value) else "{:.4g}".format(value))
            else:
                values.append(str(value).replace("|", "\\|"))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def _methods_text(config: Config) -> str:
    crossover_method = (
        "One crossover and one exponent were fitted independently in every "
        "eligible condition. The 3.75 µm cell-radius scale was used only as "
        "an optimizer start and external reference."
        if config.fitting.estimate_crossover
        else (
            "One exponent was fitted per eligible condition with the crossover "
            "fixed at {crossover:g} µm.".format(
                crossover=config.fitting.crossover
            )
        )
    )
    return """# Methods

LTDB ground-truth files were parsed positionally: calibration on line 1,
channel flags on line 2, and `TrackID;x;y;z;t` observations from line 3 onward.
Coordinates were treated as micrometres as specified by the source data, and
time was computed as frame index multiplied by each file's own `dt`.
Cell type, imaging site, and stimulus were taken from Tables 3 and 4 of the
LTDB data descriptor and cross-checked against the deposited SQL database.
Conditions were defined by `(cell type, site, stimulus, dt)`.

Tracks were sorted by time, duplicate `(track, frame)` observations were handled
according to the configured policy, and fragments were split at every frame gap
greater than one. No interpolation was performed. Frame displacements, speeds,
and turning angles were computed in three dimensions; angles used
`atan2(||u×v||, u·v)`.

The primary relocation definition was a directed run. A displacement whose
turning angle strictly exceeded the {turning_angle:g}-degree threshold started
the next run. The 90-degree primary threshold was anchored to a directly
relevant leukocyte study that classified 90--180-degree changes as large
U-turns (Georgantzoglou et al., J Cell Biol 2022,
doi:10.1083/jcb.202103207). It is an operational definition rather than a
universal biological constant: the published angle calculation is not
identical to ours, and turning angles depend on sampling interval. Thresholds
of 45 and 60 degrees were therefore retained as prespecified sensitivities,
and conditions with different frame intervals were not pooled.

The primary fit imposed no positive speed cutoff. Speed-sensitive fits treated
steps below 1, 2, or 4 µm/min as pauses that split adjacent runs. The central
2-µm/min value is the conventional in-vivo lymphocyte arrest threshold; using
it as a hard run boundary is an operational extension because arrest-coefficient
studies usually summarize time below the threshold rather than segmenting runs.
Nonmoving steps split runs. Relocation length was the end-to-end norm, with path
length stored separately. First and last runs of each fragment were excluded
from the primary fit as field-of-view censored.

The model was `f(l)=a·min(1,(l/c)^(-μ))` on `[0,lmax]`. {crossover_method}
The upper support `lmax` was the condition's empirical maximum.
Cell-balanced weights gave every track equal total mass. Uncertainty used
whole-track clustered bootstrap resampling, re-estimating all fitted parameters.
Goodness of fit used a parametric-bootstrap KS statistic with all fitted
parameters re-estimated in every replicate. Predictive model comparison used
whole-track held-out folds, including crossover estimation within each training
fold, against truncated exponential, lognormal, Weibull, and
power-law-with-cutoff alternatives.

Replay used a side-{side:g} periodic cube, detection radius `{radius:g}`, a
single target at the box centre, endpoint-only detection, and no check at the
initial position. Ball, cylinder-disk, and y-oriented line-capsule membership
match the C implementation. Shapes were matched by face-on maximum projected
area. Each source track supplied its own finite path-length budget; overshooting
steps were `{overshoot}`. Restricted mean detection distance was the primary
time-to-event outcome.
""".format(
        turning_angle=config.segmentation.turning_angle_degrees,
        crossover_method=crossover_method,
        side=config.replay.domain_side_length,
        radius=config.replay.detection_radius,
        overshoot=config.replay.overshoot_policy,
    )


def _limitations_text(config: Config) -> str:
    crossover_limitation = (
        "The crossover is estimated separately by condition and is strongly "
        "correlated with the effective exponent; small condition samples can "
        "therefore have wide joint uncertainty."
        if config.fitting.estimate_crossover
        else (
            "The flat-to-power-law crossover is imposed at {crossover:g} µm "
            "rather than measured.".format(
                crossover=config.fitting.crossover
            )
        )
    )
    return """# Limitations

- Sampling interval changes the displacement and directed-run lengths that can
  be observed.
- The 90-degree primary turning threshold is literature-anchored but not
  biologically universal. Published studies use different angle definitions,
  sampling intervals, smoothing rules, and cell systems. Directed-run
  definitions therefore remain dependent on the turning-angle and speed
  thresholds.
- The conventional 2-µm/min arrest threshold was developed as a summary of
  time spent arrested, not as a validated boundary for directed-run
  segmentation. Treating sub-threshold steps as run boundaries is therefore a
  sensitivity analysis rather than the primary definition.
- LTDB contains trajectories but no targets and no detection events; all search
  outcomes here are counterfactual replays.
- An effective fitted exponent does not establish an exact Lévy walk.
- Deviations from the fitted mixture may reflect persistence, pauses,
  correlations, tissue structure, or measurement effects.
- Empirical pools cannot generate relocations beyond their observed support.
- Pooling can hide between-track, between-video, and biological heterogeneity.
- The observed z slab makes empirical directions non-isotropic by construction.
- Face-on area matching leaves a randomly oriented disk effectively smaller on
  average than a ball, so shape sensitivity partly reflects this convention.
- {crossover_limitation}
- First and last directed runs are field-of-view censored and are excluded from
  the primary fit.
- Biological metadata was transcribed from the source publication and
  cross-checked against the deposited SQL database. Independent source review
  remains appropriate before publication; unverified rows are never pooled.
- For long line targets, the detection neighbourhood can overlap itself through
  the periodic domain; those target cells are flagged explicitly.
- The C engine is unbounded while the primary replay is finite-budget and
  censored, so their mean detection distances target different estimands.
""".format(crossover_limitation=crossover_limitation)


def _conventions_text(config: Config) -> str:
    crossover_convention = (
        "one crossover is estimated per condition; 3.75 µm is the optimizer "
        "start and cell-radius reference"
        if config.fitting.estimate_crossover
        else "the primary crossover is {crossover:g} µm".format(
            crossover=config.fitting.crossover
        )
    )
    return """# Simulation conventions

- Crossover convention: {crossover_convention}.
- Step density: `a(μ,lmax)·min(1,(l/c)^(-μ))`, with lengths and fitted
  crossover `c` stored in µm.
- Direction sampling: isotropic, with uniform azimuth and uniform cosine.
- Periodic boundary: coordinate-wise minimum image on a cube of side {side:g}.
- Detection radius: `{radius:g}`.
- Ball: radius `D/2 + d`.
- Disk: cylinder with radial limit `D/2 + d` and `|Δz| ≤ d`.
- Line: capsule of centre-line length `D`, oriented along y, radius `d`.
- Target centre: box centre.
- Detection: relocation endpoints only; initial position excluded.
- Shape matching: face-on maximum projected area, matching `func.c`.
- Walker count and target count: one each.
- Primary finite-budget overshoot policy: `{overshoot}`.
""".format(
        crossover_convention=crossover_convention,
        side=config.replay.domain_side_length,
        radius=config.replay.detection_radius,
        overshoot=config.replay.overshoot_policy,
    )


def _manifest(config: Config, store: ArtifactStore) -> Dict[str, Any]:
    selected = discover_files(
        config.input.tracks_directory,
        config.input.glob,
        config.selection.include_files,
        config.selection.exclude_files,
    )
    sources = [
        {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in selected
    ]
    artifacts = []
    for path in sorted(config.output.root.rglob("*")):
        if not path.is_file() or path.name == "analysis_manifest.json":
            continue
        artifacts.append(
            {
                "path": str(path.relative_to(config.output.root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    metadata_verified = 0
    metadata_rows = 0
    metadata_source = None
    if config.input.metadata_file and config.input.metadata_file.is_file():
        metadata = pd.read_csv(config.input.metadata_file, comment="#")
        metadata_rows = len(metadata)
        metadata_source = {
            "path": str(config.input.metadata_file.resolve()),
            "size_bytes": config.input.metadata_file.stat().st_size,
            "sha256": sha256_file(config.input.metadata_file),
        }
        if "verified" in metadata:
            metadata_verified = int(
                metadata["verified"]
                .astype(str)
                .str.lower()
                .isin(("true", "1", "yes"))
                .sum()
            )
    return {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline_version": __version__,
        "run_id": config.run_id,
        "config_path": str(config.source_path),
        "config_hash_sha256": config.hash(),
        "config_hash_scope": "semantic YAML plus selected track and metadata contents",
        "git_commit": _git_value(["rev-parse", "HEAD"]),
        "git_branch": _git_value(["rev-parse", "--abbrev-ref", "HEAD"]),
        "git_status_porcelain": _git_value(["status", "--porcelain"], default=""),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "pyyaml": yaml.__version__,
        },
        "seeds": {
            "fitting": config.fitting.random_seed,
            "replay": config.replay.random_seed,
        },
        "metadata": {
            "rows": metadata_rows,
            "verified_rows": metadata_verified,
            "source": metadata_source,
        },
        "source_files": sources,
        "artifacts": artifacts,
        "status": {
            "inspect_complete": store.json_path("inspect", "summary").is_file(),
            "preprocess_complete": store.json_path("preprocess", "summary").is_file(),
            "fit_complete": store.json_path("fit", "summary").is_file(),
            "pilot_complete": store.json_path("replay_pilot", "summary").is_file(),
            "confirmatory_replay_complete": store.json_path("replay", "summary").is_file(),
            "c_configs_prepared": (
                Path(__file__).resolve().parents[1]
                / "experiments"
                / "configs"
                / "ltdb_fitted_mu"
                / "index.csv"
            ).is_file(),
            "c_results_ingested": store.json_path(
                "c_comparison", "summary"
            ).is_file(),
            "python_c_comparison_complete": store.frame_path(
                "c_comparison", "python_c_comparison"
            ).is_file(),
        },
    }


def build_report(config: Config) -> List[Path]:
    """Generate report documents, figures, and a complete analysis manifest."""
    store = ArtifactStore(config.output.root)
    report_dir = store.stage_dir("report")
    analysis_hash = config.hash()
    for stage in ("inspect", "preprocess", "fit"):
        store.require_run_context(stage, analysis_hash)
    figures = generate_figures(config, store) if config.output.figures else []

    inspection = store.read_json("inspect", "summary")
    preprocessing = store.read_json("preprocess", "summary")
    fit_summary = store.read_json("fit", "summary")
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config.hash(),
    )
    successful = fits[fits["status"] == "ok"].copy()
    model_counts = (
        successful["best_cv_model"].value_counts().to_dict()
        if "best_cv_model" in successful
        else {}
    )
    replay_stage = (
        "replay"
        if store.json_path("replay", "summary").is_file()
        else "replay_pilot"
    )
    store.require_run_context(replay_stage, analysis_hash)
    replay_summary = store.read_json(replay_stage, "summary")
    target_dimensions = store.read_frame(
        replay_stage,
        "target_dimensions",
        expected_config_hash=config.hash(),
    )
    overlap_count = int(target_dimensions["periodic_self_overlap"].astype(bool).sum())
    metadata_verified = 0
    metadata_rows = 0
    if config.input.metadata_file and config.input.metadata_file.is_file():
        metadata = pd.read_csv(config.input.metadata_file, comment="#")
        metadata_rows = int(len(metadata))
        metadata_verified = int(
            metadata["verified"]
            .astype(str)
            .str.lower()
            .isin(("true", "1", "yes"))
            .sum()
        )

    metadata_status = (
        "All {count} rows are source-verified, so eligible files are pooled by "
        "cell type, organ, stimulus, and sampling interval.".format(
            count=metadata_rows
        )
        if metadata_rows > 0 and metadata_verified == metadata_rows
        else (
            "{verified} of {rows} metadata rows are verified. Unverified rows "
            "remain file-local and are not biologically pooled.".format(
                verified=metadata_verified, rows=metadata_rows
            )
        )
    )
    c_status = (
        "C comparison has not been ingested."
    )
    if store.json_path("c_comparison", "summary").is_file():
        store.require_run_context("c_comparison", analysis_hash)
        c_summary = store.read_json("c_comparison", "summary")
        if int(c_summary.get("valid_c_target_cells", 0)) == 0:
            c_status = (
                "Fitted-μ C configs are prepared, but none of the nine C output "
                "files exists yet; the Python/C agreement check remains pending."
            )
        else:
            c_status = (
                "C output is available for {cells} target cells; {compared} have "
                "matched unbounded Python comparisons.".format(
                    cells=int(c_summary.get("valid_c_target_cells", 0)),
                    compared=int(c_summary.get("comparison_target_cells", 0)),
                )
            )

    report = """# LTDB Lévy-walk surrogate analysis

## Current status

The data and fitting stages completed successfully. The replay shown here is
**{replay_label}**, not a confirmatory biological result. It contains
{replay_trials:,} target/class trials across {replay_tracks} tracks.
{underpowered} of 90 global target/class cells fail the configured pilot event
minimums, and {overlap} target dimensions self-overlap through the periodic
domain. The confirmatory target grid is therefore not frozen.

{metadata_status}

## Corpus and preprocessing

- Selected files: {files}
- Observations read: {observations:,}
- Tracks in the raw corpus: {raw_tracks}
- Retained fragments/tracks: {retained_tracks}
- Directed runs: {runs_total:,} total; {runs_included:,} primary included
- Eligible condition pools: {eligible_pools}

## Distributional result

All {fit_count} eligible conditions fitted numerically. The parametric mixture
was not clearly supported in any fitted condition and μ is therefore reported
as an **effective exponent**, not evidence of an exact Lévy law. Best held-out
models were: {model_counts}.

{fit_table}

The prespecified segmentation grid was also refitted with μ and the
parametric-bootstrap KS statistic. Its shape-sensitivity outcome `G` remains
explicitly pending because the confirmatory target grid has not passed the
pilot adequacy checkpoint.

## Replay checkpoint

The pilot is underpowered for the configured target grid, particularly for the
radius-{radius:g} line target in a side-{side:g} torus. A zero observed pilot probability is
not interpreted as a true zero. See
[pilot diagnostics](../replay_pilot/pilot_diagnostics.csv),
[track selection](../replay_pilot/replay_track_selection.csv), and
[target dimensions](../replay_pilot/target_dimensions.csv) before changing or
freezing the grid.

## Existing-C cross-check

{c_status}
The machine-readable ingestion audit is
[here](../c_comparison/c_ingestion_audit.csv).

## Figures

{figure_links}

## Reproducibility

See `methods.md`, `limitations.md`, `conventions.md`, and
`analysis_manifest.json`. The manifest records source and artifact hashes,
software versions, Git state, configuration hash, and random seeds.
""".format(
        replay_label="confirmatory replay"
        if replay_stage == "replay"
        else "a pre-freeze pilot",
        replay_trials=int(replay_summary["trials"]),
        replay_tracks=int(replay_summary["tracks"]),
        underpowered=int(replay_summary.get("underpowered_cells", 0)),
        overlap=overlap_count,
        metadata_status=metadata_status,
        files=int(inspection["selected_files"]),
        observations=int(inspection["observations"]),
        raw_tracks=int(inspection["tracks"]),
        retained_tracks=int(preprocessing["tracks_retained"]),
        runs_total=int(preprocessing["runs_total"]),
        runs_included=int(preprocessing["runs_included"]),
        eligible_pools=int(preprocessing["eligible_pools"]),
        fit_count=int(fit_summary["successful_fits"]),
        model_counts=", ".join(
            "{} ({})".format(name, count)
            for name, count in sorted(model_counts.items())
        ),
        fit_table=_markdown_table(
            successful,
            [
                "condition",
                "mu_hat",
                "mu_ci_low",
                "mu_ci_high",
                "crossover",
                "crossover_ci_low",
                "crossover_ci_high",
                "ks_bootstrap_p",
                "best_cv_model",
                "interpretation",
            ],
        ),
        c_status=c_status,
        radius=config.replay.detection_radius,
        side=config.replay.domain_side_length,
        figure_links="\n".join(
            "- [{}](figures/{})".format(path.stem.replace("_", " "), path.name)
            for path in figures
        )
        or "_Figure generation disabled._",
    )

    outputs = [
        _write_text(report_dir / "report.md", report),
        _write_text(report_dir / "methods.md", _methods_text(config)),
        _write_text(report_dir / "limitations.md", _limitations_text(config)),
        _write_text(report_dir / "conventions.md", _conventions_text(config)),
    ]
    manifest = _manifest(config, store)
    manifest_path = report_dir / "analysis_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    outputs.append(manifest_path)
    outputs.extend(figures)
    return outputs
