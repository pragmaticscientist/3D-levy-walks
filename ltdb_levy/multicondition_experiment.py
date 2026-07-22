"""Equal-run-weight replay orchestration across all eligible conditions."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactError, ArtifactStore
from .dist.mixture import cdf as mixture_cdf
from .logging_utils import get_logger
from .nk_experiment import (
    FINITE_CLASSES,
    FINITE_CONTRASTS,
    FINITE_PLOT_ORDER,
    UNBOUNDED_CLASSES,
    FocusedReplaySettings,
    run_focused_replay,
)
from .targets import dimension_from_projected_area

logger = get_logger(__name__)

REPLAY_SIGNATURE_VERSION = 1
REPLAY_ALGORITHM_VERSION = "endpoint_replay_v2_ordered_length_rotation"


@dataclass(frozen=True)
class AllConditionsReplaySettings:
    """Runtime controls for the all-condition replay."""

    projected_areas: Tuple[float, ...] = (4.0, 8.0, 16.0, 24.0)
    finite_replays_per_track: int = 10_000
    unbounded_trials_per_cell: int = 2_000
    workers: int = 32
    finite_batch_size: int = 1_000
    unbounded_chunk_trials: int = 25
    unbounded_step_block: int = 256
    random_seed: int = 24680
    resume: bool = True

    def validate(self) -> None:
        if not self.projected_areas:
            raise ValueError("At least one projected area is required")
        if any(area <= np.pi for area in self.projected_areas):
            raise ValueError("Every projected area must exceed pi for d=1")
        if self.finite_replays_per_track < 1:
            raise ValueError("finite_replays_per_track must be positive")
        if self.unbounded_trials_per_cell < 1:
            raise ValueError("unbounded_trials_per_cell must be positive")
        if self.workers < 1:
            raise ValueError("workers must be positive")


def _parse_condition(condition: str) -> Dict[str, object]:
    values: Dict[str, object] = {}
    for part in condition.split(" | "):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key.strip()] = value.strip()
    if "dt_seconds" in values:
        values["dt_seconds"] = float(str(values["dt_seconds"]))
    return values


def _short_token(value: str) -> str:
    replacements = {
        "Natural Killer Cells": "nk",
        "Neutrophils": "neutrophils",
        "T Cells": "t",
        "B Cells": "b",
        "popliteal lymph node": "pln",
        "Influenza Vaccine": "influenza",
        "Ovalbumin": "ova",
        "Vaccinia Virus": "vaccinia",
        "HIV-infected humanized T cell": "hiv",
        "Steady State": "steady",
    }
    return replacements.get(value, value)


def _condition_id(condition: str) -> str:
    parsed = _parse_condition(condition)
    tokens = [
        _short_token(str(parsed.get("cell_type", "condition"))),
        _short_token(str(parsed.get("organ", ""))),
        _short_token(str(parsed.get("stimulus", ""))),
        "dt{:g}".format(float(parsed.get("dt_seconds", 0.0))),
    ]
    slug = re.sub(r"[^a-z0-9]+", "_", "_".join(tokens).lower()).strip("_")
    digest = hashlib.sha256(condition.encode("utf-8")).hexdigest()[:8]
    return "{}_{}".format(slug or "condition", digest)


def _condition_label(condition: str) -> str:
    parsed = _parse_condition(condition)
    cell = str(parsed.get("cell_type", ""))
    cell = {
        "Natural Killer Cells": "NK",
        "Neutrophils": "Neutrophils",
        "T Cells": "T",
        "B Cells": "B",
    }.get(cell, cell)
    organ = {
        "popliteal lymph node": "PLN",
        "spleen": "spleen",
    }.get(str(parsed.get("organ", "")), str(parsed.get("organ", "")))
    stimulus = {
        "Influenza Vaccine": "influenza",
        "Ovalbumin": "OVA",
        "Vaccinia Virus": "vaccinia",
        "HIV-infected humanized T cell": "HIV model",
        "Steady State": "steady state",
    }.get(str(parsed.get("stimulus", "")), str(parsed.get("stimulus", "")))
    return "{} · {} · {} · {:g} s".format(
        cell,
        organ,
        stimulus,
        float(parsed.get("dt_seconds", 0.0)),
    )


def _condition_fields(condition: str) -> Dict[str, object]:
    parsed = _parse_condition(condition)
    return {
        "condition": condition,
        "condition_id": _condition_id(condition),
        "condition_label": _condition_label(condition),
        "cell_type": parsed.get("cell_type", ""),
        "organ": parsed.get("organ", ""),
        "stimulus": parsed.get("stimulus", ""),
        "dt_seconds": parsed.get("dt_seconds", np.nan),
    }


def _condition_seed(master_seed: int, condition: str) -> int:
    digest = hashlib.sha256(
        "{}|{}".format(master_seed, condition).encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "little") % (2**32 - 1)


def _canonical_signature(payload: Mapping[str, object]) -> str:
    """Return a deterministic SHA-256 signature for semantic replay inputs."""
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _stored_signature_matches(
    stored_settings: object,
    expected_signature: str,
) -> bool:
    """Reject unsigned, obsolete, or mismatched condition outputs."""
    if not isinstance(stored_settings, dict):
        return False
    return bool(
        stored_settings.get("replay_signature_version")
        == REPLAY_SIGNATURE_VERSION
        and stored_settings.get("replay_signature") == expected_signature
    )


def _artifact_sha256(
    store: ArtifactStore,
    stage: str,
    name: str,
) -> str:
    """Read one already-validated artifact digest from its sidecar."""
    meta_path = store.frame_path(stage, name).with_suffix(".meta.json")
    if not meta_path.is_file():
        raise ArtifactError("Artifact sidecar is missing: {}".format(meta_path))
    with meta_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    digest = metadata.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64:
        raise ArtifactError(
            "Artifact sidecar has no valid SHA-256 digest: {}".format(
                meta_path
            )
        )
    return digest


def _replay_signature_payload(
    config: Config,
    settings: AllConditionsReplaySettings,
    condition: str,
    fit: Mapping[str, object],
    source_artifact_hashes: Mapping[str, str],
) -> Dict[str, object]:
    """Collect every input that can change a condition replay result."""
    length_unit_um = float(fit["crossover"])
    return {
        "signature_version": REPLAY_SIGNATURE_VERSION,
        "algorithm_version": REPLAY_ALGORITHM_VERSION,
        "config_hash": config.hash(),
        "source_artifact_hashes": dict(source_artifact_hashes),
        "condition": condition,
        "weighting": config.empirical_pool.primary_weighting,
        "normalization": "fitted_crossover",
        "length_unit_um": length_unit_um,
        "fitted_mu": float(fit["mu_hat"]),
        "fitted_crossover_model_units": 1.0,
        "fitted_lmax_model_units": float(fit["lmax"]) / length_unit_um,
        "torus_width_rule": "2*fitted_lmax_model_units",
        "detection_radius_model_units": float(
            config.replay.detection_radius
        ),
        "shapes": list(config.replay.shapes),
        "projected_areas": [float(value) for value in settings.projected_areas],
        "finite_classes": list(FINITE_CLASSES),
        "unbounded_classes": list(UNBOUNDED_CLASSES),
        "finite_replays_per_track": settings.finite_replays_per_track,
        "unbounded_trials_per_cell": settings.unbounded_trials_per_cell,
        # These chunking controls change how a fixed random stream is assigned
        # to trials, so they are semantic rather than performance-only.
        "finite_batch_size": settings.finite_batch_size,
        "unbounded_chunk_trials": settings.unbounded_chunk_trials,
        "unbounded_step_block": settings.unbounded_step_block,
        "random_seed": _condition_seed(settings.random_seed, condition),
        "finite_bootstraps": config.fitting.clustered_bootstraps,
        "overshoot_policy": "discard",
        "check_initial_position": False,
        "figure_format": "pdf",
    }


def _expected_replay_row_counts(
    n_tracks: int,
    n_pool_vectors: int,
    n_shapes: int,
    settings: AllConditionsReplaySettings,
) -> Dict[str, int]:
    """Return exact per-condition row counts implied by the runtime settings."""
    n_areas = len(settings.projected_areas)
    finite_cells = len(FINITE_CLASSES) * n_shapes * n_areas
    unbounded_cells = len(UNBOUNDED_CLASSES) * n_shapes * n_areas
    return {
        "track_inputs": int(n_tracks),
        "empirical_pool_inputs": int(n_pool_vectors),
        "finite_track_summaries": int(n_tracks * finite_cells),
        "finite_summary": int(finite_cells),
        "finite_contrasts": int(
            len(FINITE_CONTRASTS) * n_shapes * n_areas
        ),
        "unbounded_trials": int(
            unbounded_cells * settings.unbounded_trials_per_cell
        ),
        "unbounded_summary": int(unbounded_cells),
        "unbounded_contrasts": int(n_shapes * n_areas),
        "targets": int(n_shapes * n_areas),
    }


def _ordered_shape_area_pairs(
    frame: pd.DataFrame,
) -> List[Tuple[str, float]]:
    """Order plot cells by canonical shape and then numeric target area."""
    present = {
        (str(shape), float(area))
        for shape, area in frame[["shape", "projected_area"]].itertuples(
            index=False, name=None
        )
    }
    canonical = ("Ball", "Disk", "Line")
    extra_shapes = sorted(
        {shape for shape, _ in present}.difference(canonical)
    )
    shape_order = [
        shape
        for shape in tuple(canonical) + tuple(extra_shapes)
        if any(pair[0] == shape for pair in present)
    ]
    return [
        (shape, area)
        for shape in shape_order
        for area in sorted(
            value for candidate, value in present if candidate == shape
        )
    ]


def _validate_replay_frames(
    frames: Mapping[str, pd.DataFrame],
    summary: Mapping[str, object],
    settings: AllConditionsReplaySettings,
    n_tracks: int,
    n_pool_vectors: int,
    shapes: Sequence[str],
) -> None:
    """Fail closed when a purportedly complete condition replay is inconsistent."""
    expected = _expected_replay_row_counts(
        n_tracks,
        n_pool_vectors,
        len(shapes),
        settings,
    )
    missing = sorted(set(expected).difference(frames))
    if missing:
        raise ValueError(
            "condition replay is missing frames: {}".format(
                ", ".join(missing)
            )
        )
    for name, expected_rows in expected.items():
        actual_rows = len(frames[name])
        if actual_rows != expected_rows:
            raise ValueError(
                "{} has {} rows; expected {}".format(
                    name, actual_rows, expected_rows
                )
            )

    expected_areas = sorted(float(value) for value in settings.projected_areas)
    expected_shapes = set(str(value) for value in shapes)

    def validate_grid(
        frame: pd.DataFrame,
        classes: Optional[Sequence[str]] = None,
    ) -> None:
        actual_shapes = set(frame["shape"].astype(str))
        actual_areas = sorted(
            frame["projected_area"].astype(float).unique()
        )
        if (
            actual_shapes != expected_shapes
            or len(actual_areas) != len(expected_areas)
            or not np.allclose(
                actual_areas, expected_areas, rtol=0.0, atol=1e-12
            )
        ):
            raise ValueError("condition replay has an unexpected target grid")
        if classes is not None and set(
            frame["trajectory_class"].astype(str)
        ) != set(classes):
            raise ValueError(
                "condition replay has unexpected trajectory classes"
            )

    track_inputs = frames["track_inputs"]
    if track_inputs["track_uid"].astype(str).nunique() != n_tracks:
        raise ValueError("track_inputs does not contain one row per track")

    pool_inputs = frames["empirical_pool_inputs"]
    if "uniform_run_sampling_probability" in pool_inputs:
        probabilities = pool_inputs[
            "uniform_run_sampling_probability"
        ].to_numpy(dtype=float)
        if (
            np.any(~np.isfinite(probabilities))
            or not np.allclose(
                probabilities,
                np.full(n_pool_vectors, 1.0 / n_pool_vectors),
                rtol=0.0,
                atol=1e-12,
            )
        ):
            raise ValueError(
                "empirical_pool_inputs is not uniformly run-weighted"
            )

    finite_tracks = frames["finite_track_summaries"]
    validate_grid(finite_tracks, FINITE_CLASSES)
    finite_keys = [
        "track_uid",
        "trajectory_class",
        "shape",
        "projected_area",
    ]
    if finite_tracks.duplicated(finite_keys).any():
        raise ValueError("finite_track_summaries has duplicate cells")
    if not np.all(
        finite_tracks["n_trials"].to_numpy(dtype=int)
        == settings.finite_replays_per_track
    ):
        raise ValueError("finite_track_summaries has an incorrect trial count")

    finite_summary = frames["finite_summary"]
    validate_grid(finite_summary, FINITE_CLASSES)
    if finite_summary.duplicated(
        ["trajectory_class", "shape", "projected_area"]
    ).any():
        raise ValueError("finite_summary has duplicate cells")
    if not np.all(
        finite_summary["n_tracks"].to_numpy(dtype=int) == n_tracks
    ) or not np.all(
        finite_summary["n_trials"].to_numpy(dtype=int)
        == n_tracks * settings.finite_replays_per_track
    ):
        raise ValueError("finite_summary has inconsistent sample sizes")

    finite_contrasts = frames["finite_contrasts"]
    validate_grid(finite_contrasts)
    expected_contrasts = {
        (left, right) for left, right in FINITE_CONTRASTS
    }
    actual_contrasts = {
        (str(left), str(right))
        for left, right in finite_contrasts[
            ["left_class", "right_class"]
        ].itertuples(index=False, name=None)
    }
    if actual_contrasts != expected_contrasts or finite_contrasts.duplicated(
        ["left_class", "right_class", "shape", "projected_area"]
    ).any():
        raise ValueError("finite_contrasts has an unexpected contrast grid")

    unbounded_trials = frames["unbounded_trials"]
    validate_grid(unbounded_trials, UNBOUNDED_CLASSES)
    unbounded_keys = [
        "trajectory_class",
        "shape",
        "projected_area",
        "trial_index",
    ]
    if unbounded_trials.duplicated(unbounded_keys).any():
        raise ValueError("unbounded_trials has duplicate trial identifiers")
    grouped_trials = unbounded_trials.groupby(
        ["trajectory_class", "shape", "projected_area"],
        sort=False,
    )
    if not np.all(
        grouped_trials.size().to_numpy(dtype=int)
        == settings.unbounded_trials_per_cell
    ):
        raise ValueError("unbounded_trials has an incorrect cell size")
    trial_min = grouped_trials["trial_index"].min().to_numpy(dtype=int)
    trial_max = grouped_trials["trial_index"].max().to_numpy(dtype=int)
    if not np.all(trial_min == 0) or not np.all(
        trial_max == settings.unbounded_trials_per_cell - 1
    ):
        raise ValueError("unbounded_trials has an incomplete trial index")

    unbounded_summary = frames["unbounded_summary"]
    validate_grid(unbounded_summary, UNBOUNDED_CLASSES)
    if unbounded_summary.duplicated(
        ["trajectory_class", "shape", "projected_area"]
    ).any() or not np.all(
        unbounded_summary["n_trials"].to_numpy(dtype=int)
        == settings.unbounded_trials_per_cell
    ):
        raise ValueError("unbounded_summary has inconsistent cells")

    unbounded_contrasts = frames["unbounded_contrasts"]
    validate_grid(unbounded_contrasts)
    if unbounded_contrasts.duplicated(["shape", "projected_area"]).any():
        raise ValueError("unbounded_contrasts has duplicate cells")

    targets = frames["targets"]
    validate_grid(targets)
    if targets.duplicated(["shape", "projected_area"]).any():
        raise ValueError("targets has duplicate cells")
    if "fits_without_periodic_self_overlap" in targets and not targets[
        "fits_without_periodic_self_overlap"
    ].astype(bool).all():
        raise ValueError("one or more replay targets self-overlaps")

    expected_finite_trials = (
        n_tracks
        * settings.finite_replays_per_track
        * len(FINITE_CLASSES)
        * len(shapes)
        * len(settings.projected_areas)
    )
    expected_unbounded_trials = (
        settings.unbounded_trials_per_cell
        * len(UNBOUNDED_CLASSES)
        * len(shapes)
        * len(settings.projected_areas)
    )
    if (
        int(summary.get("tracks", -1)) != n_tracks
        or int(summary.get("pool_vectors", -1)) != n_pool_vectors
        or int(summary.get("finite_trial_count", -1))
        != expected_finite_trials
        or int(summary.get("unbounded_trial_count", -1))
        != expected_unbounded_trials
    ):
        raise ValueError("condition summary has inconsistent totals")


def _load_resumable_condition(
    store: ArtifactStore,
    condition_stage: str,
    frame_names: Sequence[str],
    expected_signature: str,
    expected_config_hash: str,
    settings: AllConditionsReplaySettings,
    n_tracks: int,
    n_pool_vectors: int,
    shapes: Sequence[str],
) -> Optional[Tuple[Dict[str, object], Dict[str, pd.DataFrame]]]:
    """Load a complete signed checkpoint, otherwise require a clean rerun."""
    if not store.json_path(condition_stage, "settings").is_file():
        return None
    try:
        stored_settings = store.read_json(condition_stage, "settings")
        if not _stored_signature_matches(
            stored_settings, expected_signature
        ):
            logger.warning(
                "Condition replay checkpoint %s is unsigned, obsolete, or "
                "does not match the requested settings; rerunning it",
                condition_stage,
            )
            return None
        summary = store.read_json(condition_stage, "summary")
        if not isinstance(summary, dict):
            raise ValueError("condition summary is not a JSON object")
        if (
            summary.get("replay_signature_version")
            != REPLAY_SIGNATURE_VERSION
            or summary.get("replay_signature") != expected_signature
        ):
            raise ValueError("condition summary signature does not match")
        frames = {
            name: store.read_frame(
                condition_stage,
                name,
                expected_config_hash=expected_config_hash,
            )
            for name in frame_names
        }
        _validate_replay_frames(
            frames,
            summary,
            settings,
            n_tracks,
            n_pool_vectors,
            shapes,
        )
        figure_directory = (
            store.root / condition_stage / "figures"
        )
        for filename in (
            "finite_detection_probability.pdf",
            "unbounded_detection_distance.pdf",
        ):
            if not (figure_directory / filename).is_file():
                raise ValueError(
                    "condition replay figure is missing: {}".format(filename)
                )
        return summary, frames
    except (ArtifactError, OSError, TypeError, ValueError, KeyError) as exc:
        logger.warning(
            "Condition replay checkpoint %s failed validation (%s); rerunning it",
            condition_stage,
            exc,
        )
        return None


def _attach_condition(frame: pd.DataFrame, condition: str) -> pd.DataFrame:
    out = frame.copy()
    fields = _condition_fields(condition)
    for key, value in fields.items():
        out[key] = value
    ordered = list(fields) + [
        column for column in out.columns if column not in fields
    ]
    return out[ordered]


def write_runtime_plan(
    config: Config, settings: AllConditionsReplaySettings
) -> Dict[str, object]:
    """Persist the pre-run workload and wall-time estimate."""
    settings.validate()
    store = ArtifactStore(config.output.root)
    pool_audit = store.read_frame(
        "preprocess",
        "pool_audit",
        schema_name="pool_audit",
        expected_config_hash=config.hash(),
    )
    manifest_rows: List[Dict[str, object]] = []
    for _, row in pool_audit.sort_values("condition").iterrows():
        condition = str(row["condition"])
        values = _condition_fields(condition)
        values.update(
            {
                "n_tracks": int(row["n_tracks"]),
                "n_runs": int(row["n_runs"]),
                "n_videos": int(row["n_videos"]),
                "eligible": bool(row["eligible"]),
                "ineligible_reason": (
                    ""
                    if pd.isna(row["ineligible_reason"])
                    else str(row["ineligible_reason"])
                ),
                "replay_selected": bool(row["eligible"]),
            }
        )
        manifest_rows.append(values)
    manifest = pd.DataFrame(manifest_rows)
    selected = manifest[manifest["replay_selected"]]
    n_conditions = int(len(selected))
    n_tracks = int(selected["n_tracks"].sum())
    n_areas = len(settings.projected_areas)
    n_shapes = len(config.replay.shapes)
    finite_evaluations = (
        n_tracks
        * settings.finite_replays_per_track
        * len(FINITE_CLASSES)
        * n_shapes
        * n_areas
    )
    unbounded_trials = (
        n_conditions
        * settings.unbounded_trials_per_cell
        * len(UNBOUNDED_CLASSES)
        * n_shapes
        * n_areas
    )
    condition_scale = n_conditions / 8.0
    finite_scale = finite_evaluations / 295_680_000.0
    unbounded_scale = unbounded_trials / 384_000.0
    replay_scale = 0.65 * finite_scale + 0.35 * unbounded_scale
    runtime = pd.DataFrame(
        [
            {
                "stage": "inspect_and_preprocess",
                "estimated_minutes": 1.0,
                "lower_minutes": 0.1,
                "upper_minutes": 2.0,
                "basis": "44-file preprocessing benchmark",
            },
            {
                "stage": "equal_run_fit_and_sensitivities",
                "estimated_minutes": 18.0 * condition_scale,
                "lower_minutes": 15.0 * condition_scale,
                "upper_minutes": 25.0 * condition_scale,
                "basis": "measured eight-condition fit, scaled by conditions",
            },
            {
                "stage": "finite_and_unbounded_replay",
                "estimated_minutes": 12.0 * replay_scale,
                "lower_minutes": 8.0 * replay_scale,
                "upper_minutes": 15.0 * replay_scale,
                "basis": "measured 32-core replay, scaled by finite and unbounded workloads",
            },
            {
                "stage": "combined_tables_pdfs_validation",
                "estimated_minutes": 3.0,
                "lower_minutes": 1.0,
                "upper_minutes": 5.0,
                "basis": "CSV and vector-PDF output size",
            },
        ]
    )
    stage = "00_plan"
    metadata = {
        "config_hash": config.hash(),
        "weighting": config.empirical_pool.primary_weighting,
    }
    store.write_frame(
        stage,
        "condition_manifest",
        manifest,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "runtime_estimate",
        runtime,
        extra_metadata=metadata,
    )
    plan = {
        "weighting": config.empirical_pool.primary_weighting,
        "eligible_conditions": n_conditions,
        "eligible_tracks": n_tracks,
        "eligible_runs": int(selected["n_runs"].sum()),
        "projected_areas": list(settings.projected_areas),
        "finite_class_target_evaluations": int(finite_evaluations),
        "unbounded_trials": int(unbounded_trials),
        "workers": settings.workers,
        "estimated_total_minutes": float(runtime["estimated_minutes"].sum()),
        "estimated_range_minutes": [
            float(runtime["lower_minutes"].sum()),
            float(runtime["upper_minutes"].sum()),
        ],
        # CSVs are intentionally duplicated at condition and combined levels.
        # The interval is calibrated to the measured 208 MB eight-condition
        # analysis and scaled by both replay workloads.
        "estimated_storage_mb": [
            int(round(30.0 + 75.0 * finite_scale + 115.0 * unbounded_scale)),
            int(round(40.0 + 95.0 * finite_scale + 145.0 * unbounded_scale)),
        ],
    }
    store.write_json(stage, "plan", plan)
    return plan


def _target_feasibility(
    config: Config,
    fits: pd.DataFrame,
    projected_areas: Sequence[float],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    d = float(config.replay.detection_radius)
    for _, fit in fits.sort_values("condition").iterrows():
        condition = str(fit["condition"])
        c = float(fit["crossover"])
        lmax = float(fit["lmax"])
        side = 2.0 * lmax / c
        fields = _condition_fields(condition)
        for area in projected_areas:
            for shape in config.replay.shapes:
                dimension = dimension_from_projected_area(area, shape, d)
                row = dict(fields)
                row.update(
                    {
                        "shape": shape,
                        "projected_area_model_units2": float(area),
                        "projected_area_um2": float(area * c * c),
                        "dimension_model_units": dimension,
                        "detection_radius_model_units": d,
                        "crossover_um": c,
                        "lmax_model_units": lmax / c,
                        "torus_side_model_units": side,
                        "fits_without_periodic_self_overlap": bool(
                            dimension + 2.0 * d <= side
                        ),
                        "clearance_to_torus_half_width": (
                            side / 2.0 - dimension / 2.0 - d
                        ),
                    }
                )
                rows.append(row)
    return pd.DataFrame(rows)


def _plot_all_conditions(
    output_directory: Path,
    fits: pd.DataFrame,
    pools: pd.DataFrame,
    track_inputs: pd.DataFrame,
    finite: pd.DataFrame,
    finite_contrasts: pd.DataFrame,
    unbounded: pd.DataFrame,
    unbounded_contrasts: pd.DataFrame,
    old_fits: Optional[pd.DataFrame] = None,
) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    figures = output_directory / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []

    fit_plot = fits.sort_values(["cell_type", "condition_label"]).reset_index(
        drop=True
    )
    cell_types = list(dict.fromkeys(fit_plot["cell_type"].astype(str)))
    cell_colors = dict(
        zip(cell_types, plt.get_cmap("tab10").colors[: len(cell_types)])
    )
    fig, axis = plt.subplots(figsize=(8.4, 5.4))
    y = np.arange(len(fit_plot))
    for index, row in fit_plot.iterrows():
        color = cell_colors[str(row["cell_type"])]
        value = float(row["mu_hat"])
        low = float(row["mu_ci_low"])
        high = float(row["mu_ci_high"])
        supported = float(row["ks_bootstrap_p"]) >= 0.05
        axis.errorbar(
            value,
            index,
            xerr=[[value - low], [high - value]],
            marker="o",
            markerfacecolor=color if supported else "white",
            markeredgecolor=color,
            color=color,
            capsize=2,
            linewidth=1.2,
        )
    axis.axvline(2.0, color="0.35", linestyle="--", linewidth=1.0)
    axis.set_yticks(y)
    axis.set_yticklabels(fit_plot["condition_label"], fontsize=8)
    axis.invert_yaxis()
    axis.set_xlabel(r"Equal-run fitted effective exponent $\mu$")
    axis.set_title("Effective exponents by eligible condition")
    axis.grid(axis="x", alpha=0.25)
    axis.legend(
        handles=[
            Patch(facecolor=color, label=cell)
            for cell, color in cell_colors.items()
        ]
        + [
            Line2D(
                [0],
                [0],
                marker="o",
                color="0.3",
                markerfacecolor="0.3",
                linestyle="",
                label=r"KS $p\geq0.05$",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="0.3",
                markerfacecolor="white",
                linestyle="",
                label=r"KS $p<0.05$",
            ),
        ],
        fontsize=7,
        frameon=False,
        loc="best",
    )
    fig.tight_layout()
    path = figures / "fitted_exponents_by_condition.pdf"
    fig.savefig(path)
    plt.close(fig)
    paths.append(path)

    conditions = list(fits.sort_values("condition_label")["condition"])
    condition_colors = dict(
        zip(conditions, plt.get_cmap("tab20").colors[: len(conditions)])
    )
    class_styles = {
        "exact_orientation": ("-", "o", "Exact orientation"),
        "exact_global_rotation": ("--", "s", "Exact + global rotation"),
        "ordered_length_rotated": (
            (0, (1, 1)),
            "X",
            "Exact lengths/order + independent rotations",
        ),
        "vector_uniform_rotated": (
            "-.",
            "^",
            "Uniform empirical run + independent rotation",
        ),
        "fitted_mu": ((0, (5, 1, 1, 1)), "D", r"Fitted $\mu$"),
    }
    row_cells = [
        value
        for value in ("B Cells", "Natural Killer Cells", "Neutrophils", "T Cells")
        if value in set(finite["cell_type"])
    ]
    shapes = list(dict.fromkeys(finite["shape"]))
    finite_areas = sorted(finite["projected_area"].astype(float).unique())
    fig, axes = plt.subplots(
        len(row_cells),
        len(shapes),
        figsize=(12.0, 2.7 * len(row_cells)),
        sharex=True,
        squeeze=False,
    )
    for row_index, cell_type in enumerate(row_cells):
        for column_index, shape in enumerate(shapes):
            axis = axes[row_index, column_index]
            selected = finite[
                (finite["cell_type"] == cell_type)
                & (finite["shape"] == shape)
            ]
            for condition in sorted(selected["condition"].unique()):
                for class_name in FINITE_PLOT_ORDER:
                    group = selected[
                        (selected["condition"] == condition)
                        & (selected["trajectory_class"] == class_name)
                    ].sort_values("projected_area")
                    if group.empty:
                        continue
                    linestyle, marker, _ = class_styles[class_name]
                    axis.plot(
                        group["projected_area"],
                        100.0 * group["detection_probability"],
                        color=condition_colors[condition],
                        linestyle=linestyle,
                        marker=marker,
                        markersize=3.2,
                        linewidth=1.0,
                    )
            axis.set_xscale("log", base=2)
            axis.set_xticks(finite_areas)
            axis.set_xticklabels(["{:g}".format(value) for value in finite_areas])
            axis.grid(alpha=0.2)
            if row_index == 0:
                axis.set_title(shape)
            if column_index == 0:
                axis.set_ylabel(
                    "{}\nDetection probability (%)".format(
                        {
                            "Natural Killer Cells": "NK",
                            "Neutrophils": "Neutrophils",
                            "T Cells": "T",
                            "B Cells": "B",
                        }[cell_type]
                    )
                )
            if row_index == len(row_cells) - 1:
                axis.set_xlabel(r"Projected area $A/c^2$")
    condition_handles = [
        Line2D(
            [0],
            [0],
            color=condition_colors[condition],
            linewidth=2,
            label=_condition_label(condition),
        )
        for condition in conditions
    ]
    class_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=style[0],
            marker=style[1],
            linewidth=1.2,
            label=style[2],
        )
        for style in class_styles.values()
    ]
    fig.legend(
        handles=condition_handles,
        loc="upper center",
        ncol=5,
        fontsize=6.5,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.legend(
        handles=class_handles,
        loc="lower center",
        ncol=4,
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )
    fig.suptitle(
        "Finite-budget detection: colors are conditions; line styles are trajectory classes",
        y=1.055,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.94))
    path = figures / "finite_detection_probability_by_condition.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    primary_contrasts = [
        (
            "exact_global_rotation",
            "exact_orientation",
            "Global rotation − exact orientation\n(laboratory-frame orientation)",
        ),
        (
            "ordered_length_rotated",
            "exact_global_rotation",
            "Independent vector rotations − global rotation\n"
            "(directional structure)",
        ),
        (
            "vector_uniform_rotated",
            "ordered_length_rotated",
            "Condition pool − exact ordered lengths\n"
            "(track/pool length structure; composite)",
        ),
        (
            "fitted_mu",
            "vector_uniform_rotated",
            r"Fitted $\mu$ − condition pool"
            "\n(parametric marginal law)",
        ),
    ]
    primary_pairs = {(left, right) for left, right, _ in primary_contrasts}
    selected_contrasts = finite_contrasts[
        [
            (str(left), str(right)) in primary_pairs
            for left, right in finite_contrasts[
                ["left_class", "right_class"]
            ].itertuples(index=False, name=None)
        ]
    ].copy()
    if not selected_contrasts.empty:
        contrast_pairs = _ordered_shape_area_pairs(selected_contrasts)
        contrast_columns = pd.MultiIndex.from_tuples(
            contrast_pairs,
            names=["shape", "projected_area"],
        )
        contrast_labels = list(
            fits.sort_values("condition_label")["condition_label"].astype(str)
        )
        maximum = float(
            np.nanmax(
                np.abs(
                    selected_contrasts[
                        "detection_probability_difference"
                    ].to_numpy(dtype=float)
                )
            )
            * 100.0
        )
        maximum = max(maximum, 1e-6)
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(14.0, 9.2),
            sharex=True,
            sharey=True,
            squeeze=False,
        )
        image = None
        for axis, (left, right, title) in zip(
            axes.ravel(), primary_contrasts
        ):
            cell = selected_contrasts[
                (selected_contrasts["left_class"] == left)
                & (selected_contrasts["right_class"] == right)
            ]
            values = (
                cell.pivot(
                    index="condition_label",
                    columns=["shape", "projected_area"],
                    values="detection_probability_difference",
                )
                .reindex(index=contrast_labels, columns=contrast_columns)
                .to_numpy(dtype=float)
                * 100.0
            )
            image = axis.imshow(
                values,
                aspect="auto",
                cmap="coolwarm",
                vmin=-maximum,
                vmax=maximum,
            )
            significant = (
                (
                    cell["difference_ci_low"].to_numpy(dtype=float) > 0.0
                )
                | (
                    cell["difference_ci_high"].to_numpy(dtype=float) < 0.0
                )
            )
            markers = cell.assign(significant=significant).pivot(
                index="condition_label",
                columns=["shape", "projected_area"],
                values="significant",
            ).reindex(index=contrast_labels, columns=contrast_columns)
            for row_index, column_index in np.argwhere(
                markers.fillna(False).to_numpy(dtype=bool)
            ):
                axis.text(
                    column_index,
                    row_index,
                    "•",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=8,
                )
            axis.set_title(title, fontsize=9)
            axis.set_yticks(np.arange(len(contrast_labels)))
            axis.set_yticklabels(contrast_labels, fontsize=7)
            axis.set_xticks(np.arange(len(contrast_pairs)))
            axis.set_xticklabels(
                [
                    "{}\nA={:g}".format(shape, area)
                    for shape, area in contrast_pairs
                ],
                rotation=45,
                ha="right",
                fontsize=7,
            )
        if image is not None:
            colorbar_axis = fig.add_axes([0.91, 0.16, 0.018, 0.68])
            colorbar = fig.colorbar(
                image,
                cax=colorbar_axis,
            )
            colorbar.set_label("Detection-probability difference (percentage points)")
        fig.suptitle(
            "Finite-budget randomization contrasts; • marks a 95% interval excluding zero",
            y=0.995,
        )
        fig.subplots_adjust(
            left=0.22,
            right=0.88,
            bottom=0.12,
            top=0.91,
            wspace=0.08,
            hspace=0.24,
        )
        path = figures / "finite_randomization_contrasts_by_condition.pdf"
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    unbounded_styles = {
        "vector_uniform_rotated": (
            "-",
            "o",
            "Uniform empirical run + independent rotation",
        ),
        "fitted_mu": ("--", "s", r"Fitted $\mu$"),
    }
    unbounded_areas = sorted(
        unbounded["projected_area"].astype(float).unique()
    )
    fig, axes = plt.subplots(
        len(row_cells),
        len(shapes),
        figsize=(12.0, 2.7 * len(row_cells)),
        sharex=True,
        squeeze=False,
    )
    for row_index, cell_type in enumerate(row_cells):
        for column_index, shape in enumerate(shapes):
            axis = axes[row_index, column_index]
            selected = unbounded[
                (unbounded["cell_type"] == cell_type)
                & (unbounded["shape"] == shape)
            ]
            for condition in sorted(selected["condition"].unique()):
                for class_name in UNBOUNDED_CLASSES:
                    group = selected[
                        (selected["condition"] == condition)
                        & (selected["trajectory_class"] == class_name)
                    ].sort_values("projected_area")
                    if group.empty:
                        continue
                    linestyle, marker, _ = unbounded_styles[class_name]
                    axis.plot(
                        group["projected_area"],
                        group["mean_detection_distance_model_units"],
                        color=condition_colors[condition],
                        linestyle=linestyle,
                        marker=marker,
                        markersize=3.2,
                        linewidth=1.0,
                    )
            axis.set_xscale("log", base=2)
            axis.set_xticks(unbounded_areas)
            axis.set_xticklabels(
                ["{:g}".format(value) for value in unbounded_areas]
            )
            axis.set_yscale("log")
            axis.grid(alpha=0.2)
            if row_index == 0:
                axis.set_title(shape)
            if column_index == 0:
                axis.set_ylabel(
                    "{}\nMean detection distance ($c$ units)".format(
                        {
                            "Natural Killer Cells": "NK",
                            "Neutrophils": "Neutrophils",
                            "T Cells": "T",
                            "B Cells": "B",
                        }[cell_type]
                    )
                )
            if row_index == len(row_cells) - 1:
                axis.set_xlabel(r"Projected area $A/c^2$")
    class_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=style[0],
            marker=style[1],
            linewidth=1.2,
            label=style[2],
        )
        for style in unbounded_styles.values()
    ]
    fitted_mu_by_condition = (
        fits.set_index("condition")["mu_hat"].astype(float).to_dict()
    )
    unbounded_condition_handles = [
        Line2D(
            [0],
            [0],
            color=condition_colors[condition],
            linewidth=2,
            label=r"{} · $\mu={:.3f}$".format(
                _condition_label(condition),
                fitted_mu_by_condition[condition],
            ),
        )
        for condition in conditions
    ]
    fig.legend(
        handles=unbounded_condition_handles,
        loc="upper center",
        ncol=4,
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.legend(
        handles=class_handles,
        loc="lower center",
        ncol=2,
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )
    fig.suptitle(
        "Unbounded first-passage distance: renewable trajectory classes only",
        y=1.055,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.94))
    path = figures / "unbounded_detection_distance_by_condition.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    n_fit_panels = len(fits)
    n_fit_columns = min(4, max(1, n_fit_panels))
    n_fit_rows = int(math.ceil(n_fit_panels / n_fit_columns))
    fig, axes = plt.subplots(
        n_fit_rows,
        n_fit_columns,
        figsize=(3.0 * n_fit_columns, 3.1 * n_fit_rows),
        squeeze=False,
    )
    for axis, (_, fit) in zip(
        axes.ravel(),
        fits.sort_values("condition_label").iterrows(),
    ):
        condition = str(fit["condition"])
        lengths = pools.loc[
            pools["condition"] == condition, "run_length_um"
        ].to_numpy(dtype=float)
        lengths = np.sort(lengths[lengths > 0.0])
        # At the i-th ordered observation, this is the empirical
        # non-strict survival P(L >= lengths[i]); in particular the largest
        # observation has mass 1 / n rather than zero.
        empirical = (len(lengths) - np.arange(len(lengths))) / len(lengths)
        grid = np.geomspace(max(lengths.min(), 1e-6), float(fit["lmax"]), 400)
        fitted_survival = 1.0 - mixture_cdf(
            grid,
            float(fit["mu_hat"]),
            float(fit["lmax"]),
            float(fit["crossover"]),
        )
        axis.step(lengths, empirical, where="post")
        axis.plot(grid, np.maximum(fitted_survival, 1e-8), color="#d62728")
        axis.axvline(float(fit["crossover"]), color="0.4", linestyle=":", linewidth=0.8)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(_condition_label(condition), fontsize=8)
        axis.grid(alpha=0.2)
    for axis in axes.ravel()[n_fit_panels:]:
        axis.set_visible(False)
    bottom_row_start = (n_fit_rows - 1) * n_fit_columns
    for axis in axes.ravel()[bottom_row_start:n_fit_panels]:
        axis.set_xlabel(r"Complete run length ($\mu$m)")
    for axis in axes[:, 0]:
        if axis.get_visible():
            axis.set_ylabel("CCDF")
    fig.legend(
        handles=[
            Line2D([0], [0], color="#1f77b4", label="Equal-run empirical"),
            Line2D([0], [0], color="#d62728", label="Fitted mixture"),
            Line2D([0], [0], color="0.4", linestyle=":", label="Fitted crossover"),
        ],
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.955),
    )
    fig.suptitle("Equal-run empirical distributions and fitted models", y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    path = figures / "fit_ccdf_diagnostics.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    contrast = unbounded_contrasts.copy()
    pivot = contrast.pivot(
        index="condition_label",
        columns=["shape", "projected_area"],
        values="fitted_to_empirical_distance_ratio",
    )
    ordered_pairs = _ordered_shape_area_pairs(contrast)
    ordered_columns = pd.MultiIndex.from_tuples(
        ordered_pairs,
        names=["shape", "projected_area"],
    )
    pivot = pivot.reindex(columns=ordered_columns)
    fig, axis = plt.subplots(figsize=(11.0, 5.2))
    image = axis.imshow(
        np.log2(pivot.to_numpy(dtype=float)),
        aspect="auto",
        cmap="coolwarm",
        vmin=-0.5,
        vmax=0.5,
    )
    axis.set_yticks(np.arange(len(pivot)))
    axis.set_yticklabels(pivot.index, fontsize=7)
    axis.set_xticks(np.arange(len(pivot.columns)))
    axis.set_xticklabels(
        [
            "{}\nA={:g}".format(shape, area)
            for shape, area in ordered_pairs
        ],
        rotation=45,
        ha="right",
        fontsize=7,
    )
    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label(r"$\log_2(D_{\rm fitted}/D_{\rm empirical})$")
    axis.set_title("Transfer of the fitted model to unbounded detection distance")
    fig.tight_layout()
    path = figures / "unbounded_fitted_to_empirical_ratio.pdf"
    fig.savefig(path)
    plt.close(fig)
    paths.append(path)

    fig, axis = plt.subplots(figsize=(9.0, 4.8))
    budget_groups = [
        group["budget_model_units"].to_numpy(dtype=float)
        for _, group in track_inputs.groupby("condition_label", sort=True)
    ]
    budget_labels = [
        name for name, _ in track_inputs.groupby("condition_label", sort=True)
    ]
    axis.boxplot(budget_groups, labels=budget_labels, vert=False, showfliers=False)
    axis.set_xlabel("Observed finite-path budget (fitted-crossover units)")
    axis.set_title("Recorded path-budget distributions by condition")
    axis.grid(axis="x", alpha=0.2)
    axis.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    path = figures / "finite_budget_distributions.pdf"
    fig.savefig(path)
    plt.close(fig)
    paths.append(path)

    if old_fits is not None and not old_fits.empty:
        old = old_fits[["condition", "mu_hat"]].rename(
            columns={"mu_hat": "cell_balanced_mu"}
        )
        comparison = fits.merge(old, on="condition", how="inner").sort_values(
            "condition_label"
        )
        if not comparison.empty:
            fig, axis = plt.subplots(figsize=(8.5, 5.0))
            for index, row in comparison.reset_index(drop=True).iterrows():
                axis.plot(
                    [row["cell_balanced_mu"], row["mu_hat"]],
                    [index, index],
                    color="0.65",
                    linewidth=1.2,
                )
                axis.scatter(row["cell_balanced_mu"], index, color="#1f77b4")
                axis.scatter(row["mu_hat"], index, color="#d62728")
            axis.set_yticks(np.arange(len(comparison)))
            axis.set_yticklabels(comparison["condition_label"], fontsize=7)
            axis.invert_yaxis()
            axis.set_xlabel(r"Effective fitted exponent $\mu$")
            axis.set_title("Effect of likelihood weighting on the fitted exponent")
            axis.grid(axis="x", alpha=0.2)
            axis.legend(
                handles=[
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="#1f77b4",
                        linestyle="",
                        label="Cell-balanced",
                    ),
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="#d62728",
                        linestyle="",
                        label="Equal-run",
                    ),
                ],
                frameon=False,
            )
            fig.tight_layout()
            path = figures / "weighting_effect_on_mu.pdf"
            fig.savefig(path)
            plt.close(fig)
            paths.append(path)

    return paths


def run_all_conditions_replay(
    config: Config, settings: AllConditionsReplaySettings
) -> Dict[str, object]:
    """Run every successful eligible condition and write combined artifacts."""
    settings.validate()
    if config.empirical_pool.primary_weighting != "step_weighted":
        raise ValueError(
            "The all-condition equal-run experiment requires "
            "empirical_pool.primary_weighting=step_weighted"
        )
    store = ArtifactStore(config.output.root)
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config.hash(),
    )
    fits = fits[fits["status"] == "ok"].sort_values("condition").reset_index(
        drop=True
    )
    if fits.empty:
        raise ValueError("No successful eligible fits are available")
    pools = store.read_frame(
        "preprocess",
        "pools",
        expected_config_hash=config.hash(),
    )
    # The exact-orientation replay consumes the segmented run artifact even
    # though the orchestration layer itself only needs the empirical pool.
    # Verify it here so a signed checkpoint cannot mask a corrupted input.
    store.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config.hash(),
    )
    source_artifact_hashes = {
        "runs": _artifact_sha256(store, "preprocess", "runs"),
        "pools": _artifact_sha256(store, "preprocess", "pools"),
        "mixture_fits": _artifact_sha256(store, "fit", "mixture_fits"),
    }

    feasibility = _target_feasibility(
        config, fits, settings.projected_areas
    )
    stage = "all_conditions_replay"
    metadata = {
        "config_hash": config.hash(),
        "weighting": config.empirical_pool.primary_weighting,
    }
    store.write_frame(
        stage,
        "target_feasibility",
        feasibility,
        extra_metadata=metadata,
    )
    invalid = feasibility[
        ~feasibility["fits_without_periodic_self_overlap"]
    ]
    if not invalid.empty:
        raise ValueError(
            "The common target grid contains {} invalid condition/shape cells; "
            "see target_feasibility.csv".format(len(invalid))
        )

    combined: Dict[str, List[pd.DataFrame]] = {
        "track_inputs": [],
        "empirical_pool_inputs": [],
        "finite_track_summaries": [],
        "finite_summary": [],
        "finite_contrasts": [],
        "unbounded_trials": [],
        "unbounded_summary": [],
        "unbounded_contrasts": [],
        "targets": [],
    }
    condition_summaries: List[Dict[str, object]] = []
    for index, fit in fits.iterrows():
        condition = str(fit["condition"])
        condition_id = _condition_id(condition)
        condition_stage = "condition_replays/{}".format(condition_id)
        condition_pool = pools[
            (pools["condition"] == condition)
            & np.isfinite(
                pools["run_length_um"].to_numpy(dtype=float)
            )
            & (pools["run_length_um"].to_numpy(dtype=float) > 0.0)
        ]
        n_pool_vectors = int(len(condition_pool))
        n_source_tracks = int(condition_pool["track_uid"].nunique())
        signature_payload = _replay_signature_payload(
            config,
            settings,
            condition,
            fit,
            source_artifact_hashes,
        )
        replay_signature = _canonical_signature(signature_payload)
        checkpoint = (
            _load_resumable_condition(
                store,
                condition_stage,
                tuple(combined),
                replay_signature,
                config.hash(),
                settings,
                n_source_tracks,
                n_pool_vectors,
                config.replay.shapes,
            )
            if settings.resume
            else None
        )
        if checkpoint is not None:
            logger.info(
                "Reusing completed condition replay %d/%d: %s",
                index + 1,
                len(fits),
                condition,
            )
            condition_summary, replay_frames = checkpoint
        else:
            logger.info(
                "Running condition replay %d/%d: %s",
                index + 1,
                len(fits),
                condition,
            )
            replay_settings = FocusedReplaySettings(
                condition=condition,
                length_unit_um=None,
                projected_areas=tuple(settings.projected_areas),
                finite_replays_per_track=settings.finite_replays_per_track,
                unbounded_trials_per_cell=settings.unbounded_trials_per_cell,
                workers=settings.workers,
                finite_batch_size=settings.finite_batch_size,
                unbounded_chunk_trials=settings.unbounded_chunk_trials,
                unbounded_step_block=settings.unbounded_step_block,
                random_seed=_condition_seed(settings.random_seed, condition),
                output_stage=condition_stage,
                figure_format="pdf",
                replay_signature=replay_signature,
                replay_signature_version=REPLAY_SIGNATURE_VERSION,
            )
            condition_summary = run_focused_replay(config, replay_settings)
            replay_frames = {
                name: store.read_frame(
                    condition_stage,
                    name,
                    expected_config_hash=config.hash(),
                )
                for name in combined
            }
            _validate_replay_frames(
                replay_frames,
                condition_summary,
                settings,
                n_source_tracks,
                n_pool_vectors,
                config.replay.shapes,
            )
        summary_row = dict(_condition_fields(condition))
        summary_row.update(condition_summary)
        condition_summaries.append(summary_row)
        for name in combined:
            combined[name].append(
                _attach_condition(replay_frames[name], condition)
            )

    for name, frames in combined.items():
        store.write_frame(
            stage,
            name,
            pd.concat(frames, ignore_index=True),
            extra_metadata=metadata,
        )
    condition_summary_frame = pd.DataFrame(condition_summaries)
    store.write_frame(
        stage,
        "condition_summaries",
        condition_summary_frame,
        extra_metadata=metadata,
    )

    fit_plot = fits.copy()
    for key in (
        "condition_id",
        "condition_label",
        "cell_type",
        "organ",
        "stimulus",
        "dt_seconds",
    ):
        fit_plot[key] = [
            _condition_fields(str(condition))[key]
            for condition in fit_plot["condition"]
        ]
    store.write_frame(
        stage,
        "fit_summary",
        fit_plot,
        extra_metadata=metadata,
    )

    old_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "primary"
        / "fit"
        / "mixture_fits.csv"
    )
    old_fits = pd.read_csv(old_path) if old_path.is_file() else None
    output_directory = store.stage_dir(stage)
    figure_paths = _plot_all_conditions(
        output_directory,
        fit_plot,
        pools,
        pd.concat(combined["track_inputs"], ignore_index=True),
        pd.concat(combined["finite_summary"], ignore_index=True),
        pd.concat(combined["finite_contrasts"], ignore_index=True),
        pd.concat(combined["unbounded_summary"], ignore_index=True),
        pd.concat(combined["unbounded_contrasts"], ignore_index=True),
        old_fits=old_fits,
    )
    summary = {
        "weighting": config.empirical_pool.primary_weighting,
        "conditions": int(len(fits)),
        "tracks": int(condition_summary_frame["tracks"].sum()),
        "pool_vectors": int(condition_summary_frame["pool_vectors"].sum()),
        "finite_trial_count": int(
            condition_summary_frame["finite_trial_count"].sum()
        ),
        "unbounded_trial_count": int(
            condition_summary_frame["unbounded_trial_count"].sum()
        ),
        "projected_areas": list(settings.projected_areas),
        "figures": [str(path) for path in figure_paths],
    }
    store.write_json(stage, "summary", summary)
    store.write_json(
        stage,
        "settings",
        {
            "projected_areas": list(settings.projected_areas),
            "finite_replays_per_track": settings.finite_replays_per_track,
            "unbounded_trials_per_cell": settings.unbounded_trials_per_cell,
            "workers": settings.workers,
            "finite_batch_size": settings.finite_batch_size,
            "unbounded_chunk_trials": settings.unbounded_chunk_trials,
            "unbounded_step_block": settings.unbounded_step_block,
            "random_seed": settings.random_seed,
            "resume": settings.resume,
            "replay_algorithm_version": REPLAY_ALGORITHM_VERSION,
            "finite_classes": list(FINITE_CLASSES),
            "finite_plot_order": list(FINITE_PLOT_ORDER),
            "finite_contrasts": [list(pair) for pair in FINITE_CONTRASTS],
            "normalization": "each condition divided by its fitted crossover",
            "detection_radius_model_units": config.replay.detection_radius,
            "torus_width": "2 * fitted Lmax in normalized units",
        },
    )
    return summary


__all__ = [
    "AllConditionsReplaySettings",
    "run_all_conditions_replay",
    "write_runtime_plan",
]
