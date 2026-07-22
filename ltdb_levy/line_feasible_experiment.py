"""Short, condition-specific replay on a line-feasible target-area grid.

This additive experiment extends the original common-area replay without
modifying it.  Each condition uses the same projected areas for Ball, Disk,
and Line, but retains only candidate areas for which the Line detection
neighbourhood fits in that condition's ``2*Lmax/c`` torus.  Since Line is the
most restrictive configured shape at large projected area, all three targets
then fit and remain directly comparable at equal ``A/c^2``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactError, ArtifactStore
from .dataio.artifacts import sha256_file
from .logging_utils import get_logger
from .multicondition_experiment import (
    _condition_fields,
    _condition_id,
    _condition_seed,
)
from .nk_experiment import (
    FINITE_CLASSES,
    FINITE_CONTRASTS,
    UNBOUNDED_CLASSES,
    FocusedReplaySettings,
    run_focused_replay,
)
from .publication_figures import (
    SELECTED_CONDITIONS,
    _unbounded_ratio_table,
)

logger = get_logger(__name__)

ALGORITHM_VERSION = "line_feasible_short_replay_v1"
DEFAULT_OUTPUT_STAGE = "line_feasible_short_replay"
DEFAULT_CANDIDATE_AREAS: Tuple[float, ...] = (
    4.0,
    6.0,
    8.0,
    12.0,
    16.0,
    24.0,
    32.0,
    48.0,
    64.0,
    96.0,
    128.0,
)

COMBINED_FRAME_NAMES: Tuple[str, ...] = (
    "track_inputs",
    "empirical_pool_inputs",
    "finite_track_summaries",
    "finite_summary",
    "finite_contrasts",
    "unbounded_trials",
    "unbounded_summary",
    "unbounded_contrasts",
    "targets",
)


@dataclass(frozen=True)
class LineFeasibleSweepSettings:
    """Runtime and output controls for the short line-feasible sweep."""

    candidate_areas: Tuple[float, ...] = DEFAULT_CANDIDATE_AREAS
    finite_replays_per_track: int = 2_000
    unbounded_trials_per_cell: int = 500
    workers: int = 64
    finite_batch_size: int = 500
    unbounded_chunk_trials: int = 50
    unbounded_step_block: int = 256
    random_seed: int = 731_291
    output_stage: str = DEFAULT_OUTPUT_STAGE
    resume: bool = True

    def validate(self, detection_radius: float = 1.0) -> None:
        if not self.candidate_areas:
            raise ValueError("At least one candidate area is required")
        if tuple(sorted(set(self.candidate_areas))) != self.candidate_areas:
            raise ValueError("candidate_areas must be unique and increasing")
        minimum = math.pi * detection_radius**2
        if any(not np.isfinite(area) or area <= minimum for area in self.candidate_areas):
            raise ValueError("Every candidate area must exceed pi*d^2")
        if self.finite_replays_per_track < 1:
            raise ValueError("finite_replays_per_track must be positive")
        if self.unbounded_trials_per_cell < 2:
            raise ValueError("unbounded_trials_per_cell must be at least two")
        if self.workers < 1:
            raise ValueError("workers must be positive")
        if self.finite_batch_size < 1 or self.unbounded_chunk_trials < 1:
            raise ValueError("batch and chunk sizes must be positive")
        if self.unbounded_step_block < 1:
            raise ValueError("unbounded_step_block must be positive")
        if not self.output_stage.strip():
            raise ValueError("output_stage must be non-empty")


def line_projected_area_limit(
    torus_side_model_units: float,
    detection_radius_model_units: float = 1.0,
) -> float:
    """Largest projected area whose line capsule fits without self-overlap."""
    side = float(torus_side_model_units)
    radius = float(detection_radius_model_units)
    if side <= 0.0 or radius <= 0.0:
        raise ValueError("torus side and detection radius must be positive")
    # Line dimension is (A - pi*d^2)/(2*d), while the full capsule extent is
    # dimension + 2*d.  Solving that extent <= side gives this bound.
    return float(2.0 * radius * side + (math.pi - 4.0) * radius**2)


def line_feasible_areas(
    candidate_areas: Sequence[float],
    torus_side_model_units: float,
    detection_radius_model_units: float = 1.0,
) -> Tuple[float, ...]:
    """Return candidate areas admitted by the line-neighbourhood constraint."""
    limit = line_projected_area_limit(
        torus_side_model_units,
        detection_radius_model_units,
    )
    tolerance = 1e-12 * max(1.0, abs(limit))
    result = tuple(float(area) for area in candidate_areas if area <= limit + tolerance)
    if not result:
        raise ValueError(
            "No candidate target area fits the side-{:.6g} torus".format(
                torus_side_model_units
            )
        )
    return result


def _experiment_signature(
    config: Config,
    settings: LineFeasibleSweepSettings,
    condition: str,
    fit: Mapping[str, object],
    areas: Sequence[float],
) -> str:
    payload = {
        "algorithm_version": ALGORITHM_VERSION,
        "config_hash": config.hash(),
        "condition": condition,
        "mu": float(fit["mu_hat"]),
        "crossover_um": float(fit["crossover"]),
        "lmax_um": float(fit["lmax"]),
        "areas": [float(value) for value in areas],
        "finite_replays_per_track": settings.finite_replays_per_track,
        "unbounded_trials_per_cell": settings.unbounded_trials_per_cell,
        "finite_batch_size": settings.finite_batch_size,
        "unbounded_chunk_trials": settings.unbounded_chunk_trials,
        "unbounded_step_block": settings.unbounded_step_block,
        "random_seed": _condition_seed(settings.random_seed, condition),
        "shapes": list(config.replay.shapes),
        "detection_radius": float(config.replay.detection_radius),
        "normalization": "fitted_crossover",
        "torus_rule": "2*fitted_lmax/c",
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _fit_table(config: Config, store: ArtifactStore) -> pd.DataFrame:
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
        raise ValueError("No successful fitted conditions are available")
    return fits


def _manifest(
    config: Config,
    fits: pd.DataFrame,
    settings: LineFeasibleSweepSettings,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    detection_radius = float(config.replay.detection_radius)
    n_shapes = len(config.replay.shapes)
    for _, fit in fits.iterrows():
        condition = str(fit["condition"])
        crossover = float(fit["crossover"])
        lmax = float(fit["lmax"])
        side = 2.0 * lmax / crossover
        areas = line_feasible_areas(
            settings.candidate_areas,
            side,
            detection_radius,
        )
        row = dict(_condition_fields(condition))
        row.update(
            {
                "mu_hat": float(fit["mu_hat"]),
                "crossover_um": crossover,
                "lmax_um": lmax,
                "torus_side_model_units": side,
                "torus_side_um": side * crossover,
                "line_projected_area_limit_model_units2": (
                    line_projected_area_limit(side, detection_radius)
                ),
                "projected_areas_model_units2": ",".join(
                    "{:g}".format(value) for value in areas
                ),
                "n_areas": len(areas),
                "n_target_cells": n_shapes * len(areas),
                "n_tracks": int(fit["n_tracks"]),
                "n_runs": int(fit["n_runs"]),
                "finite_evaluations": (
                    int(fit["n_tracks"])
                    * n_shapes
                    * len(areas)
                    * len(FINITE_CLASSES)
                    * settings.finite_replays_per_track
                ),
                "unbounded_trials": (
                    n_shapes
                    * len(areas)
                    * len(UNBOUNDED_CLASSES)
                    * settings.unbounded_trials_per_cell
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def write_runtime_plan(
    config: Config,
    settings: LineFeasibleSweepSettings,
) -> Dict[str, object]:
    """Persist the exact approved grid and a conservative short-run estimate."""
    settings.validate(config.replay.detection_radius)
    store = ArtifactStore(config.output.root)
    fits = _fit_table(config, store)
    manifest = _manifest(config, fits, settings)
    stage = "{}/00_plan".format(settings.output_stage)
    metadata = {
        "config_hash": config.hash(),
        "algorithm_version": ALGORITHM_VERSION,
    }
    store.write_frame(
        stage,
        "condition_target_manifest",
        manifest,
        extra_metadata=metadata,
    )
    runtime = pd.DataFrame(
        [
            {
                "stage": "condition_replays",
                "central_minutes": 4.0,
                "lower_minutes": 2.5,
                "upper_minutes": 7.0,
                "basis": (
                    "measured 32-worker replay scaled to 64 workers, 500 "
                    "unbounded trials/cell, and the line-feasible grid"
                ),
            },
            {
                "stage": "validation_tables_figures",
                "central_minutes": 1.0,
                "lower_minutes": 0.5,
                "upper_minutes": 2.0,
                "basis": "combined hashes, ragged-grid validation, PDF and PNG output",
            },
        ]
    )
    store.write_frame(
        stage,
        "runtime_estimate",
        runtime,
        extra_metadata=metadata,
    )
    plan = {
        "algorithm_version": ALGORITHM_VERSION,
        "conditions": int(len(manifest)),
        "target_cells": int(manifest["n_target_cells"].sum()),
        "finite_evaluations": int(manifest["finite_evaluations"].sum()),
        "unbounded_trials": int(manifest["unbounded_trials"].sum()),
        "candidate_areas_model_units2": list(settings.candidate_areas),
        "grid_rule": (
            "retain each candidate area when the Line detection neighbourhood "
            "fits; use that same area for Ball, Disk, and Line"
        ),
        "normalization": "length/c and area/c^2",
        "finite_replays_per_track": settings.finite_replays_per_track,
        "unbounded_trials_per_class_target_cell": (
            settings.unbounded_trials_per_cell
        ),
        "workers": settings.workers,
        "estimated_minutes": float(runtime["central_minutes"].sum()),
        "estimated_range_minutes": [
            float(runtime["lower_minutes"].sum()),
            float(runtime["upper_minutes"].sum()),
        ],
        "estimated_storage_mib": [150, 250],
        "output_stage": settings.output_stage,
    }
    store.write_json(stage, "plan", plan)
    return plan


def _expected_shape_area_pairs(
    shapes: Sequence[str], areas: Sequence[float]
) -> set:
    return {(str(shape), float(area)) for shape in shapes for area in areas}


def _validate_condition_frames(
    frames: Mapping[str, pd.DataFrame],
    summary: Mapping[str, object],
    settings_payload: Mapping[str, object],
    shapes: Sequence[str],
    areas: Sequence[float],
    finite_replays: int,
    unbounded_trials: int,
    expected_signature: str,
) -> None:
    """Fail closed on incomplete or stale focused-replay checkpoints."""
    missing = sorted(set(COMBINED_FRAME_NAMES).difference(frames))
    if missing:
        raise ValueError("Missing condition frames: {}".format(", ".join(missing)))
    if settings_payload.get("replay_signature") != expected_signature:
        raise ValueError("Condition checkpoint signature does not match")
    if [float(value) for value in settings_payload.get("projected_areas", [])] != [
        float(value) for value in areas
    ]:
        raise ValueError("Condition checkpoint uses a different area grid")
    if int(settings_payload.get("finite_replays_per_track", -1)) != finite_replays:
        raise ValueError("Condition checkpoint uses a different finite replication")
    if int(settings_payload.get("unbounded_trials_per_cell", -1)) != unbounded_trials:
        raise ValueError("Condition checkpoint uses a different unbounded replication")

    target_pairs = _expected_shape_area_pairs(shapes, areas)
    targets = frames["targets"]
    actual_target_pairs = set(
        targets[["shape", "projected_area"]].itertuples(index=False, name=None)
    )
    if actual_target_pairs != target_pairs or len(targets) != len(target_pairs):
        raise ValueError("Condition target grid is incomplete or duplicated")
    if not targets["fits_without_periodic_self_overlap"].astype(bool).all():
        raise ValueError("Condition target grid contains periodic self-overlap")
    by_shape = targets.groupby("shape")["projected_area"].apply(
        lambda values: tuple(sorted(float(value) for value in values))
    )
    if any(values != tuple(float(value) for value in areas) for values in by_shape):
        raise ValueError("Shapes do not share the line-feasible area grid")

    n_tracks = int(frames["track_inputs"]["track_uid"].nunique())
    finite_tracks = frames["finite_track_summaries"]
    expected_finite_rows = n_tracks * len(FINITE_CLASSES) * len(target_pairs)
    finite_key = ["track_uid", "trajectory_class", "shape", "projected_area"]
    if len(finite_tracks) != expected_finite_rows or finite_tracks.duplicated(finite_key).any():
        raise ValueError("Condition finite track grid has an incorrect row count")
    if set(finite_tracks["trajectory_class"].astype(str)) != set(FINITE_CLASSES):
        raise ValueError("Condition finite track grid has the wrong classes")
    if not np.all(finite_tracks["n_trials"].to_numpy(dtype=int) == finite_replays):
        raise ValueError("Condition finite track grid has the wrong trial count")
    probabilities = finite_tracks["detection_probability"].to_numpy(dtype=float)
    detected = finite_tracks["n_detected"].to_numpy(dtype=int)
    trials = finite_tracks["n_trials"].to_numpy(dtype=int)
    if np.any(detected < 0) or np.any(detected > trials) or not np.allclose(
        probabilities, detected / trials, rtol=0.0, atol=1e-12
    ):
        raise ValueError("Condition finite detections are internally inconsistent")

    finite_summary = frames["finite_summary"]
    if len(finite_summary) != len(FINITE_CLASSES) * len(target_pairs):
        raise ValueError("Condition finite summary has an incorrect row count")
    if len(frames["finite_contrasts"]) != len(FINITE_CONTRASTS) * len(target_pairs):
        raise ValueError("Condition finite contrasts have an incorrect row count")

    unbounded = frames["unbounded_trials"]
    expected_unbounded_rows = len(UNBOUNDED_CLASSES) * len(target_pairs) * unbounded_trials
    unbounded_key = ["trajectory_class", "shape", "projected_area", "trial_index"]
    if len(unbounded) != expected_unbounded_rows or unbounded.duplicated(unbounded_key).any():
        raise ValueError("Condition unbounded grid has an incorrect row count")
    counts = unbounded.groupby(["trajectory_class", "shape", "projected_area"]).size()
    if not np.all(counts.to_numpy(dtype=int) == unbounded_trials):
        raise ValueError("Condition unbounded cells have the wrong trial count")
    if set(unbounded["trajectory_class"].astype(str)) != set(UNBOUNDED_CLASSES):
        raise ValueError("Condition unbounded grid has the wrong classes")
    if len(frames["unbounded_summary"]) != len(UNBOUNDED_CLASSES) * len(target_pairs):
        raise ValueError("Condition unbounded summary has an incorrect row count")
    if len(frames["unbounded_contrasts"]) != len(target_pairs):
        raise ValueError("Condition unbounded contrasts have an incorrect row count")

    expected_finite_trials = expected_finite_rows * finite_replays
    if int(summary.get("finite_trial_count", -1)) != expected_finite_trials:
        raise ValueError("Condition summary has the wrong finite trial count")
    if int(summary.get("unbounded_trial_count", -1)) != expected_unbounded_rows:
        raise ValueError("Condition summary has the wrong unbounded trial count")


def _load_checkpoint(
    store: ArtifactStore,
    stage: str,
    config: Config,
    settings: LineFeasibleSweepSettings,
    areas: Sequence[float],
    signature: str,
) -> Tuple[Dict[str, object], Dict[str, pd.DataFrame]] | None:
    try:
        payload = store.read_json(stage, "settings")
        summary = store.read_json(stage, "summary")
        frames = {
            name: store.read_frame(
                stage,
                name,
                expected_config_hash=config.hash(),
            )
            for name in COMBINED_FRAME_NAMES
        }
        _validate_condition_frames(
            frames,
            summary,
            payload,
            config.replay.shapes,
            areas,
            settings.finite_replays_per_track,
            settings.unbounded_trials_per_cell,
            signature,
        )
    except (ArtifactError, KeyError, TypeError, ValueError):
        return None
    return dict(summary), frames


def _attach_condition_fields(frame: pd.DataFrame, condition: str) -> pd.DataFrame:
    result = frame.copy()
    fields = _condition_fields(condition)
    for key, value in fields.items():
        result[key] = value
    preferred = list(fields)
    return result[preferred + [column for column in result if column not in preferred]]


def _enhance_targets(
    targets: pd.DataFrame,
    crossover_um: float,
) -> pd.DataFrame:
    result = targets.copy()
    result["projected_area_model_units2"] = result["projected_area"].astype(float)
    result["projected_area_um2"] = (
        result["projected_area_model_units2"] * crossover_um**2
    )
    result["dimension_model_units"] = result["dimension"].astype(float)
    result["dimension_um"] = result["dimension_model_units"] * crossover_um
    result["torus_side_model_units"] = result["side_length"].astype(float)
    result["torus_side_um"] = result["torus_side_model_units"] * crossover_um
    result["crossover_um"] = crossover_um
    return result


def _combined_grid_validation(
    manifest: pd.DataFrame,
    frames: Mapping[str, pd.DataFrame],
    settings: LineFeasibleSweepSettings,
    shapes: Sequence[str],
) -> None:
    expected_target_cells = int(manifest["n_target_cells"].sum())
    if len(frames["targets"]) != expected_target_cells:
        raise ValueError("Combined target table has the wrong row count")
    if not frames["targets"]["fits_without_periodic_self_overlap"].astype(bool).all():
        raise ValueError("Combined target table contains an invalid target")
    expected_unbounded_rows = expected_target_cells * len(UNBOUNDED_CLASSES) * settings.unbounded_trials_per_cell
    if len(frames["unbounded_trials"]) != expected_unbounded_rows:
        raise ValueError("Combined unbounded table has the wrong row count")
    if len(frames["unbounded_summary"]) != expected_target_cells * len(UNBOUNDED_CLASSES):
        raise ValueError("Combined unbounded summary has the wrong row count")
    if len(frames["unbounded_contrasts"]) != expected_target_cells:
        raise ValueError("Combined unbounded contrast table has the wrong row count")
    if len(frames["finite_summary"]) != expected_target_cells * len(FINITE_CLASSES):
        raise ValueError("Combined finite summary has the wrong row count")
    if len(frames["finite_contrasts"]) != expected_target_cells * len(FINITE_CONTRASTS):
        raise ValueError("Combined finite contrast table has the wrong row count")
    expected_finite_rows = int(
        (manifest["n_tracks"] * manifest["n_target_cells"]).sum()
    ) * len(FINITE_CLASSES)
    if len(frames["finite_track_summaries"]) != expected_finite_rows:
        raise ValueError("Combined finite track table has the wrong row count")
    if set(frames["targets"]["shape"].astype(str)) != set(shapes):
        raise ValueError("Combined target table has the wrong shapes")


def _write_manifest(
    store: ArtifactStore,
    stage: str,
    config: Config,
    settings: LineFeasibleSweepSettings,
    manifest: pd.DataFrame,
    frames: Mapping[str, pd.DataFrame],
    figures: Sequence[Path],
) -> Path:
    entries = []
    for name in tuple(frames) + ("fit_summary", "unbounded_ratios"):
        path = store.frame_path(stage, name)
        entries.append(
            {
                "artifact": name,
                "path": str(path.resolve()),
                "rows": int(pd.read_json(path.with_suffix(".meta.json"), typ="series")["rows"]),
                "sha256": sha256_file(path),
            }
        )
    payload = {
        "algorithm_version": ALGORITHM_VERSION,
        "config_hash": config.hash(),
        "normalization": "length/c and area/c^2",
        "figure_annotation": "condition label and fitted c in micrometres only",
        "grid_rule": (
            "same projected areas for Ball, Disk, and Line; areas retained "
            "condition-wise using Line feasibility"
        ),
        "runtime_settings": {
            "finite_replays_per_track": settings.finite_replays_per_track,
            "unbounded_trials_per_cell": settings.unbounded_trials_per_cell,
            "workers": settings.workers,
            "random_seed": settings.random_seed,
        },
        "conditions": int(len(manifest)),
        "target_cells": int(manifest["n_target_cells"].sum()),
        "tables": entries,
        "figures": [
            {
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in figures
        ],
    }
    return store.write_json(stage, "manifest", payload)


def run_line_feasible_sweep(
    config: Config,
    settings: LineFeasibleSweepSettings,
) -> Dict[str, object]:
    """Run, validate, combine, and plot the approved short target sweep."""
    settings.validate(config.replay.detection_radius)
    store = ArtifactStore(config.output.root)
    fits = _fit_table(config, store)
    manifest = _manifest(config, fits, settings)
    write_runtime_plan(config, settings)

    combined: Dict[str, List[pd.DataFrame]] = {
        name: [] for name in COMBINED_FRAME_NAMES
    }
    condition_rows: List[Dict[str, object]] = []
    for index, fit in fits.iterrows():
        condition = str(fit["condition"])
        manifest_row = manifest[manifest["condition"] == condition].iloc[0]
        areas = tuple(
            float(value)
            for value in str(manifest_row["projected_areas_model_units2"]).split(",")
        )
        condition_id = _condition_id(condition)
        condition_stage = "{}/condition_replays/{}".format(
            settings.output_stage,
            condition_id,
        )
        signature = _experiment_signature(
            config,
            settings,
            condition,
            fit,
            areas,
        )
        checkpoint = (
            _load_checkpoint(
                store,
                condition_stage,
                config,
                settings,
                areas,
                signature,
            )
            if settings.resume
            else None
        )
        if checkpoint is None:
            logger.info(
                "Running line-feasible condition %d/%d: %s",
                index + 1,
                len(fits),
                condition,
            )
            focused = FocusedReplaySettings(
                condition=condition,
                length_unit_um=None,
                projected_areas=areas,
                finite_replays_per_track=settings.finite_replays_per_track,
                unbounded_trials_per_cell=settings.unbounded_trials_per_cell,
                workers=settings.workers,
                finite_batch_size=settings.finite_batch_size,
                unbounded_chunk_trials=settings.unbounded_chunk_trials,
                unbounded_step_block=settings.unbounded_step_block,
                random_seed=_condition_seed(settings.random_seed, condition),
                output_stage=condition_stage,
                figure_format="pdf",
                replay_signature=signature,
                replay_signature_version=1,
            )
            summary = run_focused_replay(config, focused)
            frames = {
                name: store.read_frame(
                    condition_stage,
                    name,
                    expected_config_hash=config.hash(),
                )
                for name in COMBINED_FRAME_NAMES
            }
            payload = store.read_json(condition_stage, "settings")
            _validate_condition_frames(
                frames,
                summary,
                payload,
                config.replay.shapes,
                areas,
                settings.finite_replays_per_track,
                settings.unbounded_trials_per_cell,
                signature,
            )
        else:
            logger.info(
                "Reusing validated line-feasible condition %d/%d: %s",
                index + 1,
                len(fits),
                condition,
            )
            summary, frames = checkpoint

        crossover = float(fit["crossover"])
        for name, frame in frames.items():
            attached = _attach_condition_fields(frame, condition)
            if name == "targets":
                attached = _enhance_targets(attached, crossover)
            combined[name].append(attached)
        condition_row = dict(_condition_fields(condition))
        condition_row.update(summary)
        condition_row.update(
            {
                "crossover_um": crossover,
                "torus_side_model_units": float(manifest_row["torus_side_model_units"]),
                "torus_side_um": float(manifest_row["torus_side_um"]),
                "n_areas": int(manifest_row["n_areas"]),
                "experiment_signature": signature,
            }
        )
        condition_rows.append(condition_row)

    frames = {
        name: pd.concat(parts, ignore_index=True)
        for name, parts in combined.items()
    }
    _combined_grid_validation(
        manifest,
        frames,
        settings,
        config.replay.shapes,
    )

    stage = settings.output_stage
    metadata = {
        "config_hash": config.hash(),
        "algorithm_version": ALGORITHM_VERSION,
        "normalization": "fitted_crossover",
    }
    for name, frame in frames.items():
        store.write_frame(stage, name, frame, extra_metadata=metadata)

    fit_summary = fits.copy()
    for key in (
        "condition_id",
        "condition_label",
        "cell_type",
        "organ",
        "stimulus",
        "dt_seconds",
    ):
        fit_summary[key] = [
            _condition_fields(str(condition))[key]
            for condition in fit_summary["condition"]
        ]
    fit_summary = fit_summary.merge(
        manifest[
            [
                "condition",
                "torus_side_model_units",
                "torus_side_um",
                "line_projected_area_limit_model_units2",
                "projected_areas_model_units2",
                "n_areas",
                "n_target_cells",
            ]
        ],
        on="condition",
        how="left",
        validate="one_to_one",
    )
    ratios = _unbounded_ratio_table(frames["unbounded_summary"])
    store.write_frame(stage, "fit_summary", fit_summary, extra_metadata=metadata)
    store.write_frame(stage, "unbounded_ratios", ratios, extra_metadata=metadata)
    store.write_frame(
        stage,
        "condition_summaries",
        pd.DataFrame(condition_rows),
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "condition_target_manifest",
        manifest,
        extra_metadata=metadata,
    )

    from .line_feasible_figures import plot_line_feasible_figures

    figures = plot_line_feasible_figures(
        store.stage_dir(stage),
        fit_summary,
        frames["finite_summary"],
        frames["unbounded_summary"],
        ratios,
    )
    _write_manifest(
        store,
        stage,
        config,
        settings,
        manifest,
        frames,
        figures,
    )
    summary = {
        "algorithm_version": ALGORITHM_VERSION,
        "conditions": int(len(fits)),
        "selected_figure2_conditions": int(
            fit_summary["condition"].isin(SELECTED_CONDITIONS).sum()
        ),
        "target_cells": int(len(frames["targets"])),
        "finite_track_summary_rows": int(len(frames["finite_track_summaries"])),
        "finite_trial_count": int(frames["finite_track_summaries"]["n_trials"].sum()),
        "unbounded_trial_count": int(len(frames["unbounded_trials"])),
        "figures": [str(path.resolve()) for path in figures],
        "output_stage": stage,
        "validated": True,
    }
    store.write_json(stage, "summary", summary)
    return summary


__all__ = [
    "ALGORITHM_VERSION",
    "DEFAULT_CANDIDATE_AREAS",
    "DEFAULT_OUTPUT_STAGE",
    "LineFeasibleSweepSettings",
    "line_feasible_areas",
    "line_projected_area_limit",
    "run_line_feasible_sweep",
    "write_runtime_plan",
]
