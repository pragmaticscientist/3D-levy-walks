"""Additive finite replay with independently rotated vectors from each track.

This module deliberately does not modify or replace the baseline
``all_conditions_replay`` stage.  It verifies the signed per-condition
checkpoints behind that stage, runs only ``ordered_length_rotated``, and writes
both the new results and baseline-plus-new augmented tables to a separate
stage.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactStore
from .logging_utils import get_logger
from .multicondition_experiment import (
    _artifact_sha256,
    _attach_condition,
    _canonical_signature,
    _condition_fields,
    _condition_seed,
    _plot_all_conditions,
)
from .nk_experiment import _finite_track_worker, _track_bootstrap_summary
from .targets import Target

logger = get_logger(__name__)

NEW_CLASS = "ordered_length_rotated"
BASELINE_STAGE = "all_conditions_replay"
OUTPUT_STAGE = "ordered_length_rotation_replay"
BASELINE_REPLAY_SIGNATURE_VERSION = 1
BASELINE_REPLAY_ALGORITHM_VERSION = "endpoint_replay_v1"
EXPERIMENT_SIGNATURE_VERSION = 1
EXPERIMENT_ALGORITHM_VERSION = "ordered_length_rotation_replay_v1"

BASELINE_FINITE_CLASSES: Tuple[str, ...] = (
    "exact_orientation",
    "exact_global_rotation",
    "vector_uniform_rotated",
    "fitted_mu",
)
BASELINE_UNBOUNDED_CLASSES: Tuple[str, ...] = (
    "vector_uniform_rotated",
    "fitted_mu",
)
BASELINE_FINITE_CONTRASTS: Tuple[Tuple[str, str], ...] = (
    ("exact_global_rotation", "exact_orientation"),
    ("vector_uniform_rotated", "exact_global_rotation"),
    ("fitted_mu", "vector_uniform_rotated"),
)
NEW_FINITE_CONTRASTS: Tuple[Tuple[str, str], ...] = (
    (NEW_CLASS, "exact_global_rotation"),
    ("vector_uniform_rotated", NEW_CLASS),
)
NEW_FINITE_CONTRAST_ROLES = {
    (NEW_CLASS, "exact_global_rotation"): "directional_structure",
    (
        "vector_uniform_rotated",
        NEW_CLASS,
    ): "condition_pool_and_length_sequence_composite",
}

_SIGNED_FRAME_NAMES: Tuple[str, ...] = (
    "track_inputs",
    "targets",
    "finite_track_summaries",
    "finite_summary",
    "finite_contrasts",
)


@dataclass(frozen=True)
class OrderedRotationReplaySettings:
    """Runtime controls for the isolated additive experiment."""

    finite_replays_per_track: int = 10_000
    workers: int = 32
    finite_batch_size: int = 1_000
    random_seed: int = 913_579
    baseline_stage: str = BASELINE_STAGE
    output_stage: str = OUTPUT_STAGE

    def validate(self) -> None:
        if self.finite_replays_per_track < 1:
            raise ValueError("finite_replays_per_track must be positive")
        if self.workers < 0:
            raise ValueError("workers must be non-negative")
        if self.finite_batch_size < 1:
            raise ValueError("finite_batch_size must be positive")
        if self.random_seed < 0:
            raise ValueError("random_seed must be non-negative")
        if not self.baseline_stage.strip() or not self.output_stage.strip():
            raise ValueError("stage names must be non-empty")
        if self.baseline_stage == self.output_stage:
            raise ValueError("output_stage must differ from baseline_stage")

    @property
    def resolved_workers(self) -> int:
        available = os.cpu_count() or 1
        return min(available, self.workers) if self.workers > 0 else available


@dataclass(frozen=True)
class _BaselineInputs:
    """Validated source tables and provenance for the additive replay."""

    settings: Mapping[str, object]
    signatures: Mapping[str, str]
    source_hashes: Mapping[str, str]
    runs: pd.DataFrame
    pools: pd.DataFrame
    fits: pd.DataFrame
    track_inputs: pd.DataFrame
    targets: pd.DataFrame
    finite_track_summaries: pd.DataFrame
    finite_summary: pd.DataFrame
    finite_contrasts: pd.DataFrame


def _read_frame_metadata(
    store: ArtifactStore, stage: str, name: str
) -> Mapping[str, object]:
    path = store.frame_path(stage, name).with_suffix(".meta.json")
    if not path.is_file():
        raise ValueError("Missing artifact sidecar: {}".format(path))
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("Artifact sidecar is not a JSON object: {}".format(path))
    return value


def _baseline_signature_payload(
    config: Config,
    baseline_settings: Mapping[str, object],
    condition: str,
    fit: Mapping[str, object],
    source_hashes: Mapping[str, str],
) -> Dict[str, object]:
    """Reconstruct the semantic payload used by the immutable baseline."""
    crossover = float(fit["crossover"])
    return {
        "signature_version": BASELINE_REPLAY_SIGNATURE_VERSION,
        "algorithm_version": BASELINE_REPLAY_ALGORITHM_VERSION,
        "config_hash": config.hash(),
        "source_artifact_hashes": dict(source_hashes),
        "condition": condition,
        "weighting": config.empirical_pool.primary_weighting,
        "normalization": "fitted_crossover",
        "length_unit_um": crossover,
        "fitted_mu": float(fit["mu_hat"]),
        "fitted_crossover_model_units": 1.0,
        "fitted_lmax_model_units": float(fit["lmax"]) / crossover,
        "torus_width_rule": "2*fitted_lmax_model_units",
        "detection_radius_model_units": float(config.replay.detection_radius),
        "shapes": list(config.replay.shapes),
        "projected_areas": [
            float(value) for value in baseline_settings["projected_areas"]
        ],
        "finite_classes": list(BASELINE_FINITE_CLASSES),
        "unbounded_classes": list(BASELINE_UNBOUNDED_CLASSES),
        "finite_replays_per_track": int(
            baseline_settings["finite_replays_per_track"]
        ),
        "unbounded_trials_per_cell": int(
            baseline_settings["unbounded_trials_per_cell"]
        ),
        "finite_batch_size": int(baseline_settings["finite_batch_size"]),
        "unbounded_chunk_trials": int(
            baseline_settings["unbounded_chunk_trials"]
        ),
        "unbounded_step_block": int(
            baseline_settings["unbounded_step_block"]
        ),
        "random_seed": _condition_seed(
            int(baseline_settings["random_seed"]), condition
        ),
        "finite_bootstraps": int(config.fitting.clustered_bootstraps),
        "overshoot_policy": "discard",
        "check_initial_position": False,
        "figure_format": "pdf",
    }


def _sort_frame(frame: pd.DataFrame, keys: Sequence[str]) -> pd.DataFrame:
    return frame.sort_values(list(keys), kind="mergesort").reset_index(drop=True)


def _assert_frames_equivalent(
    combined: pd.DataFrame,
    signed: pd.DataFrame,
    keys: Sequence[str],
    label: str,
) -> None:
    """Require a combined baseline table to equal its signed source tables."""
    if set(combined.columns) != set(signed.columns):
        raise ValueError("{} column set does not match signed sources".format(label))
    columns = list(combined.columns)
    left = _sort_frame(combined[columns], keys)
    right = _sort_frame(signed[columns], keys)
    try:
        pd.testing.assert_frame_equal(
            left,
            right,
            check_dtype=False,
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
    except AssertionError as exc:
        raise ValueError(
            "{} does not match its signed per-condition sources".format(label)
        ) from exc


def _validate_class_grid(
    frame: pd.DataFrame,
    conditions: Sequence[str],
    classes: Sequence[str],
    shapes: Sequence[str],
    areas: Sequence[float],
    key_prefix: Sequence[str],
    label: str,
) -> None:
    expected_classes = set(classes)
    expected_shapes = set(shapes)
    expected_areas = sorted(float(value) for value in areas)
    if set(frame["condition"].astype(str)) != set(conditions):
        raise ValueError("{} has an unexpected condition set".format(label))
    if set(frame["trajectory_class"].astype(str)) != expected_classes:
        raise ValueError("{} has an unexpected trajectory class set".format(label))
    if set(frame["shape"].astype(str)) != expected_shapes:
        raise ValueError("{} has an unexpected target shape set".format(label))
    actual_areas = sorted(frame["projected_area"].astype(float).unique())
    if len(actual_areas) != len(expected_areas) or not np.allclose(
        actual_areas, expected_areas, rtol=0.0, atol=1e-12
    ):
        raise ValueError("{} has an unexpected projected-area grid".format(label))
    keys = list(key_prefix) + [
        "trajectory_class",
        "shape",
        "projected_area",
    ]
    if frame.duplicated(keys).any():
        raise ValueError("{} contains duplicate grid cells".format(label))


def _validate_baseline_grids(
    settings: OrderedRotationReplaySettings,
    baseline_settings: Mapping[str, object],
    condition_summaries: pd.DataFrame,
    track_inputs: pd.DataFrame,
    targets: pd.DataFrame,
    finite_tracks: pd.DataFrame,
    finite_summary: pd.DataFrame,
    finite_contrasts: pd.DataFrame,
    shapes: Sequence[str],
) -> None:
    conditions = sorted(condition_summaries["condition"].astype(str).unique())
    areas = [float(value) for value in baseline_settings["projected_areas"]]
    n_cells = len(shapes) * len(areas)
    if int(baseline_settings["finite_replays_per_track"]) != (
        settings.finite_replays_per_track
    ):
        raise ValueError(
            "The additive replay count must equal the signed baseline count"
        )
    if track_inputs.duplicated(["condition", "track_uid"]).any():
        raise ValueError("Baseline track_inputs contains duplicate tracks")
    if set(track_inputs["condition"].astype(str)) != set(conditions):
        raise ValueError("Baseline track_inputs has an unexpected condition set")
    if targets.duplicated(["condition", "shape", "projected_area"]).any():
        raise ValueError("Baseline targets contains duplicate cells")
    expected_target_rows = len(conditions) * n_cells
    if len(targets) != expected_target_rows:
        raise ValueError(
            "Baseline targets has {} rows; expected {}".format(
                len(targets), expected_target_rows
            )
        )
    if set(targets["shape"].astype(str)) != set(shapes):
        raise ValueError("Baseline targets has an unexpected shape set")

    _validate_class_grid(
        finite_tracks,
        conditions,
        BASELINE_FINITE_CLASSES,
        shapes,
        areas,
        ("condition", "track_uid"),
        "baseline finite_track_summaries",
    )
    expected_track_rows = len(track_inputs) * len(BASELINE_FINITE_CLASSES) * n_cells
    if len(finite_tracks) != expected_track_rows:
        raise ValueError(
            "Baseline finite_track_summaries has {} rows; expected {}".format(
                len(finite_tracks), expected_track_rows
            )
        )
    if not np.all(
        finite_tracks["n_trials"].to_numpy(dtype=int)
        == settings.finite_replays_per_track
    ):
        raise ValueError("Baseline finite_track_summaries has bad trial counts")

    _validate_class_grid(
        finite_summary,
        conditions,
        BASELINE_FINITE_CLASSES,
        shapes,
        areas,
        ("condition",),
        "baseline finite_summary",
    )
    if len(finite_summary) != (
        len(conditions) * len(BASELINE_FINITE_CLASSES) * n_cells
    ):
        raise ValueError("Baseline finite_summary has an incorrect row count")

    expected_pairs = set(BASELINE_FINITE_CONTRASTS)
    actual_pairs = {
        (str(left), str(right))
        for left, right in finite_contrasts[
            ["left_class", "right_class"]
        ].itertuples(index=False, name=None)
    }
    if actual_pairs != expected_pairs:
        raise ValueError("Baseline finite_contrasts has unexpected class pairs")
    contrast_keys = [
        "condition",
        "left_class",
        "right_class",
        "shape",
        "projected_area",
    ]
    if finite_contrasts.duplicated(contrast_keys).any() or len(
        finite_contrasts
    ) != len(conditions) * len(expected_pairs) * n_cells:
        raise ValueError("Baseline finite_contrasts has an invalid grid")


def _load_signed_baseline(
    config: Config, settings: OrderedRotationReplaySettings
) -> _BaselineInputs:
    """Load, hash-check, signature-check, and reconcile the baseline."""
    store = ArtifactStore(config.output.root)
    config_hash = config.hash()
    baseline_settings = store.read_json(settings.baseline_stage, "settings")
    if not isinstance(baseline_settings, dict):
        raise ValueError("Baseline settings must be a JSON object")
    required_settings = {
        "projected_areas",
        "finite_replays_per_track",
        "unbounded_trials_per_cell",
        "finite_batch_size",
        "unbounded_chunk_trials",
        "unbounded_step_block",
        "random_seed",
    }
    missing_settings = sorted(required_settings.difference(baseline_settings))
    if missing_settings:
        raise ValueError(
            "Baseline settings is missing: {}".format(
                ", ".join(missing_settings)
            )
        )

    runs = store.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config_hash,
    )
    pools = store.read_frame(
        "preprocess", "pools", expected_config_hash=config_hash
    )
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config_hash,
    )
    fits = fits[fits["status"] == "ok"].sort_values("condition").reset_index(
        drop=True
    )
    source_hashes = {
        "runs": _artifact_sha256(store, "preprocess", "runs"),
        "pools": _artifact_sha256(store, "preprocess", "pools"),
        "mixture_fits": _artifact_sha256(store, "fit", "mixture_fits"),
    }

    combined = {
        name: store.read_frame(
            settings.baseline_stage,
            name,
            expected_config_hash=config_hash,
        )
        for name in _SIGNED_FRAME_NAMES
    }
    condition_summaries = store.read_frame(
        settings.baseline_stage,
        "condition_summaries",
        expected_config_hash=config_hash,
    )
    baseline_fit_summary = store.read_frame(
        settings.baseline_stage,
        "fit_summary",
        expected_config_hash=config_hash,
    )
    fit_columns = ["condition", "status", "mu_hat", "crossover", "lmax"]
    _assert_frames_equivalent(
        baseline_fit_summary[fit_columns],
        fits[fit_columns],
        ("condition",),
        "baseline fit_summary",
    )
    if condition_summaries.duplicated(["condition"]).any():
        raise ValueError("Baseline condition_summaries contains duplicates")
    if set(condition_summaries["condition"].astype(str)) != set(
        fits["condition"].astype(str)
    ):
        raise ValueError("Baseline conditions do not match successful fits")

    fit_map = fits.set_index("condition")
    signatures: Dict[str, str] = {}
    signed_parts: Dict[str, List[pd.DataFrame]] = {
        name: [] for name in _SIGNED_FRAME_NAMES
    }
    for _, condition_row in condition_summaries.sort_values(
        "condition"
    ).iterrows():
        condition = str(condition_row["condition"])
        condition_id = str(condition_row["condition_id"])
        condition_stage = "condition_replays/{}".format(condition_id)
        fit = fit_map.loc[condition]
        expected_signature = _canonical_signature(
            _baseline_signature_payload(
                config,
                baseline_settings,
                condition,
                fit,
                source_hashes,
            )
        )
        signatures[condition] = expected_signature
        if (
            str(condition_row["replay_signature"]) != expected_signature
            or int(condition_row["replay_signature_version"])
            != BASELINE_REPLAY_SIGNATURE_VERSION
        ):
            raise ValueError(
                "Baseline condition summary signature mismatch for {}".format(
                    condition
                )
            )
        signed_settings = store.read_json(condition_stage, "settings")
        signed_summary = store.read_json(condition_stage, "summary")
        for value, label in (
            (signed_settings, "settings"),
            (signed_summary, "summary"),
        ):
            if not isinstance(value, dict):
                raise ValueError(
                    "Signed condition {} is not an object: {}".format(
                        label, condition
                    )
                )
            if (
                value.get("replay_signature") != expected_signature
                or value.get("replay_signature_version")
                != BASELINE_REPLAY_SIGNATURE_VERSION
            ):
                raise ValueError(
                    "Signed condition {} signature mismatch for {}".format(
                        label, condition
                    )
                )
        if tuple(signed_settings.get("finite_classes", ())) != (
            BASELINE_FINITE_CLASSES
        ):
            raise ValueError(
                "Signed baseline class ladder is unexpected for {}".format(
                    condition
                )
            )
        for name in _SIGNED_FRAME_NAMES:
            frame = store.read_frame(
                condition_stage,
                name,
                expected_config_hash=config_hash,
            )
            metadata = _read_frame_metadata(store, condition_stage, name)
            if (
                metadata.get("replay_signature") != expected_signature
                or metadata.get("replay_signature_version")
                != BASELINE_REPLAY_SIGNATURE_VERSION
            ):
                raise ValueError(
                    "Signed {} metadata mismatch for {}".format(name, condition)
                )
            signed_parts[name].append(_attach_condition(frame, condition))

    signed_combined = {
        name: pd.concat(parts, ignore_index=True)
        for name, parts in signed_parts.items()
    }
    comparison_keys = {
        "track_inputs": ("condition", "track_uid"),
        "targets": ("condition", "shape", "projected_area"),
        "finite_track_summaries": (
            "condition",
            "track_uid",
            "trajectory_class",
            "shape",
            "projected_area",
        ),
        "finite_summary": (
            "condition",
            "trajectory_class",
            "shape",
            "projected_area",
        ),
        "finite_contrasts": (
            "condition",
            "left_class",
            "right_class",
            "shape",
            "projected_area",
        ),
    }
    for name in _SIGNED_FRAME_NAMES:
        _assert_frames_equivalent(
            combined[name],
            signed_combined[name],
            comparison_keys[name],
            "baseline {}".format(name),
        )

    _validate_baseline_grids(
        settings,
        baseline_settings,
        condition_summaries,
        combined["track_inputs"],
        combined["targets"],
        combined["finite_track_summaries"],
        combined["finite_summary"],
        combined["finite_contrasts"],
        config.replay.shapes,
    )
    return _BaselineInputs(
        settings=baseline_settings,
        signatures=signatures,
        source_hashes=source_hashes,
        runs=runs,
        pools=pools,
        fits=fits,
        track_inputs=combined["track_inputs"],
        targets=combined["targets"],
        finite_track_summaries=combined["finite_track_summaries"],
        finite_summary=combined["finite_summary"],
        finite_contrasts=combined["finite_contrasts"],
    )


def _targets_for_condition(
    target_frame: pd.DataFrame, shape_order: Sequence[str]
) -> List[Tuple[float, Target]]:
    rows: List[Tuple[float, Target]] = []
    by_key = target_frame.set_index(["shape", "projected_area"])
    for area in sorted(target_frame["projected_area"].astype(float).unique()):
        for shape in shape_order:
            try:
                target = by_key.loc[(shape, area)]
            except KeyError as exc:
                raise ValueError(
                    "Missing target {} at area {}".format(shape, area)
                ) from exc
            if isinstance(target, pd.DataFrame):
                raise ValueError("Target grid contains duplicate cells")
            rows.append(
                (
                    float(area),
                    Target(
                        shape=str(shape),
                        dimension=float(target["dimension"]),
                        side_length=float(target["side_length"]),
                        detection_radius=float(target["detection_radius"]),
                    ),
                )
            )
    return rows


def _build_condition_tasks(
    config: Config,
    settings: OrderedRotationReplaySettings,
    condition: str,
    fit: Mapping[str, object],
    condition_runs: pd.DataFrame,
    condition_pool: pd.DataFrame,
    baseline_tracks: pd.DataFrame,
    baseline_targets: pd.DataFrame,
) -> List[Dict[str, object]]:
    """Rebuild the exact baseline track sequences and one-class worker tasks."""
    vector_columns = ["run_vector_x", "run_vector_y", "run_vector_z"]
    crossover = float(fit["crossover"])
    fitted_mu = float(fit["mu_hat"])
    fitted_crossover = 1.0
    fitted_lmax = float(fit["lmax"]) / crossover
    pool_vectors = (
        condition_pool[vector_columns].to_numpy(dtype=float) / crossover
    )
    if len(pool_vectors) == 0 or np.any(~np.isfinite(pool_vectors)):
        raise ValueError("Condition has no finite empirical pool: {}".format(condition))
    targets = _targets_for_condition(
        baseline_targets, config.replay.shapes
    )
    side_lengths = {float(target.side_length) for _, target in targets}
    if len(side_lengths) != 1:
        raise ValueError("Condition targets do not share one torus width")
    side_length = next(iter(side_lengths))
    if not np.isclose(
        side_length, 2.0 * fitted_lmax, rtol=1e-12, atol=1e-12
    ):
        raise ValueError("Condition target torus does not equal 2 Lmax/c")

    expected_tracks = baseline_tracks.set_index("track_uid")
    if expected_tracks.index.has_duplicates:
        raise ValueError("Baseline condition track_inputs contains duplicates")
    tasks_without_seed: List[Dict[str, object]] = []
    for track_uid, group in condition_runs.groupby("track_uid", sort=True):
        track_uid_string = str(track_uid)
        if track_uid_string not in expected_tracks.index:
            continue
        ordered = group.sort_values(
            ["frame_from", "frame_to", "fragment_uid", "run_id"],
            kind="mergesort",
        )
        vectors = (
            ordered[vector_columns].to_numpy(dtype=float) / crossover
        )
        lengths = np.linalg.norm(vectors, axis=1)
        if len(vectors) == 0 or np.any(~np.isfinite(lengths)) or np.any(
            lengths <= 0.0
        ):
            raise ValueError("Track contains invalid vectors: {}".format(track_uid))
        expected = expected_tracks.loc[track_uid_string]
        budget = float(np.sum(lengths))
        if int(expected["n_runs"]) != len(vectors) or not np.isclose(
            float(expected["budget_model_units"]),
            budget,
            rtol=1e-11,
            atol=1e-11,
        ):
            raise ValueError(
                "Rebuilt exact track differs from baseline: {}".format(
                    track_uid_string
                )
            )
        video_id = str(ordered["video_id"].iloc[0])
        if video_id != str(expected["video_id"]):
            raise ValueError(
                "Rebuilt track video differs from baseline: {}".format(
                    track_uid_string
                )
            )
        tasks_without_seed.append(
            {
                "condition": condition,
                "video_id": video_id,
                "track_uid": track_uid_string,
                "source_vectors": vectors,
                "pool_vectors": pool_vectors,
                "targets": targets,
                "n_replays": settings.finite_replays_per_track,
                "batch_size": settings.finite_batch_size,
                "side_length": side_length,
                "fitted_mu": fitted_mu,
                "fitted_crossover": fitted_crossover,
                "fitted_lmax": fitted_lmax,
                "finite_classes": (NEW_CLASS,),
            }
        )
    rebuilt_ids = {str(task["track_uid"]) for task in tasks_without_seed}
    expected_ids = set(expected_tracks.index.astype(str))
    if rebuilt_ids != expected_ids:
        missing = sorted(expected_ids.difference(rebuilt_ids))
        raise ValueError(
            "Could not rebuild {} baseline tracks for {}: {}".format(
                len(missing), condition, ", ".join(missing[:5])
            )
        )
    seed_sequence = np.random.SeedSequence(
        _condition_seed(settings.random_seed, condition)
    )
    seeds = seed_sequence.spawn(len(tasks_without_seed))
    tasks: List[Dict[str, object]] = []
    for task, child_seed in zip(tasks_without_seed, seeds):
        complete = dict(task)
        complete["seed"] = int(
            child_seed.generate_state(1, dtype=np.uint64)[0]
        )
        tasks.append(complete)
    return tasks


def _attach_condition_fields(frame: pd.DataFrame) -> pd.DataFrame:
    parts = [
        _attach_condition(group.copy(), str(condition))
        for condition, group in frame.groupby("condition", sort=True)
    ]
    if not parts:
        return frame.copy()
    return pd.concat(parts, ignore_index=True)


def _paired_track_contrasts(
    tracks: pd.DataFrame,
    n_bootstraps: int,
    random_seed: int,
) -> pd.DataFrame:
    """Compute the two additive paired, whole-track probability contrasts."""
    if n_bootstraps < 1:
        raise ValueError("n_bootstraps must be positive")
    rows: List[Dict[str, object]] = []
    for condition, condition_frame in tracks.groupby("condition", sort=True):
        rng = np.random.default_rng(
            _condition_seed(random_seed, str(condition))
        )
        for keys, cell in condition_frame.groupby(
            ["shape", "projected_area"], sort=True
        ):
            pivot = cell.pivot(
                index="track_uid",
                columns="trajectory_class",
                values="detection_probability",
            )
            for left, right in NEW_FINITE_CONTRASTS:
                if left not in pivot or right not in pivot:
                    raise ValueError(
                        "Missing {} or {} for {}".format(left, right, condition)
                    )
                paired = (pivot[left] - pivot[right]).dropna().to_numpy(
                    dtype=float
                )
                if len(paired) != len(pivot):
                    raise ValueError(
                        "Contrast does not contain every track for {}".format(
                            condition
                        )
                    )
                sampled = rng.integers(
                    0, len(paired), size=(n_bootstraps, len(paired))
                )
                draws = np.mean(paired[sampled], axis=1)
                low, high = np.quantile(draws, [0.025, 0.975])
                rows.append(
                    {
                        "condition": str(condition),
                        "left_class": left,
                        "right_class": right,
                        "contrast": "{} minus {}".format(left, right),
                        "contrast_role": NEW_FINITE_CONTRAST_ROLES[
                            (left, right)
                        ],
                        "shape": keys[0],
                        "projected_area": float(keys[1]),
                        "n_tracks": int(len(paired)),
                        "detection_probability_difference": float(
                            np.mean(paired)
                        ),
                        "difference_ci_low": float(low),
                        "difference_ci_high": float(high),
                        "n_bootstraps": int(n_bootstraps),
                    }
                )
    return _attach_condition_fields(pd.DataFrame(rows))


def _summarise_new_tracks(
    tracks: pd.DataFrame,
    n_bootstraps: int,
    random_seed: int,
) -> pd.DataFrame:
    parts: List[pd.DataFrame] = []
    for condition, condition_frame in tracks.groupby("condition", sort=True):
        summary = _track_bootstrap_summary(
            condition_frame,
            n_bootstraps,
            np.random.default_rng(
                _condition_seed(random_seed, str(condition))
            ),
        )
        summary["condition"] = str(condition)
        parts.append(_attach_condition(summary, str(condition)))
    return pd.concat(parts, ignore_index=True)


def _validate_new_outputs(
    settings: OrderedRotationReplaySettings,
    shapes: Sequence[str],
    areas: Sequence[float],
    baseline_tracks: pd.DataFrame,
    baseline_summary: pd.DataFrame,
    baseline_contrasts: pd.DataFrame,
    new_tracks: pd.DataFrame,
    new_summary: pd.DataFrame,
    new_contrasts: pd.DataFrame,
    augmented_summary: pd.DataFrame,
    augmented_contrasts: pd.DataFrame,
    condition_summaries: pd.DataFrame,
) -> None:
    """Validate exact additive and augmented row grids before any write."""
    conditions = sorted(baseline_tracks["condition"].astype(str).unique())
    n_cells = len(shapes) * len(areas)
    n_tracks = baseline_tracks[
        ["condition", "track_uid"]
    ].drop_duplicates()
    _validate_class_grid(
        new_tracks,
        conditions,
        (NEW_CLASS,),
        shapes,
        areas,
        ("condition", "track_uid"),
        "new finite_track_summaries",
    )
    if len(new_tracks) != len(n_tracks) * n_cells:
        raise ValueError("New finite_track_summaries has an incorrect row count")
    if not np.all(
        new_tracks["n_trials"].to_numpy(dtype=int)
        == settings.finite_replays_per_track
    ):
        raise ValueError("New finite_track_summaries has bad trial counts")
    n_trials = new_tracks["n_trials"].to_numpy(dtype=int)
    n_detected = new_tracks["n_detected"].to_numpy(dtype=int)
    probabilities = new_tracks["detection_probability"].to_numpy(dtype=float)
    if (
        np.any(n_detected < 0)
        or np.any(n_detected > n_trials)
        or not np.allclose(
            probabilities,
            n_detected / n_trials,
            rtol=0.0,
            atol=1e-12,
        )
    ):
        raise ValueError("New finite_track_summaries has invalid detections")

    _validate_class_grid(
        new_summary,
        conditions,
        (NEW_CLASS,),
        shapes,
        areas,
        ("condition",),
        "new finite_summary",
    )
    if len(new_summary) != len(conditions) * n_cells:
        raise ValueError("New finite_summary has an incorrect row count")

    expected_new_pairs = set(NEW_FINITE_CONTRASTS)
    actual_new_pairs = {
        (str(left), str(right))
        for left, right in new_contrasts[
            ["left_class", "right_class"]
        ].itertuples(index=False, name=None)
    }
    contrast_keys = [
        "condition",
        "left_class",
        "right_class",
        "shape",
        "projected_area",
    ]
    if (
        actual_new_pairs != expected_new_pairs
        or new_contrasts.duplicated(contrast_keys).any()
        or len(new_contrasts)
        != len(conditions) * len(expected_new_pairs) * n_cells
    ):
        raise ValueError("New finite_contrasts has an invalid row grid")

    summary_keys = [
        "condition",
        "trajectory_class",
        "shape",
        "projected_area",
    ]
    expected_augmented_summary_rows = len(baseline_summary) + len(new_summary)
    if (
        len(augmented_summary) != expected_augmented_summary_rows
        or augmented_summary.duplicated(summary_keys).any()
        or set(augmented_summary["trajectory_class"].astype(str))
        != set(BASELINE_FINITE_CLASSES).union({NEW_CLASS})
    ):
        raise ValueError("Augmented finite_summary has an invalid row grid")
    expected_augmented_contrast_rows = len(baseline_contrasts) + len(
        new_contrasts
    )
    if (
        len(augmented_contrasts) != expected_augmented_contrast_rows
        or augmented_contrasts.duplicated(contrast_keys).any()
    ):
        raise ValueError("Augmented finite_contrasts has an invalid row grid")
    _assert_frames_equivalent(
        augmented_summary[
            augmented_summary["trajectory_class"].isin(BASELINE_FINITE_CLASSES)
        ][baseline_summary.columns],
        baseline_summary,
        summary_keys,
        "augmented baseline finite_summary rows",
    )
    baseline_pair_set = set(BASELINE_FINITE_CONTRASTS)
    augmented_baseline_contrasts = augmented_contrasts[
        [
            (str(left), str(right)) in baseline_pair_set
            for left, right in augmented_contrasts[
                ["left_class", "right_class"]
            ].itertuples(index=False, name=None)
        ]
    ]
    _assert_frames_equivalent(
        augmented_baseline_contrasts[baseline_contrasts.columns],
        baseline_contrasts,
        contrast_keys,
        "augmented baseline finite_contrast rows",
    )
    if (
        condition_summaries.duplicated(["condition"]).any()
        or set(condition_summaries["condition"].astype(str)) != set(conditions)
        or len(condition_summaries) != len(conditions)
    ):
        raise ValueError("Condition summaries has an invalid row grid")


def _experiment_signature(
    config: Config,
    settings: OrderedRotationReplaySettings,
    baseline: _BaselineInputs,
) -> str:
    payload = {
        "signature_version": EXPERIMENT_SIGNATURE_VERSION,
        "algorithm_version": EXPERIMENT_ALGORITHM_VERSION,
        "config_hash": config.hash(),
        "source_artifact_hashes": dict(baseline.source_hashes),
        "baseline_stage": settings.baseline_stage,
        "baseline_condition_signatures": dict(baseline.signatures),
        "trajectory_class": NEW_CLASS,
        "contrasts": [list(pair) for pair in NEW_FINITE_CONTRASTS],
        "finite_replays_per_track": settings.finite_replays_per_track,
        "finite_batch_size": settings.finite_batch_size,
        "random_seed": settings.random_seed,
        "clustered_bootstraps": config.fitting.clustered_bootstraps,
        "projected_areas": [
            float(value) for value in baseline.settings["projected_areas"]
        ],
        "shapes": list(config.replay.shapes),
        "normalization": "each condition divided by its fitted crossover",
        "torus_width_rule": "2*fitted_lmax_model_units",
        "overshoot_policy": "discard",
        "check_initial_position": False,
    }
    return _canonical_signature(payload)


def run_ordered_length_rotation_replay(
    config: Config,
    settings: OrderedRotationReplaySettings,
) -> Dict[str, object]:
    """Run and persist only the ordered-length, independently rotated class."""
    settings.validate()
    baseline = _load_signed_baseline(config, settings)
    fit_map = baseline.fits.set_index("condition")
    all_tasks: List[Dict[str, object]] = []
    condition_task_counts: Dict[str, int] = {}
    for condition in sorted(baseline.signatures):
        pool = baseline.pools[
            (baseline.pools["condition"].astype(str) == condition)
            & np.isfinite(
                baseline.pools["run_length_um"].to_numpy(dtype=float)
            )
            & (baseline.pools["run_length_um"].to_numpy(dtype=float) > 0.0)
        ].copy()
        track_ids = set(
            baseline.track_inputs.loc[
                baseline.track_inputs["condition"].astype(str) == condition,
                "track_uid",
            ].astype(str)
        )
        runs = baseline.runs[
            (baseline.runs["condition"].astype(str) == condition)
            & baseline.runs["track_uid"].astype(str).isin(track_ids)
            & np.isfinite(
                baseline.runs["run_length_um"].to_numpy(dtype=float)
            )
            & (baseline.runs["run_length_um"].to_numpy(dtype=float) > 0.0)
        ].copy()
        tasks = _build_condition_tasks(
            config,
            settings,
            condition,
            fit_map.loc[condition],
            runs,
            pool,
            baseline.track_inputs[
                baseline.track_inputs["condition"].astype(str) == condition
            ],
            baseline.targets[
                baseline.targets["condition"].astype(str) == condition
            ],
        )
        condition_task_counts[condition] = len(tasks)
        all_tasks.extend(tasks)

    logger.info(
        "Ordered-length rotation replay: %d tracks x %d replicates on %d workers",
        len(all_tasks),
        settings.finite_replays_per_track,
        settings.resolved_workers,
    )
    rows: List[Dict[str, object]] = []
    with ProcessPoolExecutor(
        max_workers=settings.resolved_workers
    ) as executor:
        futures = [
            executor.submit(_finite_track_worker, task) for task in all_tasks
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            rows.extend(future.result())
            if completed % max(1, len(futures) // 10) == 0:
                logger.info(
                    "Ordered-length rotation completed %d/%d tracks",
                    completed,
                    len(futures),
                )
    new_tracks = _attach_condition_fields(pd.DataFrame(rows))
    new_tracks = _sort_frame(
        new_tracks,
        (
            "condition",
            "video_id",
            "track_uid",
            "trajectory_class",
            "shape",
            "projected_area",
        ),
    )
    new_summary = _summarise_new_tracks(
        new_tracks,
        int(config.fitting.clustered_bootstraps),
        settings.random_seed + 1,
    )
    contrast_input = pd.concat(
        [
            baseline.finite_track_summaries[
                baseline.finite_track_summaries["trajectory_class"].isin(
                    ("exact_global_rotation", "vector_uniform_rotated")
                )
            ],
            new_tracks,
        ],
        ignore_index=True,
    )
    new_contrasts = _paired_track_contrasts(
        contrast_input,
        int(config.fitting.clustered_bootstraps),
        settings.random_seed + 2,
    )
    augmented_summary = pd.concat(
        [baseline.finite_summary, new_summary], ignore_index=True
    )
    augmented_summary = _sort_frame(
        augmented_summary,
        (
            "condition",
            "trajectory_class",
            "shape",
            "projected_area",
        ),
    )
    augmented_contrasts = pd.concat(
        [baseline.finite_contrasts, new_contrasts], ignore_index=True
    )
    augmented_contrasts = _sort_frame(
        augmented_contrasts,
        (
            "condition",
            "left_class",
            "right_class",
            "shape",
            "projected_area",
        ),
    )
    signature = _experiment_signature(config, settings, baseline)
    condition_rows: List[Dict[str, object]] = []
    n_target_cells = len(config.replay.shapes) * len(
        baseline.settings["projected_areas"]
    )
    for condition in sorted(baseline.signatures):
        row = _condition_fields(condition)
        tracks = int(condition_task_counts[condition])
        row.update(
            {
                "trajectory_class": NEW_CLASS,
                "tracks": tracks,
                "target_cells": n_target_cells,
                "finite_replays_per_track": settings.finite_replays_per_track,
                "finite_trial_count": (
                    tracks
                    * n_target_cells
                    * settings.finite_replays_per_track
                ),
                "baseline_replay_signature": baseline.signatures[condition],
                "experiment_signature": signature,
                "experiment_signature_version": (
                    EXPERIMENT_SIGNATURE_VERSION
                ),
            }
        )
        condition_rows.append(row)
    condition_summaries = pd.DataFrame(condition_rows)

    _validate_new_outputs(
        settings,
        config.replay.shapes,
        [float(value) for value in baseline.settings["projected_areas"]],
        baseline.finite_track_summaries,
        baseline.finite_summary,
        baseline.finite_contrasts,
        new_tracks,
        new_summary,
        new_contrasts,
        augmented_summary,
        augmented_contrasts,
        condition_summaries,
    )

    store = ArtifactStore(config.output.root)
    metadata = {
        "config_hash": config.hash(),
        "weighting": config.empirical_pool.primary_weighting,
        "baseline_stage": settings.baseline_stage,
        "experiment_signature": signature,
        "experiment_signature_version": EXPERIMENT_SIGNATURE_VERSION,
        "trajectory_class": NEW_CLASS,
    }
    output_frames = {
        "finite_track_summaries": new_tracks,
        "finite_summary": new_summary,
        "finite_contrasts": new_contrasts,
        "condition_summaries": condition_summaries,
        "augmented_finite_summary": augmented_summary,
        "augmented_finite_contrasts": augmented_contrasts,
    }
    for name, frame in output_frames.items():
        store.write_frame(
            settings.output_stage,
            name,
            frame,
            extra_metadata=metadata,
        )
    settings_payload = {
        "algorithm_version": EXPERIMENT_ALGORITHM_VERSION,
        "experiment_signature": signature,
        "experiment_signature_version": EXPERIMENT_SIGNATURE_VERSION,
        "config_hash": config.hash(),
        "baseline_stage": settings.baseline_stage,
        "baseline_condition_signatures": dict(baseline.signatures),
        "output_stage": settings.output_stage,
        "trajectory_class": NEW_CLASS,
        "construction": (
            "each exact track vector is independently Haar-rotated on every "
            "replay; vector lengths, length order, count, and track budget are "
            "preserved"
        ),
        "finite_contrasts": [list(pair) for pair in NEW_FINITE_CONTRASTS],
        "finite_replays_per_track": settings.finite_replays_per_track,
        "finite_batch_size": settings.finite_batch_size,
        "workers": settings.resolved_workers,
        "random_seed": settings.random_seed,
        "clustered_bootstraps": config.fitting.clustered_bootstraps,
        "projected_areas": [
            float(value) for value in baseline.settings["projected_areas"]
        ],
        "shapes": list(config.replay.shapes),
        "normalization": "each condition divided by its fitted crossover",
        "torus_width": "2 * fitted Lmax in normalized units",
        "overshoot_policy": "discard",
        "check_initial_position": False,
        "source_artifact_hashes": dict(baseline.source_hashes),
    }
    store.write_json(settings.output_stage, "settings", settings_payload)
    fit_plot = store.read_frame(
        settings.baseline_stage,
        "fit_summary",
        expected_config_hash=config.hash(),
    )
    unbounded_summary = store.read_frame(
        settings.baseline_stage,
        "unbounded_summary",
        expected_config_hash=config.hash(),
    )
    unbounded_contrasts = store.read_frame(
        settings.baseline_stage,
        "unbounded_contrasts",
        expected_config_hash=config.hash(),
    )
    old_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "primary"
        / "fit"
        / "mixture_fits.csv"
    )
    old_fits = pd.read_csv(old_path) if old_path.is_file() else None
    figure_paths = _plot_all_conditions(
        store.stage_dir(settings.output_stage),
        fit_plot,
        baseline.pools,
        baseline.track_inputs,
        augmented_summary,
        augmented_contrasts,
        unbounded_summary,
        unbounded_contrasts,
        old_fits=old_fits,
    )
    summary: Dict[str, object] = {
        "conditions": int(len(condition_summaries)),
        "tracks": int(len(all_tasks)),
        "trajectory_class": NEW_CLASS,
        "finite_trial_count": int(new_tracks["n_trials"].sum()),
        "finite_track_summary_rows": int(len(new_tracks)),
        "finite_summary_rows": int(len(new_summary)),
        "finite_contrast_rows": int(len(new_contrasts)),
        "augmented_finite_summary_rows": int(len(augmented_summary)),
        "augmented_finite_contrast_rows": int(len(augmented_contrasts)),
        "experiment_signature": signature,
        "experiment_signature_version": EXPERIMENT_SIGNATURE_VERSION,
        "output_artifacts": sorted(output_frames),
        "figures": [str(path) for path in figure_paths],
    }
    store.write_json(settings.output_stage, "summary", summary)
    return summary


__all__ = [
    "NEW_CLASS",
    "NEW_FINITE_CONTRASTS",
    "OrderedRotationReplaySettings",
    "run_ordered_length_rotation_replay",
]
