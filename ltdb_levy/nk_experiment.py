"""Focused empirical-versus-synthetic replay for one fitted condition.

The experiment has two complementary parts.

Finite-budget replay compares five trajectory constructions:

``exact_orientation``
    The observed ordered directed-run vectors, translated to a uniform start.
``exact_global_rotation``
    The same trajectory after one Haar-uniform rotation of the whole path.
``ordered_length_rotated``
    The same ordered run-length sequence after an independent Haar-uniform
    rotation of every vector on every replay.
``vector_uniform_rotated``
    Directed-run vectors drawn uniformly with replacement from the condition
    pool, with an independent Haar-uniform rotation applied to every draw.
``fitted_mu``
    Isotropic relocations drawn from the jointly fitted crossover mixture.

Unbounded replay compares the two renewable constructions
(``vector_uniform_rotated`` and ``fitted_mu``) by simulating until first
detection.  Exact recorded trajectories are deliberately not recycled.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import math
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactStore
from .directions import isotropic_directions, random_quaternions
from .dist.mixture import ppf as mixture_ppf
from .logging_utils import get_logger
from .replay.engine import prepare_budgeted_paths, replay_prepared
from .targets import Target, dimension_from_projected_area

logger = get_logger(__name__)

DEFAULT_NK_CONDITION = (
    "cell_type=Natural Killer Cells | organ=popliteal lymph node | "
    "stimulus=Influenza Vaccine | dt_seconds=30.0"
)

FINITE_CLASSES: Tuple[str, ...] = (
    "exact_orientation",
    "exact_global_rotation",
    "vector_uniform_rotated",
    "fitted_mu",
    # Append new constructions so established classes retain their previous
    # random-number streams for a fixed seed and batch size.
    "ordered_length_rotated",
)
FINITE_PLOT_ORDER: Tuple[str, ...] = (
    "exact_orientation",
    "exact_global_rotation",
    "ordered_length_rotated",
    "vector_uniform_rotated",
    "fitted_mu",
)
UNBOUNDED_CLASSES: Tuple[str, ...] = (
    "vector_uniform_rotated",
    "fitted_mu",
)
FINITE_CONTRASTS: Tuple[Tuple[str, str], ...] = (
    ("exact_global_rotation", "exact_orientation"),
    ("ordered_length_rotated", "exact_global_rotation"),
    ("vector_uniform_rotated", "ordered_length_rotated"),
    # Retain the original aggregate contrast for backward comparability.
    ("vector_uniform_rotated", "exact_global_rotation"),
    ("fitted_mu", "vector_uniform_rotated"),
)
FINITE_CONTRAST_ROLES = {
    ("exact_global_rotation", "exact_orientation"): "laboratory_frame_orientation",
    ("ordered_length_rotated", "exact_global_rotation"): "directional_structure",
    (
        "vector_uniform_rotated",
        "ordered_length_rotated",
    ): "condition_pool_and_length_sequence_composite",
    (
        "vector_uniform_rotated",
        "exact_global_rotation",
    ): "total_randomization_secondary",
    ("fitted_mu", "vector_uniform_rotated"): "fitted_marginal_law",
}


@dataclass(frozen=True)
class FocusedReplaySettings:
    """Runtime settings for the focused replay."""

    condition: str = DEFAULT_NK_CONDITION
    length_unit_um: Optional[float] = None
    projected_areas: Tuple[float, ...] = (4.0, 8.0, 16.0, 32.0)
    finite_replays_per_track: int = 10_000
    unbounded_trials_per_cell: int = 400
    workers: int = 0
    finite_batch_size: int = 2_000
    unbounded_chunk_trials: int = 25
    unbounded_step_block: int = 512
    random_seed: int = 24680
    output_stage: str = "focused_replay"
    figure_format: str = "png"
    replay_signature: Optional[str] = None
    replay_signature_version: Optional[int] = None

    def validate(self) -> None:
        if self.length_unit_um is not None and self.length_unit_um <= 0.0:
            raise ValueError("length_unit_um must be positive")
        if not self.projected_areas:
            raise ValueError("At least one projected area is required")
        if any(value <= np.pi for value in self.projected_areas):
            raise ValueError("Projected areas must exceed pi for detection radius 1")
        if self.finite_replays_per_track < 1:
            raise ValueError("finite_replays_per_track must be positive")
        if self.unbounded_trials_per_cell < 1:
            raise ValueError("unbounded_trials_per_cell must be positive")
        if self.finite_batch_size < 1 or self.unbounded_chunk_trials < 1:
            raise ValueError("Batch and chunk sizes must be positive")
        if self.unbounded_step_block < 1:
            raise ValueError("unbounded_step_block must be positive")
        if not self.output_stage.strip():
            raise ValueError("output_stage must be non-empty")
        if self.figure_format not in ("png", "pdf"):
            raise ValueError("figure_format must be 'png' or 'pdf'")
        if (self.replay_signature is None) != (
            self.replay_signature_version is None
        ):
            raise ValueError(
                "replay_signature and replay_signature_version must be "
                "provided together"
            )
        if self.replay_signature is not None:
            if (
                len(self.replay_signature) != 64
                or any(
                    character not in "0123456789abcdef"
                    for character in self.replay_signature
                )
            ):
                raise ValueError(
                    "replay_signature must be a lowercase SHA-256 hex digest"
                )
            if int(self.replay_signature_version) < 1:
                raise ValueError("replay_signature_version must be positive")

    @property
    def resolved_workers(self) -> int:
        available = os.cpu_count() or 1
        return min(available, self.workers) if self.workers > 0 else available


def _rotate_vector_batch(
    vectors: np.ndarray, quaternions: np.ndarray
) -> np.ndarray:
    """Rotate ``(batch, steps, 3)`` vectors by batch or per-vector quaternions."""
    values = np.asarray(vectors, dtype=float)
    rotations = np.asarray(quaternions, dtype=float)
    if values.ndim != 3 or values.shape[2] != 3:
        raise ValueError("vectors must have shape (batch, steps, 3)")
    if rotations.shape == (len(values), 4):
        q_xyz = rotations[:, None, :3]
        q_w = rotations[:, None, 3:4]
    elif rotations.shape == values.shape[:2] + (4,):
        q_xyz = rotations[..., :3]
        q_w = rotations[..., 3:4]
    else:
        raise ValueError(
            "quaternions must have shape (batch, 4) or (batch, steps, 4)"
        )
    cross = np.cross(q_xyz, values)
    t = 2.0 * cross
    return values + q_w * t + np.cross(q_xyz, t)


def _candidate_count(budget: float, mean_length: float) -> int:
    """Return a conservative draw count for a finite path-length budget."""
    if budget <= 0.0 or mean_length <= 0.0:
        raise ValueError("budget and mean_length must be positive")
    return max(64, int(np.ceil(6.0 * budget / mean_length)) + 32)


def _draw_renewable_batch(
    class_name: str,
    n_trials: int,
    n_steps: int,
    rng: np.random.Generator,
    pool_vectors: np.ndarray,
    fitted_mu: float,
    fitted_crossover: float,
    fitted_lmax: float,
) -> np.ndarray:
    """Draw a rectangular batch for either renewable trajectory class."""
    if class_name == "vector_uniform_rotated":
        indices = rng.integers(0, len(pool_vectors), size=(n_trials, n_steps))
        sampled = pool_vectors[indices]
        quaternions = random_quaternions(n_trials * n_steps, rng).reshape(
            n_trials, n_steps, 4
        )
        return _rotate_vector_batch(sampled, quaternions)
    if class_name == "fitted_mu":
        uniforms = rng.random((n_trials, n_steps))
        lengths = mixture_ppf(
            uniforms.ravel(),
            fitted_mu,
            fitted_lmax,
            fitted_crossover,
        ).reshape(n_trials, n_steps)
        directions = isotropic_directions(n_trials * n_steps, rng).reshape(
            n_trials, n_steps, 3
        )
        return lengths[..., None] * directions
    raise ValueError("Unknown renewable class {!r}".format(class_name))


def _finite_sequences(
    class_name: str,
    source_vectors: np.ndarray,
    pool_vectors: np.ndarray,
    budget: float,
    n_trials: int,
    rng: np.random.Generator,
    fitted_mu: float,
    fitted_crossover: float,
    fitted_lmax: float,
) -> np.ndarray:
    """Construct one finite-budget rectangular trajectory batch."""
    if class_name == "exact_orientation":
        return np.broadcast_to(
            source_vectors, (n_trials,) + source_vectors.shape
        ).copy()
    if class_name == "exact_global_rotation":
        repeated = np.broadcast_to(
            source_vectors, (n_trials,) + source_vectors.shape
        )
        return _rotate_vector_batch(
            repeated, random_quaternions(n_trials, rng)
        )
    if class_name == "ordered_length_rotated":
        repeated = np.broadcast_to(
            source_vectors, (n_trials,) + source_vectors.shape
        )
        quaternions = random_quaternions(
            n_trials * len(source_vectors), rng
        ).reshape(n_trials, len(source_vectors), 4)
        return _rotate_vector_batch(repeated, quaternions)

    mean_length = (
        float(np.mean(np.linalg.norm(pool_vectors, axis=1)))
        if class_name == "vector_uniform_rotated"
        else _mixture_mean(fitted_mu, fitted_lmax, fitted_crossover)
    )
    n_steps = _candidate_count(budget, mean_length)
    while True:
        vectors = _draw_renewable_batch(
            class_name,
            n_trials,
            n_steps,
            rng,
            pool_vectors,
            fitted_mu,
            fitted_crossover,
            fitted_lmax,
        )
        coverage = np.sum(np.linalg.norm(vectors, axis=2), axis=1)
        if np.all(coverage + 1e-12 >= budget):
            return vectors
        n_steps *= 2


def _mixture_mean(mu: float, lmax: float, crossover: float) -> float:
    """Analytic mean of the flat-plus-power-law mixture."""
    if lmax <= crossover:
        return lmax / 2.0
    q0 = 1.0 - mu
    q1 = 2.0 - mu
    ratio = lmax / crossover
    tail_mass = (
        np.log(ratio)
        if abs(q0) < 1e-12
        else np.expm1(q0 * np.log(ratio)) / q0
    )
    tail_first = (
        np.log(ratio)
        if abs(q1) < 1e-12
        else np.expm1(q1 * np.log(ratio)) / q1
    )
    normalizer = crossover * (1.0 + tail_mass)
    first_moment = crossover ** 2 * (0.5 + tail_first)
    return float(first_moment / normalizer)


def _finite_track_worker(task: Dict[str, object]) -> List[Dict[str, object]]:
    """Replay one source track; used as a process-pool worker."""
    source_vectors = np.asarray(task["source_vectors"], dtype=float)
    pool_vectors = np.asarray(task["pool_vectors"], dtype=float)
    targets: Sequence[Tuple[float, Target]] = task["targets"]  # type: ignore[assignment]
    n_replays = int(task["n_replays"])
    batch_size = int(task["batch_size"])
    rng = np.random.default_rng(int(task["seed"]))
    budget = float(np.cumsum(np.linalg.norm(source_vectors, axis=1))[-1])
    safe_budget = budget + 1e-12 * max(1.0, budget)
    fitted_mu = float(task["fitted_mu"])
    fitted_crossover = float(task["fitted_crossover"])
    fitted_lmax = float(task["fitted_lmax"])
    finite_classes = tuple(
        str(value) for value in task.get("finite_classes", FINITE_CLASSES)
    )
    unknown_classes = set(finite_classes).difference(FINITE_CLASSES)
    if unknown_classes:
        raise ValueError(
            "Unknown finite classes: {}".format(
                ", ".join(sorted(unknown_classes))
            )
        )

    aggregates: Dict[Tuple[str, str, float], Dict[str, float]] = {}
    for class_name in finite_classes:
        for area, target in targets:
            aggregates[(class_name, target.shape, float(area))] = {
                "n_trials": 0.0,
                "n_detected": 0.0,
                "detection_distance_sum": 0.0,
                "restricted_distance_sum": 0.0,
                "hit_step_sum": 0.0,
            }

    for start_index in range(0, n_replays, batch_size):
        count = min(batch_size, n_replays - start_index)
        starts = rng.uniform(
            0.0,
            float(task["side_length"]),
            size=(count, 3),
        )
        budgets = np.full(count, safe_budget, dtype=float)
        for class_name in finite_classes:
            # The additive ordered-length control uses an independent,
            # deterministic stream. It therefore cannot advance the legacy
            # stream used by established classes and later batch starts.
            class_rng = (
                np.random.default_rng(
                    np.random.SeedSequence(
                        [
                            int(task["seed"]),
                            int(start_index),
                            0x4F4C52,
                        ]
                    )
                )
                if class_name == "ordered_length_rotated"
                else rng
            )
            vectors = _finite_sequences(
                class_name,
                source_vectors,
                pool_vectors,
                safe_budget,
                count,
                class_rng,
                fitted_mu,
                fitted_crossover,
                fitted_lmax,
            )
            paths = prepare_budgeted_paths(
                vectors,
                starts,
                budgets,
                float(task["side_length"]),
                overshoot_policy="discard",
            )
            for area, target in targets:
                result = replay_prepared(paths, target, check_initial_position=False)
                detected = result.detected
                key = (class_name, target.shape, float(area))
                cell = aggregates[key]
                cell["n_trials"] += count
                cell["n_detected"] += int(np.sum(detected))
                cell["detection_distance_sum"] += float(
                    np.nansum(result.detection_distance)
                )
                cell["restricted_distance_sum"] += float(
                    np.sum(result.restricted_detection_distance)
                )
                if np.any(detected):
                    cell["hit_step_sum"] += float(
                        np.sum(result.first_hit_step[detected] + 1)
                    )

    rows: List[Dict[str, object]] = []
    for (class_name, shape, area), cell in aggregates.items():
        n_trials = int(cell["n_trials"])
        n_detected = int(cell["n_detected"])
        rows.append(
            {
                "condition": str(task["condition"]),
                "video_id": str(task["video_id"]),
                "track_uid": str(task["track_uid"]),
                "trajectory_class": class_name,
                "shape": shape,
                "projected_area": area,
                "budget_model_units": budget,
                "n_runs": int(len(source_vectors)),
                "n_trials": n_trials,
                "n_detected": n_detected,
                "detection_probability": n_detected / n_trials,
                "conditional_mean_detection_distance": (
                    cell["detection_distance_sum"] / n_detected
                    if n_detected
                    else np.nan
                ),
                "restricted_mean_detection_distance": (
                    cell["restricted_distance_sum"] / n_trials
                ),
                "conditional_mean_hit_step": (
                    cell["hit_step_sum"] / n_detected
                    if n_detected
                    else np.nan
                ),
            }
        )
    return rows


def _unbounded_worker(task: Dict[str, object]) -> List[Dict[str, object]]:
    """Simulate an independent unbounded trial chunk until first detection."""
    class_name = str(task["trajectory_class"])
    pool_vectors = np.asarray(task["pool_vectors"], dtype=float)
    target: Target = task["target"]  # type: ignore[assignment]
    n_trials = int(task["n_trials"])
    block = int(task["step_block"])
    rng = np.random.default_rng(int(task["seed"]))
    positions = rng.uniform(0.0, target.side_length, size=(n_trials, 3))
    active = np.ones(n_trials, dtype=bool)
    distances = np.zeros(n_trials, dtype=float)
    steps = np.zeros(n_trials, dtype=np.int64)
    detection_distances = np.full(n_trials, np.nan, dtype=float)
    detection_steps = np.full(n_trials, -1, dtype=np.int64)
    fitted_mu = float(task["fitted_mu"])
    fitted_crossover = float(task["fitted_crossover"])
    fitted_lmax = float(task["fitted_lmax"])

    while np.any(active):
        active_indices = np.flatnonzero(active)
        count = len(active_indices)
        vectors = _draw_renewable_batch(
            class_name,
            count,
            block,
            rng,
            pool_vectors,
            fitted_mu,
            fitted_crossover,
            fitted_lmax,
        )
        lengths = np.linalg.norm(vectors, axis=2)
        relative = np.cumsum(vectors, axis=1)
        endpoints = np.mod(
            positions[active_indices, None, :] + relative,
            target.side_length,
        )
        hits = target.contains(endpoints)
        detected_now = np.any(hits, axis=1)

        if np.any(detected_now):
            local_hit_rows = np.flatnonzero(detected_now)
            first = np.argmax(hits[local_hit_rows], axis=1)
            global_hit_rows = active_indices[local_hit_rows]
            cumulative_lengths = np.cumsum(lengths[local_hit_rows], axis=1)
            detection_distances[global_hit_rows] = (
                distances[global_hit_rows]
                + cumulative_lengths[np.arange(len(first)), first]
            )
            detection_steps[global_hit_rows] = (
                steps[global_hit_rows] + first + 1
            )
            active[global_hit_rows] = False

        surviving_local = np.flatnonzero(~detected_now)
        if len(surviving_local):
            global_survivors = active_indices[surviving_local]
            positions[global_survivors] = endpoints[surviving_local, -1, :]
            distances[global_survivors] += np.sum(
                lengths[surviving_local], axis=1
            )
            steps[global_survivors] += block

    return [
        {
            "condition": str(task["condition"]),
            "trajectory_class": class_name,
            "shape": target.shape,
            "projected_area": float(task["projected_area"]),
            "trial_index": int(task["trial_offset"]) + index,
            "detection_distance_model_units": float(detection_distances[index]),
            "first_hit_step": int(detection_steps[index]),
        }
        for index in range(n_trials)
    ]


def _track_bootstrap_summary(
    track_summaries: pd.DataFrame,
    n_bootstraps: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Combine equal-replay track estimates and bootstrap whole tracks."""
    fields = ["trajectory_class", "shape", "projected_area"]
    rows: List[Dict[str, object]] = []
    for keys, group in track_summaries.groupby(fields, sort=True):
        probabilities = group["detection_probability"].to_numpy(dtype=float)
        restricted = group["restricted_mean_detection_distance"].to_numpy(
            dtype=float
        )
        sampled = rng.integers(
            0, len(group), size=(n_bootstraps, len(group))
        )
        p_draws = np.mean(probabilities[sampled], axis=1)
        r_draws = np.mean(restricted[sampled], axis=1)
        p_low, p_high = np.quantile(p_draws, [0.025, 0.975])
        r_low, r_high = np.quantile(r_draws, [0.025, 0.975])
        n_trials = int(group["n_trials"].sum())
        n_detected = int(group["n_detected"].sum())
        p_hat = float(np.mean(probabilities))
        rows.append(
            {
                "trajectory_class": keys[0],
                "shape": keys[1],
                "projected_area": keys[2],
                "n_tracks": int(group["track_uid"].nunique()),
                "n_trials": n_trials,
                "n_detected": n_detected,
                "detection_probability": p_hat,
                "detection_probability_mc_se": math.sqrt(
                    max(0.0, p_hat * (1.0 - p_hat) / n_trials)
                ),
                "detection_probability_ci_low": float(p_low),
                "detection_probability_ci_high": float(p_high),
                "restricted_mean_detection_distance": float(
                    np.mean(restricted)
                ),
                "rmdd_ci_low": float(r_low),
                "rmdd_ci_high": float(r_high),
            }
        )
    return pd.DataFrame(rows)


def _unbounded_summary(
    trials: pd.DataFrame, length_unit_um: float
) -> pd.DataFrame:
    """Summarise independent unbounded Monte-Carlo trials."""
    fields = ["trajectory_class", "shape", "projected_area"]
    rows: List[Dict[str, object]] = []
    for keys, group in trials.groupby(fields, sort=True):
        distances = group["detection_distance_model_units"].to_numpy(dtype=float)
        steps = group["first_hit_step"].to_numpy(dtype=float)
        n = len(group)
        distance_se = (
            float(np.std(distances, ddof=1) / np.sqrt(n))
            if n > 1
            else np.nan
        )
        step_se = (
            float(np.std(steps, ddof=1) / np.sqrt(n))
            if n > 1
            else np.nan
        )
        rows.append(
            {
                "trajectory_class": keys[0],
                "shape": keys[1],
                "projected_area": keys[2],
                "n_trials": n,
                "mean_detection_distance_model_units": float(np.mean(distances)),
                "detection_distance_mc_se_model_units": distance_se,
                "median_detection_distance_model_units": float(np.median(distances)),
                "q25_detection_distance_model_units": float(
                    np.quantile(distances, 0.25)
                ),
                "q75_detection_distance_model_units": float(
                    np.quantile(distances, 0.75)
                ),
                "mean_detection_distance_um": float(
                    np.mean(distances) * length_unit_um
                ),
                "detection_distance_mc_se_um": distance_se * length_unit_um,
                "mean_first_hit_step": float(np.mean(steps)),
                "first_hit_step_mc_se": step_se,
            }
        )
    return pd.DataFrame(rows)


def _finite_contrasts(
    track_summaries: pd.DataFrame,
    n_bootstraps: int,
    rng: np.random.Generator,
    contrasts: Sequence[Tuple[str, str]] = FINITE_CONTRASTS,
) -> pd.DataFrame:
    """Compute paired, whole-track finite detection-probability contrasts."""
    rows: List[Dict[str, object]] = []
    cell_fields = ["shape", "projected_area"]
    for keys, cell in track_summaries.groupby(cell_fields, sort=True):
        pivot = cell.pivot(
            index="track_uid",
            columns="trajectory_class",
            values="detection_probability",
        )
        for left, right in contrasts:
            if left not in pivot or right not in pivot:
                continue
            paired = (pivot[left] - pivot[right]).dropna().to_numpy(dtype=float)
            sampled = rng.integers(
                0, len(paired), size=(n_bootstraps, len(paired))
            )
            draws = np.mean(paired[sampled], axis=1)
            low, high = np.quantile(draws, [0.025, 0.975])
            rows.append(
                {
                    "left_class": left,
                    "right_class": right,
                    "contrast": "{} minus {}".format(left, right),
                    "contrast_role": FINITE_CONTRAST_ROLES.get(
                        (left, right), "custom"
                    ),
                    "shape": keys[0],
                    "projected_area": keys[1],
                    "n_tracks": len(paired),
                    "detection_probability_difference": float(
                        np.mean(paired)
                    ),
                    "difference_ci_low": float(low),
                    "difference_ci_high": float(high),
                    "n_bootstraps": n_bootstraps,
                }
            )
    return pd.DataFrame(rows)


def _unbounded_contrasts(summary: pd.DataFrame) -> pd.DataFrame:
    """Compare fitted-mu and empirical-vector mean first-detection distance."""
    rows: List[Dict[str, object]] = []
    for keys, cell in summary.groupby(["shape", "projected_area"], sort=True):
        indexed = cell.set_index("trajectory_class")
        if not set(UNBOUNDED_CLASSES).issubset(indexed.index):
            continue
        fitted = indexed.loc["fitted_mu"]
        empirical = indexed.loc["vector_uniform_rotated"]
        fitted_mean = float(fitted["mean_detection_distance_model_units"])
        empirical_mean = float(empirical["mean_detection_distance_model_units"])
        fitted_se = float(fitted["detection_distance_mc_se_model_units"])
        empirical_se = float(
            empirical["detection_distance_mc_se_model_units"]
        )
        difference = fitted_mean - empirical_mean
        difference_se = math.sqrt(fitted_se ** 2 + empirical_se ** 2)
        ratio = fitted_mean / empirical_mean
        log_ratio_se = math.sqrt(
            (fitted_se / fitted_mean) ** 2
            + (empirical_se / empirical_mean) ** 2
        )
        rows.append(
            {
                "shape": keys[0],
                "projected_area": keys[1],
                "fitted_minus_empirical_distance": difference,
                "difference_mc_se": difference_se,
                "difference_ci_low": difference - 1.96 * difference_se,
                "difference_ci_high": difference + 1.96 * difference_se,
                "fitted_to_empirical_distance_ratio": ratio,
                "ratio_ci_low": float(
                    np.exp(np.log(ratio) - 1.96 * log_ratio_se)
                ),
                "ratio_ci_high": float(
                    np.exp(np.log(ratio) + 1.96 * log_ratio_se)
                ),
                "relative_distance_change_percent": 100.0 * (ratio - 1.0),
                "interval_method": "independent Monte Carlo normal/delta approximation",
            }
        )
    return pd.DataFrame(rows)


def _make_targets(
    config: Config,
    projected_areas: Sequence[float],
    side_length: float,
) -> List[Tuple[float, Target]]:
    """Construct target cells valid for every configured shape."""
    targets: List[Tuple[float, Target]] = []
    detection_radius = config.replay.detection_radius
    for area in projected_areas:
        dimensions = {
            shape: dimension_from_projected_area(
                area, shape, detection_radius
            )
            for shape in config.replay.shapes
        }
        invalid = [
            shape
            for shape, dimension in dimensions.items()
            if dimension + 2.0 * detection_radius > side_length
        ]
        if invalid:
            raise ValueError(
                "Projected area {} has a detection neighbourhood that does not "
                "fit without periodic self-overlap in a side-{} torus for: {}".format(
                    area,
                    side_length,
                    ", ".join(invalid),
                )
            )
        for shape in config.replay.shapes:
            targets.append(
                (
                    float(area),
                    Target(
                        shape=shape,
                        dimension=dimensions[shape],
                        side_length=side_length,
                        detection_radius=detection_radius,
                    ),
                )
            )
    return targets


def _plot_results(
    output_directory: Path,
    finite: pd.DataFrame,
    unbounded: pd.DataFrame,
    figure_format: str = "png",
) -> List[Path]:
    """Write the two reviewer-facing comparison figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = output_directory / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    shapes = sorted(finite["shape"].unique())
    class_labels = {
        "exact_orientation": "Exact orientation",
        "exact_global_rotation": "Exact + global rotation",
        "ordered_length_rotated": "Exact lengths/order + independent rotations",
        "vector_uniform_rotated": "Uniform empirical vector + rotation",
        "fitted_mu": r"Fitted $\mu$",
    }
    colors = {
        "exact_orientation": "#1f77b4",
        "exact_global_rotation": "#17becf",
        "ordered_length_rotated": "#9467bd",
        "vector_uniform_rotated": "#ff7f0e",
        "fitted_mu": "#2ca02c",
    }

    fig, axes = plt.subplots(1, len(shapes), figsize=(4.3 * len(shapes), 3.8))
    axes = np.atleast_1d(axes)
    for axis, shape in zip(axes, shapes):
        selected = finite[finite["shape"] == shape]
        for class_name in FINITE_PLOT_ORDER:
            group = selected[selected["trajectory_class"] == class_name].sort_values(
                "projected_area"
            )
            if group.empty:
                continue
            y = group["detection_probability"].to_numpy(dtype=float)
            lower = group["detection_probability_ci_low"].to_numpy(dtype=float)
            upper = group["detection_probability_ci_high"].to_numpy(dtype=float)
            axis.errorbar(
                group["projected_area"],
                y,
                yerr=np.vstack((y - lower, upper - y)),
                marker="o",
                capsize=2,
                color=colors[class_name],
                label=class_labels[class_name],
            )
        axis.set_xscale("log", base=2)
        axis.set_title(shape)
        axis.set_xlabel("Projected target area (model units²)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Finite-budget detection probability")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.84))
    finite_path = figures / "finite_detection_probability.{}".format(
        figure_format
    )
    fig.savefig(finite_path, dpi=220 if figure_format == "png" else None)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(shapes), figsize=(4.3 * len(shapes), 3.8))
    axes = np.atleast_1d(axes)
    for axis, shape in zip(axes, shapes):
        selected = unbounded[unbounded["shape"] == shape]
        for class_name in UNBOUNDED_CLASSES:
            group = selected[selected["trajectory_class"] == class_name].sort_values(
                "projected_area"
            )
            if group.empty:
                continue
            y = group["mean_detection_distance_model_units"].to_numpy(dtype=float)
            se = group[
                "detection_distance_mc_se_model_units"
            ].to_numpy(dtype=float)
            axis.errorbar(
                group["projected_area"],
                y,
                yerr=1.96 * se,
                marker="o",
                capsize=2,
                color=colors[class_name],
                label=class_labels[class_name],
            )
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_title(shape)
        axis.set_xlabel("Projected target area (model units²)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Mean detection path length (model units)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    unbounded_path = figures / "unbounded_detection_distance.{}".format(
        figure_format
    )
    fig.savefig(unbounded_path, dpi=220 if figure_format == "png" else None)
    plt.close(fig)
    return [finite_path, unbounded_path]


def run_focused_replay(
    config: Config, settings: FocusedReplaySettings
) -> Dict[str, object]:
    """Run and persist the focused finite and unbounded experiments."""
    settings.validate()
    store = ArtifactStore(config.output.root)
    runs = store.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config.hash(),
    )
    pools = store.read_frame(
        "preprocess",
        "pools",
        expected_config_hash=config.hash(),
    )
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config.hash(),
    )
    fit = fits[
        (fits["condition"] == settings.condition) & (fits["status"] == "ok")
    ]
    if len(fit) != 1:
        raise ValueError(
            "Expected one successful fit for {!r}; found {}".format(
                settings.condition, len(fit)
            )
        )

    condition_runs = runs[
        (runs["condition"] == settings.condition)
        & np.isfinite(runs["run_length_um"].to_numpy(dtype=float))
        & (runs["run_length_um"].to_numpy(dtype=float) > 0.0)
    ].copy()
    condition_pool = pools[
        (pools["condition"] == settings.condition)
        & np.isfinite(pools["run_length_um"].to_numpy(dtype=float))
        & (pools["run_length_um"].to_numpy(dtype=float) > 0.0)
    ].copy()
    if condition_runs.empty or condition_pool.empty:
        raise ValueError("Selected condition has no replayable directed runs")
    # Keep the exact-path comparison on the same fitted track population. Two
    # NK tracks have only boundary-censored runs and therefore did not
    # contribute to the 27-track fit or its 517-vector empirical pool.
    fitted_track_ids = set(condition_pool["track_uid"].astype(str))
    condition_runs = condition_runs[
        condition_runs["track_uid"].astype(str).isin(fitted_track_ids)
    ].copy()

    fit_row = fit.iloc[0]
    length_unit_um = (
        float(fit_row["crossover"])
        if settings.length_unit_um is None
        else float(settings.length_unit_um)
    )
    vector_columns = ["run_vector_x", "run_vector_y", "run_vector_z"]
    pool_vectors = (
        condition_pool[vector_columns].to_numpy(dtype=float)
        / length_unit_um
    )
    pool_input_columns = [
        column
        for column in (
            "condition",
            "video_id",
            "track_uid",
            "fragment_uid",
            "run_uid",
            "run_length_um",
            "run_duration_s",
        )
        if column in condition_pool.columns
    ]
    pool_inputs = condition_pool[pool_input_columns].copy()
    for index, column in enumerate(vector_columns):
        pool_inputs["{}_model_units".format(column)] = pool_vectors[:, index]
    pool_inputs["run_length_model_units"] = np.linalg.norm(
        pool_vectors, axis=1
    )
    pool_inputs["uniform_run_sampling_probability"] = 1.0 / len(pool_inputs)
    fitted_mu = float(fit_row["mu_hat"])
    fitted_crossover = float(fit_row["crossover"]) / length_unit_um
    fitted_lmax = float(fit_row["lmax"]) / length_unit_um
    side_length = 2.0 * fitted_lmax
    targets = _make_targets(config, settings.projected_areas, side_length)

    sequence_rows = []
    sequences: List[Tuple[str, str, np.ndarray]] = []
    for track_uid, group in condition_runs.groupby("track_uid", sort=True):
        ordered = group.sort_values(
            ["frame_from", "frame_to", "fragment_uid", "run_id"],
            kind="mergesort",
        )
        vectors = (
            ordered[vector_columns].to_numpy(dtype=float)
            / length_unit_um
        )
        if len(vectors) == 0:
            continue
        video_id = str(ordered["video_id"].iloc[0])
        sequences.append((str(track_uid), video_id, vectors))
        sequence_rows.append(
            {
                "condition": settings.condition,
                "video_id": video_id,
                "track_uid": str(track_uid),
                "n_runs": len(vectors),
                "budget_model_units": float(
                    np.cumsum(np.linalg.norm(vectors, axis=1))[-1]
                ),
            }
        )

    workers = settings.resolved_workers
    seed_sequence = np.random.SeedSequence(settings.random_seed)
    finite_seeds = seed_sequence.spawn(len(sequences))
    finite_tasks = []
    for (track_uid, video_id, vectors), child_seed in zip(
        sequences, finite_seeds
    ):
        finite_tasks.append(
            {
                "condition": settings.condition,
                "video_id": video_id,
                "track_uid": track_uid,
                "source_vectors": vectors,
                "pool_vectors": pool_vectors,
                "targets": targets,
                "n_replays": settings.finite_replays_per_track,
                "batch_size": settings.finite_batch_size,
                "side_length": side_length,
                "fitted_mu": fitted_mu,
                "fitted_crossover": fitted_crossover,
                "fitted_lmax": fitted_lmax,
                "seed": int(child_seed.generate_state(1, dtype=np.uint64)[0]),
            }
        )

    logger.info(
        "Focused finite replay: %d tracks x %d replicates on %d workers",
        len(finite_tasks),
        settings.finite_replays_per_track,
        workers,
    )
    finite_rows: List[Dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_finite_track_worker, task) for task in finite_tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            finite_rows.extend(future.result())
            if completed % max(1, len(futures) // 5) == 0:
                logger.info("Finite replay completed %d/%d tracks", completed, len(futures))
    finite_tracks = (
        pd.DataFrame(finite_rows)
        .sort_values(
            [
                "condition",
                "video_id",
                "track_uid",
                "trajectory_class",
                "shape",
                "projected_area",
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )
    finite_summary = _track_bootstrap_summary(
        finite_tracks,
        config.fitting.clustered_bootstraps,
        np.random.default_rng(settings.random_seed + 1),
    )
    finite_contrasts = _finite_contrasts(
        finite_tracks,
        config.fitting.clustered_bootstraps,
        np.random.default_rng(settings.random_seed + 2),
    )

    unbounded_tasks: List[Dict[str, object]] = []
    task_specs = []
    for class_name in UNBOUNDED_CLASSES:
        for area, target in targets:
            for offset in range(
                0,
                settings.unbounded_trials_per_cell,
                settings.unbounded_chunk_trials,
            ):
                task_specs.append(
                    (
                        class_name,
                        area,
                        target,
                        offset,
                        min(
                            settings.unbounded_chunk_trials,
                            settings.unbounded_trials_per_cell - offset,
                        ),
                    )
                )
    unbounded_seeds = seed_sequence.spawn(len(task_specs))
    for spec, child_seed in zip(task_specs, unbounded_seeds):
        class_name, area, target, offset, count = spec
        unbounded_tasks.append(
            {
                "condition": settings.condition,
                "trajectory_class": class_name,
                "projected_area": area,
                "target": target,
                "trial_offset": offset,
                "n_trials": count,
                "step_block": settings.unbounded_step_block,
                "pool_vectors": pool_vectors,
                "fitted_mu": fitted_mu,
                "fitted_crossover": fitted_crossover,
                "fitted_lmax": fitted_lmax,
                "seed": int(child_seed.generate_state(1, dtype=np.uint64)[0]),
            }
        )

    logger.info(
        "Focused unbounded replay: %d class/target trials in %d worker tasks",
        len(UNBOUNDED_CLASSES)
        * len(targets)
        * settings.unbounded_trials_per_cell,
        len(unbounded_tasks),
    )
    unbounded_rows: List[Dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(_unbounded_worker, task) for task in unbounded_tasks
        ]
        for completed, future in enumerate(as_completed(futures), start=1):
            unbounded_rows.extend(future.result())
            if completed % max(1, len(futures) // 10) == 0:
                logger.info(
                    "Unbounded replay completed %d/%d chunks",
                    completed,
                    len(futures),
                )
    unbounded_trials = (
        pd.DataFrame(unbounded_rows)
        .sort_values(
            [
                "trajectory_class",
                "shape",
                "projected_area",
                "trial_index",
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )
    unbounded_summary = _unbounded_summary(unbounded_trials, length_unit_um)
    unbounded_contrasts = _unbounded_contrasts(unbounded_summary)

    stage = settings.output_stage
    metadata = {
        "config_hash": config.hash(),
        "condition": settings.condition,
        "length_unit_um": length_unit_um,
    }
    if settings.replay_signature is not None:
        metadata.update(
            {
                "replay_signature": settings.replay_signature,
                "replay_signature_version": settings.replay_signature_version,
            }
        )
    store.write_frame(
        stage,
        "track_inputs",
        pd.DataFrame(sequence_rows),
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "empirical_pool_inputs",
        pool_inputs,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "finite_track_summaries",
        finite_tracks,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "finite_summary",
        finite_summary,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "finite_contrasts",
        finite_contrasts,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "unbounded_trials",
        unbounded_trials,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "unbounded_summary",
        unbounded_summary,
        extra_metadata=metadata,
    )
    store.write_frame(
        stage,
        "unbounded_contrasts",
        unbounded_contrasts,
        extra_metadata=metadata,
    )
    output_directory = store.stage_dir(stage)
    figure_paths = _plot_results(
        output_directory,
        finite_summary,
        unbounded_summary,
        figure_format=settings.figure_format,
    )
    target_rows = [
        {
            "shape": target.shape,
            "projected_area": area,
            "dimension": target.dimension,
            "side_length": target.side_length,
            "detection_radius": target.detection_radius,
            "maximum_centered_extent": (
                target.dimension / 2.0 + target.detection_radius
            ),
            "clearance_to_torus_half_width": (
                target.side_length / 2.0
                - target.dimension / 2.0
                - target.detection_radius
            ),
            "fits_without_periodic_self_overlap": bool(
                target.dimension + 2.0 * target.detection_radius
                <= target.side_length
            ),
        }
        for area, target in targets
    ]
    store.write_frame(
        stage,
        "targets",
        pd.DataFrame(target_rows),
        extra_metadata=metadata,
    )
    settings_dict = {
        "condition": settings.condition,
        "normalization": (
            "fitted_crossover"
            if settings.length_unit_um is None
            else "user_length_unit"
        ),
        "length_unit_um": length_unit_um,
        "detection_radius_model_units": config.replay.detection_radius,
        "detection_radius_um": config.replay.detection_radius * length_unit_um,
        "projected_areas": list(settings.projected_areas),
        "finite_replays_per_track": settings.finite_replays_per_track,
        "unbounded_trials_per_cell": settings.unbounded_trials_per_cell,
        "workers": workers,
        "finite_batch_size": settings.finite_batch_size,
        "unbounded_chunk_trials": settings.unbounded_chunk_trials,
        "unbounded_step_block": settings.unbounded_step_block,
        "random_seed": settings.random_seed,
        "output_stage": settings.output_stage,
        "figure_format": settings.figure_format,
        "replay_signature": settings.replay_signature,
        "replay_signature_version": settings.replay_signature_version,
        "finite_classes": list(FINITE_CLASSES),
        "unbounded_classes": list(UNBOUNDED_CLASSES),
        "fitted_mu": fitted_mu,
        "fitted_crossover_model_units": fitted_crossover,
        "fitted_lmax_model_units": fitted_lmax,
        "uniform_pool_mean_length_model_units": float(
            np.mean(np.linalg.norm(pool_vectors, axis=1))
        ),
        "uniform_pool_mean_length_um": float(
            np.mean(np.linalg.norm(pool_vectors, axis=1))
            * length_unit_um
        ),
        "fitted_mean_length_model_units": _mixture_mean(
            fitted_mu, fitted_lmax, fitted_crossover
        ),
        "fitted_mean_length_um": _mixture_mean(
            fitted_mu, fitted_lmax, fitted_crossover
        )
        * length_unit_um,
        "torus_side_model_units": side_length,
        "torus_side_um": side_length * length_unit_um,
        "pool_vectors": int(len(pool_vectors)),
        "source_tracks": int(len(sequences)),
        "engine": "Python process pool with C-compatible target geometry",
        "c_engine_note": (
            "The normalization matches the C sampler's crossover=1 and d=1 "
            "convention. The existing C sampler was not used because it truncates "
            "lmax and side to integers and cannot consume empirical vector pools "
            "or finite recorded paths."
        ),
    }
    store.write_json(stage, "settings", settings_dict)
    summary: Dict[str, object] = {
        "condition": settings.condition,
        "tracks": len(sequences),
        "pool_vectors": len(pool_vectors),
        "finite_trial_count": int(finite_tracks["n_trials"].sum()),
        "unbounded_trial_count": int(len(unbounded_trials)),
        "figures": [str(path) for path in figure_paths],
        "replay_signature": settings.replay_signature,
        "replay_signature_version": settings.replay_signature_version,
    }
    store.write_json(stage, "summary", summary)
    return summary


__all__ = [
    "DEFAULT_NK_CONDITION",
    "FINITE_CLASSES",
    "FINITE_CONTRASTS",
    "FINITE_CONTRAST_ROLES",
    "FINITE_PLOT_ORDER",
    "FocusedReplaySettings",
    "UNBOUNDED_CLASSES",
    "run_focused_replay",
]
