"""Fixed-arrest-scale target-shape replay for the three selected conditions.

This module implements the publication sensitivity requested after the primary
joint-crossover analysis.  It fixes one operational arrest displacement per
condition, refits the flat-plus-power-law exponent, and compares four finite
trajectory constructions against targets matched either by projected
detection area or by detection-neighbourhood volume.

The legacy C fixed-volume program is deliberately not used: it cannot consume
recorded paths and its historical "volume" mapping is not the volume tested by
the target-membership code.  The geometry here is the tested and C-compatible
``Target`` geometry used by the empirical replay pipeline.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .bootstrap import clustered_mu_bootstrap
from .config import Config, load_config
from .dataio import ArtifactStore
from .dist.gof import parametric_bootstrap_ks
from .dist.mle import fit_mixture
from .logging_utils import get_logger
from .nk_experiment import _finite_sequences
from .replay.engine import prepare_budgeted_paths, replay_prepared
from .targets import (
    SUPPORTED_SHAPES,
    Target,
    dimension_from_neighbourhood_volume,
    dimension_from_projected_area,
    projected_area,
)

logger = get_logger(__name__)

ALGORITHM_VERSION = "fixed_arrest_scale_area_volume_v1"
FINITE_CLASSES: Tuple[str, ...] = (
    "exact_orientation",
    "exact_global_rotation",
    "ordered_length_rotated",
    "fitted_mu",
)
AREA_VALUES: Tuple[float, ...] = (
    4.0,
    8.0,
    16.0,
    32.0,
    64.0,
    125.0,
    250.0,
    375.0,
    500.0,
)
VOLUME_VALUES: Tuple[float, ...] = (
    7.0,
    10.0,
    20.0,
    40.0,
    80.0,
    160.0,
    320.0,
    500.0,
    750.0,
    925.0,
)


@dataclass(frozen=True)
class SelectedCondition:
    condition: str
    condition_id: str
    label: str
    dt_seconds: float
    c_arrest_um: float
    expected_mu: float


SELECTED_CONDITIONS: Tuple[SelectedCondition, ...] = (
    SelectedCondition(
        condition=(
            "cell_type=Natural Killer Cells | organ=popliteal lymph node | "
            "stimulus=Influenza Vaccine | dt_seconds=30.0"
        ),
        condition_id="nk_pln_influenza_dt30",
        label="NK, influenza, 30 s",
        dt_seconds=30.0,
        c_arrest_um=1.0,
        expected_mu=1.6730531986544044,
    ),
    SelectedCondition(
        condition=(
            "cell_type=Neutrophils | organ=popliteal lymph node | "
            "stimulus=Influenza Vaccine | dt_seconds=15.0"
        ),
        condition_id="neutrophils_pln_influenza_dt15",
        label="Neutrophils, influenza, 15 s",
        dt_seconds=15.0,
        c_arrest_um=0.5,
        expected_mu=1.2741683698252115,
    ),
    SelectedCondition(
        condition=(
            "cell_type=T Cells | organ=popliteal lymph node | "
            "stimulus=HIV-infected humanized T cell | dt_seconds=15.0"
        ),
        condition_id="t_pln_hiv_dt15",
        label="T cells, HIV model, 15 s",
        dt_seconds=15.0,
        c_arrest_um=0.5,
        expected_mu=1.0724623114743956,
    ),
)


@dataclass(frozen=True)
class ArrestScaleExperimentSettings:
    source_config: Path = Path(
        "configs/equal_run_weight_run_eligible_conditions.yaml"
    )
    output_root: Path = Path(
        "outputs/arrest_scale_selected_area_volume"
    )
    area_values: Tuple[float, ...] = AREA_VALUES
    volume_values: Tuple[float, ...] = VOLUME_VALUES
    finite_replays_per_track: int = 10_000
    finite_batch_size: int = 1_000
    workers: int = 64
    clustered_bootstraps: int = 1_000
    goodness_of_fit_bootstraps: int = 500
    random_seed: int = 731_905
    resume: bool = True

    @property
    def resolved_workers(self) -> int:
        return min(max(1, self.workers), os.cpu_count() or 1)

    def validate(self) -> None:
        if not self.area_values or not self.volume_values:
            raise ValueError("Both target grids must be non-empty")
        if any(value <= math.pi for value in self.area_values):
            raise ValueError("Projected detection areas must exceed pi at d=1")
        if any(value < 2.0 * math.pi for value in self.volume_values):
            raise ValueError(
                "Common detection volumes must be at least 2*pi at d=1"
            )
        if self.finite_replays_per_track < 1 or self.finite_batch_size < 1:
            raise ValueError("Replay and batch counts must be positive")
        if self.clustered_bootstraps < 1:
            raise ValueError("clustered_bootstraps must be positive")
        if self.goodness_of_fit_bootstraps < 1:
            raise ValueError("goodness_of_fit_bootstraps must be positive")


def line_projected_area_limit(side_length: float, d: float = 1.0) -> float:
    """Largest common projected area allowed by the line capsule."""
    return float(2.0 * d * side_length + (math.pi - 4.0) * d * d)


def line_neighbourhood_volume_limit(
    side_length: float, d: float = 1.0
) -> float:
    """Largest common detection volume allowed by the line capsule."""
    return float(
        math.pi * d * d * (side_length - 2.0 * d)
        + 4.0 * math.pi * d ** 3 / 3.0
    )


def _stable_seed(master_seed: int, *parts: object) -> int:
    payload = "|".join([str(master_seed)] + [str(part) for part in parts])
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little")


def _signature(payload: Mapping[str, object]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _fixed_arrest_fits(
    config: Config,
    pools: pd.DataFrame,
    settings: ArrestScaleExperimentSettings,
) -> pd.DataFrame:
    """Refit one exponent per selected condition with fixed c_arrest."""
    rows: List[Dict[str, object]] = []
    for index, spec in enumerate(SELECTED_CONDITIONS):
        group = pools[
            (pools["condition"] == spec.condition)
            & np.isfinite(pools["run_length_um"].to_numpy(dtype=float))
            & (pools["run_length_um"].to_numpy(dtype=float) > 0.0)
        ].copy()
        if group.empty:
            raise ValueError("No eligible runs for {}".format(spec.label))
        lengths = group["run_length_um"].to_numpy(dtype=float)
        weights = group["step_weighted_weight"].to_numpy(dtype=float)
        lmax = float(np.max(lengths))
        fitted = fit_mixture(
            lengths,
            lmax=lmax,
            crossover=spec.c_arrest_um,
            mu_bounds=config.fitting.mu_bounds,
            weights=weights,
            estimate_crossover=False,
        )
        if not np.isclose(fitted.mu, spec.expected_mu, rtol=0.0, atol=5e-7):
            raise ValueError(
                "Recomputed mu for {} changed from the approved value: "
                "{:.9g} versus {:.9g}".format(
                    spec.label, fitted.mu, spec.expected_mu
                )
            )
        rng = np.random.default_rng(
            _stable_seed(settings.random_seed, spec.condition, "fit")
        )
        gof = parametric_bootstrap_ks(
            lengths,
            fitted,
            settings.goodness_of_fit_bootstraps,
            rng,
            mu_bounds=config.fitting.mu_bounds,
            weights=weights,
            estimate_crossover=False,
        )
        clustered = clustered_mu_bootstrap(
            group,
            settings.clustered_bootstraps,
            rng,
            crossover=spec.c_arrest_um,
            mu_bounds=config.fitting.mu_bounds,
            weighting="step_weighted",
            estimate_crossover=False,
            lmax_method="empirical_max",
            user_lmax=None,
            lmax_quantile=0.99,
        )
        lmax_normalized = lmax / spec.c_arrest_um
        side = 2.0 * lmax_normalized
        rows.append(
            {
                "condition": spec.condition,
                "condition_id": spec.condition_id,
                "condition_label": spec.label,
                "dt_seconds": spec.dt_seconds,
                "c_arrest_um": spec.c_arrest_um,
                "c_arrest_definition": (
                    "2 um/min conventional arrest speed multiplied by dt"
                ),
                "mu_hat": fitted.mu,
                "mu_ci_low": clustered.ci_low,
                "mu_ci_high": clustered.ci_high,
                "lmax_um": lmax,
                "lmax_model_units": lmax_normalized,
                "torus_side_model_units": side,
                "torus_side_um": side * spec.c_arrest_um,
                "line_area_limit_model_units2": line_projected_area_limit(side),
                "line_volume_limit_model_units3": (
                    line_neighbourhood_volume_limit(side)
                ),
                "n_runs": len(group),
                "n_tracks": int(group["track_uid"].nunique()),
                "ks_statistic": gof.statistic,
                "ks_bootstrap_p": gof.p_value,
                "gof_bootstraps": gof.n_bootstraps,
                "clustered_bootstraps_successful": clustered.n_successful,
                "clustered_bootstraps_failed": clustered.n_failed,
                "fit_at_mu_bound": fitted.at_bound,
                "interpretation": "fixed_arrest_scale_sensitivity",
            }
        )
    return pd.DataFrame(rows)


def _target_manifest_for_fit(
    fit: Mapping[str, object],
    settings: ArrestScaleExperimentSettings,
) -> pd.DataFrame:
    """Build and validate the common-three-shape target cells."""
    side = float(fit["torus_side_model_units"])
    c_um = float(fit["c_arrest_um"])
    d = 1.0
    limits = {
        "projected_area": line_projected_area_limit(side, d),
        "neighbourhood_volume": line_neighbourhood_volume_limit(side, d),
    }
    candidate_values = {
        "projected_area": settings.area_values,
        "neighbourhood_volume": settings.volume_values,
    }
    rows: List[Dict[str, object]] = []
    for matching_rule, candidates in candidate_values.items():
        retained = [float(value) for value in candidates if value <= limits[matching_rule]]
        if not retained:
            raise ValueError(
                "No {} target fits {}".format(
                    matching_rule, fit["condition_label"]
                )
            )
        for matched_value in retained:
            for shape in SUPPORTED_SHAPES:
                dimension = (
                    dimension_from_projected_area(matched_value, shape, d)
                    if matching_rule == "projected_area"
                    else dimension_from_neighbourhood_volume(
                        matched_value, shape, d
                    )
                )
                target = Target(
                    shape=shape,
                    dimension=dimension,
                    side_length=side,
                    detection_radius=d,
                )
                extent = dimension + 2.0 * d
                if extent > side + 1e-10:
                    raise ValueError(
                        "Target does not fit: {} {} {}".format(
                            fit["condition_label"], matching_rule, matched_value
                        )
                    )
                actual_area = projected_area(dimension, shape, d)
                actual_volume = target.neighbourhood_volume()
                expected = (
                    actual_area
                    if matching_rule == "projected_area"
                    else actual_volume
                )
                if not np.isclose(expected, matched_value, rtol=1e-11, atol=1e-11):
                    raise ValueError("Target geometry failed its round-trip check")
                rows.append(
                    {
                        "condition": fit["condition"],
                        "condition_id": fit["condition_id"],
                        "condition_label": fit["condition_label"],
                        "matching_rule": matching_rule,
                        "matched_value_model_units": matched_value,
                        "matched_value_physical": matched_value
                        * (c_um ** (2 if matching_rule == "projected_area" else 3)),
                        "matched_value_physical_units": (
                            "um^2"
                            if matching_rule == "projected_area"
                            else "um^3"
                        ),
                        "shape": shape,
                        "dimension_model_units": dimension,
                        "dimension_um": dimension * c_um,
                        "projected_area_model_units2": actual_area,
                        "projected_area_um2": actual_area * c_um * c_um,
                        "neighbourhood_volume_model_units3": actual_volume,
                        "neighbourhood_volume_um3": actual_volume * c_um ** 3,
                        "detection_radius_model_units": d,
                        "detection_radius_um": d * c_um,
                        "torus_side_model_units": side,
                        "torus_side_um": side * c_um,
                        "full_target_extent_model_units": extent,
                        "clearance_model_units": side - extent,
                        "fits_without_periodic_self_overlap": True,
                        "common_shape_limit_model_units": limits[matching_rule],
                    }
                )
    return pd.DataFrame(rows)


def _target_objects(manifest: pd.DataFrame) -> List[Tuple[str, float, Target]]:
    result: List[Tuple[str, float, Target]] = []
    for row in manifest.itertuples(index=False):
        result.append(
            (
                str(row.matching_rule),
                float(row.matched_value_model_units),
                Target(
                    shape=str(row.shape),
                    dimension=float(row.dimension_model_units),
                    side_length=float(row.torus_side_model_units),
                    detection_radius=float(row.detection_radius_model_units),
                ),
            )
        )
    return result


def _finite_track_worker(task: Dict[str, object]) -> List[Dict[str, object]]:
    """Replay one recorded track against both target-matching families."""
    source_vectors = np.asarray(task["source_vectors"], dtype=float)
    pool_vectors = np.asarray(task["pool_vectors"], dtype=float)
    targets: Sequence[Tuple[str, float, Target]] = task["targets"]  # type: ignore[assignment]
    n_replays = int(task["n_replays"])
    batch_size = int(task["batch_size"])
    seed = int(task["seed"])
    rng = np.random.default_rng(seed)
    budget = float(np.sum(np.linalg.norm(source_vectors, axis=1)))
    safe_budget = budget + 1e-12 * max(1.0, budget)
    side = float(task["side_length"])
    fitted_mu = float(task["fitted_mu"])
    fitted_lmax = float(task["fitted_lmax"])

    aggregates: Dict[Tuple[str, str, str, float], Dict[str, float]] = {}
    for class_name in FINITE_CLASSES:
        for matching_rule, matched_value, target in targets:
            aggregates[
                (class_name, matching_rule, target.shape, matched_value)
            ] = {
                "n_trials": 0.0,
                "n_detected": 0.0,
                "detection_distance_sum": 0.0,
                "restricted_distance_sum": 0.0,
                "hit_step_sum": 0.0,
            }

    for start_index in range(0, n_replays, batch_size):
        count = min(batch_size, n_replays - start_index)
        starts = rng.uniform(0.0, side, size=(count, 3))
        budgets = np.full(count, safe_budget, dtype=float)
        for class_name in FINITE_CLASSES:
            class_rng = (
                np.random.default_rng(
                    np.random.SeedSequence(
                        [seed, start_index, 0x41525245]
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
                1.0,
                fitted_lmax,
            )
            paths = prepare_budgeted_paths(
                vectors,
                starts,
                budgets,
                side,
                overshoot_policy="discard",
            )
            for matching_rule, matched_value, target in targets:
                replay = replay_prepared(
                    paths, target, check_initial_position=False
                )
                detected = replay.detected
                cell = aggregates[
                    (class_name, matching_rule, target.shape, matched_value)
                ]
                cell["n_trials"] += count
                cell["n_detected"] += int(np.sum(detected))
                cell["detection_distance_sum"] += float(
                    np.nansum(replay.detection_distance)
                )
                cell["restricted_distance_sum"] += float(
                    np.sum(replay.restricted_detection_distance)
                )
                if np.any(detected):
                    cell["hit_step_sum"] += float(
                        np.sum(replay.first_hit_step[detected] + 1)
                    )

    rows: List[Dict[str, object]] = []
    for key, cell in aggregates.items():
        class_name, matching_rule, shape, matched_value = key
        n_trials = int(cell["n_trials"])
        n_detected = int(cell["n_detected"])
        rows.append(
            {
                "condition": str(task["condition"]),
                "condition_id": str(task["condition_id"]),
                "condition_label": str(task["condition_label"]),
                "video_id": str(task["video_id"]),
                "track_uid": str(task["track_uid"]),
                "trajectory_class": class_name,
                "matching_rule": matching_rule,
                "shape": shape,
                "matched_value_model_units": matched_value,
                "budget_model_units": budget,
                "n_runs": len(source_vectors),
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


def _summarize_finite(
    track_rows: pd.DataFrame,
    n_bootstraps: int,
    seed: int,
) -> pd.DataFrame:
    fields = [
        "condition",
        "condition_id",
        "condition_label",
        "trajectory_class",
        "matching_rule",
        "shape",
        "matched_value_model_units",
    ]
    rows: List[Dict[str, object]] = []
    for keys, group in track_rows.groupby(fields, sort=True):
        probabilities = group["detection_probability"].to_numpy(dtype=float)
        restricted = group["restricted_mean_detection_distance"].to_numpy(
            dtype=float
        )
        rng = np.random.default_rng(_stable_seed(seed, *keys))
        indices = rng.integers(
            0, len(group), size=(n_bootstraps, len(group))
        )
        p_draws = np.mean(probabilities[indices], axis=1)
        r_draws = np.mean(restricted[indices], axis=1)
        p_hat = float(np.mean(probabilities))
        n_trials = int(group["n_trials"].sum())
        row = dict(zip(fields, keys))
        row.update(
            {
                "n_tracks": int(group["track_uid"].nunique()),
                "n_trials": n_trials,
                "n_detected": int(group["n_detected"].sum()),
                "detection_probability": p_hat,
                "detection_probability_mc_se": math.sqrt(
                    max(0.0, p_hat * (1.0 - p_hat) / n_trials)
                ),
                "detection_probability_ci_low": float(
                    np.quantile(p_draws, 0.025)
                ),
                "detection_probability_ci_high": float(
                    np.quantile(p_draws, 0.975)
                ),
                "restricted_mean_detection_distance": float(
                    np.mean(restricted)
                ),
                "rmdd_ci_low": float(np.quantile(r_draws, 0.025)),
                "rmdd_ci_high": float(np.quantile(r_draws, 0.975)),
                "track_weighting": "equal weight per source track",
                "interval_method": "whole-track percentile bootstrap, 95%",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _finite_contrasts(
    track_rows: pd.DataFrame,
    n_bootstraps: int,
    seed: int,
) -> pd.DataFrame:
    contrasts = (
        ("exact_global_rotation", "exact_orientation", "laboratory orientation"),
        (
            "ordered_length_rotated",
            "exact_global_rotation",
            "within-track directional structure",
        ),
        (
            "fitted_mu",
            "ordered_length_rotated",
            "fitted run-length law",
        ),
        ("fitted_mu", "exact_orientation", "total fitted-model difference"),
    )
    cell_fields = [
        "condition",
        "condition_id",
        "condition_label",
        "matching_rule",
        "shape",
        "matched_value_model_units",
    ]
    rows: List[Dict[str, object]] = []
    for keys, group in track_rows.groupby(cell_fields, sort=True):
        pivot = group.pivot(
            index="track_uid",
            columns="trajectory_class",
            values="detection_probability",
        )
        for left, right, role in contrasts:
            paired = (pivot[left] - pivot[right]).dropna().to_numpy(dtype=float)
            rng = np.random.default_rng(
                _stable_seed(seed, *keys, left, right)
            )
            indices = rng.integers(
                0, len(paired), size=(n_bootstraps, len(paired))
            )
            draws = np.mean(paired[indices], axis=1)
            row = dict(zip(cell_fields, keys))
            row.update(
                {
                    "left_class": left,
                    "right_class": right,
                    "contrast": "{} minus {}".format(left, right),
                    "contrast_role": role,
                    "n_tracks": len(paired),
                    "detection_probability_difference": float(
                        np.mean(paired)
                    ),
                    "difference_ci_low": float(np.quantile(draws, 0.025)),
                    "difference_ci_high": float(np.quantile(draws, 0.975)),
                    "n_bootstraps": n_bootstraps,
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def _make_figures(
    output_root: Path,
    summary: pd.DataFrame,
    fits: pd.DataFrame,
) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    import matplotlib.pyplot as plt

    figures = output_root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    from matplotlib.lines import Line2D

    # This is the compact paper-facing comparison.  Shape is encoded by
    # colour (and a redundant marker), while trajectory construction is
    # encoded only by line style.  The observed and whole-track-rotation
    # controls remain in the saved tables but are deliberately omitted here.
    class_style = {
        "ordered_length_rotated": ("Empirical", "-", True),
        "fitted_mu": ("Fitted $\\mu$", "--", False),
    }
    shape_style = {
        "Ball": ("#0072B2", "o"),
        "Disk": ("#D55E00", "s"),
        "Line": ("#009E73", "^"),
    }
    shape_handles = [
        Line2D(
            [0],
            [0],
            color=color,
            marker=marker,
            linestyle="none",
            markersize=4.0,
            label=shape,
        )
        for shape, (color, marker) in shape_style.items()
    ]
    class_handles = [
        Line2D(
            [0],
            [0],
            color="#222222",
            linestyle=linestyle,
            linewidth=1.35,
            marker="o",
            markerfacecolor="#222222" if filled else "white",
            markeredgecolor="#222222",
            markeredgewidth=0.6,
            markersize=3.0,
            label=label,
        )
        for label, linestyle, filled in class_style.values()
    ]
    panel_titles = {
        "nk_pln_influenza_dt30": "NK",
        "neutrophils_pln_influenza_dt15": "Neutrophils",
        "t_pln_hiv_dt15": "T cells",
    }
    compact_limits = {
        "nk_pln_influenza_dt30": {
            "xlim": (0.0, 130.0),
            "xticks": (0.0, 60.0, 120.0),
            "ylim": (0.0, 0.012),
            "yticks": (0.0, 0.005, 0.010),
            "y_exponent": -2,
        },
        "neutrophils_pln_influenza_dt15": {
            "xlim": (0.0, 515.0),
            "xticks": (0.0, 250.0, 500.0),
            "ylim": (0.0, 0.0005),
            "yticks": (0.0, 0.0002, 0.0004),
            "y_exponent": -4,
        },
        "t_pln_hiv_dt15": {
            "xlim": (0.0, 515.0),
            "xticks": (0.0, 250.0, 500.0),
            "ylim": (0.0, 0.0035),
            "yticks": (0.0, 0.0015, 0.0030),
            "y_exponent": -3,
        },
    }
    expected_ids = {spec.condition_id for spec in SELECTED_CONDITIONS}
    if not expected_ids.issubset(set(fits["condition_id"].astype(str))):
        raise ValueError("The figure fit table lacks a selected condition")

    def draw(
        shared_y: bool,
        stem: str,
        one_column: bool = False,
        compact_row: bool = False,
    ) -> Path:
        if one_column and compact_row:
            raise ValueError("A figure cannot use both compact layouts")
        if compact_row:
            # Rerender the three panels, rather than scaling the full-width
            # figure, so text remains at least 6 pt at the 8.7-cm PNAS width.
            fig, axes = plt.subplots(
                1,
                len(SELECTED_CONDITIONS),
                figsize=(8.7 / 2.54, 2.05),
                sharex=False,
                sharey=False,
            )
        elif one_column:
            # PNAS final one-column width is 8.7 cm.  Re-rendering as a
            # vertical stack keeps every text element above the journal's
            # 6-point minimum; shrinking the horizontal layout would not.
            fig, axes = plt.subplots(
                len(SELECTED_CONDITIONS),
                1,
                figsize=(8.7 / 2.54, 6.2),
                sharex=False,
                sharey=False,
            )
        else:
            fig, axes = plt.subplots(
                1,
                len(SELECTED_CONDITIONS),
                figsize=(7.2, 2.45),
                sharex=False,
                sharey=False,
            )
        selected = summary[
            (summary["matching_rule"] == "projected_area")
            & summary["trajectory_class"].isin(class_style)
            & summary["condition_id"].isin(expected_ids)
        ]
        common_y_high = float(
            1.08
            * selected["detection_probability_ci_high"].to_numpy(dtype=float).max()
        )
        plot_specs = (
            tuple(reversed(SELECTED_CONDITIONS))
            if compact_row
            else SELECTED_CONDITIONS
        )
        for column_index, spec in enumerate(plot_specs):
            axis = axes[column_index]
            condition = selected[selected["condition_id"] == spec.condition_id]
            for shape, (color, marker) in shape_style.items():
                cell = condition[condition["shape"] == shape]
                for class_name, (_, linestyle, filled) in class_style.items():
                    group = cell[
                        cell["trajectory_class"] == class_name
                    ].sort_values("matched_value_model_units")
                    axis.plot(
                        group["matched_value_model_units"],
                        group["detection_probability"],
                        color=color,
                        marker=marker,
                        markerfacecolor=color if filled else "white",
                        markeredgecolor=color,
                        markeredgewidth=0.7,
                        linestyle=linestyle,
                        linewidth=1.1 if compact_row else 1.25,
                        markersize=2.9 if compact_row else 3.2,
                        solid_capstyle="round",
                        dash_capstyle="round",
                    )
                    axis.fill_between(
                        group["matched_value_model_units"].to_numpy(dtype=float),
                        group["detection_probability_ci_low"].to_numpy(dtype=float),
                        group["detection_probability_ci_high"].to_numpy(dtype=float),
                        color=color,
                        alpha=0.075 if compact_row else 0.065,
                        linewidth=0.0,
                    )
            condition_y_high = float(
                1.08
                * condition["detection_probability_ci_high"]
                .to_numpy(dtype=float)
                .max()
            )
            if compact_row:
                from matplotlib.ticker import FuncFormatter

                limits = compact_limits[spec.condition_id]
                exponent = int(limits["y_exponent"])
                y_scale = 10.0 ** exponent
                axis.set_xlim(*limits["xlim"])
                axis.set_xticks(limits["xticks"])
                axis.set_ylim(*limits["ylim"])
                axis.set_yticks(limits["yticks"])
                axis.yaxis.set_major_formatter(
                    FuncFormatter(
                        lambda value, _position, scale=y_scale: (
                            "{:g}".format(value / scale)
                        )
                    )
                )
                axis.set_xlabel(r"$\Delta_P$", fontsize=6.7, labelpad=1.0)
                axis.set_ylabel("")
                axis.grid(
                    axis="y",
                    color="#AAB4C3",
                    alpha=0.32,
                    linewidth=0.35,
                )
                axis.set_axisbelow(True)
                axis.set_title(
                    panel_titles[spec.condition_id],
                    fontsize=7.2,
                    fontweight="semibold",
                    pad=2.0,
                )
                axis.text(
                    0.035,
                    0.965,
                    r"$\times 10^{{{}}}$".format(exponent),
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=5.9,
                    color="#4B5563",
                )
                axis.tick_params(
                    axis="both",
                    direction="out",
                    labelsize=6.1,
                    length=2.0,
                    width=0.5,
                    pad=1.4,
                )
                for spine_name in ("top", "right"):
                    axis.spines[spine_name].set_visible(False)
                for spine_name in ("left", "bottom"):
                    axis.spines[spine_name].set_color("#30343B")
                    axis.spines[spine_name].set_linewidth(0.55)
            else:
                axis.set_xlim(left=0.0)
                axis.set_ylim(
                    0.0, common_y_high if shared_y else condition_y_high
                )
                axis.set_xlabel(
                    r"$\Delta_P$",
                    fontsize=8.0,
                    labelpad=2.0 if one_column else None,
                )
                axis.set_ylabel(
                    r"$\Pr(\mathrm{Detection})$",
                    fontsize=7.5,
                    labelpad=2.0 if one_column else None,
                )
                axis.grid(alpha=0.20, linewidth=0.45)
                axis.set_title(
                    panel_titles[spec.condition_id],
                    fontsize=8.5,
                    pad=2.5 if one_column else 4,
                )
                axis.tick_params(
                    labelsize=7.0 if one_column else 6.7,
                    pad=1.5 if one_column else 3.5,
                )
            if one_column:
                axis.locator_params(axis="x", nbins=5)
                axis.locator_params(axis="y", nbins=4)
            if not compact_row:
                axis.ticklabel_format(
                    axis="y", style="sci", scilimits=(0, 0)
                )
                axis.yaxis.get_offset_text().set_fontsize(
                    7.0 if one_column else 6.5
                )

        if compact_row:
            fig.text(
                0.050,
                0.565,
                r"$\Pr(\mathrm{Detection})$",
                ha="center",
                va="center",
                rotation="vertical",
                fontsize=6.7,
            )
            fig.legend(
                shape_handles + class_handles,
                [
                    handle.get_label()
                    for handle in shape_handles + class_handles
                ],
                loc="lower center",
                bbox_to_anchor=(0.5, 0.035),
                ncol=5,
                columnspacing=0.45,
                handlelength=1.25,
                handletextpad=0.22,
                markerscale=0.82,
                borderaxespad=0.0,
                frameon=False,
                fontsize=6.0,
            )
            fig.subplots_adjust(
                left=0.12,
                right=0.88,
                top=0.86,
                bottom=0.23,
                wspace=0.33,
            )
        elif one_column:
            fig.legend(
                shape_handles,
                [handle.get_label() for handle in shape_handles],
                loc="lower center",
                bbox_to_anchor=(0.5, 0.052),
                ncol=3,
                columnspacing=1.15,
                handletextpad=0.35,
                markerscale=1.0,
                frameon=False,
                fontsize=6.8,
            )
            fig.legend(
                class_handles,
                [handle.get_label() for handle in class_handles],
                loc="lower center",
                bbox_to_anchor=(0.5, 0.012),
                ncol=2,
                columnspacing=1.4,
                handlelength=2.2,
                handletextpad=0.45,
                frameon=False,
                fontsize=6.8,
            )
            fig.subplots_adjust(
                left=0.22,
                right=0.98,
                top=0.98,
                bottom=0.14,
                hspace=0.58,
            )
        else:
            handles = shape_handles + class_handles
            fig.legend(
                handles,
                [handle.get_label() for handle in handles],
                loc="lower center",
                bbox_to_anchor=(0.5, 0.005),
                ncol=5,
                columnspacing=1.3,
                handlelength=2.1,
                frameon=False,
                fontsize=7.0,
            )
            fig.subplots_adjust(
                left=0.075,
                right=0.995,
                top=0.86,
                bottom=0.29,
                wspace=0.42,
            )
        path = figures / "{}.pdf".format(stem)
        save_kwargs = (
            {} if one_column or compact_row else {"bbox_inches": "tight"}
        )
        fig.savefig(path, **save_kwargs)
        fig.savefig(figures / "{}.png".format(stem), dpi=300, **save_kwargs)
        plt.close(fig)
        return path

    return [
        draw(False, "finite_detection_probability_equal_area"),
        draw(True, "finite_detection_probability_equal_area_shared_y"),
        draw(
            False,
            "finite_detection_probability_equal_area_one_column",
            one_column=True,
        ),
        draw(
            True,
            "finite_detection_probability_equal_area_one_column_shared_y",
            one_column=True,
        ),
        draw(
            False,
            "finite_detection_probability_equal_area_compact_row",
            compact_row=True,
        ),
    ]


def make_arrest_scale_figures(
    output_root: Path = Path("outputs/arrest_scale_selected_area_volume"),
) -> List[Path]:
    """Regenerate the paper-facing figures from completed replay tables only."""
    root = Path(output_root)
    summary = pd.read_csv(root / "finite_replay" / "finite_summary.csv")
    fits = pd.read_csv(root / "fit" / "fixed_arrest_fits.csv")
    return _make_figures(root, summary, fits)


def _validate_results(
    track_rows: pd.DataFrame,
    summary: pd.DataFrame,
    manifests: pd.DataFrame,
    settings: ArrestScaleExperimentSettings,
) -> Dict[str, object]:
    key = [
        "condition_id",
        "track_uid",
        "trajectory_class",
        "matching_rule",
        "shape",
        "matched_value_model_units",
    ]
    if track_rows.duplicated(key).any():
        raise ValueError("Finite track output contains duplicate cells")
    if not (track_rows["n_trials"] == settings.finite_replays_per_track).all():
        raise ValueError("One or more finite cells has an incomplete trial count")
    if set(track_rows["trajectory_class"]) != set(FINITE_CLASSES):
        raise ValueError("Finite output is missing a requested trajectory class")
    if not manifests["fits_without_periodic_self_overlap"].astype(bool).all():
        raise ValueError("Target manifest contains a self-overlapping target")
    expected_summary = len(manifests) * len(FINITE_CLASSES)
    if len(summary) != expected_summary:
        raise ValueError(
            "Expected {} summary rows, found {}".format(
                expected_summary, len(summary)
            )
        )
    return {
        "status": "passed",
        "finite_track_rows": len(track_rows),
        "finite_summary_rows": len(summary),
        "target_cells": len(manifests),
        "trajectory_classes": list(FINITE_CLASSES),
        "all_trials_complete": True,
        "all_targets_fit": True,
        "duplicate_key_rows": 0,
    }


def _condition_sequences(
    runs: pd.DataFrame,
    pools: pd.DataFrame,
    fit: Mapping[str, object],
) -> Tuple[List[Tuple[str, str, np.ndarray]], np.ndarray, pd.DataFrame]:
    condition = str(fit["condition"])
    c_um = float(fit["c_arrest_um"])
    condition_pool = pools[
        (pools["condition"] == condition)
        & np.isfinite(pools["run_length_um"].to_numpy(dtype=float))
        & (pools["run_length_um"].to_numpy(dtype=float) > 0.0)
    ].copy()
    fitted_tracks = set(condition_pool["track_uid"].astype(str))
    condition_runs = runs[
        (runs["condition"] == condition)
        & runs["track_uid"].astype(str).isin(fitted_tracks)
        & np.isfinite(runs["run_length_um"].to_numpy(dtype=float))
        & (runs["run_length_um"].to_numpy(dtype=float) > 0.0)
    ].copy()
    vector_columns = ["run_vector_x", "run_vector_y", "run_vector_z"]
    pool_vectors = condition_pool[vector_columns].to_numpy(dtype=float) / c_um
    sequences: List[Tuple[str, str, np.ndarray]] = []
    rows: List[Dict[str, object]] = []
    for track_uid, group in condition_runs.groupby("track_uid", sort=True):
        ordered = group.sort_values(
            ["frame_from", "frame_to", "fragment_uid", "run_id"],
            kind="mergesort",
        )
        vectors = ordered[vector_columns].to_numpy(dtype=float) / c_um
        if not len(vectors):
            continue
        video_id = str(ordered["video_id"].iloc[0])
        sequences.append((str(track_uid), video_id, vectors))
        rows.append(
            {
                "condition": condition,
                "condition_id": fit["condition_id"],
                "condition_label": fit["condition_label"],
                "video_id": video_id,
                "track_uid": str(track_uid),
                "n_runs": len(vectors),
                "budget_model_units": float(
                    np.sum(np.linalg.norm(vectors, axis=1))
                ),
                "budget_um": float(
                    np.sum(np.linalg.norm(vectors, axis=1)) * c_um
                ),
            }
        )
    return sequences, pool_vectors, pd.DataFrame(rows)


def run_arrest_scale_experiment(
    settings: ArrestScaleExperimentSettings = ArrestScaleExperimentSettings(),
) -> Dict[str, object]:
    """Run the full fixed-arrest-scale finite replay and write PDF figures."""
    settings.validate()
    config = load_config(settings.source_config)
    source = ArtifactStore(config.output.root)
    output = ArtifactStore(settings.output_root)
    runs = source.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config.hash(),
    )
    pools = source.read_frame(
        "preprocess", "pools", expected_config_hash=config.hash()
    )

    fit_payload = {
        "algorithm_version": ALGORITHM_VERSION,
        "source_config_hash": config.hash(),
        "conditions": [asdict(spec) for spec in SELECTED_CONDITIONS],
        "area_values": list(settings.area_values),
        "volume_values": list(settings.volume_values),
        "finite_classes": list(FINITE_CLASSES),
        "finite_replays_per_track": settings.finite_replays_per_track,
        "finite_batch_size": settings.finite_batch_size,
        "clustered_bootstraps": settings.clustered_bootstraps,
        "goodness_of_fit_bootstraps": settings.goodness_of_fit_bootstraps,
        "random_seed": settings.random_seed,
        "detection_radius_model_units": 1.0,
        "torus_rule": "W'=2*Lmax/c_arrest",
        "detection_semantics": "endpoint-only; initial position unchecked",
    }
    signature = _signature(fit_payload)
    completion = output.stage_dir("validation") / "complete.json"
    if settings.resume and completion.is_file():
        prior = json.loads(completion.read_text(encoding="utf-8"))
        if prior.get("signature") == signature:
            logger.info("Reusing completed fixed-arrest-scale experiment")
            return prior

    fits = _fixed_arrest_fits(config, pools, settings)
    manifests = pd.concat(
        [
            _target_manifest_for_fit(row, settings)
            for row in fits.to_dict("records")
        ],
        ignore_index=True,
    )

    all_tasks: List[Dict[str, object]] = []
    input_frames: List[pd.DataFrame] = []
    plan_rows: List[Dict[str, object]] = []
    for fit in fits.to_dict("records"):
        condition_manifest = manifests[
            manifests["condition_id"] == fit["condition_id"]
        ]
        targets = _target_objects(condition_manifest)
        sequences, pool_vectors, inputs = _condition_sequences(
            runs, pools, fit
        )
        input_frames.append(inputs)
        for track_uid, video_id, vectors in sequences:
            all_tasks.append(
                {
                    "condition": fit["condition"],
                    "condition_id": fit["condition_id"],
                    "condition_label": fit["condition_label"],
                    "video_id": video_id,
                    "track_uid": track_uid,
                    "source_vectors": vectors,
                    "pool_vectors": pool_vectors,
                    "targets": targets,
                    "n_replays": settings.finite_replays_per_track,
                    "batch_size": settings.finite_batch_size,
                    "side_length": fit["torus_side_model_units"],
                    "fitted_mu": fit["mu_hat"],
                    "fitted_lmax": fit["lmax_model_units"],
                    "seed": _stable_seed(
                        settings.random_seed, fit["condition"], track_uid
                    ),
                }
            )
        plan_rows.append(
            {
                "condition": fit["condition"],
                "condition_id": fit["condition_id"],
                "condition_label": fit["condition_label"],
                "tracks": len(sequences),
                "target_cells": len(targets),
                "area_values_retained": ",".join(
                    "{:g}".format(value)
                    for value in sorted(
                        condition_manifest.loc[
                            condition_manifest["matching_rule"]
                            == "projected_area",
                            "matched_value_model_units",
                        ].unique()
                    )
                ),
                "volume_values_retained": ",".join(
                    "{:g}".format(value)
                    for value in sorted(
                        condition_manifest.loc[
                            condition_manifest["matching_rule"]
                            == "neighbourhood_volume",
                            "matched_value_model_units",
                        ].unique()
                    )
                ),
                "finite_target_evaluations": (
                    len(sequences)
                    * len(targets)
                    * len(FINITE_CLASSES)
                    * settings.finite_replays_per_track
                ),
            }
        )

    plan = pd.DataFrame(plan_rows)
    metadata = {
        "signature": signature,
        "algorithm_version": ALGORITHM_VERSION,
        "source_config_hash": config.hash(),
    }
    output.write_frame("00_plan", "condition_plan", plan, extra_metadata=metadata)
    output.write_frame("fit", "fixed_arrest_fits", fits, extra_metadata=metadata)
    output.write_frame("targets", "target_manifest", manifests, extra_metadata=metadata)
    output.write_frame(
        "inputs",
        "track_inputs",
        pd.concat(input_frames, ignore_index=True),
        extra_metadata=metadata,
    )
    output.write_json(
        "00_plan",
        "settings",
        dict(
            fit_payload,
            signature=signature,
            workers=settings.resolved_workers,
            output_root=str(settings.output_root),
        ),
    )

    logger.info(
        "Fixed-arrest finite replay: %d tracks, %d target evaluations, %d workers",
        len(all_tasks),
        int(plan["finite_target_evaluations"].sum()),
        settings.resolved_workers,
    )
    rows: List[Dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=settings.resolved_workers) as executor:
        futures = [
            executor.submit(_finite_track_worker, task) for task in all_tasks
        ]
        for completed_count, future in enumerate(as_completed(futures), start=1):
            rows.extend(future.result())
            if completed_count % max(1, len(futures) // 20) == 0:
                logger.info(
                    "Finite replay completed %d/%d tracks",
                    completed_count,
                    len(futures),
                )

    track_rows = pd.DataFrame(rows).sort_values(
        [
            "condition_id",
            "video_id",
            "track_uid",
            "matching_rule",
            "shape",
            "matched_value_model_units",
            "trajectory_class",
        ],
        kind="mergesort",
    ).reset_index(drop=True)
    summary = _summarize_finite(
        track_rows,
        settings.clustered_bootstraps,
        settings.random_seed + 1,
    )
    contrasts = _finite_contrasts(
        track_rows,
        settings.clustered_bootstraps,
        settings.random_seed + 2,
    )
    summary = summary.merge(
        manifests[
            [
                "condition_id",
                "matching_rule",
                "shape",
                "matched_value_model_units",
                "matched_value_physical",
                "matched_value_physical_units",
                "dimension_model_units",
                "dimension_um",
                "projected_area_model_units2",
                "projected_area_um2",
                "neighbourhood_volume_model_units3",
                "neighbourhood_volume_um3",
            ]
        ],
        on=[
            "condition_id",
            "matching_rule",
            "shape",
            "matched_value_model_units",
        ],
        how="left",
        validate="many_to_one",
    )
    validation = _validate_results(track_rows, summary, manifests, settings)
    output.write_frame(
        "finite_replay",
        "finite_track_summaries",
        track_rows,
        extra_metadata=metadata,
    )
    output.write_frame(
        "finite_replay", "finite_summary", summary, extra_metadata=metadata
    )
    output.write_frame(
        "finite_replay", "finite_contrasts", contrasts, extra_metadata=metadata
    )
    output.write_json("validation", "checks", validation)
    figures = _make_figures(settings.output_root, summary, fits)

    result = {
        "status": "complete",
        "signature": signature,
        "algorithm_version": ALGORITHM_VERSION,
        "conditions": len(fits),
        "tracks": len(all_tasks),
        "target_cells": len(manifests),
        "finite_track_rows": len(track_rows),
        "finite_summary_rows": len(summary),
        "finite_target_evaluations": int(
            plan["finite_target_evaluations"].sum()
        ),
        "finite_classes": list(FINITE_CLASSES),
        "figures": [str(path) for path in figures],
        "output_root": str(settings.output_root.resolve()),
    }
    output.write_json("validation", "complete", result)
    return result


__all__ = [
    "ALGORITHM_VERSION",
    "AREA_VALUES",
    "ArrestScaleExperimentSettings",
    "FINITE_CLASSES",
    "SELECTED_CONDITIONS",
    "VOLUME_VALUES",
    "line_neighbourhood_volume_limit",
    "line_projected_area_limit",
    "make_arrest_scale_figures",
    "run_arrest_scale_experiment",
]
