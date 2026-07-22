"""Command-line entry points for the LTDB Lévy analysis pipeline."""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
warnings.filterwarnings(
    "ignore",
    message="Pandas requires version .* of 'bottleneck'",
    category=UserWarning,
)
import pandas as pd

from .bootstrap import clustered_mu_bootstrap
from .classes import (
    TRAJECTORY_CLASSES,
    build_block_pool,
    build_trajectory_class,
    make_common_random_numbers,
    mixture_mean,
)
from .conditions import assign_conditions, load_metadata, write_metadata_template
from .config import Config, load_config
from .dataio import ArtifactStore, discover_files, inspect_file, read_dataset
from .dist.compare import grouped_model_comparison
from .dist.families import fit_family
from .dist.gof import ks_power_curve, parametric_bootstrap_ks
from .dist.mle import fit_mixture, select_lmax
from .errors import FitError, LtdbLevyError
from .logging_utils import configure_logging, get_logger
from .outcomes import combine_trial_summaries, shape_sensitivity, summarise_trials
from .pools import build_pools
from .preprocess import preprocess
from .replay import replay_targets_batch
from .segment import SegmentParams, segment_runs, segmentation_sweep
from .stats import clustered_outcome_intervals
from .targets import Target, dimension_from_projected_area, projected_area

logger = get_logger(__name__)


def _store_run_context(store: ArtifactStore, config: Config, stage: str) -> None:
    store.write_json(
        stage,
        "run_context",
        {
            "run_id": config.run_id,
            "config_hash": config.hash(),
            "config_path": str(config.source_path) if config.source_path else None,
            "config": config.to_dict(),
        },
    )


def run_inspect(config: Config) -> Dict[str, int]:
    """Inspect every selected source file and write a metadata template if needed."""
    store = ArtifactStore(config.output.root)
    candidates = sorted(config.input.tracks_directory.glob(config.input.glob))
    selected = discover_files(
        config.input.tracks_directory,
        config.input.glob,
        config.selection.include_files,
        config.selection.exclude_files,
    )
    selected_names = {path.name for path in selected}
    selection_rows: List[Dict[str, object]] = []
    include_names = set(config.selection.include_files)
    exclude_names = set(config.selection.exclude_files)
    for path in candidates:
        if include_names and path.name not in include_names:
            reason = "not_in_include_files"
        elif path.name in exclude_names:
            reason = "listed_in_exclude_files"
        else:
            reason = ""
        selection_rows.append(
            {"file": path.name, "selected": path.name in selected_names, "exclusion_reason": reason}
        )
    inspection = pd.DataFrame(
        [
            inspect_file(path, config.input.coordinates_are_physical)
            for path in selected
        ]
    )
    store.write_frame("inspect", "file_inspection", inspection)
    store.write_frame("inspect", "file_selection_audit", pd.DataFrame(selection_rows))

    ok = inspection[inspection["ok"].astype(bool)] if not inspection.empty else inspection
    if config.input.metadata_file is not None and not config.input.metadata_file.exists():
        template_columns = [
            "video_id",
            "video",
            "population",
            "cohort",
            "dt_s",
        ]
        write_metadata_template(ok[template_columns], config.input.metadata_file)

    summary = {
        "candidate_files": int(len(candidates)),
        "selected_files": int(len(selected)),
        "ok_files": int(inspection["ok"].sum()) if not inspection.empty else 0,
        "parse_failures": int((~inspection["ok"].astype(bool)).sum())
        if not inspection.empty
        else 0,
        "observations": int(ok["n_observations"].sum()) if not ok.empty else 0,
        "tracks": int(ok["n_tracks"].sum()) if not ok.empty else 0,
        "id_mismatches": int(ok["id_mismatch"].sum()) if not ok.empty else 0,
        "frame_gaps": int(ok["n_frame_gaps"].sum()) if not ok.empty else 0,
    }
    store.write_json("inspect", "summary", summary)
    _store_run_context(store, config, "inspect")
    logger.info("Inspection summary: %s", json.dumps(summary, sort_keys=True))
    return summary


def run_preprocess(config: Config) -> Dict[str, int]:
    """Parse, clean, segment, assign conditions, and build empirical pools."""
    store = ArtifactStore(config.output.root)
    observations, file_table = read_dataset(
        config.input.tracks_directory,
        config.input.glob,
        config.selection.include_files,
        config.selection.exclude_files,
        config.input.coordinates_are_physical,
    )
    cleaned, displacements, pre_audit = preprocess(
        observations,
        duplicate_policy=config.preprocessing.duplicate_policy,
        split_at_missing_frames=config.preprocessing.split_at_missing_frames,
        minimum_observations_per_track=config.selection.minimum_observations_per_track,
    )
    segment_params = SegmentParams(
        turning_angle_degrees=config.segmentation.turning_angle_degrees,
        speed_threshold=config.segmentation.speed_threshold,
        minimum_run_frames=config.segmentation.minimum_run_frames,
        minimum_run_length=config.segmentation.minimum_run_length,
        exclude_first_and_last_run=config.segmentation.exclude_first_and_last_run,
    )
    runs, segment_audit = segment_runs(displacements, segment_params)
    metadata = load_metadata(config.input.metadata_file)
    runs = assign_conditions(
        runs,
        metadata,
        config.conditions.grouping_fields,
        file_table=file_table,
        require_verified=config.conditions.require_verified_metadata,
        fallback_to_file=config.conditions.fallback_to_file,
    )
    pools, pool_audit = build_pools(
        runs,
        minimum_tracks=config.empirical_pool.minimum_tracks_per_pool,
        minimum_runs=config.empirical_pool.minimum_runs_per_pool,
        complete_only=config.empirical_pool.use_complete_runs_only,
    )
    angles = sorted(
        set(
            [config.segmentation.turning_angle_degrees]
            + list(config.segmentation.alternative_turning_angles)
        )
    )
    speeds = list(config.segmentation.alternative_speed_thresholds)
    if config.segmentation.speed_threshold not in speeds:
        speeds.append(config.segmentation.speed_threshold)
    sweep = segmentation_sweep(displacements, angles, speeds, segment_params)

    common_meta = {"config_hash": config.hash()}
    store.write_frame(
        "preprocess", "file_table", file_table, extra_metadata=common_meta
    )
    store.write_frame(
        "preprocess",
        "observations",
        cleaned,
        schema_name="observations",
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess",
        "displacements",
        displacements,
        schema_name="displacements",
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess", "runs", runs, schema_name="runs", extra_metadata=common_meta
    )
    store.write_frame("preprocess", "pools", pools, extra_metadata=common_meta)
    store.write_frame(
        "preprocess",
        "pool_audit",
        pool_audit,
        schema_name="pool_audit",
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess",
        "preprocess_audit",
        pre_audit.to_frame(),
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess",
        "preprocess_exclusions",
        pd.DataFrame(pre_audit.rows),
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess",
        "segment_audit",
        segment_audit.to_frame(),
        extra_metadata=common_meta,
    )
    store.write_frame(
        "preprocess",
        "segmentation_sweep",
        sweep,
        extra_metadata=common_meta,
    )
    summary = {
        "observations_read": int(len(observations)),
        "observations_retained": int(len(cleaned)),
        "tracks_retained": int(cleaned["track_uid"].nunique()),
        "fragments_retained": int(cleaned["fragment_uid"].nunique()),
        "displacements": int(len(displacements)),
        "runs_total": int(len(runs)),
        "runs_included": int(runs["included"].sum()) if len(runs) else 0,
        "conditions": int(runs["condition"].nunique()) if len(runs) else 0,
        "eligible_pools": int(pool_audit["eligible"].sum()) if len(pool_audit) else 0,
    }
    store.write_json("preprocess", "summary", summary)
    _store_run_context(store, config, "preprocess")
    logger.info("Preprocessing summary: %s", json.dumps(summary, sort_keys=True))
    return summary


def _weight_column(config: Config) -> str:
    return "{}_weight".format(config.empirical_pool.primary_weighting)


def _fit_segmentation_sensitivity(
    config: Config,
    displacements: pd.DataFrame,
    file_table: pd.DataFrame,
) -> pd.DataFrame:
    """Refit mixture parameters and bootstrap KS across the segmentation grid."""
    metadata = load_metadata(config.input.metadata_file)
    angles = sorted(
        set(
            [config.segmentation.turning_angle_degrees]
            + list(config.segmentation.alternative_turning_angles)
        )
    )
    speeds = list(config.segmentation.alternative_speed_thresholds)
    if config.segmentation.speed_threshold not in speeds:
        speeds.append(config.segmentation.speed_threshold)
    # Keep sensitivity random numbers independent from the primary-fit stream.
    rng = np.random.default_rng(config.fitting.random_seed + 104729)
    rows: List[Dict[str, object]] = []
    weight_column = _weight_column(config)
    for angle in angles:
        for speed in speeds:
            params = SegmentParams(
                turning_angle_degrees=angle,
                speed_threshold=speed,
                minimum_run_frames=config.segmentation.minimum_run_frames,
                minimum_run_length=config.segmentation.minimum_run_length,
                exclude_first_and_last_run=config.segmentation.exclude_first_and_last_run,
            )
            runs, _ = segment_runs(displacements, params)
            runs = assign_conditions(
                runs,
                metadata,
                config.conditions.grouping_fields,
                file_table=file_table,
                require_verified=config.conditions.require_verified_metadata,
                fallback_to_file=config.conditions.fallback_to_file,
            )
            pools, pool_audit = build_pools(
                runs,
                minimum_tracks=config.empirical_pool.minimum_tracks_per_pool,
                minimum_runs=config.empirical_pool.minimum_runs_per_pool,
                complete_only=config.empirical_pool.use_complete_runs_only,
            )
            for _, audit in pool_audit.iterrows():
                condition = str(audit["condition"])
                base: Dict[str, object] = {
                    "turning_angle_degrees": float(angle),
                    "speed_threshold": (
                        float(speed) if speed is not None else np.nan
                    ),
                    "condition": condition,
                    "n_tracks": int(audit["n_tracks"]),
                    "n_runs": int(audit["n_runs"]),
                    "pool_eligible": bool(audit["eligible"]),
                    "status": "ineligible",
                    "mu_hat": np.nan,
                    "crossover": np.nan,
                    "lmax": np.nan,
                    "ks_statistic": np.nan,
                    "ks_bootstrap_p": np.nan,
                    "gof_bootstraps": 0,
                    "note": str(audit["ineligible_reason"])
                    if not pd.isna(audit["ineligible_reason"])
                    else "",
                }
                if not bool(audit["eligible"]):
                    rows.append(base)
                    continue
                group = pools[pools["condition"] == condition].copy()
                lengths_all = group["run_length_um"].to_numpy(dtype=float)
                lmax = select_lmax(
                    lengths_all,
                    method=config.fitting.lmax_method,
                    user_lmax=config.fitting.user_lmax,
                    quantile=config.fitting.lmax_quantile,
                )
                group = group[lengths_all <= lmax]
                lengths = group["run_length_um"].to_numpy(dtype=float)
                weights = group[weight_column].to_numpy(dtype=float)
                try:
                    fitted = fit_mixture(
                        lengths,
                        lmax=lmax,
                        crossover=config.fitting.crossover,
                        mu_bounds=config.fitting.mu_bounds,
                        weights=weights,
                        estimate_crossover=config.fitting.estimate_crossover,
                        crossover_bounds=config.fitting.crossover_bounds,
                    )
                    gof = parametric_bootstrap_ks(
                        lengths,
                        fitted,
                        config.fitting.segmentation_goodness_of_fit_bootstraps,
                        rng,
                        mu_bounds=config.fitting.mu_bounds,
                        weights=weights,
                        estimate_crossover=config.fitting.estimate_crossover,
                        crossover_bounds=config.fitting.crossover_bounds,
                    )
                    base.update(
                        {
                            "status": "ok",
                            "mu_hat": fitted.mu,
                            "crossover": fitted.crossover,
                            "lmax": fitted.lmax,
                            "ks_statistic": gof.statistic,
                            "ks_bootstrap_p": gof.p_value,
                            "gof_bootstraps": gof.n_bootstraps,
                            "note": "",
                        }
                    )
                except FitError as exc:
                    base.update({"status": "failed", "note": str(exc)})
                rows.append(base)
    return pd.DataFrame(rows)


def run_fit(config: Config) -> Dict[str, int]:
    """Fit one mixture, including one crossover, per eligible condition."""
    store = ArtifactStore(config.output.root)
    pools = store.read_frame(
        "preprocess", "pools", expected_config_hash=config.hash()
    )
    pool_audit = store.read_frame(
        "preprocess",
        "pool_audit",
        schema_name="pool_audit",
        expected_config_hash=config.hash(),
    )
    eligible = set(pool_audit.loc[pool_audit["eligible"].astype(bool), "condition"])
    rng = np.random.default_rng(config.fitting.random_seed)
    fit_rows: List[Dict[str, object]] = []
    sensitivity_rows: List[Dict[str, object]] = []
    bootstrap_rows: List[Dict[str, object]] = []
    comparison_rows: List[Dict[str, object]] = []
    power_rows: List[Dict[str, object]] = []
    weight_column = _weight_column(config)

    for condition in sorted(eligible):
        group = pools[pools["condition"] == condition].copy()
        lengths_all = group["run_length_um"].to_numpy(dtype=float)
        lmax = select_lmax(
            lengths_all,
            method=config.fitting.lmax_method,
            user_lmax=config.fitting.user_lmax,
            quantile=config.fitting.lmax_quantile,
        )
        keep = lengths_all <= lmax
        excluded_above_lmax = int((~keep).sum())
        group_fit = group.loc[keep].copy()
        lengths = group_fit["run_length_um"].to_numpy(dtype=float)
        weights = group_fit[weight_column].to_numpy(dtype=float)
        try:
            fitted = fit_mixture(
                lengths,
                lmax=lmax,
                crossover=config.fitting.crossover,
                mu_bounds=config.fitting.mu_bounds,
                weights=weights,
                estimate_crossover=config.fitting.estimate_crossover,
                crossover_bounds=config.fitting.crossover_bounds,
            )
            gof = parametric_bootstrap_ks(
                lengths,
                fitted,
                config.fitting.goodness_of_fit_bootstraps,
                rng,
                mu_bounds=config.fitting.mu_bounds,
                weights=weights,
                estimate_crossover=config.fitting.estimate_crossover,
                crossover_bounds=config.fitting.crossover_bounds,
            )
            clustered = clustered_mu_bootstrap(
                group_fit,
                config.fitting.clustered_bootstraps,
                rng,
                crossover=config.fitting.crossover,
                mu_bounds=config.fitting.mu_bounds,
                weighting=config.empirical_pool.primary_weighting,
                estimate_crossover=config.fitting.estimate_crossover,
                crossover_bounds=config.fitting.crossover_bounds,
                lmax_method=config.fitting.lmax_method,
                user_lmax=config.fitting.user_lmax,
                lmax_quantile=config.fitting.lmax_quantile,
            )
            comparison = grouped_model_comparison(
                group_fit,
                lmax=lmax,
                crossover=config.fitting.crossover,
                mu_bounds=config.fitting.mu_bounds,
                n_folds=config.fitting.cv_folds,
                rng=rng,
                weighting=config.empirical_pool.primary_weighting,
                estimate_crossover=config.fitting.estimate_crossover,
                crossover_bounds=config.fitting.crossover_bounds,
            )
            full_family_fits = {}
            for _, comparison_row in comparison.iterrows():
                model = str(comparison_row["model"])
                row = comparison_row.to_dict()
                row["condition"] = condition
                if model == "c_mixture":
                    row["full_log_likelihood"] = fitted.log_likelihood
                    row["n_parameters"] = fitted.n_parameters
                    row["parameters_json"] = json.dumps(
                        {
                            "mu": fitted.mu,
                            "crossover": fitted.crossover,
                        },
                        sort_keys=True,
                    )
                else:
                    try:
                        alternative_fit = fit_family(
                            model,
                            lengths,
                            lmax=lmax,
                            weights=weights,
                            crossover=fitted.crossover,
                        )
                        full_family_fits[model] = alternative_fit
                        row["full_log_likelihood"] = alternative_fit.log_likelihood
                        row["n_parameters"] = (
                            alternative_fit.n_parameters
                            + int(
                                config.fitting.estimate_crossover
                                and model == "power_law_cutoff"
                            )
                        )
                        row["parameters_json"] = json.dumps(
                            dict(
                                alternative_fit.parameters,
                                crossover=alternative_fit.crossover,
                            ),
                            sort_keys=True,
                        )
                    except FitError as exc:
                        row["full_log_likelihood"] = np.nan
                        row["n_parameters"] = np.nan
                        row["parameters_json"] = ""
                        row["full_fit_error"] = str(exc)
                comparison_rows.append(row)

            mixture_cv = comparison[comparison["model"] == "c_mixture"].iloc[0]
            alternatives_cv = comparison[comparison["model"] != "c_mixture"]
            best_alternative_cv = alternatives_cv.sort_values(
                "mean_test_log_density", ascending=False
            ).iloc[0]
            score_difference = float(
                mixture_cv["mean_test_log_density"]
                - best_alternative_cv["mean_test_log_density"]
            )
            se_values = np.asarray(
                [
                    mixture_cv["se_test_log_density"],
                    best_alternative_cv["se_test_log_density"],
                ],
                dtype=float,
            )
            score_se = (
                float(np.sqrt(np.sum(np.square(se_values))))
                if np.all(np.isfinite(se_values))
                else np.nan
            )
            mixture_clearly_supported = bool(
                gof.p_value >= 0.05
                and score_difference > (1.96 * score_se if np.isfinite(score_se) else 0.0)
            )
            interpretation = (
                "fitted_mixture_exponent"
                if mixture_clearly_supported
                else "effective_exponent"
            )

            best_alternative_name = str(best_alternative_cv["model"])
            if (
                gof.n_bootstraps > 0
                and best_alternative_name in full_family_fits
                and gof.bootstrap_statistics.size
            ):
                critical = float(np.quantile(gof.bootstrap_statistics, 0.95))
                power = ks_power_curve(
                    full_family_fits[best_alternative_name],
                    sample_sizes=[
                        max(25, len(lengths) // 4),
                        max(25, len(lengths) // 2),
                        len(lengths),
                    ],
                    critical_value=critical,
                    n_replicates=min(200, max(50, gof.n_bootstraps)),
                    rng=rng,
                    crossover=fitted.crossover,
                    mu_bounds=config.fitting.mu_bounds,
                    estimate_crossover=config.fitting.estimate_crossover,
                    crossover_bounds=config.fitting.crossover_bounds,
                )
                power["condition"] = condition
                power_rows.extend(power.to_dict("records"))
            fit_rows.append(
                {
                    "condition": condition,
                    "status": "ok",
                    "weighting": config.empirical_pool.primary_weighting,
                    "interpretation": interpretation,
                    "mu_hat": fitted.mu,
                    "mu_ci_low": clustered.ci_low,
                    "mu_ci_high": clustered.ci_high,
                    "lmax": fitted.lmax,
                    "crossover": fitted.crossover,
                    "crossover_ci_low": clustered.crossover_ci_low,
                    "crossover_ci_high": clustered.crossover_ci_high,
                    "crossover_estimated": fitted.crossover_estimated,
                    "crossover_bound_low": fitted.crossover_bound_low,
                    "crossover_bound_high": fitted.crossover_bound_high,
                    "n_runs": fitted.n_observations,
                    "n_tracks": int(group_fit["track_uid"].nunique()),
                    "effective_sample_size": fitted.effective_sample_size,
                    "log_likelihood": fitted.log_likelihood,
                    "ks_statistic": gof.statistic,
                    "ks_bootstrap_p": gof.p_value,
                    "gof_bootstraps": gof.n_bootstraps,
                    "clustered_bootstraps_successful": clustered.n_successful,
                    "clustered_bootstraps_failed": clustered.n_failed,
                    "fit_at_mu_bound": fitted.at_bound,
                    "fit_at_crossover_bound": fitted.crossover_at_bound,
                    "excluded_above_lmax": excluded_above_lmax,
                    "best_cv_model": str(comparison.iloc[0]["model"]),
                    "best_alternative_model": best_alternative_name,
                    "mixture_minus_best_alternative_log_score": score_difference,
                    "mixture_clearly_supported": mixture_clearly_supported,
                    "note": (
                        ""
                        if mixture_clearly_supported
                        else "Mixture not clearly supported; report mu as an effective exponent."
                    ),
                }
            )
            for index, draw in enumerate(clustered.draws):
                bootstrap_rows.append(
                    {
                        "condition": condition,
                        "replicate": index,
                        "mu_hat": draw,
                        "crossover": clustered.crossover_draws[index],
                        "successful": bool(np.isfinite(draw)),
                    }
                )
        except FitError as exc:
            fit_rows.append(
                {
                    "condition": condition,
                    "status": "failed",
                    "weighting": config.empirical_pool.primary_weighting,
                    "interpretation": "not_fitted",
                    "note": str(exc),
                }
            )

        for crossover in config.fitting.crossover_sensitivity:
            try:
                sensitivity_fit = fit_mixture(
                    lengths,
                    lmax=lmax,
                    crossover=crossover,
                    mu_bounds=config.fitting.mu_bounds,
                    weights=weights,
                )
                sensitivity_rows.append(
                    {
                        "condition": condition,
                        "crossover": crossover,
                        "mu_hat": sensitivity_fit.mu,
                        "status": "ok",
                        "note": "",
                    }
                )
            except FitError as exc:
                sensitivity_rows.append(
                    {
                        "condition": condition,
                        "crossover": crossover,
                        "mu_hat": np.nan,
                        "status": "failed",
                        "note": str(exc),
                    }
                )

    fits = pd.DataFrame(fit_rows)
    sensitivity = pd.DataFrame(sensitivity_rows)
    bootstrap_draws = pd.DataFrame(bootstrap_rows)
    model_comparison = pd.DataFrame(comparison_rows)
    power_curve = pd.DataFrame(power_rows)
    displacements = store.read_frame(
        "preprocess",
        "displacements",
        schema_name="displacements",
        expected_config_hash=config.hash(),
    )
    file_table = store.read_frame(
        "preprocess", "file_table", expected_config_hash=config.hash()
    )
    segmentation_sensitivity = _fit_segmentation_sensitivity(
        config, displacements, file_table
    )
    grid = (
        segmentation_sensitivity[
            ["turning_angle_degrees", "speed_threshold"]
        ]
        .drop_duplicates()
        .sort_values(
            ["turning_angle_degrees", "speed_threshold"],
            na_position="first",
        )
    )
    segmentation_shape_sensitivity = grid.copy()
    segmentation_shape_sensitivity["headline_G"] = np.nan
    segmentation_shape_sensitivity["status"] = "pending_target_grid"
    segmentation_shape_sensitivity["note"] = (
        "Confirmatory target grid is not frozen because the replay pilot failed "
        "its event-count adequacy checks."
    )
    common_meta = {"config_hash": config.hash()}
    store.write_frame(
        "fit", "mixture_fits", fits, schema_name="fits", extra_metadata=common_meta
    )
    store.write_frame(
        "fit", "crossover_sensitivity", sensitivity, extra_metadata=common_meta
    )
    store.write_frame(
        "fit", "clustered_bootstrap_draws", bootstrap_draws, extra_metadata=common_meta
    )
    store.write_frame(
        "fit", "model_comparison", model_comparison, extra_metadata=common_meta
    )
    store.write_frame("fit", "gof_power_curve", power_curve, extra_metadata=common_meta)
    store.write_frame(
        "fit",
        "segmentation_fit_sensitivity",
        segmentation_sensitivity,
        extra_metadata=common_meta,
    )
    store.write_frame(
        "fit",
        "segmentation_shape_sensitivity",
        segmentation_shape_sensitivity,
        extra_metadata=common_meta,
    )
    summary = {
        "eligible_conditions": int(len(eligible)),
        "successful_fits": int((fits["status"] == "ok").sum()) if not fits.empty else 0,
        "failed_fits": int((fits["status"] == "failed").sum()) if not fits.empty else 0,
    }
    store.write_json("fit", "summary", summary)
    _store_run_context(store, config, "fit")
    logger.info("Fit summary: %s", json.dumps(summary, sort_keys=True))
    return summary


def _track_vector_sequence(runs: pd.DataFrame) -> np.ndarray:
    ordered = runs.sort_values(
        ["frame_from", "frame_to", "fragment_uid", "run_id"], kind="mergesort"
    )
    return ordered[["run_vector_x", "run_vector_y", "run_vector_z"]].to_numpy(
        dtype=float
    )


def _pad_vector_sequences(
    sequences: Sequence[np.ndarray],
) -> tuple:
    """Pad variable-length vector sequences for one vectorized replay batch."""
    if not sequences:
        raise ValueError("Cannot pad an empty sequence list")
    maximum = max(len(sequence) for sequence in sequences)
    vectors = np.zeros((len(sequences), maximum, 3), dtype=float)
    valid = np.zeros((len(sequences), maximum), dtype=bool)
    for index, sequence in enumerate(sequences):
        length = len(sequence)
        vectors[index, :length] = sequence
        valid[index, :length] = True
    return vectors, valid


def _select_budget_quantile_tracks(
    budget_by_track: Dict[str, float], maximum_tracks: int
) -> List[str]:
    """Select deterministic, budget-representative pilot tracks."""
    if maximum_tracks < 1:
        raise ValueError("maximum_tracks must be positive")
    ordered = sorted(
        budget_by_track,
        key=lambda track: (budget_by_track[track], str(track)),
    )
    if len(ordered) <= maximum_tracks:
        return ordered
    quantiles = (np.arange(maximum_tracks, dtype=float) + 0.5) / maximum_tracks
    indices = np.minimum(
        (quantiles * len(ordered)).astype(int), len(ordered) - 1
    )
    return [ordered[index] for index in indices]


def run_replay(config: Config, pilot: bool = False) -> Dict[str, int]:
    """Replay the six-class control ladder against every configured target."""
    store = ArtifactStore(config.output.root)
    runs = store.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config.hash(),
    )
    pools = store.read_frame(
        "preprocess", "pools", expected_config_hash=config.hash()
    )
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config.hash(),
    )
    fits = fits[fits["status"] == "ok"]
    if fits.empty:
        raise FitError("No successful fitted conditions are available for replay")

    stage = "replay_pilot" if pilot else "replay"
    rng = np.random.default_rng(config.replay.random_seed)
    n_replays = (
        min(
            config.replay.pilot_replays_per_track,
            config.replay.replays_per_track_target,
        )
        if pilot
        else config.replay.replays_per_track_target
    )
    max_tracks = config.replay.pilot_tracks_per_condition if pilot else None
    weight_column = _weight_column(config)

    target_rows: List[Dict[str, object]] = []
    targets = []
    for area in config.replay.projected_areas:
        for shape in config.replay.shapes:
            dimension = dimension_from_projected_area(
                area, shape, config.replay.detection_radius
            )
            target = Target(
                shape=shape,
                dimension=dimension,
                side_length=config.replay.domain_side_length,
                detection_radius=config.replay.detection_radius,
            )
            targets.append((area, target))
            target_rows.append(
                {
                    "shape": shape,
                    "projected_area": area,
                    "dimension": dimension,
                    "detection_radius": config.replay.detection_radius,
                    "reconstructed_projected_area": projected_area(
                        dimension, shape, config.replay.detection_radius
                    ),
                    "periodic_self_overlap": bool(
                        (
                            shape in ("Ball", "Disk")
                            and dimension / 2.0 + config.replay.detection_radius
                            > config.replay.domain_side_length / 2.0
                        )
                        or (
                            shape == "Line"
                            and dimension + 2.0 * config.replay.detection_radius
                            > config.replay.domain_side_length
                        )
                    ),
                }
            )

    common_meta = {"config_hash": config.hash(), "pilot": pilot}
    trial_rows: List[Dict[str, object]] = []
    selection_rows: List[Dict[str, object]] = []
    track_summary_frames: List[pd.DataFrame] = []
    n_trial_rows = 0
    trial_writer = (
        None
        if pilot
        else store.incremental_frame_writer(
            stage,
            "replay_trials",
            schema_name="replay_trials",
            extra_metadata=common_meta,
        )
    )
    fit_map = fits.set_index("condition")
    for condition in sorted(fit_map.index):
        condition_runs = runs[
            (runs["condition"] == condition)
            & np.isfinite(runs["run_length_um"].to_numpy(dtype=float))
            & (runs["run_length_um"].to_numpy(dtype=float) > 0.0)
        ]
        condition_pool = pools[pools["condition"] == condition]
        if condition_runs.empty or condition_pool.empty:
            logger.warning("%s: no replayable runs or pool; skipping", condition)
            continue
        sequence_by_track = {
            str(track_uid): _track_vector_sequence(group)
            for track_uid, group in condition_runs.groupby("track_uid", sort=True)
        }
        all_track_sequences = list(sequence_by_track.values())
        condition_block_pool = build_block_pool(
            all_track_sequences, config.replay.block_lengths[0]
        )
        budget_by_track = {
            track_uid: float(
                np.cumsum(np.linalg.norm(sequence, axis=1))[-1]
            )
            for track_uid, sequence in sequence_by_track.items()
        }
        track_ids = sorted(budget_by_track)
        if max_tracks is not None:
            track_ids = _select_budget_quantile_tracks(
                budget_by_track, max_tracks
            )
        selected_tracks = set(track_ids)
        for track_uid in sorted(budget_by_track):
            selection_rows.append(
                {
                    "condition": condition,
                    "track_uid": track_uid,
                    "budget": budget_by_track[track_uid],
                    "n_runs": int(len(sequence_by_track[track_uid])),
                    "selected": track_uid in selected_tracks,
                    "selection_method": (
                        "budget_quantiles" if pilot else "all_tracks"
                    ),
                }
            )
        pool_lengths = condition_pool["run_length_um"].to_numpy(dtype=float)
        pool_weights = condition_pool[weight_column].to_numpy(dtype=float)
        pool_probabilities = pool_weights / pool_weights.sum()
        pool_mean = float(np.dot(pool_probabilities, pool_lengths))
        fit_row = fit_map.loc[condition]
        fitted_mu = float(fit_row["mu_hat"])
        lmax = float(fit_row["lmax"])
        crossover = float(fit_row["crossover"])

        for track_uid in track_ids:
            track_trial_rows: List[Dict[str, object]] = []
            track_runs = condition_runs[condition_runs["track_uid"] == track_uid]
            source_vectors = sequence_by_track[str(track_uid)]
            budget = budget_by_track[str(track_uid)]
            if budget <= 0.0:
                continue
            expected_counts = [
                len(source_vectors),
                int(np.ceil(4.0 * budget / pool_mean)) + 16,
                int(np.ceil(4.0 * budget / mixture_mean(fitted_mu, lmax, crossover)))
                + 16,
                int(np.ceil(4.0 * budget / mixture_mean(2.0, lmax, crossover))) + 16,
            ]
            common_size = max(64, max(expected_counts))

            starts = np.empty((n_replays, 3), dtype=float)
            sequences_by_class = {
                class_name: [] for class_name in TRAJECTORY_CLASSES
            }
            for replay_index in range(n_replays):
                common = make_common_random_numbers(common_size, rng)
                starts[replay_index] = rng.uniform(
                    0.0, config.replay.domain_side_length, size=3
                )
                for class_name in TRAJECTORY_CLASSES:
                    sequences_by_class[class_name].append(
                        build_trajectory_class(
                            class_name,
                            source_vectors=source_vectors,
                            budget=budget,
                            rng=rng,
                            pool_lengths=pool_lengths,
                            pool_weights=pool_weights,
                            fitted_mu=fitted_mu,
                            lmax=lmax,
                            crossover=crossover,
                            track_sequences=all_track_sequences,
                            block_pool=condition_block_pool,
                            block_length=config.replay.block_lengths[0],
                            common=common,
                        )
                    )
            target_objects = [target for _, target in targets]
            budgets = np.full(n_replays, budget, dtype=float)
            for class_name, sequences in sequences_by_class.items():
                padded_vectors, valid_steps = _pad_vector_sequences(sequences)
                target_results = replay_targets_batch(
                    padded_vectors,
                    starts,
                    target_objects,
                    budgets,
                    valid_steps=valid_steps,
                    overshoot_policy=config.replay.overshoot_policy,
                    check_initial_position=config.replay.check_initial_position,
                )
                for (area, target), result in zip(targets, target_results):
                    for replay_index in range(n_replays):
                        track_trial_rows.append(
                            {
                                "condition": condition,
                                "video_id": str(track_runs["video_id"].iloc[0]),
                                "track_uid": track_uid,
                                "replay_index": replay_index,
                                "trajectory_class": class_name,
                                "shape": target.shape,
                                "projected_area": area,
                                "budget": budget,
                                "detected": bool(result.detected[replay_index]),
                                "detection_distance": float(
                                    result.detection_distance[replay_index]
                                ),
                                "restricted_detection_distance": (
                                    float(
                                        result.restricted_detection_distance[
                                            replay_index
                                        ]
                                    )
                                ),
                                "first_hit_step": int(
                                    result.first_hit_step[replay_index]
                                ),
                                "endpoints_checked": int(
                                    result.endpoints_checked[replay_index]
                                ),
                                "overshoot_policy": config.replay.overshoot_policy,
                            }
                        )
            track_trials = pd.DataFrame(track_trial_rows)
            if track_trials.empty:
                continue
            n_trial_rows += len(track_trials)
            track_summary_frames.append(
                summarise_trials(
                    track_trials,
                    [
                        "condition",
                        "video_id",
                        "track_uid",
                        "trajectory_class",
                        "shape",
                        "projected_area",
                    ],
                )
            )
            if pilot:
                trial_rows.extend(track_trial_rows)
            else:
                assert trial_writer is not None
                trial_writer.append(track_trials)
        logger.info(
            "%s: replayed %d tracks x %d paired replicates",
            condition,
            len(track_ids),
            n_replays,
        )

    if not track_summary_frames:
        raise RuntimeError("Replay produced no trials")
    if trial_writer is not None:
        trial_writer.close()
    trials = pd.DataFrame(trial_rows)
    summary_fields = [
        "condition",
        "trajectory_class",
        "shape",
        "projected_area",
    ]
    track_summaries = pd.concat(track_summary_frames, ignore_index=True)
    summaries = combine_trial_summaries(track_summaries, summary_fields)
    global_summaries = combine_trial_summaries(
        track_summaries, ["trajectory_class", "shape", "projected_area"]
    )
    sensitivity = shape_sensitivity(
        summaries,
        ["condition", "trajectory_class", "projected_area"],
    )
    boundary_cells = (
        (global_summaries["detection_probability"] <= 0.0)
        | (global_summaries["detection_probability"] >= 1.0)
    )
    diagnostics = global_summaries.loc[
        (
            (global_summaries["n_detected"] < config.replay.pilot_min_detected)
            | (
                global_summaries["n_trials"] - global_summaries["n_detected"]
                < config.replay.pilot_min_undetected
            )
        ),
        [
            "trajectory_class",
            "shape",
            "projected_area",
            "n_trials",
            "detection_probability",
        ],
    ].copy()
    diagnostics["n_detected"] = global_summaries.loc[diagnostics.index, "n_detected"]
    diagnostics["n_undetected"] = (
        diagnostics["n_trials"] - diagnostics["n_detected"]
    )
    diagnostics["minimum_detected_required"] = config.replay.pilot_min_detected
    diagnostics["minimum_undetected_required"] = config.replay.pilot_min_undetected
    diagnostics["issue"] = "insufficient_pilot_events"
    underpowered_cells = len(diagnostics)

    if pilot:
        store.write_frame(
            stage,
            "replay_trials",
            trials,
            schema_name="replay_trials",
            extra_metadata=common_meta,
        )
    store.write_frame(stage, "track_summaries", track_summaries, extra_metadata=common_meta)
    store.write_frame(stage, "outcome_summaries", summaries, extra_metadata=common_meta)
    store.write_frame(
        stage, "global_outcome_summaries", global_summaries, extra_metadata=common_meta
    )
    store.write_frame(stage, "shape_sensitivity", sensitivity, extra_metadata=common_meta)
    store.write_frame(
        stage,
        "replay_track_selection",
        pd.DataFrame(selection_rows),
        extra_metadata=common_meta,
    )
    store.write_frame(
        stage, "target_dimensions", pd.DataFrame(target_rows), extra_metadata=common_meta
    )
    store.write_frame(
        stage, "pilot_diagnostics", diagnostics, extra_metadata=common_meta
    )
    if not pilot:
        outcome_intervals, adjacent_contrasts, shape_intervals = (
            clustered_outcome_intervals(
                track_summaries,
                config.fitting.clustered_bootstraps,
                rng,
            )
        )
        store.write_frame(
            stage, "clustered_outcome_intervals", outcome_intervals, extra_metadata=common_meta
        )
        store.write_frame(
            stage, "adjacent_class_contrasts", adjacent_contrasts, extra_metadata=common_meta
        )
        store.write_frame(
            stage, "clustered_shape_intervals", shape_intervals, extra_metadata=common_meta
        )
    result_summary = {
        "pilot": int(pilot),
        "conditions": int(track_summaries["condition"].nunique()),
        "tracks": int(track_summaries["track_uid"].nunique()),
        "trials": int(n_trial_rows),
        "boundary_probability_cells": int(boundary_cells.sum()),
        "underpowered_cells": int(underpowered_cells),
    }
    store.write_json(stage, "summary", result_summary)
    _store_run_context(store, config, stage)
    logger.info("Replay summary: %s", json.dumps(result_summary, sort_keys=True))
    if pilot and underpowered_cells:
        logger.warning(
            "Pilot has %d target/class cells below the configured event-count minimums; "
            "review pilot_diagnostics.csv before freezing projected areas",
            int(underpowered_cells),
        )
    return result_summary


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ltdb-levy",
        description="LTDB directed-run and Levy-surrogate analysis pipeline",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("inspect", "Inspect source files and surface anomalies"),
        ("preprocess", "Clean, segment, assign conditions, and build pools"),
        ("fit", "Fit the C-mixture exponent for eligible pools"),
        ("report", "Generate figures, report documents, and analysis manifest"),
        ("prepare-c", "Generate fitted-mu configs for the existing C simulator"),
        ("compare-c", "Ingest C output and optionally run matched Python trials"),
    ):
        command = subparsers.add_parser(name, help=help_text)
        command.add_argument(
            "--config", type=Path, default=Path("configs/primary.yaml")
        )
    replay = subparsers.add_parser(
        "replay", help="Replay trajectory controls against shape-matched targets"
    )
    replay.add_argument("--config", type=Path, default=Path("configs/primary.yaml"))
    replay.add_argument(
        "--pilot",
        action="store_true",
        help="Use configured pilot replicates and budget-quantile tracks",
    )
    focused = subparsers.add_parser(
        "focused-replay",
        help="Run focused finite and unbounded empirical-versus-fitted replay",
    )
    focused.add_argument(
        "--config", type=Path, default=Path("configs/primary.yaml")
    )
    focused.add_argument(
        "--condition",
        default=None,
        help="Exact fitted condition label (defaults to the eligible NK condition)",
    )
    focused.add_argument(
        "--length-unit-um",
        type=float,
        default=None,
        help="Physical size of one model unit; default uses the fitted crossover c",
    )
    focused.add_argument(
        "--projected-areas",
        type=float,
        nargs="+",
        default=[4.0, 8.0, 16.0, 32.0],
    )
    focused.add_argument("--finite-replays", type=int, default=10_000)
    focused.add_argument("--unbounded-trials", type=int, default=400)
    focused.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Process count; 0 uses all available CPUs",
    )
    focused.add_argument("--finite-batch-size", type=int, default=2_000)
    focused.add_argument("--unbounded-chunk-trials", type=int, default=25)
    focused.add_argument("--unbounded-step-block", type=int, default=512)
    focused.add_argument("--seed", type=int, default=24680)
    all_conditions = subparsers.add_parser(
        "all-conditions-replay",
        help="Run equal-run finite and unbounded replay for every eligible condition",
    )
    all_conditions.add_argument(
        "--config",
        type=Path,
        default=Path("configs/equal_run_weight_all_conditions.yaml"),
    )
    all_conditions.add_argument(
        "--estimate-only",
        action="store_true",
        help="Write the condition manifest and runtime estimate without replay",
    )
    all_conditions.add_argument(
        "--projected-areas",
        type=float,
        nargs="+",
        default=None,
        help="Common dimensionless target areas; defaults to replay.projected_areas",
    )
    all_conditions.add_argument("--finite-replays", type=int, default=10_000)
    all_conditions.add_argument("--unbounded-trials", type=int, default=2_000)
    all_conditions.add_argument("--workers", type=int, default=32)
    all_conditions.add_argument("--finite-batch-size", type=int, default=1_000)
    all_conditions.add_argument("--unbounded-chunk-trials", type=int, default=25)
    all_conditions.add_argument("--unbounded-step-block", type=int, default=256)
    all_conditions.add_argument("--seed", type=int, default=24680)
    all_conditions.add_argument(
        "--no-resume",
        action="store_true",
        help="Rerun completed condition directories instead of reusing them",
    )
    ordered_rotation = subparsers.add_parser(
        "ordered-rotation-replay",
        help=(
            "Add exact-track, independently rotated-vector finite replay "
            "without replacing the signed baseline"
        ),
    )
    ordered_rotation.add_argument(
        "--config",
        type=Path,
        default=Path("configs/equal_run_weight_run_eligible_conditions.yaml"),
    )
    ordered_rotation.add_argument("--finite-replays", type=int, default=10_000)
    ordered_rotation.add_argument("--workers", type=int, default=32)
    ordered_rotation.add_argument("--finite-batch-size", type=int, default=1_000)
    ordered_rotation.add_argument("--seed", type=int, default=913_579)
    ordered_rotation.add_argument(
        "--baseline-stage",
        default="all_conditions_replay",
    )
    ordered_rotation.add_argument(
        "--output-stage",
        default="ordered_length_rotation_replay",
    )
    selected_figures = subparsers.add_parser(
        "selected-publication-figures",
        help=(
            "Build validated publication figures for the selected NK, "
            "neutrophil, and T-cell conditions"
        ),
    )
    selected_figures.add_argument(
        "--config",
        type=Path,
        default=Path("configs/equal_run_weight_run_eligible_conditions.yaml"),
    )
    selected_figures.add_argument(
        "--baseline-stage",
        default="all_conditions_replay",
    )
    selected_figures.add_argument(
        "--augmented-stage",
        default="ordered_length_rotation_replay",
    )
    selected_figures.add_argument(
        "--output-stage",
        default="publication_selected_conditions",
    )
    line_feasible = subparsers.add_parser(
        "line-feasible-short-replay",
        help=(
            "Run the short condition-specific target sweep using the Line "
            "feasibility limit for all three shapes"
        ),
    )
    line_feasible.add_argument(
        "--config",
        type=Path,
        default=Path("configs/equal_run_weight_run_eligible_conditions.yaml"),
    )
    line_feasible.add_argument(
        "--estimate-only",
        action="store_true",
        help="Write the exact target manifest and runtime estimate without replay",
    )
    line_feasible.add_argument(
        "--candidate-areas",
        type=float,
        nargs="+",
        default=None,
        help="Candidate A/c^2 values; each condition retains its Line-feasible subset",
    )
    line_feasible.add_argument("--finite-replays", type=int, default=2_000)
    line_feasible.add_argument("--unbounded-trials", type=int, default=500)
    line_feasible.add_argument("--workers", type=int, default=64)
    line_feasible.add_argument("--finite-batch-size", type=int, default=500)
    line_feasible.add_argument("--unbounded-chunk-trials", type=int, default=50)
    line_feasible.add_argument("--unbounded-step-block", type=int, default=256)
    line_feasible.add_argument("--seed", type=int, default=731_291)
    line_feasible.add_argument(
        "--output-stage",
        default="line_feasible_short_replay",
    )
    line_feasible.add_argument(
        "--no-resume",
        action="store_true",
        help="Rerun completed condition checkpoints",
    )
    compare_c = subparsers.choices["compare-c"]
    compare_c.add_argument(
        "--python-trials",
        type=int,
        default=0,
        help="Matched unbounded Python trials per available C target cell",
    )
    compare_c.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Optional cap on target cells for an incremental cross-check",
    )
    return parser


def main(argv: Sequence[str] = None) -> int:
    """CLI main function; returns a process exit code for easy testing."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    try:
        config = load_config(args.config)
        if args.command == "inspect":
            run_inspect(config)
        elif args.command == "preprocess":
            run_preprocess(config)
        elif args.command == "fit":
            run_fit(config)
        elif args.command == "replay":
            run_replay(config, pilot=args.pilot)
        elif args.command == "focused-replay":
            from .nk_experiment import (
                DEFAULT_NK_CONDITION,
                FocusedReplaySettings,
                run_focused_replay,
            )

            settings = FocusedReplaySettings(
                condition=args.condition or DEFAULT_NK_CONDITION,
                length_unit_um=args.length_unit_um,
                projected_areas=tuple(args.projected_areas),
                finite_replays_per_track=args.finite_replays,
                unbounded_trials_per_cell=args.unbounded_trials,
                workers=args.workers,
                finite_batch_size=args.finite_batch_size,
                unbounded_chunk_trials=args.unbounded_chunk_trials,
                unbounded_step_block=args.unbounded_step_block,
                random_seed=args.seed,
            )
            summary = run_focused_replay(config, settings)
            logger.info(
                "Focused replay summary: %s", json.dumps(summary, sort_keys=True)
            )
        elif args.command == "all-conditions-replay":
            from .multicondition_experiment import (
                AllConditionsReplaySettings,
                run_all_conditions_replay,
                write_runtime_plan,
            )

            settings = AllConditionsReplaySettings(
                projected_areas=tuple(
                    args.projected_areas
                    if args.projected_areas is not None
                    else config.replay.projected_areas
                ),
                finite_replays_per_track=args.finite_replays,
                unbounded_trials_per_cell=args.unbounded_trials,
                workers=args.workers,
                finite_batch_size=args.finite_batch_size,
                unbounded_chunk_trials=args.unbounded_chunk_trials,
                unbounded_step_block=args.unbounded_step_block,
                random_seed=args.seed,
                resume=not args.no_resume,
            )
            plan = write_runtime_plan(config, settings)
            logger.info(
                "All-condition replay plan: %s",
                json.dumps(plan, sort_keys=True),
            )
            if not args.estimate_only:
                summary = run_all_conditions_replay(config, settings)
                logger.info(
                    "All-condition replay summary: %s",
                    json.dumps(summary, sort_keys=True),
                )
        elif args.command == "ordered-rotation-replay":
            from .ordered_rotation_experiment import (
                OrderedRotationReplaySettings,
                run_ordered_length_rotation_replay,
            )

            settings = OrderedRotationReplaySettings(
                finite_replays_per_track=args.finite_replays,
                workers=args.workers,
                finite_batch_size=args.finite_batch_size,
                random_seed=args.seed,
                baseline_stage=args.baseline_stage,
                output_stage=args.output_stage,
            )
            summary = run_ordered_length_rotation_replay(config, settings)
            logger.info(
                "Ordered-rotation replay summary: %s",
                json.dumps(summary, sort_keys=True),
            )
        elif args.command == "selected-publication-figures":
            from .publication_figures import run_selected_publication_figures

            summary = run_selected_publication_figures(
                config,
                baseline_stage=args.baseline_stage,
                augmented_stage=args.augmented_stage,
                output_stage=args.output_stage,
            )
            logger.info(
                "Selected publication figure summary: %s",
                json.dumps(summary, sort_keys=True),
            )
        elif args.command == "line-feasible-short-replay":
            from .line_feasible_experiment import (
                DEFAULT_CANDIDATE_AREAS,
                LineFeasibleSweepSettings,
                run_line_feasible_sweep,
                write_runtime_plan,
            )

            settings = LineFeasibleSweepSettings(
                candidate_areas=tuple(
                    args.candidate_areas
                    if args.candidate_areas is not None
                    else DEFAULT_CANDIDATE_AREAS
                ),
                finite_replays_per_track=args.finite_replays,
                unbounded_trials_per_cell=args.unbounded_trials,
                workers=args.workers,
                finite_batch_size=args.finite_batch_size,
                unbounded_chunk_trials=args.unbounded_chunk_trials,
                unbounded_step_block=args.unbounded_step_block,
                random_seed=args.seed,
                output_stage=args.output_stage,
                resume=not args.no_resume,
            )
            plan = write_runtime_plan(config, settings)
            logger.info(
                "Line-feasible short replay plan: %s",
                json.dumps(plan, sort_keys=True),
            )
            if not args.estimate_only:
                summary = run_line_feasible_sweep(config, settings)
                logger.info(
                    "Line-feasible short replay summary: %s",
                    json.dumps(summary, sort_keys=True),
                )
        elif args.command == "report":
            from .report import build_report

            paths = build_report(config)
            logger.info("Generated %d report artifacts", len(paths))
        elif args.command == "prepare-c":
            from .c_experiment import prepare_c_experiments

            paths = prepare_c_experiments(config)
            logger.info("Generated %d C experiment files", len(paths))
        elif args.command == "compare-c":
            from .c_compare import run_c_comparison

            summary = run_c_comparison(
                config,
                python_trials=args.python_trials,
                max_cells=args.max_cells,
            )
            logger.info("C comparison summary: %s", json.dumps(summary, sort_keys=True))
        else:  # pragma: no cover - argparse enforces the choices
            parser.error("unknown command")
    except (LtdbLevyError, FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
        logger.error("%s", exc)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
