"""Unbounded Python replay and ingestion of the existing C simulator output."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactStore
from .directions import isotropic_directions
from .dist.mixture import sample
from .targets import Target, dimension_from_projected_area, wrap_positions


@dataclass(frozen=True)
class UnboundedReplayResult:
    """First-passage distances and step counts for independent walkers."""

    detection_distance: np.ndarray
    first_hit_step: np.ndarray


def simulate_unbounded_mixture(
    target: Target,
    mu: float,
    lmax: float,
    n_trials: int,
    rng: np.random.Generator,
    crossover: float = 1.0,
    chunk_steps: int = 256,
    max_chunks: int = 100_000,
) -> UnboundedReplayResult:
    """Simulate C-compatible endpoint first passage in bounded memory.

    Walkers are advanced in vectorized chunks. The finite ``max_chunks`` guard
    turns an unexpectedly expensive target into a clear failure instead of an
    infinite run.
    """
    if n_trials < 1:
        raise ValueError("n_trials must be positive")
    if chunk_steps < 1 or max_chunks < 1:
        raise ValueError("chunk_steps and max_chunks must be positive")

    positions = rng.uniform(0.0, target.side_length, size=(n_trials, 3))
    travelled = np.zeros(n_trials, dtype=float)
    steps = np.zeros(n_trials, dtype=np.int64)
    detection_distance = np.full(n_trials, np.nan, dtype=float)
    first_hit_step = np.full(n_trials, -1, dtype=np.int64)
    detected = np.zeros(n_trials, dtype=bool)

    for _ in range(max_chunks):
        active = np.flatnonzero(~detected)
        if active.size == 0:
            break
        n_draws = int(active.size * chunk_steps)
        lengths = sample(
            mu,
            lmax,
            n_draws,
            rng,
            crossover=crossover,
        ).reshape(active.size, chunk_steps)
        directions = isotropic_directions(n_draws, rng).reshape(
            active.size, chunk_steps, 3
        )
        vectors = directions * lengths[..., None]
        cumulative_vectors = np.cumsum(vectors, axis=1)
        endpoints = wrap_positions(
            positions[active, None, :] + cumulative_vectors,
            target.side_length,
        )
        hits = target.contains(endpoints.reshape(-1, 3)).reshape(
            active.size, chunk_steps
        )
        hit_in_chunk = np.any(hits, axis=1)
        first = np.argmax(hits, axis=1)
        cumulative_distance = np.cumsum(lengths, axis=1)

        if np.any(hit_in_chunk):
            rows = np.flatnonzero(hit_in_chunk)
            trial_indices = active[rows]
            hit_steps = first[rows]
            detection_distance[trial_indices] = (
                travelled[trial_indices]
                + cumulative_distance[rows, hit_steps]
            )
            first_hit_step[trial_indices] = (
                steps[trial_indices] + hit_steps + 1
            )
            detected[trial_indices] = True

        missed_rows = np.flatnonzero(~hit_in_chunk)
        if missed_rows.size:
            trial_indices = active[missed_rows]
            positions[trial_indices] = endpoints[missed_rows, -1]
            travelled[trial_indices] += cumulative_distance[missed_rows, -1]
            steps[trial_indices] += chunk_steps
    else:
        unresolved = int((~detected).sum())
        raise RuntimeError(
            "Unbounded replay reached the {}-step guard with {} of {} trials "
            "still unresolved".format(
                chunk_steps * max_chunks, unresolved, n_trials
            )
        )

    return UnboundedReplayResult(
        detection_distance=detection_distance,
        first_hit_step=first_hit_step,
    )


def _resolve_repo_path(repository_root: Path, value: str) -> Path:
    path = (repository_root / value).resolve()
    root = repository_root.resolve()
    if root != path and root not in path.parents:
        raise ValueError("C output path escapes the repository: {}".format(value))
    return path


def ingest_c_outputs(
    index_path: Path,
    repository_root: Path,
    expected_analysis_hash: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read generated C outputs and return target summaries plus an audit."""
    index = pd.read_csv(index_path)
    required_index = {
        "condition",
        "mu_hat",
        "c_lmax",
        "num_runs",
        "output_path",
    }
    missing_index = sorted(required_index - set(index.columns))
    if missing_index:
        raise ValueError(
            "C experiment index is missing columns: {}".format(
                ", ".join(missing_index)
            )
        )
    if expected_analysis_hash is not None:
        if "analysis_hash" not in index:
            raise ValueError(
                "C experiment index has no analysis_hash; rerun prepare-c"
            )
        hashes = set(index["analysis_hash"].astype(str))
        if hashes != {expected_analysis_hash}:
            raise ValueError(
                "C experiment index is stale for the current analysis hash"
            )

    summaries: List[Dict[str, object]] = []
    audit: List[Dict[str, object]] = []
    required_output = {"TargetShape", "surface", "detection_time"}
    for _, entry in index.iterrows():
        output_path = _resolve_repo_path(
            repository_root, str(entry["output_path"])
        )
        base_audit: Dict[str, object] = {
            "condition": str(entry["condition"]),
            "output_path": str(output_path),
            "expected_runs_per_target": int(entry["num_runs"]),
            "expected_target_cells": (
                int(entry["expected_target_cells"])
                if "expected_target_cells" in entry
                and not pd.isna(entry["expected_target_cells"])
                else np.nan
            ),
        }
        if not output_path.is_file():
            audit.append(
                dict(
                    base_audit,
                    status="missing",
                    rows=0,
                    target_cells=0,
                    note="Run experiments/detection_time_fitted_mu.sh",
                )
            )
            continue
        try:
            frame = pd.read_csv(output_path)
        except (OSError, pd.errors.ParserError) as exc:
            audit.append(
                dict(
                    base_audit,
                    status="unreadable",
                    rows=0,
                    target_cells=0,
                    note=str(exc),
                )
            )
            continue
        missing_output = sorted(required_output - set(frame.columns))
        if missing_output:
            audit.append(
                dict(
                    base_audit,
                    status="invalid_schema",
                    rows=int(len(frame)),
                    target_cells=0,
                    note="missing {}".format(", ".join(missing_output)),
                )
            )
            continue
        frame["detection_time"] = pd.to_numeric(
            frame["detection_time"], errors="coerce"
        )
        frame["surface"] = pd.to_numeric(frame["surface"], errors="coerce")
        valid = frame[
            np.isfinite(frame["detection_time"])
            & np.isfinite(frame["surface"])
            & frame["TargetShape"].isin(("Ball", "Disk", "Line"))
        ].copy()
        groups = valid.groupby(["TargetShape", "surface"], sort=True)
        incomplete = 0
        for (shape, area), group in groups:
            values = group["detection_time"].to_numpy(dtype=float)
            n = int(len(values))
            if n != int(entry["num_runs"]):
                incomplete += 1
            standard_deviation = (
                float(np.std(values, ddof=1)) if n > 1 else np.nan
            )
            summaries.append(
                {
                    "condition": str(entry["condition"]),
                    "shape": str(shape),
                    "projected_area": float(area),
                    "mu_hat": float(entry["mu_hat"]),
                    "c_lmax": int(entry["c_lmax"]),
                    "n_trials": n,
                    "mean_detection_distance": float(np.mean(values)),
                    "sd_detection_distance": standard_deviation,
                    "se_detection_distance": (
                        standard_deviation / np.sqrt(n)
                        if n > 1
                        else np.nan
                    ),
                    "source_path": str(output_path),
                }
            )
        target_count_mismatch = (
            "expected_target_cells" in entry
            and not pd.isna(entry["expected_target_cells"])
            and groups.ngroups != int(entry["expected_target_cells"])
        )
        status = (
            "ok"
            if (
                len(valid) == len(frame)
                and incomplete == 0
                and not target_count_mismatch
            )
            else "incomplete"
        )
        notes = []
        if len(valid) != len(frame):
            notes.append("{} invalid rows".format(len(frame) - len(valid)))
        if incomplete:
            notes.append(
                "{} target cells differ from configured run count".format(
                    incomplete
                )
            )
        if target_count_mismatch:
            notes.append(
                "{} target cells found; {} expected".format(
                    groups.ngroups, int(entry["expected_target_cells"])
                )
            )
        audit.append(
            dict(
                base_audit,
                status=status,
                rows=int(len(frame)),
                target_cells=int(groups.ngroups),
                note="; ".join(notes),
            )
        )

    summary_columns = [
        "condition",
        "shape",
        "projected_area",
        "mu_hat",
        "c_lmax",
        "n_trials",
        "mean_detection_distance",
        "sd_detection_distance",
        "se_detection_distance",
        "source_path",
    ]
    audit_columns = [
        "condition",
        "output_path",
        "expected_runs_per_target",
        "expected_target_cells",
        "status",
        "rows",
        "target_cells",
        "note",
    ]
    return (
        pd.DataFrame(summaries, columns=summary_columns),
        pd.DataFrame(audit, columns=audit_columns),
    )


def run_c_comparison(
    config: Config,
    python_trials: int = 0,
    max_cells: Optional[int] = None,
) -> Dict[str, int]:
    """Ingest C results and optionally run matched unbounded Python trials."""
    if python_trials < 0:
        raise ValueError("python_trials cannot be negative")
    if max_cells is not None and max_cells < 1:
        raise ValueError("max_cells must be positive")

    repository_root = Path(__file__).resolve().parents[1]
    index_path = (
        repository_root
        / "experiments"
        / "configs"
        / "ltdb_fitted_mu"
        / "index.csv"
    )
    if not index_path.is_file():
        raise FileNotFoundError(
            "C experiment index is missing; run `ltdb-levy prepare-c` first"
        )
    analysis_hash = config.hash()
    c_summary, audit = ingest_c_outputs(
        index_path,
        repository_root,
        expected_analysis_hash=analysis_hash,
    )
    store = ArtifactStore(config.output.root)
    common_meta = {"config_hash": analysis_hash}
    store.write_frame(
        "c_comparison", "c_results_summary", c_summary, extra_metadata=common_meta
    )
    store.write_frame(
        "c_comparison", "c_ingestion_audit", audit, extra_metadata=common_meta
    )

    python_rows: List[Dict[str, object]] = []
    comparison = pd.DataFrame()
    if python_trials > 0 and not c_summary.empty:
        cells = c_summary.sort_values(
            ["condition", "shape", "projected_area"]
        )
        if max_cells is not None:
            cells = cells.head(max_cells)
        rng = np.random.default_rng(config.replay.random_seed)
        for _, cell in cells.iterrows():
            shape = str(cell["shape"])
            area = float(cell["projected_area"])
            target = Target(
                shape=shape,
                dimension=dimension_from_projected_area(
                    area, shape, config.replay.detection_radius
                ),
                side_length=config.replay.domain_side_length,
                detection_radius=config.replay.detection_radius,
            )
            result = simulate_unbounded_mixture(
                target,
                mu=float(cell["mu_hat"]),
                lmax=float(cell["c_lmax"]),
                n_trials=python_trials,
                rng=rng,
                crossover=config.fitting.crossover,
            )
            values = result.detection_distance
            standard_deviation = (
                float(np.std(values, ddof=1))
                if len(values) > 1
                else np.nan
            )
            python_rows.append(
                {
                    "condition": str(cell["condition"]),
                    "shape": shape,
                    "projected_area": area,
                    "n_trials": int(len(values)),
                    "mean_detection_distance": float(np.mean(values)),
                    "sd_detection_distance": standard_deviation,
                    "se_detection_distance": (
                        standard_deviation / np.sqrt(len(values))
                        if len(values) > 1
                        else np.nan
                    ),
                }
            )
        python_summary = pd.DataFrame(python_rows)
        store.write_frame(
            "c_comparison",
            "python_unbounded_summary",
            python_summary,
            extra_metadata=common_meta,
        )
        comparison = c_summary.merge(
            python_summary,
            on=["condition", "shape", "projected_area"],
            suffixes=("_c", "_python"),
        )
        comparison["mean_difference_python_minus_c"] = (
            comparison["mean_detection_distance_python"]
            - comparison["mean_detection_distance_c"]
        )
        comparison["combined_se"] = np.sqrt(
            comparison["se_detection_distance_python"] ** 2
            + comparison["se_detection_distance_c"] ** 2
        )
        comparison["z_difference"] = (
            comparison["mean_difference_python_minus_c"]
            / comparison["combined_se"]
        )
        comparison["agree_within_3se"] = (
            np.abs(comparison["z_difference"]) <= 3.0
        )
        store.write_frame(
            "c_comparison",
            "python_c_comparison",
            comparison,
            extra_metadata=common_meta,
        )

    result_summary = {
        "configured_conditions": int(len(audit)),
        "available_c_conditions": int((audit["status"] != "missing").sum()),
        "valid_c_target_cells": int(len(c_summary)),
        "python_target_cells": int(len(python_rows)),
        "comparison_target_cells": int(len(comparison)),
    }
    store.write_json("c_comparison", "summary", result_summary)
    store.write_json(
        "c_comparison",
        "run_context",
        {
            "run_id": config.run_id,
            "config_hash": analysis_hash,
            "config_path": str(config.source_path),
            "config": config.to_dict(),
        },
    )
    return result_summary
