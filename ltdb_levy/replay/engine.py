"""Vectorized finite-budget endpoint-only replay engine.

There is no trajectory ``while`` loop: budget inclusion is derived from
left-to-right cumulative path lengths, and every endpoint is tested in one
vectorized target-membership call.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..targets import Target, wrap_positions


@dataclass(frozen=True)
class ReplayBatchResult:
    """Results for a batch of padded trajectory sequences."""

    detected: np.ndarray
    first_hit_step: np.ndarray
    detection_distance: np.ndarray
    restricted_detection_distance: np.ndarray
    endpoints_checked: np.ndarray
    budget: np.ndarray


@dataclass(frozen=True)
class BudgetedPaths:
    """Target-independent endpoints and path distances for a replay batch."""

    endpoints: np.ndarray
    endpoint_mask: np.ndarray
    endpoint_distance: np.ndarray
    starts: np.ndarray
    budget: np.ndarray
    side_length: float


def prepare_budgeted_paths(
    vectors: np.ndarray,
    starts: np.ndarray,
    budgets: np.ndarray,
    side_length: float,
    valid_steps: Optional[np.ndarray] = None,
    overshoot_policy: str = "discard",
) -> BudgetedPaths:
    """Compute wrapped budget-valid endpoints once, independently of targets.

    Args:
        vectors: ``(n_trials, n_steps, 3)`` relocation vectors.
        starts: ``(n_trials, 3)`` initial positions.
        budgets: One non-negative path-length budget per trial.
        side_length: Periodic box side shared by all later targets.
        valid_steps: Optional boolean padding mask of shape ``(n_trials, n_steps)``.
        overshoot_policy: ``discard`` (primary) or labelled ``truncate`` sensitivity.
    """
    vector_array = np.asarray(vectors, dtype=float)
    start_array = np.asarray(starts, dtype=float)
    budget_array = np.asarray(budgets, dtype=float)
    if vector_array.ndim != 3 or vector_array.shape[2] != 3:
        raise ValueError("vectors must have shape (n_trials, n_steps, 3)")
    n_trials, n_steps, _ = vector_array.shape
    if start_array.shape != (n_trials, 3):
        raise ValueError("starts must have shape (n_trials, 3)")
    if budget_array.shape != (n_trials,):
        raise ValueError("budgets must have shape (n_trials,)")
    if np.any(~np.isfinite(vector_array)) or np.any(~np.isfinite(start_array)):
        raise ValueError("vectors and starts must be finite")
    if np.any(~np.isfinite(budget_array)) or np.any(budget_array < 0.0):
        raise ValueError("budgets must be finite and non-negative")
    if not np.isfinite(side_length) or side_length <= 0.0:
        raise ValueError("side_length must be finite and positive")
    if overshoot_policy not in ("discard", "truncate"):
        raise ValueError("overshoot_policy must be 'discard' or 'truncate'")
    if valid_steps is None:
        valid = np.ones((n_trials, n_steps), dtype=bool)
    else:
        valid = np.asarray(valid_steps, dtype=bool)
        if valid.shape != (n_trials, n_steps):
            raise ValueError("valid_steps must have shape (n_trials, n_steps)")

    lengths = np.linalg.norm(vector_array, axis=2)
    lengths = np.where(valid, lengths, 0.0)
    cumulative = np.cumsum(lengths, axis=1)
    previous = np.concatenate(
        (np.zeros((n_trials, 1), dtype=float), cumulative[:, :-1]), axis=1
    )
    full_endpoint = valid & (cumulative <= budget_array[:, None])
    effective_vectors = np.where(full_endpoint[..., None], vector_array, 0.0)
    endpoint_distance = np.where(full_endpoint, cumulative, np.nan)
    endpoint_mask = full_endpoint.copy()

    if overshoot_policy == "truncate":
        candidates = (
            valid
            & (lengths > 0.0)
            & (previous < budget_array[:, None])
            & (cumulative > budget_array[:, None])
        )
        first_candidate = candidates & (
            np.cumsum(candidates.astype(np.int8), axis=1) == 1
        )
        scale = np.divide(
            budget_array[:, None] - previous,
            lengths,
            out=np.zeros_like(lengths),
            where=first_candidate,
        )
        effective_vectors += np.where(
            first_candidate[..., None], vector_array * scale[..., None], 0.0
        )
        endpoint_mask |= first_candidate
        endpoint_distance = np.where(first_candidate, budget_array[:, None], endpoint_distance)

    cumulative_vectors = np.cumsum(effective_vectors, axis=1)
    endpoints = wrap_positions(
        start_array[:, None, :] + cumulative_vectors, side_length
    )
    return BudgetedPaths(
        endpoints=endpoints,
        endpoint_mask=endpoint_mask,
        endpoint_distance=endpoint_distance,
        starts=start_array,
        budget=budget_array,
        side_length=float(side_length),
    )


def replay_prepared(
    paths: BudgetedPaths,
    target: Target,
    check_initial_position: bool = False,
) -> ReplayBatchResult:
    """Apply one target membership rule to already prepared endpoints."""
    n_trials, n_steps, _ = paths.endpoints.shape
    if not np.isclose(target.side_length, paths.side_length):
        raise ValueError("target and prepared paths must share one side_length")
    hits = target.contains(paths.endpoints.reshape(-1, 3)).reshape(n_trials, n_steps)
    hits &= paths.endpoint_mask
    detected = np.any(hits, axis=1)
    first = np.argmax(hits, axis=1) if n_steps else np.zeros(n_trials, dtype=int)
    first = np.where(detected, first, -1)
    detection_distance = np.full(n_trials, np.nan, dtype=float)
    detected_rows = np.flatnonzero(detected)
    if detected_rows.size:
        detection_distance[detected_rows] = paths.endpoint_distance[
            detected_rows, first[detected_rows]
        ]

    if check_initial_position:
        initial_hit = target.contains(wrap_positions(paths.starts, target.side_length))
        detected = detected | initial_hit
        first = np.where(initial_hit, -1, first)
        detection_distance = np.where(initial_hit, 0.0, detection_distance)

    restricted = np.where(detected, detection_distance, paths.budget)
    return ReplayBatchResult(
        detected=detected,
        first_hit_step=first.astype(np.int64),
        detection_distance=detection_distance,
        restricted_detection_distance=restricted,
        endpoints_checked=np.sum(paths.endpoint_mask, axis=1).astype(np.int64),
        budget=paths.budget,
    )


def replay_batch(
    vectors: np.ndarray,
    starts: np.ndarray,
    target: Target,
    budgets: np.ndarray,
    valid_steps: Optional[np.ndarray] = None,
    overshoot_policy: str = "discard",
    check_initial_position: bool = False,
) -> ReplayBatchResult:
    """Replay a padded batch of vector sequences against one target."""
    paths = prepare_budgeted_paths(
        vectors,
        starts,
        budgets,
        target.side_length,
        valid_steps=valid_steps,
        overshoot_policy=overshoot_policy,
    )
    return replay_prepared(
        paths, target, check_initial_position=check_initial_position
    )
