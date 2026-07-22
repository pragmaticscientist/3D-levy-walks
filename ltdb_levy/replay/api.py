"""Convenience API over the vectorized replay engine."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np

from ..targets import Target
from .engine import (
    ReplayBatchResult,
    prepare_budgeted_paths,
    replay_batch,
    replay_prepared,
)


@dataclass(frozen=True)
class ReplayTrialResult:
    """Scalar view of a one-trial replay."""

    detected: bool
    first_hit_step: int
    detection_distance: float
    restricted_detection_distance: float
    endpoints_checked: int
    budget: float


def replay_trial(
    vectors: np.ndarray,
    start: np.ndarray,
    target: Target,
    budget: float,
    overshoot_policy: str = "discard",
    check_initial_position: bool = False,
) -> ReplayTrialResult:
    """Replay one unpadded vector sequence."""
    vector_array = np.asarray(vectors, dtype=float)
    result = replay_batch(
        vector_array[None, :, :],
        np.asarray(start, dtype=float)[None, :],
        target,
        np.asarray([budget], dtype=float),
        overshoot_policy=overshoot_policy,
        check_initial_position=check_initial_position,
    )
    return ReplayTrialResult(
        detected=bool(result.detected[0]),
        first_hit_step=int(result.first_hit_step[0]),
        detection_distance=float(result.detection_distance[0]),
        restricted_detection_distance=float(result.restricted_detection_distance[0]),
        endpoints_checked=int(result.endpoints_checked[0]),
        budget=float(result.budget[0]),
    )


def replay_targets(
    vectors: np.ndarray,
    start: np.ndarray,
    targets: Sequence[Target],
    budget: float,
    overshoot_policy: str = "discard",
    check_initial_position: bool = False,
) -> List[ReplayTrialResult]:
    """Replay one trajectory against many same-domain targets with endpoint reuse."""
    if not targets:
        return []
    side_length = targets[0].side_length
    if any(not np.isclose(target.side_length, side_length) for target in targets):
        raise ValueError("All targets must share one side_length")
    vector_array = np.asarray(vectors, dtype=float)
    paths = prepare_budgeted_paths(
        vector_array[None, :, :],
        np.asarray(start, dtype=float)[None, :],
        np.asarray([budget], dtype=float),
        side_length,
        overshoot_policy=overshoot_policy,
    )
    output: List[ReplayTrialResult] = []
    for target in targets:
        result = replay_prepared(
            paths, target, check_initial_position=check_initial_position
        )
        output.append(
            ReplayTrialResult(
                detected=bool(result.detected[0]),
                first_hit_step=int(result.first_hit_step[0]),
                detection_distance=float(result.detection_distance[0]),
                restricted_detection_distance=float(
                    result.restricted_detection_distance[0]
                ),
                endpoints_checked=int(result.endpoints_checked[0]),
                budget=float(result.budget[0]),
            )
        )
    return output


def replay_targets_batch(
    vectors: np.ndarray,
    starts: np.ndarray,
    targets: Sequence[Target],
    budgets: np.ndarray,
    valid_steps: np.ndarray = None,
    overshoot_policy: str = "discard",
    check_initial_position: bool = False,
) -> List[ReplayBatchResult]:
    """Replay a padded trajectory batch against many targets with endpoint reuse."""
    if not targets:
        return []
    side_length = targets[0].side_length
    if any(not np.isclose(target.side_length, side_length) for target in targets):
        raise ValueError("All targets must share one side_length")
    paths = prepare_budgeted_paths(
        vectors,
        starts,
        budgets,
        side_length,
        valid_steps=valid_steps,
        overshoot_policy=overshoot_policy,
    )
    return [
        replay_prepared(
            paths, target, check_initial_position=check_initial_position
        )
        for target in targets
    ]


__all__ = [
    "ReplayBatchResult",
    "ReplayTrialResult",
    "replay_batch",
    "replay_targets",
    "replay_targets_batch",
    "replay_trial",
]
