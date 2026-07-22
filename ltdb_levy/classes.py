"""Trajectory-class controls used by the surrogate ladder."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np

from .directions import (
    independently_rotate_vectors,
    isotropic_directions,
    random_quaternions,
    rotate_vectors,
)
from .dist.mixture import normalization_integral, sample as sample_mixture

TRAJECTORY_CLASSES = (
    "empirical",
    "block_resampled",
    "ordered_length_rotated",
    "pool_isotropic",
    "fitted_mu",
    "mu2",
)


@dataclass(frozen=True)
class CommonRandomNumbers:
    """Random arrays reusable across class constructions for paired contrasts."""

    uniforms: np.ndarray
    directions: np.ndarray
    quaternions: np.ndarray


@dataclass(frozen=True)
class BlockPool:
    """Precomputed fixed-length empirical blocks for one condition."""

    blocks: np.ndarray
    mean_path_length: float


def make_common_random_numbers(
    size: int, rng: np.random.Generator
) -> CommonRandomNumbers:
    """Pre-draw class-independent randomness for one paired trial index."""
    return CommonRandomNumbers(
        uniforms=rng.random(size),
        directions=isotropic_directions(size, rng),
        quaternions=random_quaternions(size, rng),
    )


def mixture_mean(mu: float, lmax: float, crossover: float = 1.0) -> float:
    """Analytic mean relocation length under the C mixture."""
    z = normalization_integral(mu, lmax, crossover)
    if lmax <= crossover:
        return lmax / 2.0
    q = 2.0 - mu
    log_ratio = np.log(lmax / crossover)
    tail_factor = (
        log_ratio
        if q == 0.0
        else float(np.expm1(q * log_ratio) / q)
    )
    first_moment = crossover ** 2 * (0.5 + tail_factor)
    return float(first_moment / z)


def _draw_count(budget: float, mean_length: float) -> int:
    if budget < 0.0 or mean_length <= 0.0:
        raise ValueError("budget must be non-negative and mean length positive")
    # Four expected budgets plus a fixed margin makes exhaustion negligible for
    # the bounded distributions used here, while keeping arrays compact.
    return max(64, int(np.ceil(4.0 * budget / mean_length)) + 16)


def _ensure_budget_covered(vectors: np.ndarray, budget: float, class_name: str) -> None:
    if np.linalg.norm(vectors, axis=1).sum() + 1e-12 < budget:
        raise RuntimeError(
            "{} candidate sequence did not cover budget {}; increase the draw margin".format(
                class_name, budget
            )
        )


def empirical(source_vectors: np.ndarray) -> np.ndarray:
    """Return a copy of the observed ordered vector sequence."""
    return np.asarray(source_vectors, dtype=float).copy()


def ordered_length_rotated(
    source_vectors: np.ndarray,
    rng: np.random.Generator,
    common: Optional[CommonRandomNumbers] = None,
) -> np.ndarray:
    """Preserve length order while independently randomizing every direction."""
    vectors = np.asarray(source_vectors, dtype=float)
    quaternions = (
        common.quaternions[: len(vectors)]
        if common is not None and len(common.quaternions) >= len(vectors)
        else random_quaternions(len(vectors), rng)
    )
    return rotate_vectors(vectors, quaternions)


def pool_isotropic(
    pool_lengths: np.ndarray,
    pool_weights: np.ndarray,
    budget: float,
    rng: np.random.Generator,
    common: Optional[CommonRandomNumbers] = None,
) -> np.ndarray:
    """Destroy ordering by independently drawing pool lengths and directions."""
    lengths = np.asarray(pool_lengths, dtype=float)
    weights = np.asarray(pool_weights, dtype=float)
    probabilities = weights / weights.sum()
    mean = float(np.dot(probabilities, lengths))
    count = _draw_count(budget, mean)
    uniforms = (
        common.uniforms[:count]
        if common is not None and len(common.uniforms) >= count
        else rng.random(count)
    )
    cumulative = np.cumsum(probabilities)
    indices = np.searchsorted(cumulative, uniforms, side="right")
    indices = np.minimum(indices, len(lengths) - 1)
    sampled_lengths = lengths[indices]
    directions = (
        common.directions[:count]
        if common is not None and len(common.directions) >= count
        else isotropic_directions(count, rng)
    )
    vectors = sampled_lengths[:, None] * directions
    _ensure_budget_covered(vectors, budget, "pool_isotropic")
    return vectors


def fitted_isotropic(
    mu: float,
    lmax: float,
    budget: float,
    rng: np.random.Generator,
    crossover: float = 1.0,
    common: Optional[CommonRandomNumbers] = None,
) -> np.ndarray:
    """Draw mixture lengths and isotropic directions for fitted-μ or μ=2."""
    count = _draw_count(budget, mixture_mean(mu, lmax, crossover))
    uniforms = (
        common.uniforms[:count]
        if common is not None and len(common.uniforms) >= count
        else rng.random(count)
    )
    lengths = sample_mixture(
        mu, lmax, count, rng, crossover=crossover, uniforms=uniforms
    )
    directions = (
        common.directions[:count]
        if common is not None and len(common.directions) >= count
        else isotropic_directions(count, rng)
    )
    vectors = lengths[:, None] * directions
    _ensure_budget_covered(vectors, budget, "fitted_isotropic")
    return vectors


def build_block_pool(
    track_sequences: Sequence[np.ndarray],
    block_length: int,
) -> BlockPool:
    """Precompute every contiguous block once for repeated replay trials."""
    if block_length < 1:
        raise ValueError("block_length must be >= 1")
    candidates: List[np.ndarray] = []
    for sequence in track_sequences:
        vectors = np.asarray(sequence, dtype=float)
        if len(vectors) < block_length:
            continue
        for start in range(len(vectors) - block_length + 1):
            candidates.append(vectors[start : start + block_length])
    if not candidates:
        raise ValueError("No source sequence contains a complete block")
    blocks = np.stack(candidates, axis=0)
    block_paths = np.sum(np.linalg.norm(blocks, axis=2), axis=1)
    return BlockPool(blocks=blocks, mean_path_length=float(np.mean(block_paths)))


def block_resampled(
    track_sequences: Sequence[np.ndarray],
    block_length: int,
    budget: float,
    rng: np.random.Generator,
    block_pool: Optional[BlockPool] = None,
) -> np.ndarray:
    """Resample contiguous empirical blocks, applying one rotation per block."""
    pool = (
        block_pool
        if block_pool is not None
        else build_block_pool(track_sequences, block_length)
    )
    count = _draw_count(budget, pool.mean_path_length)
    indices = rng.integers(0, len(pool.blocks), size=count)
    quaternions = random_quaternions(count, rng)
    blocks = [
        rotate_vectors(pool.blocks[index], quaternions[i])
        for i, index in enumerate(indices)
    ]
    vectors = np.concatenate(blocks, axis=0)
    _ensure_budget_covered(vectors, budget, "block_resampled")
    return vectors


def build_trajectory_class(
    class_name: str,
    source_vectors: np.ndarray,
    budget: float,
    rng: np.random.Generator,
    pool_lengths: Optional[np.ndarray] = None,
    pool_weights: Optional[np.ndarray] = None,
    fitted_mu: Optional[float] = None,
    lmax: Optional[float] = None,
    crossover: float = 1.0,
    track_sequences: Optional[Sequence[np.ndarray]] = None,
    block_pool: Optional[BlockPool] = None,
    block_length: int = 3,
    common: Optional[CommonRandomNumbers] = None,
) -> np.ndarray:
    """Construct one rung of the control ladder with explicit dependencies."""
    if class_name == "empirical":
        return empirical(source_vectors)
    if class_name == "ordered_length_rotated":
        return ordered_length_rotated(source_vectors, rng, common)
    if class_name == "pool_isotropic":
        if pool_lengths is None or pool_weights is None:
            raise ValueError("pool_isotropic requires pool lengths and weights")
        return pool_isotropic(pool_lengths, pool_weights, budget, rng, common)
    if class_name in ("fitted_mu", "mu2"):
        if lmax is None or (class_name == "fitted_mu" and fitted_mu is None):
            raise ValueError("{} requires fitted support parameters".format(class_name))
        mu = 2.0 if class_name == "mu2" else float(fitted_mu)
        return fitted_isotropic(mu, lmax, budget, rng, crossover, common)
    if class_name == "block_resampled":
        if track_sequences is None:
            raise ValueError("block_resampled requires empirical track sequences")
        return block_resampled(
            track_sequences,
            block_length,
            budget,
            rng,
            block_pool=block_pool,
        )
    raise ValueError("Unknown trajectory class {!r}".format(class_name))
