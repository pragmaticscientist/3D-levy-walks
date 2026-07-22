"""Track-clustered bootstrap utilities.

Monte-Carlo replays and runs from one tracked cell are not independent.
Resampling therefore happens at the whole-track level throughout this package.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .dist.mle import fit_mixture, select_lmax
from .errors import FitError


@dataclass(frozen=True)
class ClusteredMuBootstrap:
    """Clustered-bootstrap draws and percentile intervals for the mixture."""

    draws: np.ndarray
    crossover_draws: np.ndarray
    ci_low: float
    ci_high: float
    crossover_ci_low: float
    crossover_ci_high: float
    n_successful: int
    n_failed: int


def clustered_mu_bootstrap(
    runs: pd.DataFrame,
    n_bootstraps: int,
    rng: np.random.Generator,
    crossover: float,
    mu_bounds: Sequence[float],
    weighting: str = "cell_balanced",
    estimate_crossover: bool = False,
    crossover_bounds: Optional[Sequence[float]] = None,
    lmax_method: str = "empirical_max",
    user_lmax: Optional[float] = None,
    lmax_quantile: float = 0.99,
) -> ClusteredMuBootstrap:
    """Resample complete tracks and refit all configured mixture parameters."""
    if n_bootstraps < 0:
        raise ValueError("n_bootstraps must be non-negative")
    if weighting not in ("cell_balanced", "step_weighted"):
        raise ValueError("Unknown weighting {!r}".format(weighting))
    track_ids = runs["track_uid"].drop_duplicates().to_numpy()
    if track_ids.size == 0:
        raise ValueError("Cannot bootstrap an empty run table")
    groups = {
        track_id: group["run_length_um"].to_numpy(dtype=float)
        for track_id, group in runs.groupby("track_uid", sort=False)
    }
    draws = np.full(n_bootstraps, np.nan, dtype=float)
    crossover_draws = np.full(n_bootstraps, np.nan, dtype=float)
    for replicate_index in range(n_bootstraps):
        chosen = rng.choice(track_ids, size=len(track_ids), replace=True)
        pieces = [groups[track_id] for track_id in chosen]
        lengths = np.concatenate(pieces)
        if weighting == "cell_balanced":
            # Each pseudo-track gets equal total mass, including repeated draws
            # of the same biological track.
            weights = np.concatenate(
                [
                    np.full(
                        len(piece),
                        1.0 / (len(chosen) * len(piece)),
                        dtype=float,
                    )
                    for piece in pieces
                ]
            )
        else:
            # Tracks remain the resampling clusters, while every run in the
            # resulting pseudo-sample contributes the same likelihood weight.
            weights = np.full(len(lengths), 1.0 / len(lengths), dtype=float)
        try:
            lmax = select_lmax(
                lengths,
                method=lmax_method,
                user_lmax=user_lmax,
                quantile=lmax_quantile,
            )
            keep = lengths <= lmax
            fitted = fit_mixture(
                lengths[keep],
                lmax=lmax,
                crossover=crossover,
                mu_bounds=mu_bounds,
                weights=weights[keep],
                estimate_crossover=estimate_crossover,
                crossover_bounds=crossover_bounds,
            )
            draws[replicate_index] = fitted.mu
            crossover_draws[replicate_index] = fitted.crossover
        except (FitError, ValueError, RuntimeError):
            # Failed replicates are reported, never silently replaced.
            continue
    successful = draws[np.isfinite(draws)]
    if successful.size:
        ci_low, ci_high = np.quantile(successful, [0.025, 0.975])
    else:
        ci_low, ci_high = np.nan, np.nan
    successful_crossovers = crossover_draws[np.isfinite(crossover_draws)]
    if successful_crossovers.size:
        crossover_ci_low, crossover_ci_high = np.quantile(
            successful_crossovers, [0.025, 0.975]
        )
    else:
        crossover_ci_low, crossover_ci_high = np.nan, np.nan
    return ClusteredMuBootstrap(
        draws=draws,
        crossover_draws=crossover_draws,
        ci_low=float(ci_low),
        ci_high=float(ci_high),
        crossover_ci_low=float(crossover_ci_low),
        crossover_ci_high=float(crossover_ci_high),
        n_successful=int(successful.size),
        n_failed=int(n_bootstraps - successful.size),
    )
