"""Goodness-of-fit diagnostics for the fitted C mixture."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..errors import FitError
from .families import FamilyFit, sample_family
from .mixture import cdf, sample
from .mle import MixtureFit, fit_mixture


@dataclass(frozen=True)
class BootstrapKS:
    """Parametric-bootstrap Kolmogorov--Smirnov result."""

    statistic: float
    p_value: float
    n_bootstraps: int
    bootstrap_statistics: np.ndarray


def ks_statistic(
    lengths: np.ndarray,
    mu: float,
    lmax: float,
    crossover: float = 1.0,
    weights: Optional[np.ndarray] = None,
) -> float:
    """Return the two-sided (optionally weighted) empirical-CDF KS distance."""
    values = np.asarray(lengths, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ordered = values[order]
    if weights is None:
        mass = np.full(len(values), 1.0 / len(values))
    else:
        mass = np.asarray(weights, dtype=float)[order]
        mass = mass / mass.sum()
    upper_empirical = np.cumsum(mass)
    lower_empirical = upper_empirical - mass
    model = cdf(ordered, mu, lmax, crossover)
    return float(
        max(
            np.max(np.abs(upper_empirical - model)),
            np.max(np.abs(lower_empirical - model)),
        )
    )


def parametric_bootstrap_ks(
    lengths: np.ndarray,
    fitted: MixtureFit,
    n_bootstraps: int,
    rng: np.random.Generator,
    mu_bounds: Sequence[float] = (1.001, 10.0),
    weights: Optional[np.ndarray] = None,
    estimate_crossover: bool = False,
    crossover_bounds: Optional[Sequence[float]] = None,
) -> BootstrapKS:
    """Bootstrap the KS null, refitting every configured model parameter."""
    values = np.asarray(lengths, dtype=float)
    observed = ks_statistic(
        values, fitted.mu, fitted.lmax, fitted.crossover, weights=weights
    )
    statistics = np.empty(n_bootstraps, dtype=float)
    for i in range(n_bootstraps):
        simulated = sample(
            fitted.mu, fitted.lmax, len(values), rng, crossover=fitted.crossover
        )
        replicate = fit_mixture(
            simulated,
            lmax=fitted.lmax,
            crossover=fitted.crossover,
            mu_bounds=mu_bounds,
            weights=weights,
            estimate_crossover=estimate_crossover,
            crossover_bounds=crossover_bounds,
        )
        statistics[i] = ks_statistic(
            simulated,
            replicate.mu,
            replicate.lmax,
            replicate.crossover,
            weights=weights,
        )
    p_value = (
        float((1 + np.count_nonzero(statistics >= observed)) / (n_bootstraps + 1))
        if n_bootstraps
        else np.nan
    )
    return BootstrapKS(
        statistic=observed,
        p_value=p_value,
        n_bootstraps=int(n_bootstraps),
        bootstrap_statistics=statistics,
    )


def ks_power_curve(
    alternative: FamilyFit,
    sample_sizes: Sequence[int],
    critical_value: float,
    n_replicates: int,
    rng: np.random.Generator,
    crossover: float,
    mu_bounds: Sequence[float],
    estimate_crossover: bool = False,
    crossover_bounds: Optional[Sequence[float]] = None,
) -> pd.DataFrame:
    """Estimate KS rejection power against one fitted alternative family.

    The null model is refitted in every simulated dataset, exactly as in the
    parametric GoF bootstrap. Biological track clustering is not represented
    here; this curve diagnoses statistical power of the length-distribution
    check itself.
    """
    rows: List[Dict[str, object]] = []
    for size in sorted(set(int(x) for x in sample_sizes if int(x) > 1)):
        rejected = 0
        successful = 0
        for _ in range(n_replicates):
            values = sample_family(alternative, size, rng)
            try:
                fitted = fit_mixture(
                    values,
                    lmax=alternative.lmax,
                    crossover=crossover,
                    mu_bounds=mu_bounds,
                    estimate_crossover=estimate_crossover,
                    crossover_bounds=crossover_bounds,
                )
            except (FitError, ValueError, RuntimeError):
                continue
            statistic = ks_statistic(
                values, fitted.mu, fitted.lmax, fitted.crossover
            )
            successful += 1
            rejected += int(statistic > critical_value)
        rows.append(
            {
                "alternative": alternative.name,
                "sample_size": size,
                "power": float(rejected / successful) if successful else np.nan,
                "replicates_successful": successful,
                "replicates_requested": int(n_replicates),
                "alpha": 0.05,
                "critical_value": float(critical_value),
            }
        )
    return pd.DataFrame(rows)
