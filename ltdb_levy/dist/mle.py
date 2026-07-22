"""Maximum-likelihood fitting of the C-mixture exponent and crossover."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from ..errors import FitError
from .mixture import log_density, normalization_integral


@dataclass(frozen=True)
class MixtureFit:
    """Result of a mixture MLE with fixed or estimated crossover."""

    mu: float
    lmax: float
    crossover: float
    log_likelihood: float
    n_observations: int
    effective_sample_size: float
    optimizer_success: bool
    at_bound: bool
    crossover_estimated: bool = False
    crossover_at_bound: bool = False
    crossover_bound_low: float = np.nan
    crossover_bound_high: float = np.nan

    @property
    def n_parameters(self) -> int:
        """Number of smoothly fitted parameters, excluding the support maximum."""
        return 2 if self.crossover_estimated else 1


def select_lmax(
    lengths: np.ndarray,
    method: str = "empirical_max",
    user_lmax: Optional[float] = None,
    quantile: float = 0.99,
) -> float:
    """Resolve the imposed upper support from data and configuration."""
    values = np.asarray(lengths, dtype=float)
    values = values[np.isfinite(values) & (values >= 0.0)]
    if values.size == 0:
        raise FitError("Cannot select lmax from an empty length sample")
    if method == "empirical_max":
        result = float(np.max(values))
    elif method == "user":
        if user_lmax is None:
            raise FitError("user_lmax is required for lmax_method='user'")
        result = float(user_lmax)
    elif method == "quantile":
        if not 0.0 < quantile <= 1.0:
            raise FitError("lmax quantile must lie in (0, 1]")
        result = float(np.quantile(values, quantile))
    else:
        raise FitError("Unknown lmax method {!r}".format(method))
    if not np.isfinite(result) or result <= 0.0:
        raise FitError("Selected lmax must be finite and positive, got {}".format(result))
    return result


def _prepare_weights(lengths: np.ndarray, weights: Optional[np.ndarray]) -> np.ndarray:
    if weights is None:
        return np.ones(len(lengths), dtype=float)
    out = np.asarray(weights, dtype=float)
    if out.shape != lengths.shape:
        raise FitError("weights must have the same shape as lengths")
    if np.any(~np.isfinite(out)) or np.any(out < 0.0) or out.sum() <= 0.0:
        raise FitError("weights must be finite, non-negative, and have positive sum")
    # Preserve the conventional n-observation log-likelihood scale.
    return out * (len(out) / out.sum())


def _resolve_crossover_bounds(
    values: np.ndarray,
    lmax: float,
    crossover_bounds: Optional[Sequence[float]],
) -> Sequence[float]:
    """Return finite data-compatible bounds for estimating the crossover."""
    positive = values[values > 0.0]
    if positive.size == 0:
        raise FitError("Crossover estimation requires at least one positive length")
    # The likelihood vanishes as c -> 0 for mu > 1. A lower guard at one tenth
    # of the smallest positive observation is therefore numerical, not a
    # biologically imposed scale. The upper guard keeps c strictly below L,
    # where mu would become unidentifiable under the uniform branch.
    lower = max(float(np.min(positive)) * 0.1, float(lmax) * 1e-12)
    upper = float(lmax) * (1.0 - 1e-9)
    if crossover_bounds is not None:
        if len(crossover_bounds) != 2:
            raise FitError("crossover_bounds must be a two-element increasing pair")
        requested_low = float(crossover_bounds[0])
        requested_high = float(crossover_bounds[1])
        if (
            not np.isfinite(requested_low)
            or not np.isfinite(requested_high)
            or requested_low <= 0.0
            or requested_low >= requested_high
        ):
            raise FitError("crossover_bounds must be finite, positive, and increasing")
        lower = max(lower, requested_low)
        upper = min(upper, requested_high)
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise FitError(
            "No identifiable crossover interval remains below lmax={}".format(lmax)
        )
    return (lower, upper)


def _fit_joint_mixture(
    values: np.ndarray,
    fit_weights: np.ndarray,
    lmax: float,
    initial_crossover: float,
    mu_bounds: Sequence[float],
    crossover_bounds: Optional[Sequence[float]],
) -> Sequence[float]:
    """Jointly maximize the weighted likelihood over ``(mu, crossover)``.

    Sorting once reduces each objective evaluation to a binary search and a
    handful of scalar operations. This is important because the same fit is
    repeated thousands of times in clustered and parametric bootstraps.
    """
    lower_c, upper_c = _resolve_crossover_bounds(
        values, lmax, crossover_bounds
    )
    lower_mu, upper_mu = float(mu_bounds[0]), float(mu_bounds[1])
    order = np.argsort(values, kind="mergesort")
    ordered = values[order]
    ordered_weights = fit_weights[order]
    safe_logs = np.zeros_like(ordered)
    positive = ordered > 0.0
    safe_logs[positive] = np.log(ordered[positive])
    tail_weights = np.cumsum(ordered_weights[::-1])[::-1]
    tail_weighted_logs = np.cumsum(
        (ordered_weights * safe_logs)[::-1]
    )[::-1]
    total_weight = float(np.sum(fit_weights))

    def objective(theta: np.ndarray) -> float:
        mu = float(theta[0])
        log_c = float(theta[1])
        if (
            not lower_mu <= mu <= upper_mu
            or not np.log(lower_c) <= log_c <= np.log(upper_c)
        ):
            return float("inf")
        crossover = float(np.exp(log_c))
        try:
            integral = normalization_integral(mu, lmax, crossover)
        except (FloatingPointError, OverflowError, ValueError):
            return float("inf")
        if not np.isfinite(integral) or integral <= 0.0:
            return float("inf")
        index = int(np.searchsorted(ordered, crossover, side="left"))
        if index < len(ordered):
            penalty = float(
                tail_weighted_logs[index] - log_c * tail_weights[index]
            )
        else:
            penalty = 0.0
        result = total_weight * np.log(integral) + mu * penalty
        return float(result) if np.isfinite(result) else float("inf")

    log_lower = float(np.log(lower_c))
    log_upper = float(np.log(upper_c))
    # A compact set of starts keeps bootstrap refitting tractable. The
    # cell-radius reference, data quartiles, and interval endpoints cover both
    # biologically plausible and boundary-seeking solutions.
    crossover_starts = [lower_c, upper_c]
    if np.isfinite(initial_crossover):
        crossover_starts.append(
            float(np.clip(initial_crossover, lower_c, upper_c))
        )
    positive_values = values[values > 0.0]
    crossover_starts.extend(
        float(np.clip(value, lower_c, upper_c))
        for value in np.quantile(positive_values, [0.25, 0.5, 0.75])
    )

    coarse = []
    for crossover in sorted(set(crossover_starts)):
        log_c = float(np.log(crossover))
        conditional = minimize_scalar(
            lambda mu: objective(np.array([mu, log_c])),
            bounds=(lower_mu, upper_mu),
            method="bounded",
            options={"xatol": 1e-8, "maxiter": 200},
        )
        if conditional.success and np.isfinite(conditional.fun):
            coarse.append(
                (float(conditional.fun), float(conditional.x), log_c)
            )
    if not coarse:
        raise FitError("Joint crossover/exponent search found no finite start")
    coarse.sort(key=lambda item: item[0])

    candidates = list(coarse)
    bounds = [(lower_mu, upper_mu), (log_lower, log_upper)]
    for _, start_mu, start_log_c in coarse[:1]:
        result = minimize(
            objective,
            np.array([start_mu, start_log_c]),
            method="Powell",
            bounds=bounds,
            options={
                "xtol": 1e-8,
                "ftol": 1e-10,
                "maxiter": 300,
            },
        )
        if np.isfinite(result.fun):
            candidates.append(
                (float(result.fun), float(result.x[0]), float(result.x[1]))
            )

    # The likelihood derivative changes where c crosses an observation. Check
    # the nearest such kinks explicitly after the continuous search.
    current_best = min(candidates, key=lambda item: item[0])
    best_c = float(np.exp(current_best[2]))
    unique_values = np.unique(
        positive_values[
            (positive_values >= lower_c) & (positive_values <= upper_c)
        ]
    )
    if unique_values.size:
        position = int(np.searchsorted(unique_values, best_c))
        nearby = unique_values[
            max(0, position - 1) : min(len(unique_values), position + 2)
        ]
        for crossover in nearby:
            log_c = float(np.log(crossover))
            conditional = minimize_scalar(
                lambda mu: objective(np.array([mu, log_c])),
                bounds=(lower_mu, upper_mu),
                method="bounded",
                options={"xatol": 1e-9, "maxiter": 200},
            )
            if conditional.success and np.isfinite(conditional.fun):
                candidates.append(
                    (float(conditional.fun), float(conditional.x), log_c)
                )

    best = min(candidates, key=lambda item: item[0])
    return (
        float(best[1]),
        float(np.exp(best[2])),
        float(-best[0]),
        float(lower_c),
        float(upper_c),
    )


def fit_mixture(
    lengths: np.ndarray,
    lmax: float,
    crossover: float = 1.0,
    mu_bounds: Sequence[float] = (1.001, 10.0),
    weights: Optional[np.ndarray] = None,
    estimate_crossover: bool = False,
    crossover_bounds: Optional[Sequence[float]] = None,
) -> MixtureFit:
    """Fit ``mu`` and optionally one crossover for the supplied condition.

    Observations above ``lmax`` are an error. A quantile-based support therefore
    requires the caller to explicitly exclude and audit values above that
    quantile; silently clipping them would create mass at the boundary.
    """
    values = np.asarray(lengths, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise FitError("Mixture fitting requires a non-empty one-dimensional sample")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise FitError("Relocation lengths must be finite and non-negative")
    if not estimate_crossover and lmax <= crossover:
        raise FitError(
            "lmax ({}) must exceed crossover ({}); otherwise mu is unidentifiable".format(
                lmax, crossover
            )
        )
    if np.any(values > lmax * (1.0 + 1e-12)):
        raise FitError("One or more relocation lengths exceed the imposed lmax")
    if len(mu_bounds) != 2 or mu_bounds[0] <= 0.0 or mu_bounds[0] >= mu_bounds[1]:
        raise FitError("mu_bounds must be a positive increasing pair")
    fit_weights = _prepare_weights(values, weights)

    if estimate_crossover:
        mu, fitted_crossover, log_likelihood, lower_c, upper_c = (
            _fit_joint_mixture(
                values,
                fit_weights,
                lmax,
                crossover,
                mu_bounds,
                crossover_bounds,
            )
        )
        tolerance = max(
            1e-6,
            1e-5 * (float(mu_bounds[1]) - float(mu_bounds[0])),
        )
        mu_at_bound = bool(
            mu - float(mu_bounds[0]) <= tolerance
            or float(mu_bounds[1]) - mu <= tolerance
        )
        log_span = float(np.log(upper_c / lower_c))
        crossover_tolerance = max(1e-7, 1e-5 * log_span)
        crossover_at_bound = bool(
            np.log(fitted_crossover / lower_c) <= crossover_tolerance
            or np.log(upper_c / fitted_crossover) <= crossover_tolerance
        )
        probability_weights = fit_weights / fit_weights.sum()
        n_eff = float(1.0 / np.sum(np.square(probability_weights)))
        return MixtureFit(
            mu=mu,
            lmax=float(lmax),
            crossover=fitted_crossover,
            log_likelihood=log_likelihood,
            n_observations=int(len(values)),
            effective_sample_size=n_eff,
            optimizer_success=True,
            at_bound=mu_at_bound,
            crossover_estimated=True,
            crossover_at_bound=crossover_at_bound,
            crossover_bound_low=lower_c,
            crossover_bound_high=upper_c,
        )

    def objective(mu: float) -> float:
        log_pdf = log_density(values, mu, lmax, crossover)
        return -float(np.dot(fit_weights, log_pdf))

    result = minimize_scalar(
        objective,
        bounds=(float(mu_bounds[0]), float(mu_bounds[1])),
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 500},
    )
    if not result.success or not np.isfinite(result.fun):
        raise FitError("Mixture MLE failed: {}".format(result.message))
    mu = float(result.x)
    tolerance = max(1e-6, 1e-5 * (float(mu_bounds[1]) - float(mu_bounds[0])))
    at_bound = bool(
        mu - float(mu_bounds[0]) <= tolerance
        or float(mu_bounds[1]) - mu <= tolerance
    )
    probability_weights = fit_weights / fit_weights.sum()
    n_eff = float(1.0 / np.sum(np.square(probability_weights)))
    return MixtureFit(
        mu=mu,
        lmax=float(lmax),
        crossover=float(crossover),
        log_likelihood=float(-result.fun),
        n_observations=int(len(values)),
        effective_sample_size=n_eff,
        optimizer_success=True,
        at_bound=at_bound,
        crossover_estimated=False,
        crossover_at_bound=False,
        crossover_bound_low=float(crossover),
        crossover_bound_high=float(crossover),
    )
