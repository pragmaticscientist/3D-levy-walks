"""Alternative truncated relocation-length families used for model comparison.

All models use the same closed support ``[0, lmax]`` as the C mixture.  The
power-law-with-cutoff alternative keeps the mixture's flat sub-crossover branch
and adds an exponential cutoff to its tail; this makes it well-defined at zero
without introducing an arbitrary positive lower cutoff.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import ndtr, ndtri

from ..errors import FitError

FAMILY_NAMES = (
    "truncated_exponential",
    "truncated_lognormal",
    "truncated_weibull",
    "power_law_cutoff",
)


@dataclass(frozen=True)
class FamilyFit:
    """A fitted alternative family."""

    name: str
    parameters: Dict[str, float]
    lmax: float
    crossover: float
    log_likelihood: float
    n_parameters: int
    optimizer_success: bool


def _normalised_weights(values: np.ndarray, weights: Optional[np.ndarray]) -> np.ndarray:
    if weights is None:
        return np.ones(len(values), dtype=float)
    result = np.asarray(weights, dtype=float)
    if result.shape != values.shape:
        raise FitError("weights must have the same shape as lengths")
    if np.any(~np.isfinite(result)) or np.any(result < 0.0) or result.sum() <= 0.0:
        raise FitError("weights must be finite, non-negative, and have positive sum")
    return result * (len(result) / result.sum())


def _validate_values(values: np.ndarray, lmax: float) -> None:
    if values.ndim != 1 or not len(values):
        raise FitError("Family fitting needs a non-empty one-dimensional sample")
    if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values > lmax):
        raise FitError("Lengths must be finite and contained in [0, lmax]")


def _cutoff_integral(alpha: float, rate: float, lmax: float, crossover: float) -> float:
    flat = min(crossover, lmax)
    if lmax <= crossover:
        return float(lmax)
    tail, _ = quad(
        lambda value: (value / crossover) ** (-alpha)
        * np.exp(-rate * (value - crossover)),
        crossover,
        lmax,
        epsabs=1e-9,
        epsrel=1e-9,
        limit=100,
    )
    return float(flat + tail)


def log_density(
    name: str,
    lengths: np.ndarray,
    parameters: Dict[str, float],
    lmax: float,
    crossover: float = 1.0,
) -> np.ndarray:
    """Evaluate a fitted alternative's log density."""
    values = np.asarray(lengths, dtype=float)
    out = np.full(values.shape, -np.inf, dtype=float)
    inside = np.isfinite(values) & (values >= 0.0) & (values <= lmax)
    x = values[inside]
    if name == "truncated_exponential":
        rate = parameters["rate"]
        log_z = np.log(-np.expm1(-rate * lmax))
        out[inside] = np.log(rate) - rate * x - log_z
    elif name == "truncated_lognormal":
        positive = inside & (values > 0.0)
        x = values[positive]
        location = parameters["log_location"]
        scale = parameters["log_scale"]
        z = ndtr((np.log(lmax) - location) / scale)
        out[positive] = (
            -np.log(x)
            - np.log(scale)
            - 0.5 * np.log(2.0 * np.pi)
            - 0.5 * ((np.log(x) - location) / scale) ** 2
            - np.log(z)
        )
    elif name == "truncated_weibull":
        positive = inside & (values > 0.0)
        x = values[positive]
        shape = parameters["shape"]
        scale = parameters["scale"]
        log_upper_power = shape * np.log(lmax / scale)
        if log_upper_power < -20.0:
            # log(1 - exp(-t)) = log(t) + O(t) for tiny t.
            log_z = log_upper_power
        elif log_upper_power > 20.0:
            log_z = 0.0
        else:
            upper_power = float(np.exp(log_upper_power))
            log_z = float(np.log(-np.expm1(-upper_power)))
        log_x_power = shape * np.log(x / scale)
        x_power = np.exp(np.minimum(log_x_power, 700.0))
        out[positive] = (
            np.log(shape)
            - np.log(scale)
            + (shape - 1.0) * np.log(x / scale)
            - x_power
            - log_z
        )
    elif name == "power_law_cutoff":
        alpha = parameters["alpha"]
        rate = parameters["cutoff_rate"]
        z = _cutoff_integral(alpha, rate, lmax, crossover)
        out[inside] = -np.log(z)
        tail = inside & (values >= crossover)
        out[tail] += (
            -alpha * np.log(values[tail] / crossover)
            - rate * (values[tail] - crossover)
        )
    else:
        raise ValueError("Unknown family {!r}".format(name))
    return out


def fit_family(
    name: str,
    lengths: np.ndarray,
    lmax: float,
    weights: Optional[np.ndarray] = None,
    crossover: float = 1.0,
) -> FamilyFit:
    """Maximum-likelihood fit one alternative family."""
    values = np.asarray(lengths, dtype=float)
    _validate_values(values, lmax)
    fit_weights = _normalised_weights(values, weights)
    positive = values[values > 0.0]
    if name in ("truncated_lognormal", "truncated_weibull") and len(positive) != len(values):
        raise FitError("{} cannot assign finite density to exact zero lengths".format(name))

    if name == "truncated_exponential":
        initial = np.array([np.log(1.0 / max(float(np.average(values, weights=fit_weights)), 1e-6))])
        bounds = [(-20.0, 20.0)]

        def unpack(theta: np.ndarray) -> Dict[str, float]:
            return {"rate": float(np.exp(theta[0]))}

    elif name == "truncated_lognormal":
        logs = np.log(values)
        location = float(np.average(logs, weights=fit_weights))
        variance = float(np.average((logs - location) ** 2, weights=fit_weights))
        initial = np.array([location, np.log(max(np.sqrt(variance), 0.1))])
        bounds = [(-30.0, 30.0), (-8.0, 5.0)]

        def unpack(theta: np.ndarray) -> Dict[str, float]:
            return {
                "log_location": float(theta[0]),
                "log_scale": float(np.exp(theta[1])),
            }

    elif name == "truncated_weibull":
        initial = np.array(
            [0.0, np.log(max(float(np.average(values, weights=fit_weights)), 1e-6))]
        )
        bounds = [(-6.0, 6.0), (-20.0, 20.0)]

        def unpack(theta: np.ndarray) -> Dict[str, float]:
            return {"shape": float(np.exp(theta[0])), "scale": float(np.exp(theta[1]))}

    elif name == "power_law_cutoff":
        if lmax <= crossover:
            raise FitError("power_law_cutoff needs lmax > crossover")
        initial = np.log(np.array([1.5, 0.01]))
        bounds = [(-8.0, np.log(50.0)), (-16.0, np.log(50.0))]

        def unpack(theta: np.ndarray) -> Dict[str, float]:
            return {
                "alpha": float(np.exp(theta[0])),
                "cutoff_rate": float(np.exp(theta[1])),
            }

    else:
        raise FitError("Unknown alternative family {!r}".format(name))

    def objective(theta: np.ndarray) -> float:
        parameters = unpack(theta)
        try:
            logs = log_density(name, values, parameters, lmax, crossover)
        except (FloatingPointError, OverflowError, ValueError):
            return float("inf")
        if np.any(~np.isfinite(logs)):
            return float("inf")
        return -float(np.dot(fit_weights, logs))

    result = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 500, "ftol": 1e-11},
    )
    if not result.success or not np.isfinite(result.fun):
        raise FitError("{} fit failed: {}".format(name, result.message))
    parameters = unpack(result.x)
    return FamilyFit(
        name=name,
        parameters=parameters,
        lmax=float(lmax),
        crossover=float(crossover),
        log_likelihood=float(-result.fun),
        n_parameters=len(result.x),
        optimizer_success=True,
    )


def sample_family(
    fitted: FamilyFit, size: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw from a fitted truncated alternative by inverse transform."""
    if size < 0:
        raise ValueError("size must be non-negative")
    u = rng.random(size)
    lmax = fitted.lmax
    p = fitted.parameters
    if fitted.name == "truncated_exponential":
        rate = p["rate"]
        return -np.log1p(-u * (-np.expm1(-rate * lmax))) / rate
    if fitted.name == "truncated_lognormal":
        location = p["log_location"]
        scale = p["log_scale"]
        upper = ndtr((np.log(lmax) - location) / scale)
        return np.exp(location + scale * ndtri(u * upper))
    if fitted.name == "truncated_weibull":
        shape = p["shape"]
        scale = p["scale"]
        log_upper_power = shape * np.log(lmax / scale)
        if log_upper_power > 20.0:
            upper = 1.0
        else:
            upper = -np.expm1(-np.exp(log_upper_power))
        return scale * (-np.log1p(-u * upper)) ** (1.0 / shape)
    if fitted.name == "power_law_cutoff":
        # A dense monotone numerical CDF is adequate here: this sampler is used
        # only for a diagnostic power curve, never for trajectory replay.
        grid = np.linspace(0.0, lmax, 20_001)
        log_pdf = log_density(
            fitted.name, grid, fitted.parameters, lmax, fitted.crossover
        )
        pdf = np.exp(log_pdf)
        increments = 0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid)
        cumulative = np.concatenate(([0.0], np.cumsum(increments)))
        cumulative /= cumulative[-1]
        return np.interp(u, cumulative, grid)
    raise ValueError("Unknown family {!r}".format(fitted.name))
