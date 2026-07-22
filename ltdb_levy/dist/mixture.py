"""The exact flat-plus-power-law relocation distribution used by ``func.c``.

For crossover ``c`` and upper support ``L`` the density is

``f(l) = a * min(1, (l / c) ** -mu),  0 <= l <= L``.

The C simulator has ``c = 1`` in its internal units. The empirical pipeline can
either impose a physical crossover or fit one crossover per condition before
normalizing a later C comparison.
"""

from __future__ import annotations

from typing import Optional

import numpy as np


def _validate(mu: float, lmax: float, crossover: float) -> None:
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu must be finite and positive")
    if not np.isfinite(lmax) or lmax <= 0.0:
        raise ValueError("lmax must be finite and positive")
    if not np.isfinite(crossover) or crossover <= 0.0:
        raise ValueError("crossover must be finite and positive")


def normalization_integral(mu: float, lmax: float, crossover: float = 1.0) -> float:
    """Return the integral of the unnormalised density over ``[0, lmax]``.

    ``expm1`` supplies a stable continuous limit at ``mu = 1``.  This remains
    accurate at ``mu = 1 +/- 1e-9``, where direct subtraction loses digits.
    """
    _validate(mu, lmax, crossover)
    if lmax <= crossover:
        return float(lmax)
    log_ratio = float(np.log(lmax / crossover))
    q = 1.0 - mu
    if abs(q) < 1e-8:
        tail_factor = log_ratio if q == 0.0 else float(np.expm1(q * log_ratio) / q)
    else:
        tail_factor = float(np.expm1(q * log_ratio) / q)
    return float(crossover * (1.0 + tail_factor))


def normalization_constant(mu: float, lmax: float, crossover: float = 1.0) -> float:
    """Return ``a = 1 / integral``, matching ``get_normalization_constant`` in C."""
    return 1.0 / normalization_integral(mu, lmax, crossover)


def log_density(
    lengths: np.ndarray, mu: float, lmax: float, crossover: float = 1.0
) -> np.ndarray:
    """Evaluate log density, returning ``-inf`` outside the closed support."""
    _validate(mu, lmax, crossover)
    values = np.asarray(lengths, dtype=float)
    out = np.full(values.shape, -np.inf, dtype=float)
    inside = np.isfinite(values) & (values >= 0.0) & (values <= lmax)
    out[inside] = -np.log(normalization_integral(mu, lmax, crossover))
    tail = inside & (values >= crossover)
    out[tail] -= mu * np.log(values[tail] / crossover)
    return out


def density(
    lengths: np.ndarray, mu: float, lmax: float, crossover: float = 1.0
) -> np.ndarray:
    """Evaluate the probability density."""
    return np.exp(log_density(lengths, mu, lmax, crossover))


def cdf(
    lengths: np.ndarray, mu: float, lmax: float, crossover: float = 1.0
) -> np.ndarray:
    """Evaluate the analytic CDF."""
    _validate(mu, lmax, crossover)
    values = np.asarray(lengths, dtype=float)
    z = normalization_integral(mu, lmax, crossover)
    out = np.zeros(values.shape, dtype=float)
    out[values >= lmax] = 1.0
    flat = (values > 0.0) & (values < min(crossover, lmax))
    out[flat] = values[flat] / z
    if lmax > crossover:
        tail = (values >= crossover) & (values < lmax)
        log_ratio = np.log(values[tail] / crossover)
        q = 1.0 - mu
        if abs(q) < 1e-8:
            tail_factor = log_ratio if q == 0.0 else np.expm1(q * log_ratio) / q
        else:
            tail_factor = np.expm1(q * log_ratio) / q
        out[tail] = crossover * (1.0 + tail_factor) / z
    return np.clip(out, 0.0, 1.0)


def ppf(
    probabilities: np.ndarray, mu: float, lmax: float, crossover: float = 1.0
) -> np.ndarray:
    """Evaluate the inverse CDF, including the analytic ``mu = 1`` limit."""
    _validate(mu, lmax, crossover)
    u = np.asarray(probabilities, dtype=float)
    if np.any(~np.isfinite(u)) or np.any((u < 0.0) | (u > 1.0)):
        raise ValueError("probabilities must be finite and lie in [0, 1]")
    z = normalization_integral(mu, lmax, crossover)
    if lmax <= crossover:
        return u * lmax

    kink_probability = crossover / z
    out = np.empty(u.shape, dtype=float)
    flat = u <= kink_probability
    out[flat] = u[flat] * z
    tail = ~flat
    q = 1.0 - mu
    tail_coordinate = u[tail] * z / crossover - 1.0
    if abs(q) < 1e-8:
        if q == 0.0:
            log_ratio = tail_coordinate
        else:
            log_ratio = np.log1p(q * tail_coordinate) / q
        out[tail] = crossover * np.exp(log_ratio)
    else:
        out[tail] = crossover * np.power(1.0 + q * tail_coordinate, 1.0 / q)
    out[u == 1.0] = lmax
    return np.clip(out, 0.0, lmax)


def sample(
    mu: float,
    lmax: float,
    size: int,
    rng: np.random.Generator,
    crossover: float = 1.0,
    uniforms: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Draw relocation lengths by inverse transform."""
    if size < 0:
        raise ValueError("size must be non-negative")
    u = rng.random(size) if uniforms is None else np.asarray(uniforms, dtype=float)
    if u.shape != (size,):
        raise ValueError("uniforms must have shape ({},)".format(size))
    return ppf(u, mu, lmax, crossover)
