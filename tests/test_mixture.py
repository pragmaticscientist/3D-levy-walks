from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad

from ltdb_levy.bootstrap import clustered_mu_bootstrap
from ltdb_levy.dist.gof import parametric_bootstrap_ks
from ltdb_levy.dist.mixture import cdf, density, normalization_constant, ppf, sample
from ltdb_levy.dist.mle import fit_mixture


@pytest.mark.parametrize("mu", [1.0 - 1e-9, 1.0, 1.0 + 1e-9, 2.3])
def test_normalization_and_inverse_cdf_are_stable_near_one(mu):
    lmax = 73.0
    integral, _ = quad(
        lambda value: float(density(np.array([value]), mu, lmax)[0]),
        0.0,
        lmax,
        points=[1.0],
    )
    assert integral == pytest.approx(1.0, abs=2e-10)
    probabilities = np.linspace(0.0, 1.0, 1001)
    values = ppf(probabilities, mu, lmax)
    assert np.max(np.abs(cdf(values, mu, lmax) - probabilities)) < 2e-12
    assert values[-1] == lmax


def test_c_normalization_constant_matches_closed_form_at_mu_two():
    lmax = 10.0
    expected = 1.0 / (1.0 + (1.0 - 1.0 / lmax))
    assert normalization_constant(2.0, lmax) == pytest.approx(expected)


def test_mle_recovers_known_exponent_from_20000_draws():
    rng = np.random.default_rng(2031)
    values = sample(2.2, 50.0, 20_000, rng)
    fitted = fit_mixture(values, lmax=50.0, mu_bounds=(1.001, 6.0))
    assert fitted.mu == pytest.approx(2.2, abs=0.06)
    assert not fitted.at_bound


def test_physical_crossover_matches_explicit_length_normalization():
    """Changing physical units must not change the fitted exponent."""
    rng = np.random.default_rng(2032)
    reference_length_um = 3.75
    values_um = sample(
        2.4,
        75.0,
        20_000,
        rng,
        crossover=reference_length_um,
    )
    physical_fit = fit_mixture(
        values_um,
        lmax=75.0,
        crossover=reference_length_um,
        mu_bounds=(1.001, 6.0),
    )
    normalized_fit = fit_mixture(
        values_um / reference_length_um,
        lmax=75.0 / reference_length_um,
        crossover=1.0,
        mu_bounds=(1.001, 6.0),
    )
    # Bounded scalar optimization can stop at slightly different floating-point
    # iterates after the objective's additive unit-conversion shift.
    assert physical_fit.mu == pytest.approx(normalized_fit.mu, abs=1e-7)


def test_joint_mle_recovers_exponent_and_crossover():
    rng = np.random.default_rng(2033)
    values = sample(2.4, 75.0, 20_000, rng, crossover=3.75)
    fitted = fit_mixture(
        values,
        lmax=75.0,
        crossover=3.75,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    assert fitted.mu == pytest.approx(2.4, abs=0.08)
    assert fitted.crossover == pytest.approx(3.75, abs=0.25)
    assert fitted.crossover_estimated
    assert fitted.n_parameters == 2
    assert not fitted.at_bound
    assert not fitted.crossover_at_bound


def test_joint_fit_is_invariant_to_length_units():
    rng = np.random.default_rng(2034)
    scale = 3.75
    values_um = sample(2.2, 60.0, 10_000, rng, crossover=scale)
    physical = fit_mixture(
        values_um,
        lmax=60.0,
        crossover=scale,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    normalized = fit_mixture(
        values_um / scale,
        lmax=60.0 / scale,
        crossover=1.0,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    assert physical.mu == pytest.approx(normalized.mu, abs=2e-6)
    assert physical.crossover / scale == pytest.approx(
        normalized.crossover, abs=2e-6
    )


def test_joint_parameters_are_refitted_in_bootstraps():
    rng = np.random.default_rng(2035)
    pieces = []
    for track_index in range(8):
        values = sample(2.3, 40.0, 80, rng, crossover=3.0)
        pieces.append(
            pd.DataFrame(
                {
                    "track_uid": "track{}".format(track_index),
                    "run_length_um": values,
                }
            )
        )
    runs = pd.concat(pieces, ignore_index=True)
    values = runs["run_length_um"].to_numpy()
    fitted = fit_mixture(
        values,
        lmax=40.0,
        crossover=3.0,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    gof = parametric_bootstrap_ks(
        values,
        fitted,
        20,
        rng,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    clustered = clustered_mu_bootstrap(
        runs,
        20,
        rng,
        crossover=3.0,
        mu_bounds=(1.001, 6.0),
        estimate_crossover=True,
    )
    assert np.isfinite(gof.p_value)
    assert clustered.n_successful == 20
    assert np.isfinite(clustered.draws).all()
    assert np.isfinite(clustered.crossover_draws).all()
    assert clustered.crossover_ci_low < clustered.crossover_ci_high

    step_weighted = clustered_mu_bootstrap(
        runs,
        10,
        rng,
        crossover=3.0,
        mu_bounds=(1.001, 6.0),
        weighting="step_weighted",
        estimate_crossover=True,
    )
    assert step_weighted.n_successful == 10
    assert np.isfinite(step_weighted.draws).all()
