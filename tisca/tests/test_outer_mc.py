"""P2-T3 acceptance tests for the outer-MC harness.

Acceptance (REVISION_PLAN.md P2-T3): on a configuration with a known closed-form
answer (bivariate normal, rho = 0, fixed J, oracle sigma), the harness recovers
the nominal Type I error 0.05 +/- 2 MCSE and the nominal power 0.80 +/- 2 MCSE.

Run:  PYTHONPATH=tisca/python python -m pytest tisca/tests/test_outer_mc.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from tisca.outermc import engine
from tisca.outermc.designs import solve_J


def _acceptance_cell(R, theta):
    c = dict(engine.DEFAULT_CONFIG)
    c.update(design="D6", family="normal", rho=0.0, sigma_a=1.0, sigma_b=1.0,
             theta=theta, delta=0.5, power_target=0.80, sigma_D=None,
             R=R, J0=50, Jmax=2000, K=1, alpha=0.05)
    return c


def test_oracle_type_I_and_power():
    """Main acceptance: within 2 MCSE of 0.05 (Type I) and 0.80 (power)."""
    delta, sd = 0.5, np.sqrt(2.0)
    J = int(np.atleast_1d(solve_J(0.80, delta, sd, 0.05)[0]).max())
    R = 20_000
    mcse0 = np.sqrt(0.05 * 0.95 / R)
    mcse1 = np.sqrt(0.80 * 0.20 / R)

    s0, _, _ = engine.run_e1(_acceptance_cell(R, 0.0))       # theta = 0
    s1, _, _ = engine.run_e1(_acceptance_cell(R, delta))     # theta = delta

    t1e, power = s0["reject_rate"], s1["reject_rate"]
    assert abs(t1e - 0.05) <= 2 * mcse0, f"Type I {t1e:.4f} off by >2 MCSE"
    assert abs(power - 0.80) <= 2 * mcse1, f"Power {power:.4f} off by >2 MCSE"
    # oracle estimate unbiased
    assert abs(s1["E_theta"] - delta) <= 5 * mcse1, "theta_hat biased"


def test_oracle_estimate_and_ci_coverage():
    """CI coverage of theta under the oracle design should be ~ 1 - alpha."""
    delta, sd = 0.5, np.sqrt(2.0)
    J = int(np.atleast_1d(solve_J(0.80, delta, sd, 0.05)[0]).max())
    c = _acceptance_cell(20_000, delta)
    c["fixed_J"] = J
    s, _, _ = engine.run_e1(c)
    mcse = np.sqrt(0.95 * 0.05 / s["R"])
    assert abs(s["ci_cover"] - 0.95) <= 3 * mcse, "CI coverage off"


def test_sigma_D_analytic():
    """sigma_D_TRUE matches sqrt(a^2 + b^2 - 2 rho a b)."""
    assert engine.sigma_D_true({"sigma_a": 1.0, "sigma_b": 1.0, "rho": 0.0}) == pytest.approx(np.sqrt(2.0))
    assert engine.sigma_D_true({"sigma_a": 2.0, "sigma_b": 1.0, "rho": 0.0}) == pytest.approx(np.sqrt(5.0))
    assert engine.sigma_D_true({"sigma_a": 1.0, "sigma_b": 1.0, "rho": 0.6}) == pytest.approx(np.sqrt(2.0 - 1.2))


def test_all_designs_run():
    """All six designs execute without error and return valid per-rep results."""
    for d in ["D1", "D2", "D3", "D4", "D5", "D6"]:
        c = dict(engine.DEFAULT_CONFIG)
        c.update(design=d, family="normal", rho=0.0, theta=0.5, delta=0.5,
                 sigma_D=None, R=300, J0=40, Jmax=500, K=1, alpha=0.05)
        s, res, meta = engine.run_e1(c)
        assert len(res) == 300
        assert np.all(np.isfinite(res.col("J"))) and np.all(res.col("J") >= 1)
        assert 0.0 <= s["reject_rate"] <= 1.0


def test_all_families_sampled():
    """Each loss family samples (n,2) with the intended D distribution scale."""
    from tisca.outermc import families
    for fam in ["normal", "lognormal", "gamma", "mix", "beta", "t3"]:
        # The 'mix' family has a 2% catastrophic tail; its mean needs a much
        # larger n to converge onto theta (it is unbiased, just high-variance).
        n = 50_000 if fam == "mix" else 2000
        mat = families.sample_pairs(fam, n=n, rho=0.0, sigma_a=1.0, sigma_b=1.0, theta=0.5,
                                    seed=0)
        assert mat.shape == (n, 2)
        D = mat[:, 0] - mat[:, 1]
        assert abs(np.mean(D) - 0.5) < 0.15, f"family {fam} effect drift"
