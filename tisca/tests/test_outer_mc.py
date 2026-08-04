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


# --------------------------------------------------------------------------- #
# Regression tests for the harness defects found in the Phase 2/3 audit.        #
# --------------------------------------------------------------------------- #

_ALL_FAMILIES = ["normal", "lognormal", "gamma", "mix", "beta", "t3"]


def test_every_family_recovers_nominal_type_I_under_the_oracle():
    """The oracle design must hit 0.05 in EVERY family, not just the normal one.

    Regression test for the per-repetition standardisation in
    ``families._standardise_scaled``. Dividing each repetition by its OWN sample
    mean and sd pinned ``D_bar`` to ``theta`` exactly, so at ``theta = 0`` the
    lognormal / gamma / beta / t3 cells returned a Type I error of exactly 0.000
    and a CI coverage of exactly 1.000 -- silently, and precisely in the four
    families the P3-T4 skewness sub-study reads.

    ``mix`` is excluded from the nominal band on purpose: its 2% catastrophic
    component leaves the paired t genuinely conservative (~0.039 at J = 200),
    which is a finding to report, not a defect to assert away. It is checked
    separately below.
    """
    R = 8000
    mcse = np.sqrt(0.05 * 0.95 / R)
    for fam in _ALL_FAMILIES:
        c = dict(engine.DEFAULT_CONFIG)
        c.update(design="D6", family=fam, rho=0.3, theta=0.0, delta=0.5,
                 sigma_D=None, R=R, Jmax=200, fixed_J=200, alpha=0.05, seed=3)
        s, _, _ = engine.run_e1(c)
        # No family may come back degenerate -- that was the bug.
        assert 0.01 < s["reject_rate"] < 0.12, (fam, s["reject_rate"])
        assert 0.88 < s["ci_cover"] < 0.99, (fam, s["ci_cover"])
        if fam != "mix":
            assert abs(s["reject_rate"] - 0.05) <= 4 * mcse, (fam, s["reject_rate"])
            assert abs(s["ci_cover"] - 0.95) <= 4 * mcse, (fam, s["ci_cover"])


def test_mix_family_is_conservative_not_degenerate():
    """The catastrophic-failure family is genuinely below nominal -- a finding."""
    R = 20000
    c = dict(engine.DEFAULT_CONFIG)
    c.update(design="D6", family="mix", rho=0.3, theta=0.0, R=R, Jmax=200,
             fixed_J=200, alpha=0.05, seed=3)
    s, _, _ = engine.run_e1(c)
    assert 0.03 < s["reject_rate"] < 0.047, s["reject_rate"]


def test_dbar_retains_its_sampling_variability():
    """sd(D_bar) across repetitions must match the i.i.d. value s_D/sqrt(J)."""
    from tisca.outermc import families
    for fam in _ALL_FAMILIES:
        b = families.draw_rep_losses(fam, 3000, 100, rho=0.3, theta=0.0, master_seed=7)
        D = b[..., 0] - b[..., 1]
        observed = D.mean(axis=1).std()
        expected = D.std(axis=1, ddof=1).mean() / np.sqrt(100)
        assert 0.85 < observed / expected < 1.25, (fam, observed, expected)


def test_skewed_families_actually_produce_a_skewed_contrast():
    """P3-T4 plots Type I error against the skewness of D, so it must be non-zero.

    With the same marginal on both methods and an exchangeable copula the pair
    is exchangeable and ``D`` is exactly symmetric, i.e. the skewness factor is
    inert. ``asym=True`` (default) breaks the exchangeability.

    ``t3`` is deliberately NOT in this list: Student-t is symmetric and its third
    moment does not exist at 3 df, so it is the heavy-TAIL family, not a skewed
    one. It is checked on kurtosis instead.
    """
    from tisca.outermc.families import contrast_skewness, draw_rep_losses
    for fam in ["lognormal", "gamma", "beta", "mix"]:
        sk = contrast_skewness(fam, rho=0.3, n=200_000, seed=1)
        assert abs(sk) > 0.15, (fam, sk)
    # the exchangeable construction is the symmetric control case
    for fam in ["lognormal", "gamma"]:
        assert abs(contrast_skewness(fam, rho=0.3, asym=False, n=200_000, seed=1)) < 0.05
    # t3 supplies tail weight rather than asymmetry
    b = draw_rep_losses("t3", 1, 200_000, rho=0.3, theta=0.0, master_seed=2)
    d = b[0, :, 0] - b[0, :, 1]
    d = (d - d.mean()) / d.std()
    assert float((d ** 4).mean()) - 3.0 > 5.0


def test_v1_unpaired_design_is_conservative_when_methods_correlate():
    """D3 must be planned from the UNPAIRED variance, or D3 == D4 and E1 is moot.

    Planning D3 from ``config['sigma_D']`` handed the v1 machinery the paired
    variance it never had, which collapsed the D3-vs-D4 contrast that the whole
    "fixing IJDA #2 buys ~2.3x" claim rests on. With the unpaired planning scale
    restored, E[J] separates with rho and the v1 test is measurably conservative
    (not liberal) exactly in the regime the case study lives in (r = 0.58-0.80).
    """
    out = {}
    for rho in (0.0, 0.6):
        for d in ("D3", "D4"):
            c = dict(engine.DEFAULT_CONFIG)
            c.update(design=d, family="normal", rho=rho, theta=0.0, delta=0.5,
                     sigma_D=None, R=3000, J0=50, Jmax=1000, alpha=0.05, seed=5)
            out[(rho, d)] = engine.run_e1(c)[0]
    # no pairing gain at rho = 0 ...
    assert out[(0.0, "D3")]["E_J"] / out[(0.0, "D4")]["E_J"] < 1.2
    # ... a substantial one at rho = 0.6
    assert out[(0.6, "D3")]["E_J"] / out[(0.6, "D4")]["E_J"] > 1.6
    # and the v1 test is conservative there, while v2 holds its level
    assert out[(0.6, "D3")]["reject_rate"] < 0.02
    assert abs(out[(0.6, "D4")]["reject_rate"] - 0.05) < 0.015


def test_solve_J_bisection_matches_brute_force():
    """The bisection solver returns the exact smallest J the dense grid would."""
    from tisca.outermc.designs import m1_power
    rng = np.random.default_rng(0)
    sd = rng.uniform(0.4, 3.0, size=40)
    J, _ = solve_J(0.80, 0.5, sd, 0.05, mode=1, Jmax=1000)
    for s, j in zip(sd, J):
        assert m1_power(j, 0.5, s, 0.05) >= 0.80
        if j > 4:
            assert m1_power(j - 1, 0.5, s, 0.05) < 0.80
