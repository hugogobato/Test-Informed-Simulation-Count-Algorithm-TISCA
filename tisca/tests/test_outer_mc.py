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


# --------------------------------------------------------------------------- #
# Regression tests for the empirical family (P3-T2 defect, found 2026-08-06).   #
# --------------------------------------------------------------------------- #

def _synthetic_loss_matrix(m=500, seed=11):
    """A correlated, right-skewed (M, 2) matrix standing in for the real one.

    Built rather than loaded so the test is hermetic, but it reproduces the
    features that matter: positive dependence, unequal marginal sds, a non-zero
    mean contrast, and a strongly skewed D.
    """
    rng = np.random.default_rng(seed)
    z = rng.multivariate_normal([0, 0], [[1, 0.6], [0.6, 1]], size=m)
    a = 9.0 + 2.2 * np.exp(0.35 * z[:, 0])
    b = 11.2 + 3.2 * np.exp(0.35 * z[:, 1])
    return np.column_stack([a, b])


def test_empirical_families_are_not_inert():
    """theta must move E[D] and the empirical null must actually be null.

    Regression for the v1 branch ``return matrix[idx]``, which ignored theta,
    rho, sigma_a and sigma_b outright. Because the real matrix carries
    E[D] = -2.20, its ``theta = 0`` cells were not nulls: Module A reported a
    Type I error of 1.000 for D1-D5 and 0.948 for the oracle D6, and its five
    rho levels were exact duplicates of each other.
    """
    from tisca.outermc import families
    mat = _synthetic_loss_matrix()
    assert abs(np.mean(mat[:, 0] - mat[:, 1])) > 1.0, "fixture must carry a real effect"

    for fam, rho in [("empirical", None), ("empirical_copula", 0.3)]:
        for theta in (0.0, 0.5):
            b = families.sample_batch(fam, 1, 200_000, rho=rho, theta=theta,
                                      matrix=mat, master_seed=1)
            d = b[0, :, 0] - b[0, :, 1]
            assert abs(np.mean(d) - theta) < 0.02, (fam, theta, np.mean(d))


def test_empirical_copula_honours_rho_and_sigma_ratio():
    """rho and the variance ratio must both change the copula variant's contrast."""
    from tisca.outermc import families
    mat = _synthetic_loss_matrix()
    sds = [families.sample_batch("empirical_copula", 1, 200_000, rho=r, theta=0.0,
                                 matrix=mat, master_seed=2)[0]
           for r in (-0.3, 0.3, 0.9)]
    widths = [np.std(x[:, 0] - x[:, 1], ddof=1) for x in sds]
    assert widths[0] > widths[1] > widths[2], widths          # rho is live

    wide = families.sample_batch("empirical_copula", 1, 200_000, rho=0.3, theta=0.0,
                                 sigma_a=2.0, matrix=mat, master_seed=2)[0]
    assert np.std(wide[:, 0] - wide[:, 1], ddof=1) > widths[1] * 1.2   # sigma_a is live


def test_empirical_row_bootstrap_preserves_the_real_joint_shape():
    """The bootstrap variant must keep the matrix's own skewness and dependence.

    This is the realism family (g) exists to supply to the P3-T4 skewness
    sub-study, so both are asserted against the source matrix rather than assumed.
    """
    from scipy import stats as _st

    from tisca.outermc import families
    mat = _synthetic_loss_matrix()
    raw = mat[:, 0] - mat[:, 1]
    raw_skew = float(_st.skew(raw))
    raw_r = float(np.corrcoef(mat[:, 0], mat[:, 1])[0, 1])
    assert abs(raw_skew) > 0.5, "fixture must be skewed for this test to bite"

    b = families.sample_batch("empirical", 1, 400_000, rho=None, theta=0.0,
                              matrix=mat, master_seed=4)[0]
    d = b[:, 0] - b[:, 1]
    assert float(_st.skew(d)) == pytest.approx(raw_skew, abs=0.06)
    assert float(np.corrcoef(b[:, 0], b[:, 1])[0, 1]) == pytest.approx(raw_r, abs=0.02)
    # ... and the real marginal variance ratio survives the common rescaling
    assert (np.std(b[:, 0]) / np.std(b[:, 1])) == pytest.approx(
        float(mat[:, 0].std() / mat[:, 1].std()), abs=0.02)


def test_empirical_rejects_an_imposed_rho():
    """A rho cannot be silently ignored: the row bootstrap refuses one."""
    from tisca.outermc import families
    mat = _synthetic_loss_matrix()
    with pytest.raises(ValueError, match="row bootstrap"):
        families.sample_batch("empirical", 1, 10, rho=0.3, matrix=mat, master_seed=0)
    with pytest.raises(ValueError, match="requires a numeric rho"):
        families.sample_batch("normal", 1, 10, rho=None, master_seed=0)
    with pytest.raises(ValueError, match="loss matrix"):
        families.sample_batch("empirical", 1, 10, rho=None, matrix=None, master_seed=0)


def test_oracle_sigma_matches_the_family_it_plans_for():
    """D6's 'true sigma' must be the family's sigma, not the normal formula.

    The closed form ``sqrt(a^2 + b^2 - 2 rho a b)`` is exact only for the
    bivariate normal. It was being used for all seven families, understating the
    true sigma_D by 37-78% for ``mix`` (whose 2% catastrophic component is pure
    added variance) and by 4-24% for the copula-transformed families at
    rho >= 0.6, because the copula attenuates the Pearson correlation relative to
    the design's rank rho. An oracle that plans from a sigma that low under-sizes
    J, so the reference design was reporting optimistic power and E[J] exactly
    where the non-normal families were meant to stress it.
    """
    from tisca.outermc import families
    mat = _synthetic_loss_matrix()
    cases = [("normal", 0.0), ("normal", 0.9), ("lognormal", 0.9), ("gamma", 0.9),
             ("beta", 0.9), ("t3", 0.6), ("mix", 0.0), ("mix", 0.9),
             ("empirical", None), ("empirical_copula", 0.6)]
    for fam, rho in cases:
        cfg = {"family": fam, "rho": rho, "sigma_a": 1.0, "sigma_b": 1.0,
               "matrix": mat if fam.startswith("empirical") else None}
        oracle = engine.sigma_D_true(cfg)
        b = families.sample_batch(fam, 1, 400_000, rho=rho, theta=0.0,
                                  matrix=cfg["matrix"], master_seed=77)
        truth = float(np.std(b[0, :, 0] - b[0, :, 1], ddof=1))
        assert oracle == pytest.approx(truth, rel=0.03), (fam, rho, oracle, truth)

    # the normal closed form is still used, and still exact, for the normal family
    assert engine.sigma_D_true({"family": "normal", "sigma_a": 1.0, "sigma_b": 1.0,
                                "rho": 0.6}) == pytest.approx(np.sqrt(0.8))


def test_empirical_families_recover_nominal_type_I_under_the_oracle():
    """Both empirical variants must sit at 0.05 at theta = 0 with a correct oracle."""
    R = 8000
    mcse = np.sqrt(0.05 * 0.95 / R)
    mat = _synthetic_loss_matrix()
    for fam, rho in [("empirical", None), ("empirical_copula", 0.6)]:
        c = dict(engine.DEFAULT_CONFIG)
        c.update(design="D6", family=fam, rho=rho, theta=0.0, delta=0.5,
                 sigma_D=None, matrix=mat, R=R, Jmax=200, fixed_J=200,
                 alpha=0.05, seed=3)
        s, _, _ = engine.run_e1(c)
        assert abs(s["reject_rate"] - 0.05) <= 4 * mcse, (fam, s["reject_rate"])
        assert abs(s["ci_cover"] - 0.95) <= 4 * mcse, (fam, s["ci_cover"])


def test_oracle_sigma_is_exact_not_sampled():
    """The D6 oracle must be a derived constant, not a Monte Carlo estimate.

    The 1e6-draw estimator it replaced claimed 0.07% accuracy from
    ``sigma/sqrt(2n)``, which presumes a finite fourth moment. ``t3`` has none
    (Student-t with 3 df has finite moments only of order < 3), so no CLT applies
    to its sample sd; measured against the closed form the estimate was off by up
    to 1.4%. These assertions pin the derivations themselves.
    """
    from tisca.outermc import sigma_d

    # normal: textbook closed form
    assert sigma_d.sigma_D_exact("normal", 0.6) == pytest.approx(np.sqrt(0.8), rel=1e-12)
    # mix: normal pair plus an independent 2% symmetric Exp(10) component whose
    # second moment is exactly 0.02 * 2 * 10^2 = 4
    assert sigma_d.sigma_D_exact("mix", 0.0) == pytest.approx(np.sqrt(6.0), rel=1e-12)
    assert sigma_d.sigma_D_exact("mix", 0.9) == pytest.approx(np.sqrt(4.2), rel=1e-12)
    # lognormal c_f has a closed form: E[e^{Z/2} Z] / sd = 0.5 e^{1/8} / sd
    sd_ln = np.sqrt((np.exp(0.25) - 1) * np.exp(0.25))
    assert sigma_d.family_c("lognormal") == pytest.approx(0.5 * np.exp(0.125) / sd_ln,
                                                          rel=1e-10)
    # a_1 of the Mehler expansion must equal c_f: two independent routes agreeing
    for fam in ("lognormal", "gamma", "beta", "t3"):
        a1 = sigma_d._mehler_coeffs(fam)[0]
        assert a1 == pytest.approx(sigma_d.family_c(fam), rel=1e-7), fam
    # determinism: no seed anywhere
    assert (sigma_d.sigma_D_exact("t3", 0.9) == sigma_d.sigma_D_exact("t3", 0.9))


def test_oracle_sigma_beats_the_million_draw_estimate_on_t3():
    """t3 is where the sampled oracle failed. Pin both halves of that claim.

    Half one: a 1e6-draw sample sd disagrees with the closed form by far more than
    the 0.07% that estimator advertised, so the sampled version cannot quietly come
    back and still pass.

    Half two: the closed form is the one that is right. That is established
    WITHOUT Monte Carlo, by recomputing ``c_f`` under a different quadrature rule
    -- fixed-order Gauss-Legendre after a tanh substitution, which places its nodes
    completely differently from the adaptive Gauss-Kronrod rule used in
    ``family_c`` and handles the endpoint singularities by a different mechanism.
    Two independent deterministic rules agreeing to 1e-9 settles it. A Monte Carlo
    cross-check cannot: the sample variance of a t3 contrast has infinite variance,
    so replicate averages converge too slowly to adjudicate anything at this scale
    -- which is precisely the defect being fixed.
    """
    from scipy import stats as _st

    from tisca.outermc import families, sigma_d

    exact = sigma_d.sigma_D_exact("t3", 0.9)
    b = families.sample_batch("t3", 1, 1_000_000, rho=0.9, theta=0.0,
                              master_seed=20260806)
    sampled = float(np.std(b[0, :, 0] - b[0, :, 1], ddof=1))
    assert abs(sampled - exact) / exact > 0.005, (exact, sampled)

    # Independent rule: u = (1 + tanh(t))/2 maps (-inf, inf) -> (0, 1), and the
    # Jacobian 1/(2 cosh^2 t) kills the endpoint singularities geometrically.
    for family in ("t3", "lognormal", "gamma", "beta"):
        mu, sd = families._MOMENTS[family]
        t, w = np.polynomial.legendre.leggauss(4000)
        t = t * 30.0                       # a wide finite window; the rest underflows
        w = w * 30.0
        u = 0.5 * (1.0 + np.tanh(t))
        jac = 0.5 / np.cosh(t) ** 2
        z = _st.norm.ppf(u)
        vals = families._family_marginal(family, u, z) * z * jac
        vals = np.where(np.isfinite(vals), vals, 0.0)
        alt = float(np.dot(w, vals)) / sd
        assert alt == pytest.approx(sigma_d.family_c(family), abs=1e-9), (
            family, alt, sigma_d.family_c(family))


def test_empirical_family_accepts_its_own_default_rho():
    """sample_batch('empirical', matrix=...) must work without naming rho.

    ``rho`` defaulted to 0.0, so the row bootstrap raised on its own default and
    every call site had to spell out ``rho=None``. Three states are now distinct:
    unnamed, explicit None, and an explicit number.
    """
    from tisca.outermc import families

    mat = _synthetic_loss_matrix()
    a = families.sample_batch("empirical", 1, 500, matrix=mat, master_seed=9)
    b = families.sample_batch("empirical", 1, 500, rho=None, matrix=mat, master_seed=9)
    assert np.array_equal(a, b)
    # an explicit numeric rho is still refused
    with pytest.raises(ValueError, match="row bootstrap"):
        families.sample_batch("empirical", 1, 10, rho=0.3, matrix=mat, master_seed=0)
    # and the copula families still default to independence
    assert np.array_equal(
        families.sample_batch("normal", 1, 500, master_seed=4),
        families.sample_batch("normal", 1, 500, rho=0.0, master_seed=4),
    )
