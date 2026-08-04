"""Tests for `estimands.py` — estimand-table rows and MCSE formulas."""

from __future__ import annotations

import numpy as np
import pytest

from tisca import estimands as E


def test_rooted_vs_unrooted_pehe():
    tau_hat = np.array([0.0, 1.0, 2.0])
    tau_true = np.array([0.0, 0.0, 0.0])
    q = E.cate_mse(tau_hat, tau_true)
    pehe = E.rooted_pehe(tau_hat, tau_true)
    assert q == pytest.approx(np.mean([0, 1, 4]))
    assert pehe == pytest.approx(np.sqrt(q))
    # Jensen: E[sqrt Q] != sqrt(E[Q]) for non-degenerate Q.
    rng = np.random.default_rng(0)
    tau_hats = rng.normal(size=(100, 20))
    taus = np.zeros((100, 20))
    qs = E.cate_mse(tau_hats, taus)
    pehes = E.rooted_pehe(tau_hats, taus)
    assert np.mean(pehes) != pytest.approx(np.sqrt(np.mean(qs)), abs=1e-3)
    assert np.mean(pehes) < np.sqrt(np.mean(qs))  # Jensen strict concavity of sqrt


def test_ate_squared_error_and_rmse():
    ate_hat = np.array([1.0, 2.0, 3.0])
    ate_true = np.array([0.0, 0.0, 0.0])
    q = E.ate_squared_error(ate_hat, ate_true)
    assert q.tolist() == [1.0, 4.0, 9.0]
    assert E.rmse_ate(q) == pytest.approx(np.sqrt(np.mean([1, 4, 9])))
    assert E.rmse_ate(q) == pytest.approx(2.1602469)


def test_rmse_mcse_delta_has_inner_root():
    """Row 3b: denominator carries sqrt(mean Q), not mean Q (P1-T1 verification)."""
    q = np.array([0.25] * 100)
    mcse_correct = E.rmse_ate_mcse_delta(q, sqrt_n=True)
    mcse_wrong = E.rmse_ate_mcse_delta(q, sqrt_n=False)
    # At E[Q]=0.25, dropping the root inflates the MCSE by exactly ~2x.
    assert mcse_correct <= mcse_wrong
    # Verify against the analytic form.
    mu, s = 0.25, float(np.std(q, ddof=1))
    assert mcse_correct == pytest.approx(s / (2 * np.sqrt(100) * np.sqrt(mu)))


def test_coverage_binary_and_cate():
    ate_hat = np.array([0.0, 5.0])
    ate_true = np.array([0.0, 5.0])
    cov_ate = E.ate_coverage(ate_hat, np.array([-1.0, 4.0]), np.array([1.0, 6.0]), ate_true)
    assert cov_ate.tolist() == [1.0, 1.0]
    cov_cate = E.cate_coverage(
        np.array([[0.0], [5.0]]),
        np.array([[-1.0], [4.0]]),
        np.array([[1.0], [6.0]]),
        np.array([[0.0], [5.0]]),
    )
    assert cov_cate.tolist() == [1.0, 1.0]


def test_interval_score_proper_scoring_rule():
    # Perfect coverage gives the width as score.
    is_p = E.interval_score_unit(0.0, 2.0, 1.0, level=0.95)
    assert is_p == pytest.approx(2.0)
    # A miss far below incurs the (2/c)(lo - x) penalty.
    is_bad_lo = E.interval_score_unit(0.0, 2.0, -10.0, level=0.95)
    assert is_bad_lo == pytest.approx(2.0 + (2 / 0.05) * (0.0 - (-10.0)))
    # Vectorised mean form equals the single-unit form.
    lo = np.array([0.0, -1.0]); hi = np.array([2.0, 3.0]); x = np.array([1.0, 0.5])
    assert E.interval_score(lo, hi, x, level=0.95) == pytest.approx(
        np.mean([E.interval_score_unit(0, 2, 1, level=0.95), E.interval_score_unit(-1, 3, 0.5, level=0.95)])
    )


def test_mean_mcse_and_ci():
    rng = np.random.default_rng(1)
    x = rng.normal(size=200)
    assert E.mean_estimate(x) == pytest.approx(np.mean(x))
    assert E.mean_mcse(x) == pytest.approx(np.std(x, ddof=1) / np.sqrt(len(x)))
    lo, hi = E.nc_ci(0.0, 1.0, level=0.95, df=199)
    assert lo < 0 < hi


def test_crps_mean():
    assert E.crps_mean(np.array([[0.1, 0.2, 0.3], [0.2, 0.2, 0.2]])) == pytest.approx([0.2, 0.2])


def test_calibration_deviation_on_estimate_not_per_rep():
    """Row 8: |mean(Cov)-0.95|, and must NOT equal mean(|Cov-0.95|)."""
    cov = np.array([0.94, 0.96, 0.95])
    dev = E.calibration_deviation(np.mean(cov), 0.95)
    assert dev == pytest.approx(0.0, abs=1e-6)
    # The forbidden object mean(|Cov-0.95|) is larger than the valid deviation.
    cov2 = np.array([0.90, 1.00])
    valid = E.calibration_deviation(np.mean(cov2), 0.95)  # |0.95-0.95|=0
    forbidden = np.mean(np.abs(cov2 - 0.95))  # 0.05
    assert valid < forbidden


def test_rmse_mcse_monotone_with_variance():
    q_low = np.array([0.25] * 100)
    q_high = np.array(list(q_low) + [0.0, 2.0, 5.0])
    assert E.rmse_ate_mcse_delta(q_high) > E.rmse_ate_mcse_delta(q_low)
