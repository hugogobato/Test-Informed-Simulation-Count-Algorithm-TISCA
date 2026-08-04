"""Tests for `inference.py` — paired t, studentized paired bootstrap, McNemar,
and listwise pairing."""

from __future__ import annotations

import numpy as np
import pytest

from tisca import inference as INF
from tisca import validate


def test_paired_t_rejects_alternative_known():
    """Paired t correctly detects a negative mean (lower-is-better win)."""
    rng = np.random.default_rng(0)
    D = rng.normal(-0.4, 1.0, size=80)
    res = INF.paired_t(D, mu0=0.0, alternative="less")
    assert res["p_value"] < 0.05
    assert res["estimate"] == pytest.approx(np.mean(D))
    assert res["df"] == 79
    assert res["se"] == pytest.approx(np.std(D, ddof=1) / np.sqrt(80))


def test_paired_t_two_sided_scale():
    li = INF.paired_t(np.linspace(0, 1, 50), alternative="two-sided")
    le = INF.paired_t(np.linspace(0, 1, 50), alternative="less")
    gr = INF.paired_t(np.linspace(0, 1, 50), alternative="greater")
    # For a mean clearly > 0: greater p is tiny, less p ~ 1, two-sided = 2*min.
    assert gr["p_value"] < 0.05
    assert le["p_value"] == pytest.approx(1.0)
    assert li["p_value"] == pytest.approx(2 * min(le["p_value"], gr["p_value"]), abs=1e-9)
    # For a mean exactly 0, p ~ 1.
    assert INF.paired_t(np.zeros(20))["p_value"] == pytest.approx(1.0)


def test_paired_t_requires_two_reps():
    with pytest.raises(validate.ValidationError):
        INF.paired_t(np.array([1.0]))


def test_studentized_bootstrap_ci_covers():
    rng = np.random.default_rng(2)
    D = rng.normal(0.2, 1.0, size=300)
    res = INF.studentized_paired_bootstrap(D, B=999, seed=3)
    assert res["estimate"] == pytest.approx(np.mean(D))
    lo, hi = res["ci"]
    assert lo < np.mean(D) < hi
    assert res["mcse"] == pytest.approx(np.std(D, ddof=1) / np.sqrt(300))


def test_studentized_bootstrap_detects_negative_mean():
    rng = np.random.default_rng(4)
    D = rng.normal(-0.3, 1.0, size=200)
    res = INF.studentized_paired_bootstrap(D, B=499, alternative="less", seed=5)
    assert res["p_value"] < 0.05


def test_studentized_bootstrap_skewed_distribution():
    """Works with a skewed contrast (heavy bootstrap case, no CLT reliance)."""
    rng = np.random.default_rng(6)
    chi = rng.chisquare(df=2, size=200)
    # Make a skewed contrast by combining a large variance component.
    D = chi - chi[0]  # centered, skewed positive
    D = D - np.median(D)
    res = INF.studentized_paired_bootstrap(D, B=499, seed=7)
    assert np.isfinite(res["ci"][0]) and np.isfinite(res["ci"][1])


def test_studentized_bootstrap_bootstrap_indices_match():
    D = np.random.default_rng(0).normal(size=50)
    idx = np.array([np.random.default_rng(b + 10).integers(0, 50, size=50) for b in range(499)]).T
    r = INF.studentized_paired_bootstrap(D, B=499, seed=1, bootstrap_indices=idx)
    assert r["B"] == 499


def test_mcnemar_exact():
    """McNemar on a paired 2x2 table gives the exact binomial p-value."""
    a = np.array([1, 1, 1, 0, 0, 1, 1, 0])  # A outcome
    b = np.array([0, 0, 1, 1, 1, 0, 0, 0])  # B outcome
    res = INF.mcnemar_exact(a, b)
    n01, n10 = res["n01"], res["n10"]
    assert res["m_discordant"] == n01 + n10
    # Recompute independently.
    from scipy.stats import binom
    k = min(n01, n10)
    assert res["p_value"] == pytest.approx(min(1.0, 2 * binom.cdf(k, n01 + n10, 0.5)))


def test_mcnemar_no_discordance_p1():
    res = INF.mcnemar_exact(np.array([1, 1, 1]), np.array([1, 1, 1]))
    assert res["p_value"] == 1.0


def test_mcnemar_rejects_non_binary():
    with pytest.raises(validate.ValidationError):
        INF.mcnemar_exact(np.array([1, 2, 3]), np.array([0, 0, 0]))


def test_contrast_summary_listwise_deletion():
    A = np.array([1.0, np.nan, 3.0, 4.0])
    B = np.array([0.0, 0.0, np.nan, 4.0])
    res = INF.contrast_summary(A, B, name="X")
    assert res["n_dropped"] == 2
    D_kept = A[[0, 3]] - B[[0, 3]]
    assert res["estimate"] == pytest.approx(np.mean(D_kept))
    assert "ci" in res and res["paired"] is True


def test_contrast_summary_mismatched_length():
    with pytest.raises(validate.ValidationError):
        INF.contrast_summary(np.array([1.0, 2.0]), np.array([1.0]))


def test_paired_t_degenerate_and_bad_alternative():
    """Degenerate contrasts (sd=0) are handled; unknown alternatives rejected."""
    D = np.zeros(20)
    res2 = INF.paired_t(D, alternative="two-sided")
    assert res2["p_value"] == pytest.approx(1.0)
    res_less = INF.paired_t(D, alternative="less")
    assert res_less["p_value"] == pytest.approx(1.0)
    res_gr = INF.paired_t(D, alternative="greater")
    assert res_gr["p_value"] == pytest.approx(1.0)
    # A non-zero degenerate mean yields a definitive decision.
    res_diff = INF.paired_t(np.zeros(20) + 1.0, mu0=0.0, alternative="two-sided")
    assert res_diff["p_value"] == 0.0
    with pytest.raises(validate.ValidationError):
        INF.paired_t(np.arange(5.0), alternative="bogus")


def test_studentized_bootstrap_degenerate_and_validation():
    """Zero-variance contrasts and bad arguments terminate gracefully."""
    res = INF.studentized_paired_bootstrap(np.zeros(20), B=99, seed=1)
    assert res["degenerate"] is True and res["estimate"] == 0.0
    with pytest.raises(validate.ValidationError):
        INF.studentized_paired_bootstrap(np.arange(5.0), alternative="bogus")
    with pytest.raises(validate.ValidationError):
        INF.studentized_paired_bootstrap(
            np.arange(20.0), B=100, bootstrap_indices=np.zeros((20, 5), dtype=int)
        )
    with pytest.raises(validate.ValidationError):
        INF.studentized_paired_bootstrap(np.array([1.0]))


def test_studentized_bootstrap_one_sided_alternatives():
    plt = INF.studentized_paired_bootstrap(np.linspace(-1, 1, 300), B=199,
                                           alternative="less", seed=1)
    pgr = INF.studentized_paired_bootstrap(np.linspace(-1, 1, 300), B=199,
                                           alternative="greater", seed=2)
    # Near-zero mean: both one-sided p-values are moderate.
    assert 0 < plt["p_value"] < 1 and 0 < pgr["p_value"] < 1
    # "less" CI is unbounded below; "greater" CI is unbounded above.
    assert plt["ci"][0] == -np.inf and np.isfinite(plt["ci"][1])
    assert np.isfinite(pgr["ci"][0]) and pgr["ci"][1] == np.inf


def test_mcnemar_with_continuity():
    a = np.array([1, 1, 1, 0, 0, 0])
    b = np.array([0, 0, 0, 1, 1, 1])
    res = INF.mcnemar_with_continuity(a, b)
    assert res["n01"] == 3 and res["n10"] == 3
    assert res["chi2"] >= 0.0 and 0 <= res["p_value"] <= 1
    # No discordance -> chi2 = 0, p = 1.
    res_none = INF.mcnemar_with_continuity(np.array([1, 1, 1]), np.array([1, 1, 1]))
    assert res_none["p_value"] == pytest.approx(1.0)
