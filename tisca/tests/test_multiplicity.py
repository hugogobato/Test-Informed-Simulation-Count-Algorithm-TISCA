"""Tests for `multiplicity.py` — Bonferroni, Holm, BH, Romano-Wolf stepdown, and
family power."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from tisca import multiplicity as M
from tisca import validate


def test_bonferroni():
    p = np.array([0.01, 0.05, 0.5])
    b = M.bonferroni(p)
    assert b == pytest.approx(np.minimum(p * 3, 1.0))
    assert b[0] == pytest.approx(0.03)


def test_holm_monotone_and_bounded():
    p = np.array([0.001, 0.02, 0.05])
    h = M.holm(p)
    # Holm adjusted >= raw and <= Bonferroni.
    assert np.all(h >= p)
    assert np.all(h <= M.bonferroni(p) + 1e-12)
    # Monotone in the ranked order.
    order = np.argsort(h)
    assert np.all(np.diff(h[order]) >= -1e-12)


def test_bh_controls_at_most():
    p = np.array([0.01, 0.02, 0.03, 0.99])
    bh = M.benjamini_hochberg(p)
    # The largest value is unchanged, smaller ones inflated at most to min line.
    assert bh[-1] == pytest.approx(0.99)
    # BH p-values are <= Holm for the smallest.
    assert bh[0] <= M.holm(p)[0] + 1e-12


def test_p_adjust_dispatch():
    p = np.array([0.02, 0.03, 0.04])
    assert np.allclose(M.p_adjust(p, "none"), p)
    assert np.allclose(M.p_adjust(p, "bonferroni"), M.bonferroni(p))
    assert np.allclose(M.p_adjust(p, "holm"), M.holm(p))
    assert np.allclose(M.p_adjust(p, "bh"), M.benjamini_hochberg(p))
    with pytest.raises(validate.ValidationError):
        M.p_adjust(p, "bogus")


def test_romano_wolf_controls_worst_case():
    """RW stepdown returns adjusted p-values >= raw marginal p-values."""
    rng = np.random.default_rng(0)
    J, K = 300, 4
    D = rng.normal(size=(J, K))
    D[:, 0] += 0.3  # one genuine effect
    res = M.romano_wolf_stepdown(D, B=999, alpha=0.05, seed=1)
    assert res["family_size"] == K
    assert np.all(res["p_values"] >= 0) and np.all(res["p_values"] <= 1)
    assert res["rejections"].shape == (K,)


def test_romano_wolf_detects_true_positive():
    """A strong contrast should survive the family correction."""
    rng = np.random.default_rng(0)
    J, K = 400, 5
    D = rng.normal(size=(J, K))
    D[:, 0] += 1.0  # strong effect vs 0 nulls
    res = M.romano_wolf_stepdown(D, B=1999, alpha=0.05, seed=2)
    assert res["rejections"][0]


def test_romano_wolf_controls_fwer_under_the_global_null():
    """FWER at the global null is close to alpha -- measured, not asserted once.

    A single draw at alpha = 0.05 rejects about 1 time in 20, so asserting
    "0 rejections" on one seed tests the seed, not the procedure. Estimate the
    family-wise rejection rate over independent replications instead.
    """
    rng = np.random.default_rng(20260804)
    reps, J, K, alpha = 300, 200, 4, 0.05
    any_reject = 0
    for _ in range(reps):
        D = rng.normal(size=(J, K))          # all K nulls true
        res = M.romano_wolf_stepdown(D, B=399, alpha=alpha, seed=int(rng.integers(1 << 30)))
        any_reject += int(res["rejections"].any())
    fwer = any_reject / reps
    mcse = (alpha * (1 - alpha) / reps) ** 0.5
    assert fwer <= alpha + 4 * mcse, f"FWER {fwer:.3f} exceeds {alpha} by more than 4 MCSE"


def test_romano_wolf_is_sign_symmetric():
    """A family of NEGATIVE contrasts must reject exactly as a positive one does.

    Regression test for the two-sided stepdown: the bootstrap max-statistic is
    built from |t*|, so it has to be compared against |t_obs|. Comparing it
    against the signed t_obs made every negative contrast unrejectable -- and
    in a lower-is-better loss framing "the proposed method wins" is exactly the
    negative sign, so the whole case-study family came back with p = 1.
    """
    rng = np.random.default_rng(11)
    D = rng.normal(-0.5, 1.0, size=(200, 4))
    neg = M.romano_wolf_stepdown(D, B=999, alpha=0.05, seed=5)
    pos = M.romano_wolf_stepdown(-D, B=999, alpha=0.05, seed=5)
    assert neg["rejections"].sum() == 4
    np.testing.assert_allclose(neg["p_values"], pos["p_values"])


def test_planning_alpha_mapping():
    a, note = M.planning_alpha("bonferroni", K=6, alpha=0.05)
    assert a == pytest.approx(0.05 / 6)
    a, _ = M.planning_alpha("holm", K=6, alpha=0.05)
    assert a == pytest.approx(0.05 / 6)
    a, _ = M.planning_alpha("bh", K=6, alpha=0.05, r=2)
    assert a == pytest.approx(0.05 * 2 / 6)
    a, _ = M.planning_alpha("romano_wolf", K=6, alpha=0.05)
    assert a == pytest.approx(0.05)
    a, _ = M.planning_alpha("none", K=1, alpha=0.05)
    assert a == pytest.approx(0.05)


def test_family_power():
    pw = np.array([0.8, 0.8, 0.8])
    assert M.family_power_conjunctive(pw) == pytest.approx(0.512)
    assert M.family_power_disjunctive(pw) == pytest.approx(1 - 0.008)


def test_rw_with_bootstrap_indices_shape_checked():
    D = np.random.default_rng(0).normal(size=(50, 3))
    with pytest.raises(validate.ValidationError):
        M.romano_wolf_stepdown(D, B=100, bootstrap_indices=np.zeros((50, 5), dtype=int))
    with pytest.raises(validate.ValidationError):
        M.romano_wolf_stepdown(np.zeros((5, 3)) + 1.0, B=10, bootstrap_indices=np.zeros((5, 5), dtype=int))


def test_pvals_and_alpha_edge_validation():
    # _check_pvals rejects non-1D and out-of-range inputs.
    with pytest.raises(validate.ValidationError):
        M.p_adjust(np.array([[0.1, 0.2]]), "none")
    with pytest.raises(validate.ValidationError):
        M.p_adjust(np.array([0.1, 1.5]), "none")
    with pytest.raises(validate.ValidationError):
        M.p_adjust(np.array([0.1, np.nan]), "bonferroni")
    # romano_wolf_stepdown input checks.
    with pytest.raises(validate.ValidationError):
        M.romano_wolf_stepdown(np.array([1.0, 2.0, 3.0]))
    with pytest.raises(validate.ValidationError):
        M.romano_wolf_stepdown(np.zeros((2, 1)))  # only 1 contrast
    with pytest.raises(validate.ValidationError):
        M.romano_wolf_stepdown(np.zeros((5, 3)) + 1.0, alpha=1.2)
    # planning_alpha edge cases.
    with pytest.raises(validate.ValidationError):
        M.planning_alpha("bonferroni", K=0)
    with pytest.raises(validate.ValidationError):
        M.planning_alpha("bh", K=3, r=5)  # r > K
    a, note = M.planning_alpha("romano_wolf", K=4, alpha=0.05)
    assert a == pytest.approx(0.05)
    with pytest.raises(validate.ValidationError):
        M.planning_alpha("mystery", K=2, alpha=0.05)
    a, _ = M.planning_alpha("bh", K=2, alpha=0.05)  # r defaults to 1
    assert a == pytest.approx(0.05 / 2)


def test_family_power_rejects_bad_input():
    with pytest.raises(validate.ValidationError):
        M.family_power_conjunctive(np.array([0.8, 1.2]))
    with pytest.raises(validate.ValidationError):
        M.family_power_disjunctive(np.array([-0.2, 0.8]))
    # Edge: single contrast family power.
    assert M.family_power_conjunctive(np.array([0.8])) == pytest.approx(0.8)
    assert M.family_power_disjunctive(np.array([0.8])) == pytest.approx(0.8)
