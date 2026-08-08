"""Tests for `planning.py` — hypothesis modes M1-M5, precision/power targets,
variance-inflation, and the scipy nct fix (REVISION_PLAN.md §1.6).

The single most important checks are the textbook power values:
* M1 two-sided, delta/sigma = 0.5, alpha = 0.05 -> smallest J for 80% is 34
  (docs/tisca_v2_spec.md §1, verified in P1_math_verification.ipynb).
* The spec's own §4.5 case-study J values (reproduced against
  legacy/Paper_Experiments/DGP1_500_results.csv).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats
from scipy.stats import nct

from tisca import planning as P
from tisca import validate


def test_m1_textbook_power_and_J():
    """delta/sigma=0.5 two-sided needs J=34 for 80% power."""
    assert P.power_M1(30, 0.5, 1.0, 0.05) < 0.80
    assert P.power_M1(34, 0.5, 1.0, 0.05) >= 0.80
    assert P.power_M1(40, 0.5, 1.0, 0.05) >= 0.80
    # Exact textbook value from the spec.
    assert P.power_M1(34, 0.5, 1.0, 0.05) == pytest.approx(0.8078, abs=2e-4)
    assert P.required_J_power("M1", 0.5, 1.0, 0.80, 0.05) == 34


def test_nct_bug_fixed():
    """v1 used scipy.stats.t.cdf with ncp in the loc slot; v2 uses nct (§1.6).

    The (wrong) t.cdf and the (right) nct must differ, and our power function must
    equal the nct-based value at the true critical value ``t_{1-alpha/2, J-1}``.
    """
    from scipy.stats import t as st

    J, df = 34, 33
    ncp = np.sqrt(J) * 0.5 / 1.0
    crit = st.ppf(1 - 0.05 / 2, df=df)  # t_{0.975, 33}
    wrong_v1 = st.cdf(crit, df, loc=ncp)  # the v1 bug (ncp in loc slot)
    right = nct.cdf(crit, df=df, nc=ncp)
    # Document that the v1 bug disagrees with the correct nct call.
    assert wrong_v1 != right
    # Our power_M1 equals the nct-based two-sided power.
    pw = P.power_M1(J, 0.5, 1.0, 0.05)
    manual = 1 - (nct.cdf(crit, df=df, nc=ncp) - nct.cdf(-crit, df=df, nc=ncp))
    assert pw == pytest.approx(manual)
    # And it is NOT equal to the value obtained with the buggy t.cdf call.
    wrong_full = 1 - (st.cdf(crit, df, loc=ncp) - st.cdf(-crit, df, loc=ncp))
    assert abs(wrong_full - manual) > 1e-6
    # The published-v1 conclusion was preserved in sign/magnitude (small J effect).
    assert abs(wrong_full - manual) <= 0.01


def test_m2_directional_power():
    """M2 lower-is-better directional; power rises as J grows; rejects in lower tail."""
    p10 = P.power_M2(10, -0.5, 1.0, 0.05)
    p50 = P.power_M2(50, -0.5, 1.0, 0.05)
    assert 0 < p10 < p50 < 1.0
    assert P.power_M2(50, -0.5, 1.0, 0.05) > P.power_M2(10, -0.5, 1.0, 0.05)
    # A positive delta in the wrong direction has essentially alpha power.
    assert P.power_M2(50, 0.5, 1.0, 0.05) < 0.06


def test_m3_m4_margins():
    """M3 min-effect and M4 non-inferiority use the margin-shifted ncp."""
    # M4: delta well inside margin -> high power.
    assert P.power_M4(50, 0.0, 0.5, 1.0, 0.05) > 0.9
    # M4: delta at the margin boundary -> converges to alpha.
    assert P.power_M4(500, 0.5, 0.5, 1.0, 0.05) == pytest.approx(0.05, abs=0.02)
    # M3: true effect past the inner boundary -Delta yields high power.
    assert P.power_M3(50, -1.0, 0.5, 1.0, 0.05) > 0.9


def test_m5_equivalence_boundaries():
    """M5 TOST power is 0 at small J (feasibility), rises with J, alpha at |delta|=Delta."""
    assert P.power_M5_approx(10, 0.2, 0.5, 1.0, 0.05) == 0.0  # boundaries cross
    p40 = P.power_M5_approx(40, 0.2, 0.5, 1.0, 0.05)
    p80 = P.power_M5_approx(80, 0.2, 0.5, 1.0, 0.05)
    assert p40 < p80
    assert p80 > 0.8
    # |delta| = Delta converges to alpha from below.
    assert P.power_M5_approx(5000, 0.5, 0.5, 1.0, 0.10) == pytest.approx(0.10, abs=0.02)


def test_m5_exact_matches_approx_in_range():
    """The exact quadrature TOST power agrees with the known-σ approximation here."""
    for J in (40, 80, 150):
        approx = P.power_M5_approx(J, 0.2, 0.5, 1.0, 0.05)
        exact = P.power_M5_exact(J, 0.2, 0.5, 1.0, 0.05)
        assert exact == pytest.approx(approx, abs=0.02), (J, exact, approx)


def test_m5_exact_at_boundary_converges_to_alpha():
    assert P.power_M5_exact(200, 0.5, 0.5, 1.0, 0.05) == pytest.approx(0.05, abs=0.01)


def test_inflation_factors_match_spec():
    """Spec §3: J0=25->1.153, J0=50->1.099, J0=100->1.067."""
    assert P.inflate_std(1.0, 25) == pytest.approx(1.153, abs=1e-3)
    assert P.inflate_std(1.0, 50) == pytest.approx(1.099, abs=1e-3)
    assert P.inflate_std(1.0, 100) == pytest.approx(1.067, abs=1e-3)


def test_required_J_precision_targets():
    """MCSE target J = ceil((sigma/m)^2); half-width solvable by scan."""
    assert P.required_J_mcse(2.0, 0.2) == 100
    assert P.required_J_mcse(1.0, 0.1) == 100
    # half-width is no less than the MCSE-equivalent and >= 2
    hw = P.required_J_halfwidth(1.0, 0.1, 0.05)
    assert hw >= 2
    crit = 1.0  # t approx; check the inequality shrinks
    j_big = P.required_J_halfwidth(5.0, 0.2, 0.05)
    assert j_big > P.required_J_halfwidth(1.0, 0.2, 0.05)


def test_solve_and_combine():
    J = P.solve_J_scan(1.0, mode="M1", delta=0.5, target_power=0.8, alpha=0.05)
    assert J >= 34  # precision is default; power drives J to >=34
    assert P.combine_J([10, 20, 30]) == 30
    assert P.combine_J([10, 20, 30], J_max=25) == 25


def test_case_study_J_values_from_dgp1_csv():
    """Reproduce the spec's §4.5 worked example from the actual case-study data."""
    import os

    p = os.path.join(
        os.path.dirname(__file__), "..", "..", "legacy", "Paper_Experiments",
        "DGP1_500_results.csv",
    )
    if not os.path.exists(p):
        pytest.skip("legacy DGP1_500_results.csv not present")
    import pandas as pd

    df = pd.read_csv(p)
    cells = [
        ("mvbcf_tau_951", "wsbcf_tau_951", 0.015, (114, 177, 213)),
        ("mvbcf_tau_952", "wsbcf_tau_952", 0.015, (97, 149, 180)),
        ("mvbcf_pehe1", "bcf_pehe1", 0.5, (205, 316, 382)),
        ("mvbcf_pehe2", "bcf_pehe2", 0.5, (185, 285, 344)),
        ("mvbcf_tau_951", "mvbart_tau_951", 0.015, (49, 75, 90)),
        ("mvbcf_tau_952", "mvbart_tau_952", 0.015, (66, 101, 122)),
    ]
    for a, b, delta, expected in cells:
        A = df[a].to_numpy(float)
        B = df[b].to_numpy(float)
        D = A - B
        sd = float(np.std(D, ddof=1))
        j_raw = P.required_J_power("M1", delta, sd, 0.80, 0.05)
        j_adj = P.required_J_power("M1", delta, sd, 0.80, 0.05 / 6)
        sig_ub = P.inflate_std(sd, 50)
        j_both = P.required_J_power("M1", delta, sig_ub, 0.80, 0.05 / 6)
        assert j_raw == expected[0], (a, j_raw, expected)
        assert j_adj == expected[1], (a, j_adj, expected)
        assert j_both == expected[2], (a, j_both, expected)


def test_degenerate_sd_zero_terminates():
    """Rule 8.5: sd=0 -> power target treated as satisfied, no infinite loop."""
    assert P.required_J_power("M1", 0.5, 0.0, 0.80, 0.05) == 2
    assert P.solve_J_scan(0.0, mode="M1", delta=0.5, target_power=0.8, alpha=0.05) == 2
    assert P.required_J_power("M2", -0.5, 0.0, 0.80, 0.05) == 2


def test_m5_infeasible_rejected():
    """Rule 8.6: |delta| >= Delta is genuinely infeasible and must be refused."""
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M5", 0.5, 1.0, 0.8, 0.05, margin=0.5)
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M5", -0.6, 1.0, 0.8, 0.05, margin=0.5)
    # But a zero-power artifact at small J is not infeasible: solver continues.
    assert P.required_J_power("M5", 0.2, 0.5, 0.8, 0.05, margin=0.5) > 10


def test_invalid_mode_and_target_rejected():
    with pytest.raises(validate.ValidationError):
        P.power_function("M9", 50, 0.5, 1.0, alpha=0.05)
    with pytest.raises(validate.ValidationError):
        P.solve_J_scan(1.0, mode="M1", delta=0.5)  # no active target
    with pytest.raises(validate.ValidationError):
        P.required_J_mcse(1.0, 0.0)  # non-positive target
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M3", 0.5, 1.0, 0.8, 0.05)  # M3 needs margin


def test_sigma_zero_branches_for_all_modes():
    """sigma=0 -> deterministic outcome for every mode (rule 8.5).

    The M1 null case is the one this test previously got wrong. It asserted power
    1.0 at ``delta = 0, sigma = 0`` -- i.e. that a contrast whose per-replication
    difference is identically zero is detected with certainty. Nothing is
    detectable there: the statistic is 0/0 and the test never rejects. M1 now
    reports ``alpha`` when the planning alternative sits inside the null, which is
    the convention M2-M4 already used on the lines below.
    """
    assert P.power_function("M1", 10, 0.5, 0.0, alpha=0.05) == 1.0
    assert P.power_function("M1", 10, 0.0, 0.0, alpha=0.05) == pytest.approx(0.05)
    assert P.power_function("M2", 10, -0.5, 0.0, alpha=0.05) == 1.0
    assert P.power_function("M2", 10, 0.5, 0.0, alpha=0.05) == pytest.approx(0.05)
    assert P.power_function("M3", 10, -1.0, 0.0, alpha=0.05, margin=0.5) == 1.0
    assert P.power_function("M4", 10, 0.0, 0.0, alpha=0.05, margin=0.5) == 1.0
    assert P.power_function("M5", 10, 0.2, 0.0, alpha=0.05, margin=0.5) == 1.0
    # A planning alternative outside the margin is infeasible and must be rejected.
    with pytest.raises(validate.ValidationError):
        P.power_function("M5", 10, 0.6, 0.0, alpha=0.05, margin=0.5)


def test_power_function_dispatch_edge_cases():
    """M3/M4/M5 dispatch, exact flag, and unknown-mode rejection."""
    assert P.power_function("M3", 50, -1.0, 1.0, margin=0.5, alpha=0.05) == pytest.approx(
        P.power_M3(50, -1.0, 0.5, 1.0, 0.05)
    )
    assert P.power_function("M4", 50, 0.0, 1.0, margin=0.5, alpha=0.05) == P.power_M4(50, 0.0, 0.5, 1.0, 0.05)
    assert P.power_function("M5", 50, 0.2, 1.0, margin=0.5, alpha=0.05, exact=True) == pytest.approx(
        P.power_M5_exact(50, 0.2, 0.5, 1.0, 0.05)
    )
    with pytest.raises(validate.ValidationError):
        P.power_function("M6", 50, 0.5, 1.0, alpha=0.05)


def test_halfwidth_solver_meets_the_target_it_returns():
    """The returned J must actually achieve the half-width, and be the smallest.

    Regression test. The solver used to short-circuit to ``J = 1`` whenever
    ``sigma < halfwidth``, which compares the wrong quantity: the achieved
    half-width is ``t_{1-a/2, J-1} sigma / sqrt(J)``, and at ``J = 2`` that
    carries a factor of ``t_{0.975,1}/sqrt(2) ~ 9``. sigma = 0.9 with a target
    of 1.0 was returned as J = 1 with a true half-width of 11.4.
    """
    def achieved(sigma, j, alpha=0.05):
        return float(stats.t.ppf(1 - alpha / 2, df=j - 1) * sigma / np.sqrt(j))

    for sigma, h in [(0.9, 1.0), (1.0, 5.0), (0.5, 0.6), (2.0, 3.0), (50.0, 0.5), (1.0, 0.1)]:
        j = P.required_J_halfwidth(sigma, h, 0.05)
        assert j >= 2, (sigma, h, j)
        assert achieved(sigma, j) <= h, (sigma, h, j, achieved(sigma, j))
        if j > 2:                              # and no smaller J would do
            assert achieved(sigma, j - 1) > h, (sigma, h, j)

    with pytest.raises(validate.ValidationError):
        P.required_J_halfwidth(1.0, 0.0, 0.05)


def test_required_J_power_edge_cases():
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M1", 0.5, 1.0, 1.5, 0.05)   # target_power > 1
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M1", 0.5, -1.0, 0.8, 0.05)  # sigma < 0
    with pytest.raises(validate.ValidationError):
        P.required_J_power("M1", None, 1.0, 0.8, 0.05)  # delta None
    # Unreachable within a tiny budget -> returns J_max (rule 8.6, not infeasible).
    j = P.required_J_power("M1", 0.01, 1.0, 0.95, 0.05, J_max=20)
    assert j == 20


def test_solve_J_scan_default_alpha_and_combine_empty():
    j_hw = P.solve_J_scan(2.0, mode="M1", delta=0.5, halfwidth=0.5)
    assert j_hw >= 2
    j_pw = P.solve_J_scan(2.0, mode="M1", delta=0.5, target_power=0.8)
    assert j_pw >= 2
    assert P.combine_J([]) == 0


def test_degenerate_sigma_never_reports_certainty_or_a_testless_J():
    """F9: the two zero-variance escapes that used to produce nonsense.

    1. ``required_J_mcse(0, m)`` returned 1. One replication satisfies an MCSE
       target vacuously but leaves ``df = 0``, so the paired t the plan is built
       on does not exist; the other two solvers already returned 2.
    2. M1 reported power 1.0 at ``sigma = 0`` for ANY delta, so a contrast with no
       detectable signal at all was published as fully powered.

    The solver must still terminate at the smallest admissible J (rule 8.5): the
    honest power is not allowed to send a degenerate cell to J_max.
    """
    assert P.required_J_mcse(0.0, 0.01) == 2
    assert P.required_J_mcse(1e-9, 0.01) == 2          # and never below 2
    assert P.required_J_halfwidth(0.0, 0.01, 0.05) == 2
    assert P.required_J_power("M1", 0.0, 0.0, 0.80, 0.05) == 2
    assert P.required_J_power("M1", 0.5, 0.0, 0.80, 0.05) == 2
    assert P.solve_J_scan(0.0, mode="M1", delta=0.0, target_mcse=0.01,
                          target_power=0.80, alpha=0.05) == 2

    assert P.power_M1(10, 0.0, 0.0, 0.05) == pytest.approx(0.05)
    assert P.power_M1(10, 0.5, 0.0, 0.05) == 1.0
