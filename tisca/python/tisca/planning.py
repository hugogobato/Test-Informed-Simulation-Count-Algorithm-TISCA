"""Design targets and hypothesis modes — `planning.py`.

Implements §1 and §2 of `docs/tisca_v2_spec.md`. All power functions use the
**noncentral t with ``J-1`` degrees of freedom** for the mean-of-paired-
contrast statistic ``T = sqrt(J)(D_bar - theta0)/s_D ~ t_{J-1, ncp}`` under a
planning alternative, via ``scipy.stats.nct``.

This is the module that fixes the v1 ``stats.t.cdf(crit, df, ncp)`` bug
(revision plan §1.6): in SciPy the third positional argument of ``t.cdf`` is
``loc``, not a noncentrality parameter; the noncentral t is ``nct``.

Design targets (§2):
* precision: smallest ``J`` with ``s_D/sqrt(J) <= m`` (MCSE) or
  ``t_{1-alpha/2, J-1} s_D/sqrt(J) <= h`` (half-width, monotone scan);
* power:     smallest ``J`` with ``Pr(reject | theta = delta) >= 1 - beta``,
  computed at the **adjusted** level ``alpha_adj`` and the same sidedness as the
  final test;
* combination: ``J_final = min( max over contrasts and active targets , J_max )``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import nct, t as _t, norm as _norm, chi2 as _chi2

from . import validate

__all__ = [
    "nct_cdf",
    "power_M1",
    "power_M2",
    "power_M3",
    "power_M4",
    "power_M5_approx",
    "power_M5_exact",
    "power_function",
    "required_J_mcse",
    "required_J_halfwidth",
    "required_J_power",
    "required_J",
    "solve_J_scan",
    "inflate_std",
    "combine_J",
]

# Problem sizes for which brute-force scanning is fine.
_DEFAULT_J_MAX = 100_000
_SCAN_STEP = 4


# --------------------------------------------------------------------------- #
# Per-mode power functions (spec §1). All return Pr(reject | theta = delta).   #
# --------------------------------------------------------------------------- #

def _ncp(J, delta, sigma):
    return float(np.sqrt(J) * delta / sigma)


def nct_cdf(x, df, nc):
    """Noncentral-t CDF, NaN-safe in the far tail.

    ``scipy.stats.nct.cdf`` returns ``NaN`` instead of an underflowed ~0 once
    ``|nc|`` gets large (e.g. ``nct.cdf(-1.96, 499, 17.9)``), which silently
    poisons two-sided power and can make a monotone ``J`` search terminate at
    the wrong point. R's ``pt(q, df, ncp)`` returns the correct value, so the
    R port never saw this; the R<->Python parity run flagged 28 such configs.

    Where SciPy returns a non-finite value, fall back to the Johnson-Welch
    normal approximation
    ``Phi( (x(1 - 1/(4 df)) - nc) / sqrt(1 + x^2/(2 df)) )``,
    which is accurate to a few 1e-4 in the region where SciPy fails and is
    exactly the ~0 / ~1 limit deep in the tail.
    """
    x = np.asarray(x, dtype=float)
    df = np.asarray(df, dtype=float)
    nc = np.asarray(nc, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        out = np.asarray(nct.cdf(x, df=df, nc=nc), dtype=float)
        bad = ~np.isfinite(out)
        if np.any(bad):
            xb, dfb, ncb = np.broadcast_arrays(x, df, nc)
            approx = _norm.cdf(
                (xb * (1.0 - 1.0 / (4.0 * dfb)) - ncb)
                / np.sqrt(1.0 + xb ** 2 / (2.0 * dfb))
            )
            out = np.where(bad, approx, out)
    return out if out.ndim else float(out)


def power_M1(J, delta, sigma, alpha) -> float:
    """M1 two-sided equality ``H0: theta = 0``.

    ``beta1(J) = 1 - [F(t_{+}) - F(t_{-})]`` with
    ``t_+ = t_{1-alpha/2, J-1}``, ``ncp = sqrt(J) delta/sigma``.
    """
    J, delta, sigma = float(J), float(delta), float(sigma)
    df = J - 1.0
    crit = _t.ppf(1.0 - alpha / 2.0, df=df)
    if sigma == 0.0:
        # `return 1.0` unconditionally was wrong at delta = 0, which is exactly the
        # case a degenerate contrast produces: when D_j is identically 0 the
        # statistic is 0/0 and nothing is detectable, yet the caller was told the
        # power was 1.0. Branch on the alternative, matching M2-M5 below, which
        # already report `alpha` when the planning alternative sits inside the null.
        return 1.0 if delta != 0.0 else alpha
    nc = _ncp(J, delta, sigma)
    return float(1.0 - (nct_cdf(crit, df, nc) - nct_cdf(-crit, df, nc)))


def power_M2(J, delta, sigma, alpha) -> float:
    """M2 directional superiority, lower is better ``H0: theta >= 0``, ``delta < 0``.

    ``beta2(J) = F(t_{alpha, J-1})`` with ``ncp = sqrt(J) delta/sigma``.
    """
    J, delta, sigma = float(J), float(delta), float(sigma)
    df = J - 1.0
    crit = _t.ppf(alpha, df=df)
    if sigma == 0.0:
        return 1.0 if delta < 0 else alpha
    nc = _ncp(J, delta, sigma)
    return float(nct_cdf(crit, df, nc))


def power_M3(J, delta, margin, sigma, alpha) -> float:
    """M3 minimum-effect superiority ``H0: theta >= -Delta``.

    ``beta3(J) = F(t_{alpha, J-1})`` with ``ncp = sqrt(J)(delta + Delta)/sigma``.
    """
    J, delta, margin, sigma = float(J), float(delta), float(margin), float(sigma)
    df = J - 1.0
    crit = _t.ppf(alpha, df=df)
    if sigma == 0.0:
        return 1.0 if delta < -margin else alpha
    nc = _ncp(J, delta + margin, sigma)
    return float(nct_cdf(crit, df, nc))


def power_M4(J, delta, margin, sigma, alpha) -> float:
    """M4 non-inferiority ``H0: theta >= Delta`` (test difference against upper edge).

    ``beta4(J) = F(t_{alpha, J-1})`` with ``ncp' = sqrt(J)(delta - Delta)/sigma``.
    """
    J, delta, margin, sigma = float(J), float(delta), float(margin), float(sigma)
    df = J - 1.0
    crit = _t.ppf(alpha, df=df)
    if sigma == 0.0:
        return 1.0 if delta < margin else alpha
    nc = _ncp(J, delta - margin, sigma)
    return float(nct_cdf(crit, df, nc))


def power_M5_approx(J, delta, margin, sigma, alpha) -> float:
    """M5 equivalence (TOST) known-sigma approximation, spec §1.

    ``beta5(J) = F(b_hi) - F(b_lo)`` where (central t, ``F = t_{J-1}`` CDF)
    ``b_lo = (-Delta + t_{1-alpha,J-1} s/sqrt(J) - delta) sqrt(J)/s``
    ``b_hi = ( Delta - t_{1-alpha,J-1} s/sqrt(J) - delta) sqrt(J)/s``.
    This is the object used for planning (accurate to <=0.009 in power over the
    planning range); the exact quadrature form is ``power_M5_exact``.
    """
    J, delta, margin, sigma = float(J), float(delta), float(margin), float(sigma)
    df = J - 1.0
    crit = _t.ppf(1.0 - alpha, df=df)  # upper tail, one-sided per TOST arm
    sJ = np.sqrt(J)
    if sigma == 0.0:
        return 1.0 if abs(delta) < margin else (alpha if abs(delta) == margin else 0.0)
    b_lo = (-margin + crit * sigma / sJ - delta) * sJ / sigma
    b_hi = (margin - crit * sigma / sJ - delta) * sJ / sigma
    power = _t.cdf(b_hi, df=df) - _t.cdf(b_lo, df=df)
    # When the two TOST acceptance boundaries cross (b_lo > b_hi), the rejection
    # triangle is empty and power is exactly 0 at this J (a small-J feasibility
    # artifact, spec §1 M5), never negative.
    return float(max(0.0, power))


def power_M5_exact(J, delta, margin, sigma, alpha, n_quad: int = 4096) -> float:
    """M5 TOST exact power via quadrature over the sample ``s_D`` distribution.

    ``beta5_exact(J) = E_{s_D}[ Phi((m(s_D)-delta)J^{1/2}/sigma)
                                   - Phi((-m(s_D)-delta)J^{1/2}/sigma) ]``
    where ``m(s) = max(Delta - t_{1-alpha,J-1} s/sqrt(J), 0)`` and
    ``s_D ~ sigma sqrt(chi^2_{J-1}/(J-1))``. The averaged value sits slightly
    below the known-sigma approximation (the random ``s_D`` enters the rejection
    boundary), and at ``|delta| = Delta`` it converges to ``alpha`` from below.
    """
    J, delta, margin, sigma = float(J), float(delta), float(margin), float(sigma)
    df = J - 1.0
    crit = _t.ppf(1.0 - alpha, df=df)
    sJ = np.sqrt(J)
    # Gauss-Legendre quadrature over s in [0, s_max].
    # s_D ~ sigma*sqrt(W/(J-1)), W ~ chi2_{J-1}, so the density of s_D is
    #   f(s) = f_W((J-1)(s/sigma)^2) * 2*(J-1)*s/sigma^2.
    s_max = sigma * np.sqrt(_chi2.ppf(1 - 1e-12, df=df) / df)
    z, w = np.polynomial.legendre.leggauss(n_quad)
    s = s_max * (z + 1.0) / 2.0
    jac = s_max / 2.0
    w_arg = df * (s / sigma) ** 2
    pdf_w = np.exp(_chi2.logpdf(w_arg, df=df))
    f_s = pdf_w * 2.0 * df * s / sigma ** 2
    m = np.maximum(margin - crit * s / sJ, 0.0)
    inner = _norm.cdf((m - delta) * sJ / sigma) - _norm.cdf((-m - delta) * sJ / sigma)
    integrand = inner * f_s * jac
    return float(np.sum(w * integrand))


# --------------------------------------------------------------------------- #
# Generic dispatcher                                                        #
# --------------------------------------------------------------------------- #

def power_function(mode, J, delta, sigma, *, margin=None, alpha, exact=False) -> float:
    """Compute power at ``(J, delta, sigma, alpha)`` for any mode M1..M5."""
    validate.validate_mode_targets(mode, None)
    mode = mode.upper()
    if mode in ("M1",):
        return power_M1(J, delta, sigma, alpha)
    if mode in ("M2",):
        return power_M2(J, delta, sigma, alpha)
    if mode in ("M3",):
        return power_M3(J, delta, margin, sigma, alpha)
    if mode in ("M4",):
        return power_M4(J, delta, margin, sigma, alpha)
    if mode in ("M5",):
        validate.validate_margin_feasibility(delta, margin, mode="M5")
        if exact:
            return power_M5_exact(J, delta, margin, sigma, alpha)
        return power_M5_approx(J, delta, margin, sigma, alpha)
    raise validate.ValidationError(f"Unknown mode {mode!r}.")


def _mode_requires_margin(mode):
    return mode.upper() in ("M3", "M4", "M5")


# --------------------------------------------------------------------------- #
# J solvers                                                                  #
# --------------------------------------------------------------------------- #

def required_J_mcse(sigma, target_mcse) -> int:
    """Smallest ``J`` with ``sigma/sqrt(J) <= target_mcse`` (MCSE target, §2.1).

    A degenerate ``sigma == 0`` returns 2, not 1. One replication satisfies an MCSE
    target vacuously, but it leaves ``df = J - 1 = 0``, so the paired t that every
    downstream target and the final test are defined on does not exist. Two is the
    smallest admissible size and is what ``required_J_halfwidth`` and
    ``required_J_power`` already return in the same situation; returning 1 here
    made ``solve_J_scan`` capable of planning a run no test could be run on.
    """
    if target_mcse <= 0:
        raise validate.ValidationError("target_mcse must be positive.")
    if sigma < 0:
        raise validate.ValidationError("sigma must be non-negative.")
    if sigma == 0:
        return 2
    return max(2, int(math.ceil((sigma / target_mcse) ** 2)))


def required_J_halfwidth(sigma, halfwidth, alpha, J_max=_DEFAULT_J_MAX) -> int:
    """Smallest ``J`` with ``t_{1-alpha/2,J-1} sigma/sqrt(J) <= halfwidth``.

    Monotone in ``J`` for fixed ``sigma``; solved by scanning.
    """
    if halfwidth <= 0:
        raise validate.ValidationError("halfwidth must be positive.")
    if sigma <= 0:
        return 2
    # The smallest admissible J is 2 (df = J-1 >= 1); at J = 1 the t quantile is
    # undefined. `sigma < halfwidth -> J = 1` was a short circuit on the WRONG
    # comparison: it is the achieved half-width t_{1-a/2,J-1} sigma/sqrt(J), not
    # sigma itself, that must fall under the target, and at J = 2 that carries a
    # factor t_{0.975,1}/sqrt(2) ~ 9. E.g. sigma = 0.9, h = 1.0 returned J = 1
    # with a true half-width of 11.4.
    #
    # t_{1-a/2,J-1} sigma/sqrt(J) is strictly decreasing in J, so bisect on the
    # first J that meets the target instead of scanning (the scan was O(J_max)
    # with a redundant second pass over the same range).
    def _hw(j):
        return float(_t.ppf(1.0 - alpha / 2.0, df=j - 1.0) * sigma / np.sqrt(j))

    if _hw(2) <= halfwidth:
        return 2
    lo, hi = 2, 4
    while hi < J_max and _hw(hi) > halfwidth:
        lo, hi = hi, min(hi * 2, J_max)
    if _hw(hi) > halfwidth:
        return J_max        # target not reachable within the budget
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if _hw(mid) <= halfwidth:
            hi = mid
        else:
            lo = mid
    return hi


def required_J_power(
    mode,
    delta,
    sigma,
    target_power,
    alpha,
    *,
    margin=None,
    J_max=_DEFAULT_J_MAX,
    exact=False,
) -> int:
    """Smallest ``J`` with ``Pr(reject | theta=delta) >= target_power``.

    Handles the degenerate ``sigma == 0`` case (rule 8.5): with zero contrast
    variance there is no estimation error, so the power target is considered
    *satisfied* at the smallest admissible ``J`` and the planner must not loop.

    Also handles the M5 small-``J`` zero-power artifact by continuing to scan
    (rule 8.6): a zero-power evaluation at one ``J`` is a ``J``-artifact, not
    global infeasibility (which is handled by ``validate_margin_feasibility``).
    """
    validate.validate_mode_targets(mode, None)
    mode = mode.upper()
    if _mode_requires_margin(mode) and margin is None:
        raise validate.ValidationError(f"Mode {mode} requires a margin.")
    if target_power <= 0 or target_power >= 1:
        raise validate.ValidationError("target_power must lie in (0, 1).")
    if sigma < 0 or delta is None:
        raise validate.ValidationError("sigma must be non-negative and delta given.")
    J_max = int(J_max)

    if sigma == 0.0:
        # Zero variance: the target is vacuous, because more replications cannot
        # reduce an estimation error that is already exactly zero. Rule 8.5 requires
        # the solver to TERMINATE here ("treat power targets as satisfied") rather
        # than escalate, so this returns the smallest admissible J unconditionally.
        #
        # It deliberately does NOT gate on the achieved power. The power functions
        # now report the honest rejection probability at sigma = 0, which is `alpha`
        # when the planning alternative lies inside the null; gating on that would
        # send every degenerate null contrast to J_max, i.e. spend the entire budget
        # on a contrast where no J can ever help. Rule 8.5's other half -- report the
        # degenerate cell -- is handled by the caller (`procedure.py` flags it in
        # `degenerate_contrasts`).
        return 2

    critical = target_power
    # Fast geometric/scan search. Power is monotone increasing in J for M1-M4 and
    # (eventually) for M5 after feasibility.
    lo, hi = 2, 2
    while power_function(mode, hi, delta, sigma, margin=margin, alpha=alpha, exact=exact) < critical and hi < J_max:
        lo = hi
        hi = min(hi * 2, J_max)
        if hi == lo:
            break
    if power_function(mode, hi, delta, sigma, margin=margin, alpha=alpha, exact=exact) < critical:
        return J_max  # target not reached within budget (rule 8.6, distinct message)
    # Binary search in [lo, hi].
    best_hi = hi
    while lo < hi:
        mid = (lo + hi) // 2
        if power_function(mode, mid, delta, sigma, margin=margin, alpha=alpha, exact=exact) >= critical:
            hi = mid
        else:
            lo = mid + 1
    return max(2, lo)


def solve_J_scan(
    sigma,
    *,
    mode,
    delta,
    target_mcse=None,
    halfwidth=None,
    alpha=None,
    target_power=None,
    margin=None,
    J_max=_DEFAULT_J_MAX,
    exact=False,
) -> int:
    """Solve for ``J`` over whichever active precision/power targets are given.

    Returns the required ``J`` (the max across active targets), per §2.3.
    ``alpha`` is the **adjusted** alpha (multiplicity-aware) for the power layer.
    """
    candidates = []
    if target_mcse is not None:
        candidates.append(required_J_mcse(sigma, target_mcse))
    if halfwidth is not None:
        if alpha is None:
            alpha = 0.05
        candidates.append(required_J_halfwidth(sigma, halfwidth, alpha, J_max=J_max))
    if target_power is not None:
        if alpha is None:
            alpha = 0.05
        candidates.append(
            required_J_power(
                mode, delta, sigma, target_power, alpha, margin=margin, J_max=J_max, exact=exact
            )
        )
    if not candidates:
        raise validate.ValidationError(
            "No active target supplied to solve_J_scan: specify at least one of "
            "target_mcse, halfwidth, or target_power."
        )
    return int(min(max(candidates), J_max))


def required_J(
    sd_pilot, J0, *, gamma=0.20, mode, delta, target_mcse=None,
    halfwidth=None, alpha=None, target_power=None, margin=None,
    J_max=_DEFAULT_J_MAX, exact=False,
) -> tuple[int, float]:
    """Plan ``J`` from a pilot SD, applying the variance-uncertainty inflation §3.

    Returns ``(J, sigma_UB)`` where ``sigma_UB = sd_pilot * sqrt((J0-1)/chi2(gamma, J0-1))``.
    """
    sigma_UB = validate.validate_pilot_samples(sd_pilot, J0, gamma=gamma)
    J = solve_J_scan(
        sigma_UB,
        mode=mode,
        delta=delta,
        target_mcse=target_mcse,
        halfwidth=halfwidth,
        alpha=alpha,
        target_power=target_power,
        margin=margin,
        J_max=J_max,
        exact=exact,
    )
    return J, sigma_UB


def inflate_std(sd, J0, gamma=0.20) -> float:
    """Variance-uncertainty inflation factor ``sd * sqrt((J0-1)/chi2(gamma, J0-1))``."""
    return validate.validate_pilot_samples(sd, J0, gamma=gamma)


def combine_J(j_values, J_max=_DEFAULT_J_MAX) -> int:
    """``J_final = min(max(J) over comparisons and active targets, J_max)`` (§2.3)."""
    if not j_values:
        return 0
    return int(min(max(j_values), J_max))
