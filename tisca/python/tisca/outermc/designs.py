"""Designs D1..D6 for the outer-MC harness (E1).

Each design wraps the *whole* procedure (pilot -> planning -> confirmatory ->
final test) for one outer repetition and is evaluated over many independent
repetitions to obtain unconditional operating characteristics (P3-T2).

The designs are implemented vectorised over the ``R`` outer repetitions (arrays
of shape ``(R, ...)``). A design receives a ``draws`` dict with two blocks of
paired losses:

  draws["pilot"]   shape (R, J0, 2)   independent-seed pilot (D1/D4/D5 use it)
  draws["conf"]    shape (R, Jmax, 2) confirmatory pool (one test run on the
                                      first ``J`` rows; extra rows are unused)

For the adaptive designs (D2, D3) the caller sets ``Jmax`` = the budgeted
maximum and the design starts its accumulation over ``conf``'s rows from ``J0``,
reusing the leading rows exactly as the internal-pilot algorithms do.

Designs (mapping to REVISION_PLAN.md P3-T2 and tisca_v2_spec.md):
  D1  fixed-J two-stage from an independent pilot (RAW pilot sigma_D, unadjusted alpha)
  D2  internal-pilot re-estimation (Algorithm 2; pilot rows reused)
  D3  TISCA v1: unpaired Welch on paired data (the v1 defect being studied).
      MEASURED, not assumed: with positively correlated methods the v1 test is
      severely CONSERVATIVE, not liberal. At delta = 0.5, J0 = 50, normal losses
      its unconditional Type I error falls 0.053 -> 0.022 -> 0.002 -> 0.000 as
      rho goes 0 -> 0.3 -> 0.6 -> 0.9, while D4 holds ~0.05 throughout. The
      case study's own paired correlations are r = 0.58-0.80 (plan 1.1), so
      that is the operative regime: v1 did not inflate the false-positive rate,
      it threw away power and oversized J.
  D4  TISCA v2: two-stage (Algorithm 1), sigma_D_UB + adjusted alpha
  D5  paired fixed-precision (MCSE / half-width precision target)
  D6  oracle fixed-J (true sigma known), one-shot test

Every design returns ``(result, meta)`` where ``result`` is a ``Result`` array
with columns ``[T, p, J, theta_hat, s_D, capped]`` per outer repetition, and
``meta`` carries design-specific diagnostics.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from ..planning import nct_cdf as _nct_cdf

# --------------------------------------------------------------------------- #
# shared planning math (vectorised)
# --------------------------------------------------------------------------- #

def _inflate(s_D, J0, gamma):
    """Variance-uncertainty inflation sigma_hat_UB (tisca_v2_spec.md section 3)."""
    q = stats.chi2.ppf(gamma, J0 - 1)
    return s_D * np.sqrt((J0 - 1) / q)


def m1_power(J, delta, sigma_D, alpha):
    """M1 two-sided noncentral-t power at planning alternative ``delta``."""
    J = np.asarray(J, dtype=float)
    ncp = np.sqrt(J) * delta / sigma_D
    crit = stats.t.ppf(1.0 - alpha / 2.0, J - 1)
    # Route through the NaN-safe wrapper: scipy's nct.cdf returns NaN instead of
    # an underflowed ~0 at large |ncp|, and a NaN power silently fails the
    # `pw >= target` test, pushing the solver to a larger J than needed.
    return 1.0 - (_nct_cdf(crit, J - 1, ncp) - _nct_cdf(-crit, J - 1, ncp))


def m2_directional_power(J, delta, sigma_D, alpha):
    """M2 one-sided power (lower tail, lower-is-better)."""
    ncp = np.sqrt(J) * delta / sigma_D
    return _nct_cdf(stats.t.ppf(alpha, J - 1), J - 1, ncp)


def _power_for_mode(mode, J, delta, sigma_D, alpha):
    if mode == 1:
        return m1_power(J, delta, sigma_D, alpha)
    if mode == 2:
        return m2_directional_power(J, delta, sigma_D, alpha)
    raise ValueError(f"planning power for mode {mode} not implemented in the harness")


def solve_J(power_target, delta, sigma_D, alpha, mode=1, Jmin=4, Jmax=10_000):
    """Smallest J with power(J) >= power_target (per repetition if sigma_D is (R,)).

    Returns ``(J, within_cap)``. Power is monotone in J for the M1/M2 modes, so
    this is a vectorised bisection: ~``log2(Jmax)`` noncentral-t evaluations per
    repetition instead of the ``Jmax`` of a dense grid. At the default
    ``R = 5000`` and ``Jmax = 1000`` that is 5e4 rather than 5e6 nct calls per
    cell -- the difference between an E1 module that runs in minutes and one
    that does not fit its notebook budget. The answer is identical (both return
    the exact smallest integer J).
    """
    aa = alpha
    sd = np.atleast_1d(np.asarray(sigma_D, dtype=float))
    Jmin, Jmax = int(Jmin), int(Jmax)

    def pw(j):
        return np.asarray(_power_for_mode(mode, np.asarray(j, dtype=float), delta, sd, aa))

    lo = np.full(sd.shape, Jmin, dtype=int)
    at_min = pw(lo) >= power_target
    hi = np.full(sd.shape, Jmax, dtype=int)
    any_ok = pw(hi) >= power_target
    # Bisect on the monotone predicate; rows already satisfied at Jmin collapse
    # to Jmin, rows never satisfied are reported at Jmax with within_cap False.
    active = any_ok & ~at_min
    lo_a = np.where(active, Jmin, Jmin)
    hi_a = np.where(active, Jmax, Jmin)
    while np.any(hi_a - lo_a > 1):
        mid = (lo_a + hi_a) // 2
        ok_mid = pw(mid) >= power_target
        hi_a = np.where(active & ok_mid, mid, hi_a)
        lo_a = np.where(active & ~ok_mid, mid, lo_a)
        if not np.any(active & (hi_a - lo_a > 1)):
            break
    J = np.where(at_min, Jmin, np.where(any_ok, hi_a, Jmax))
    return J.astype(int), any_ok


def paired_t_stats(D, J_used):
    """Vectorised one-sample t over the first ``J_used`` rows of each rep.

    ``D`` is (R, Jmax); ``J_used`` is (R,) int. Returns ``(T, Dbar, s_D)`` (R,).
    """
    R, Jmax = D.shape
    take = np.clip(J_used - 1, 0, Jmax - 1).astype(int)
    mask = np.arange(Jmax)[None, :] < J_used[:, None]
    Dm = np.where(mask, D, 0.0)
    CS = np.cumsum(Dm, axis=1)
    CSS = np.cumsum(Dm * Dm, axis=1)
    cs = CS[np.arange(R), take]
    css = CSS[np.arange(R), take]
    Dbar = cs / J_used
    var = np.clip((css - J_used * Dbar * Dbar) / (J_used - 1), 0, None)
    s_D = np.sqrt(var)
    s_D = np.where(s_D > 0, s_D, np.nan)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        T = np.sqrt(J_used) * Dbar / s_D
    return T, Dbar, s_D


def paired_pvalue(T, J_used):
    """Two-sided p-value for the paired one-sample t."""
    return 2.0 * stats.t.sf(np.abs(T), J_used - 1)


class Result(np.ndarray):
    """Per-rep output of one design (R, ncols). Columns ``T, p, J, theta_hat, s_D, capped``."""
    _COLUMNS = ["T", "p", "J", "theta_hat", "s_D", "capped"]

    def __new__(cls, arr):
        obj = np.asarray(arr).view(cls)
        return obj

    def col(self, name):
        return self[:, self._COLUMNS.index(name)]


def _emit(T, p, J, Dbar, s_D, capped):
    return Result(np.stack([T, p, J, Dbar, s_D, capped], axis=1))


# --------------------------------------------------------------------------- #
# designs
# --------------------------------------------------------------------------- #

def design_oracle(config, draws, alpha, mode=1, **k):
    """D6: oracle fixed-J, true sigma known. ``config.sigma_D`` must be set.

    ``J`` is either fixed (``config.fixed_J``) or solved from the oracle sigma and
    the planning alternative ``delta`` at ``config.power_target``.
    """
    conf = draws["conf"]
    R = conf.shape[0]
    if config.get("fixed_J") is not None:
        J = int(config["fixed_J"])
    else:
        J, _ = solve_J(config["power_target"], config["delta"], config["sigma_D"],
                       alpha, mode=mode)  # scalar sigma_D -> single J
        J = int(np.max(J))
    sub = conf[:, :J, :]
    D = sub[:, :, 0] - sub[:, :, 1]
    T, Dbar, s_D = paired_t_stats(D, np.full(R, J))
    p = paired_pvalue(T, np.full(R, J))
    return _emit(T, p, np.full(R, J, dtype=float), Dbar, s_D,
                 np.zeros(R)), dict(J_oracle=J, pilot_reused=False)


def _plan_power_J(s_p, J0, config, aa, mode):
    """Shared planning step for the power-target two-stage designs (D1/D4)."""
    Jmax = config.get("Jmax", 100_000)
    infl = config.get("gamma", 0.20)
    sigma_eff = _inflate(s_p, J0, infl)
    J, ok = solve_J(config["power_target"], config["delta"], sigma_eff,
                    aa, mode=mode, Jmax=Jmax)
    J = np.clip(J, 1, Jmax)
    return J.astype(int), sigma_eff, ok


def design_twostage(config, draws, alpha, mode=1, **k):
    """D4 (default) and D1: two-stage, independent pilot, one confirmatory run.

    ``D4`` plans from ``sigma_D_UB`` at the adjusted alpha; ``D1`` uses the raw
    pilot ``s_D`` and the unadjusted alpha (a simpler plugin two-stage).
    which is selected via ``config["inflate"]`` / ``config["alpha_adj"]``.
    """
    J0 = config.get("J0", 50)
    Jmax = config.get("Jmax", 100_000)
    inflate = config.get("inflate", True)
    aa = config.get("alpha_adj", alpha)
    pilot = draws["pilot"]
    Dp = pilot[:, :, 0] - pilot[:, :, 1]
    s_p = Dp.std(axis=1, ddof=1)
    s_p = np.where(np.isclose(s_p, 0), 1e-12, s_p)
    if inflate:
        J, sigma_eff, ok = _plan_power_J(s_p, J0, config, aa, mode)
    else:
        J, ok = solve_J(config["power_target"], config["delta"], s_p,
                        aa, mode=mode, Jmax=Jmax)
        sigma_eff = s_p
        J = np.clip(J, 1, Jmax)
    conf = draws["conf"]
    D = conf[:, :, 0] - conf[:, :, 1]
    T, Dbar, s_D = paired_t_stats(D, J)
    p = paired_pvalue(T, J)
    capped = np.where((J >= Jmax) & ~ok, 1.0, 0.0)
    return _emit(T, p, J.astype(float), Dbar, s_D, capped), dict(
        pilot_sD=s_p, sigma_eff=sigma_eff, J_planned=J, pilot_reused=False, ok=ok)


def design_internal_pilot(config, draws, alpha, mode=1, **k):
    """D2: internal-pilot re-estimation (Algorithm 2), pilot rows reused.

    Starts at ``J0`` on the confirmatory pool; looks up to ``nmax_looks`` times and
    re-estimates ``sigma_D_UB`` from the accumulated rows, growing J until the
    required closed-form size is reached; final test ONCE at ``alpha_adj`` on the
    full accumulated set INCLUDING the leading (reused) pilot rows.
    """
    J0 = config.get("J0", 50)
    Jmax = config.get("Jmax", 100_000)
    nmax = config.get("nmax_looks", 3)
    aa = config.get("alpha_adj", alpha)
    pool = draws["conf"]
    D = pool[:, :, 0] - pool[:, :, 1]   # (R, Jmax)
    R = D.shape[0]
    J = np.full(R, J0, dtype=int)
    for _ in range(nmax):
        _, _, s_D = paired_t_stats(D, J)
        sD = np.where(np.isclose(s_D, 0) | np.isnan(s_D), 1e-12, s_D)
        # Per-repetition df. After the first look J varies ACROSS repetitions, so
        # `int(J[0])` applied one repetition's df to all of them and mis-inflated
        # every other row.
        sigma_eff = _inflate(sD, J, config.get("gamma", 0.20))
        need, _ = solve_J(config["power_target"], config["delta"], sigma_eff,
                          aa, mode=mode, Jmax=Jmax)
        J_next = np.clip(np.maximum(need, J), J, Jmax)
        if np.array_equal(J_next, J):
            break
        J = J_next
    T, Dbar, s_D = paired_t_stats(D, J)
    p = paired_pvalue(T, J)
    capped = (J >= Jmax).astype(float)
    return _emit(T, p, J.astype(float), Dbar, s_D, capped), dict(
        pilot_reused=True, J_planned=J)


def design_v1_welch(config, draws, alpha, **k):
    """D3: TISCA v1's *unpaired* Welch design -- plan AND test as v1 did.

    Two things have to be faithful for the D3-vs-D4 comparison to mean anything,
    because it is the measurement behind the paper's "fixing IJDA #2 buys ~2.3x"
    claim:

    1. **The planning variance is the unpaired one.** v1 sized the study from the
       two-independent-sample formula, whose relevant scale is
       ``sqrt(s_A^2 + s_B^2)`` -- NOT ``sigma_D = sqrt(s_A^2 + s_B^2 - 2 rho
       s_A s_B)``. Planning D3 from ``config["sigma_D"]`` (as this function used
       to) hands v1 the paired variance it never had: at ``rho = 0.6`` the two
       differ by a factor of ~1.6 in sigma and ~2.5 in J, so D3 and D4 would come
       out nearly identical and the headline effect would vanish. The marginal
       sds are estimated from the leading ``J0`` rows, matching v1's plug-in
       character.
    2. **The test carries Welch-Satterthwaite df**, not ``J - 1``. Using ``J - 1``
       makes the v1 test look more conservative than it is, which understates
       exactly the defect being reported.
    """
    J0 = config.get("J0", 50)
    Jmax = config.get("Jmax", 100_000)
    pool = draws["conf"]
    R = pool.shape[0]
    LA = pool[:, :, 0]
    LB = pool[:, :, 1]

    # v1 plug-in planning scale from the leading J0 rows (the accumulating pilot).
    sA0 = LA[:, :J0].std(axis=1, ddof=1)
    sB0 = LB[:, :J0].std(axis=1, ddof=1)
    plan_sig = np.sqrt(sA0 ** 2 + sB0 ** 2)
    plan_sig = np.where(np.isclose(plan_sig, 0), 1e-12, plan_sig)

    need, ok = solve_J(config["power_target"], config["delta"], plan_sig,
                       alpha, mode=1, Jmax=Jmax)
    J = np.clip(np.maximum(need, J0), J0, Jmax).astype(int)

    sel = np.arange(LA.shape[1])[None, :] < J[:, None]
    la_sel = np.where(sel, LA, 0.0)
    lb_sel = np.where(sel, LB, 0.0)
    meanA = la_sel.sum(1) / J
    meanB = lb_sel.sum(1) / J
    s2A = np.where(sel, (LA - meanA[:, None]) ** 2, 0.0).sum(1) / (J - 1)
    s2B = np.where(sel, (LB - meanB[:, None]) ** 2, 0.0).sum(1) / (J - 1)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        vA, vB = s2A / J, s2B / J
        T = (meanA - meanB) / np.sqrt(vA + vB)
        # Welch-Satterthwaite, both groups of size J.
        df_w = (vA + vB) ** 2 / (vA ** 2 / (J - 1) + vB ** 2 / (J - 1))
    df_w = np.where(np.isfinite(df_w) & (df_w > 0), df_w, J - 1)
    p = 2.0 * stats.t.sf(np.abs(T), df_w)
    capped = np.where((J >= Jmax) & ~ok, 1.0, 0.0)
    return _emit(T, p, J.astype(float), meanA - meanB,
                 np.sqrt((s2A + s2B) / 2), capped), dict(
        pilot_reused=True, J_planned=J, plan_sigma_unpaired=plan_sig, df_welch=df_w)


def design_precision(config, draws, alpha, **k):
    """D5: paired fixed-precision (MCSE or CI-half-width target).

    Smallest J with ``s_D/sqrt(J) <= mcse`` (or ``t s_D/sqrt(J) <= half_width``),
    planned from the independent pilot; one confirmatory paired test at ``alpha``.
    """
    J0 = config.get("J0", 50)
    Jmax = config.get("Jmax", 100_000)
    pilot = draws["pilot"]
    Dp = pilot[:, :, 0] - pilot[:, :, 1]
    s_p = Dp.std(axis=1, ddof=1)
    s_p = np.where(np.isclose(s_p, 0), 1e-12, s_p)
    grid = np.arange(J0, Jmax + 1, dtype=float)
    if config.get("half_width") is not None:
        h = config["half_width"]
        v = stats.t.ppf(1 - alpha / 2, grid - 1) * s_p[:, None] / np.sqrt(grid[None, :])
    else:
        m = config.get("mcse", 0.05)
        v = s_p[:, None] / np.sqrt(grid[None, :])
        h = m
    ok = v <= h
    J = np.where(ok.any(1), grid[ok.argmax(1)].astype(int), int(Jmax))
    J = np.clip(J, 1, Jmax)
    conf = draws["conf"]
    D = conf[:, :, 0] - conf[:, :, 1]
    T, Dbar, s_D = paired_t_stats(D, J)
    p = paired_pvalue(T, J)
    return _emit(T, p, J.astype(float), Dbar, s_D, np.zeros_like(J)), dict(
        pilot_sD=s_p, J_planned=J, pilot_reused=False)


DESIGNS = {
    "D1": lambda c, d, a, **k: design_twostage({**c, "inflate": False, "alpha_adj": a}, d, a, **k),
    "D2": design_internal_pilot,
    "D3": design_v1_welch,
    "D4": design_twostage,
    "D5": design_precision,
    "D6": design_oracle,
}


def make_design(name, config):
    """Return the callable for one of the six designs with the given config."""
    if name not in DESIGNS:
        raise ValueError(f"unknown design {name!r}; choose from {sorted(DESIGNS)}")
    return DESIGNS[name]
