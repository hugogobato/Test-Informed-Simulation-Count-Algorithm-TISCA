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
  D3  TISCA v1: iterative, unpaired Welch (the v1 defect being studied)
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
    return 1.0 - (stats.nct.cdf(crit, J - 1, ncp) - stats.nct.cdf(-crit, J - 1, ncp))


def m2_directional_power(J, delta, sigma_D, alpha):
    """M2 one-sided power (lower tail, lower-is-better)."""
    ncp = np.sqrt(J) * delta / sigma_D
    return stats.nct.cdf(stats.t.ppf(alpha, J - 1), J - 1, ncp)


def _power_for_mode(mode, J, delta, sigma_D, alpha):
    if mode == 1:
        return m1_power(J, delta, sigma_D, alpha)
    if mode == 2:
        return m2_directional_power(J, delta, sigma_D, alpha)
    raise ValueError(f"planning power for mode {mode} not implemented in the harness")


def solve_J(power_target, delta, sigma_D, alpha, mode=1, Jmin=4, Jmax=10_000):
    """Smallest J with power(J) >= power_target (per repetition if sigma_D is (R,)).

    Returns ``(J, within_cap)``. Because power is monotone in J for the M1/M2
    modes, evaluating on a dense grid and taking the first index that meets the
    target is exact for M1 (up to grid resolution); for M2 a refined scan is used.
    """
    aa = alpha
    sd = np.atleast_1d(np.asarray(sigma_D, dtype=float))
    grid = np.arange(Jmin, Jmax + 1, dtype=float)
    pw = _power_for_mode(mode, grid[None, :], delta, sd[:, None], aa)
    ok = pw >= power_target
    any_ok = ok.any(axis=1)
    idx = np.where(any_ok, ok.argmax(axis=1), -1)
    J = np.where(any_ok, grid[idx].astype(int), int(Jmax))
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
        sigma_eff = _inflate(sD, int(J[0]), config.get("gamma", 0.20))
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
    """D3: TISCA v1's iterative *unpaired* Welch, stopping at 80% power.

    Reproduces the v1 defect (independent Welch on paired data): the design stops
    at the size the (mis-specified) unpaired power target requires, then runs an
    independent-samples Welch test. Its unconditional Type I error is inflated,
    which P3-T2/P3-T3 report as the finding in favour of the two-stage default.
    """
    J0 = config.get("J0", 50)
    Jmax = config.get("Jmax", 100_000)
    pool = draws["conf"]
    R = pool.shape[0]
    LA = pool[:, :, 0]
    LB = pool[:, :, 1]
    plan_sig = config.get("sigma_D", 1.0)
    need, _ = solve_J(config["power_target"], config["delta"], plan_sig,
                      alpha, mode=1, Jmax=Jmax)
    need = np.broadcast_to(np.atleast_1d(np.asarray(need, dtype=int)), (R,))
    J = np.clip(np.maximum(need, J0), J0, Jmax).astype(int)
    sel = np.arange(LA.shape[1])[None, :] < J[:, None]
    la_sel = np.where(sel, LA, 0.0)
    lb_sel = np.where(sel, LB, 0.0)
    meanA = la_sel.sum(1) / J
    meanB = lb_sel.sum(1) / J
    s2A = np.where(sel, (LA - meanA[:, None]) ** 2, 0.0).sum(1) / (J - 1)
    s2B = np.where(sel, (LB - meanB[:, None]) ** 2, 0.0).sum(1) / (J - 1)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        T = (meanA - meanB) / np.sqrt(s2A / J + s2B / J)
    p = 2.0 * stats.t.sf(np.abs(T), J - 1)
    return _emit(T, p, J.astype(float), meanA - meanB,
                 np.sqrt((s2A + s2B) / 2), np.zeros_like(J)), dict(pilot_reused=True,
                                                                    J_planned=J)


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
