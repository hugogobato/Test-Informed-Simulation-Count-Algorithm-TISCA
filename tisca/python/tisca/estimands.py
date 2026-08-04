"""Estimand table rows — `estimands.py`.

Implements the replicate-level losses, across-replication estimators, and MCSE
formulas of `docs/estimand_table.md` row by row. Every function here operates on
a **per-replication** ``L_j`` (or on per-unit quantities aggregated to one value
per replication *before* any across-replication statistic).

Conventions (estimand table §1, all verbatim):
* ``i`` = unit index within a replication, ``j`` = replication index.
* ``tau_i`` = unit-level CATE; ``tau_j`` = replication-level ATE.
* rooted PEHE is ``E[sqrt(Q)]``; unrooted CATE-MSE is ``E[Q]``; these differ by
  Jensen (``E[sqrt Q] != sqrt(E[Q])``).
* RMSE_ATE (row 3b) is ``sqrt(E[Q_ATE])`` only, never a re-rooted per-replication
  value.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "rooted_pehe",
    "cate_mse",
    "ate_squared_error",
    "rmse_ate",
    "rmse_ate_mcse_delta",
    "ate_coverage",
    "cate_coverage",
    "interval_score_unit",
    "interval_score",
    "crps_mean",
    "calibration_deviation",
    "mean_mcse",
    "mean_estimate",
    "nc_ci",
    "winkler_score",
]


def rooted_pehe(tau_hat, tau_true, *, axis=-1) -> np.ndarray:
    """Row 1: rooted PEHE, ``PEHE_j = sqrt(mean_i (tau_hat - tau_true)^2)``.

    ``tau_hat``, ``tau_true`` may be per-unit vectors (returns one scalar per cell
    along ``axis``) or already per-replication columns (returns ``sqrt(mean(..))``
    of the single squared value).
    """
    tau_hat = np.asarray(tau_hat, dtype=float)
    tau_true = np.asarray(tau_true, dtype=float)
    q = np.mean((tau_hat - tau_true) ** 2, axis=axis)
    return np.sqrt(q)


def cate_mse(tau_hat, tau_true, *, axis=-1) -> np.ndarray:
    """Row 2: unrooted CATE MSE, ``Q_j = mean_i (tau_hat - tau_true)^2``."""
    tau_hat = np.asarray(tau_hat, dtype=float)
    tau_true = np.asarray(tau_true, dtype=float)
    return np.mean((tau_hat - tau_true) ** 2, axis=axis)


def ate_squared_error(ate_hat, ate_true) -> np.ndarray:
    """Row 3a: per-replication ATE squared error ``Q_ATE,j = (tau_hat - tau)^2``."""
    ate_hat = np.asarray(ate_hat, dtype=float)
    ate_true = np.asarray(ate_true, dtype=float)
    return (ate_hat - ate_true) ** 2


def rmse_ate(sq_err) -> float:
    """Row 3b: ``RMSE_ATE = sqrt( E[Q_ATE] )``.

    Operates on the per-replication squared errors ``Q_ATE,j`` (row 3a). The
    square root is applied to the *averaged* value, never to a single
    replication. This is the estimand that exists; "per-replication RMSE" does
    not (rejected in ``validate.py``).
    """
    q = np.asarray(sq_err, dtype=float)
    return float(np.sqrt(np.mean(q)))


def rmse_ate_mcse_delta(q_sq_err, *, sqrt_n: bool = True) -> float:
    """Row 3b delta-method MCSE for RMSE_ATE.

    With ``g(x)=sqrt(x)``, ``g'(x)=1/(2*sqrt(x))`` the delta-method MCSE is
    ``MCSE ≈ s_{Q} / ( 2*sqrt(J) * sqrt(mean(Q)) )``. The denominator carries
    ``sqrt(mean(Q))``, **not** ``mean(Q)`` (see estimand table row 3b).
    Returning the delta-method cross-check; the paired bootstrap is preferred.
    """
    q = np.asarray(q_sq_err, dtype=float)
    mu = float(np.mean(q))
    s = float(np.std(q, ddof=1))
    J = float(q.size)
    denom = 2.0 * np.sqrt(J) * np.sqrt(mu)
    out = s / denom if sqrt_n else s / (2.0 * mu)
    return float(out)


def ate_coverage(ate_hat, ci_lo, ci_hi, ate_true) -> np.ndarray:
    """Row 4: ATE coverage indicator ``1( tau_j in CI_j )``, binary per replication."""
    a = np.asarray(ate_hat, dtype=float)
    lo = np.asarray(ci_lo, dtype=float)
    hi = np.asarray(ci_hi, dtype=float)
    t = np.asarray(ate_true, dtype=float)
    return ((lo <= t) & (t <= hi)).astype(float)


def cate_coverage(tau_hat, ci_lo, ci_hi, tau_true, *, axis=-1) -> np.ndarray:
    """Row 5: CATE coverage proportion ``mean_i 1(tau_i in CI_i)`` per replication."""
    tau_hat = np.asarray(tau_hat, dtype=float)
    lo = np.asarray(ci_lo, dtype=float)
    hi = np.asarray(ci_hi, dtype=float)
    tau_true = np.asarray(tau_true, dtype=float)
    return np.mean((lo <= tau_true) & (tau_true <= hi), axis=axis).astype(float)


def interval_score_unit(lo, hi, x, *, level) -> float:
    """Row 9: Gneiting-Raftery interval score for one ``(lo, hi, x)`` triplet.

    For a ``(1-c)`` central interval, ``c = 1 - level`` and
    ``IS = (hi - lo) + (2/c)(lo - x) 1[x<lo] + (2/c)(x - hi) 1[x>hi]``.
    Lower is better (a proper scoring rule jointly penalising width and
    non-coverage).
    """
    c = 1.0 - level
    if not 0.0 < level < 1.0:
        raise ValueError("level must be in (0, 1).")
    lo = float(lo)
    hi = float(hi)
    x = float(x)
    pen_lo = (2.0 / c) * (lo - x) if x < lo else 0.0
    pen_hi = (2.0 / c) * (x - hi) if x > hi else 0.0
    return (hi - lo) + pen_lo + pen_hi


def interval_score(lo, hi, x, *, level, axis=-1) -> np.ndarray:
    """Row 9: mean interval score over units ``IS_j = mean_i IS(...)``.

    Vectorised over arrays ``lo, hi, x`` (per-unit lower/upper/observed).
    """
    level = float(level)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    x = np.asarray(x, dtype=float)
    c = 1.0 - level
    width = hi - lo
    pen_lo = (2.0 / c) * (lo - x) * (x < lo)
    pen_hi = (2.0 / c) * (x - hi) * (x > hi)
    return np.mean(width + pen_lo + pen_hi, axis=axis)


#: Alias consistent with the manuscript/Winkler literature.
winkler_score = interval_score


def crps_mean(unit_crps) -> np.ndarray:
    """Row 10: per-replication CRPS ``CRPS_j = mean_i CRPS_i``.

    ``unit_crps`` holds the per-unit CRPS values; the mean over units is the
    replication-level loss (the ``scoringRules`` values already live on this
    unit scale in the case-study CSVs).
    """
    u = np.asarray(unit_crps, dtype=float)
    return np.mean(u, axis=-1)


def calibration_deviation(cov_mean, nominal: float) -> float:
    """Row 8: calibration deviation applied to the *estimate*, never per replication.

    ``|mean(Cov_j) - nominal|``. The estimand table (§3.1) forbids
    ``E[|Cov_j - nominal|]`` because it confounds miscalibration with
    replication-level dispersion. The deviation is formed **after** averaging.
    """
    return float(abs(float(cov_mean) - float(nominal)))


def _std(x, ddof: int = 1) -> float:
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        return 0.0
    return float(np.std(x, ddof=ddof))


def mean_estimate(values) -> float:
    """Across-replication estimator ``theta_hat = mean(L_j)`` (estimand table §2)."""
    return float(np.mean(np.asarray(values, dtype=float)))


def mean_mcse(values) -> float:
    """MCSE of the sample mean: ``s_L / sqrt(J)`` (estimand table §2)."""
    v = np.asarray(values, dtype=float)
    if v.size < 2:
        return float("nan")
    return _std(v, ddof=1) / np.sqrt(v.size)


def nc_ci(mean, se, *, level: float, df: int) -> tuple[float, float]:
    """Normal-approximation / t CI ``mean +/- t_{1-alpha/2, J-1} se``.

    Uses the Student-t quantile with ``df = J - 1`` degrees of freedom (the
    discrete small-``J`` correction over the plain normal CI).
    """
    from scipy.stats import t as _t

    if not 0.0 < level < 1.0:
        raise ValueError("level must be in (0, 1).")
    crit = float(_t.ppf((1.0 + level) / 2.0, df=df))
    return (float(mean) - crit * float(se), float(mean) + crit * float(se))


def _binomial_mcse(p_hat, J: int) -> float:
    """MCSE of a proportion: ``sqrt(p_hat (1 - p_hat) / J)`` (rows 4, 12)."""
    p = float(p_hat)
    return float(np.sqrt(p * (1.0 - p) / J)) if J > 0 else float("nan")
