"""Paired inference — `inference.py`.

Implements §7.1 of `docs/tisca_v2_spec.md`:
* **paired one-sample t** on ``D_j`` (df ``J-1``) as the fast default;
* **studentized paired bootstrap** (B resamples, its MCSE reported) as the
  recommended check for non-trivial skewness;
* **McNemar exact / binomial** test for paired binary indicators (rows 4, 12).

All contrasts are formed on common replications with **listwise deletion across
the pair** (convention 3); per-column NA dropping is forbidden.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import t as _t, binom as _binom

from . import validate, estimands

__all__ = [
    "paired_t",
    "paired_t_test",
    "studentized_paired_bootstrap",
    "paired_bootstrap_ci",
    "mcnemar_exact",
    "mcnemar_with_continuity",
    "contrast_summary",
]

_DEFAULT_B = 4999


def paired_t(D, mu0: float = 0.0, *, alternative: str = "two-sided") -> dict:
    """One-sample t on the paired contrasts ``D_j``.

    Returns estimate, SE, t statistic, degrees of freedom, and a ``p_value``.
    ``alternative`` is 'two-sided', 'less' (``H1: theta < 0``, lower-is-better) or
    'greater'.
    """
    D = np.asarray(D, dtype=float)
    J = D.size
    if J < 2:
        raise validate.ValidationError(f"paired_t needs >= 2 replications, got {J}.")
    mean = float(np.mean(D))
    sd = float(np.std(D, ddof=1))
    se = sd / np.sqrt(J)
    df = J - 1
    if sd == 0.0:
        # Degenerate contrasts have no sampling uncertainty on the difference.
        if alternative == "two-sided":
            p = 0.0 if abs(mean - mu0) > 0 else 1.0
        elif alternative == "less":
            p = 0.0 if mean < mu0 else 1.0
        elif alternative == "greater":
            p = 0.0 if mean > mu0 else 1.0
        else:
            raise validate.ValidationError(f"Unknown alternative {alternative!r}.")
        tstat = 0.0
    else:
        tstat = (mean - mu0) / se
        if alternative == "two-sided":
            p = 2.0 * _t.sf(abs(tstat), df=df)
        elif alternative == "less":
            p = _t.cdf(tstat, df=df)
        elif alternative == "greater":
            p = _t.sf(tstat, df=df)
        else:
            raise validate.ValidationError(f"Unknown alternative {alternative!r}.")
    return {
        "estimate": mean,
        "sd": sd,
        "se": se,
        "mcse": se,
        "t": tstat,
        "df": df,
        "J": J,
        "p_value": float(p),
        "alternative": alternative,
    }


paired_t_test = paired_t


def studentized_paired_bootstrap(
    D,
    *,
    mu0: float = 0.0,
    B: int = _DEFAULT_B,
    alternative: str = "two-sided",
    seed: int | None = None,
    bootstrap_indices: np.ndarray | None = None,
) -> dict:
    """Studentized paired bootstrap over the mean ``D_bar``.

    For each resample the studentized statistic is
    ``T_b = (D_bar_b - D_bar) / se_b``; the confidence interval is built from the
    quantiles of ``T_b`` and the p-value from their centred position. Returns an
    ``estimate``, a ``(1-level)`` ``ci``, a ``p_value``, and the bootstrap MCSE.
    """
    D = np.asarray(D, dtype=float)
    J = D.size
    if J < 2:
        raise validate.ValidationError("studentized bootstrap needs >= 2 replications.")
    if alternative not in ("two-sided", "less", "greater"):
        raise validate.ValidationError(f"Unknown alternative {alternative!r}.")
    mean = float(np.mean(D))
    sd = float(np.std(D, ddof=1))
    if sd == 0.0:
        # Degenerate: all bootstrap studentized statistics are undefined; fall back
        # to the deterministic point.
        base = dict(
            estimate=mean, sd=sd, se=0.0, mcse=0.0, p_value=None, ci=(mean, mean),
            B=B, seed=seed, degenerate=True,
        )
        return base

    se_obs = sd / np.sqrt(J)

    if bootstrap_indices is not None:
        idx = np.asarray(bootstrap_indices, dtype=int)
        if idx.ndim != 2 or idx.shape[1] != B or idx.shape[0] != J:
            raise validate.ValidationError(
                f"bootstrap_indices must have shape (J={J}, B={B}); got {idx.shape}."
            )
    else:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, J, size=(J, B))

    shape = idx.shape
    D_boot = D[idx]  # (J, B)
    D_bar_b = D_boot.mean(axis=0)
    sd_b = D_boot.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        se_b = sd_b / np.sqrt(J)
        T_b = (D_bar_b - mean) / se_b
    T_b = np.nan_to_num(T_b, nan=0.0, posinf=10.0, neginf=-10.0)

    level = 0.95
    if alternative == "two-sided":
        lo_q = np.quantile(T_b, 1 - (1 - level) / 2)
        hi_q = np.quantile(T_b, (1 - level) / 2)
        lower = mean - lo_q * se_obs
        upper = mean - hi_q * se_obs
        t_obs = (mean - mu0) / se_obs
        p = 2.0 * min(float(np.mean(T_b <= -abs(t_obs))), float(np.mean(T_b >= abs(t_obs))), 1.0)
    elif alternative == "less":
        hi_q = np.quantile(T_b, 1 - level)
        upper = mean - hi_q * se_obs
        lower = -np.inf
        t_obs = (mean - mu0) / se_obs
        p = float(np.mean(T_b <= t_obs))
    else:  # greater
        # One-sided LOWER bound [mean - q_{1-a}(T) se, inf). The upper quantile
        # is the right one here: using q_{a}(T) (which is negative) put the
        # "lower" bound ABOVE the point estimate.
        lo_q = np.quantile(T_b, level)
        lower = mean - lo_q * se_obs
        upper = np.inf
        t_obs = (mean - mu0) / se_obs
        p = float(np.mean(T_b >= t_obs))
    p = min(1.0, max(0.0, p))

    return dict(
        estimate=mean,
        sd=sd,
        se=se_obs,
        mcse=se_obs,
        p_value=float(p),
        ci=(float(lower), float(upper)),
        B=B,
        seed=seed,
        degenerate=False,
    )


paired_bootstrap_ci = studentized_paired_bootstrap


def _discordant_counts(a, b):
    """McNemar table cells: n01, n10."""
    a = np.asarray(a, dtype=int)
    b = np.asarray(b, dtype=int)
    n01 = int(np.sum((a == 0) & (b == 1)))
    n10 = int(np.sum((a == 1) & (b == 0)))
    return n01, n10


def mcnemar_exact(a, b) -> dict:
    """Exact (binomial) McNemar test on paired binary indicators.

    Under ``H0`` the discordant pairs are equally likely; conditioning on
    ``m = n01 + n10``, ``n10 ~ Binomial(m, 1/2)``. Two-sided p-value uses the
    doubled one-sided tail capped at 1.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if a.size != b.size:
        raise validate.ValidationError("McNemar requires equal-length paired indicators.")
    validate.validate_binary(a, name="a")
    validate.validate_binary(b, name="b")
    n01, n10 = _discordant_counts(a, b)
    m = n01 + n10
    if m == 0:
        p = 1.0
    else:
        k = min(n01, n10)
        p_one = _binom.cdf(k, n=m, p=0.5)
        p = min(1.0, 2.0 * p_one)
    return {
        "n": a.size,
        "n01": n01,
        "n10": n10,
        "m_discordant": m,
        "p_value": float(p),
        "exact": True,
    }


def mcnemar_with_continuity(a, b) -> dict:
    """McNemar with Yates continuity correction (chi-square approximation)."""
    from scipy.stats import chi2 as _chi2

    a = np.asarray(a)
    b = np.asarray(b)
    validate.validate_binary(a, name="a")
    validate.validate_binary(b, name="b")
    n01, n10 = _discordant_counts(a, b)
    diff = abs(n10 - n01)
    m = n01 + n10
    if m == 0:
        chi2 = 0.0
    else:
        chi2 = ((diff - 1.0) ** 2) / float(m)
    p = float(_chi2.sf(chi2, df=1))
    return {"n01": n01, "n10": n10, "m_discordant": m, "chi2": chi2, "p_value": p}


def contrast_summary(A, B, *, name: str = "contrast", paired: bool = True) -> dict:
    """Build a full per-contrast summary on the common-replication contrast ``D``.

    Applies listwise deletion across the pair (convention 3), then reports the
    paired estimate, MCSE, CI, and the number of replications dropped due to
    missing values.
    """
    a, b, n_dropped = validate.validate_contrast_pair(A, B, name=name)
    D = a - b
    res = paired_t(D)
    res["n_dropped"] = n_dropped
    res["D"] = D
    ci = estimands.nc_ci(res["estimate"], res["se"], level=0.95, df=res["df"])
    res["ci"] = ci
    res["name"] = name
    res["paired"] = paired
    return res
