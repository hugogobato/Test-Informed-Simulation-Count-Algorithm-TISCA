"""Multiplicity correction and family power — `multiplicity.py`.

Implements C5 / §4 of `docs/tisca_v2_spec.md`:
* Bonferroni, Holm, and Benjamini-Hochberg adjusted p-values;
* **Romano-Wolf stepdown** (bootstrap, FWER-controlling) for a family of paired
  contrasts — the strictly better answer to IJDA #6 because it uses the actual
  correlation structure rather than the worst case;
* the mapping from a correction to a planning level (spec §4.3);
* marginal and conjunctive/disjunctive family power (§4.2).
"""

from __future__ import annotations

import numpy as np

from . import planning as _plan, validate

__all__ = [
    "p_adjust",
    "bonferroni",
    "holm",
    "benjamini_hochberg",
    "romano_wolf_stepdown",
    "planning_alpha",
    "family_power_conjunctive",
    "family_power_disjunctive",
]


def _check_pvals(p):
    p = np.asarray(p, dtype=float)
    if p.ndim != 1:
        raise validate.ValidationError("p-values must be a 1-D array.")
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise validate.ValidationError("p-values must lie in [0, 1] and be finite.")
    return p


def bonferroni(p):
    """Bonferroni adjusted p-values: ``min(K*p, 1)``."""
    p = _check_pvals(p)
    return np.minimum(p * p.size, 1.0)


def holm(p):
    """Holm (1979) step-down adjusted p-values (FWER control)."""
    p = _check_pvals(p)
    m = p.size
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for i in range(m):
        val = min(1.0, (m - i) * p[order[i]])
        running = max(running, val)
        adj[order[i]] = running
    return adj


def benjamini_hochberg(p):
    """Benjamini-Hochberg (1995) FDR-controlling adjusted p-values."""
    p = _check_pvals(p)
    m = p.size
    order = np.argsort(p)
    adj = np.empty(m)
    running = 1.0
    for i in range(m - 1, -1, -1):
        val = min(1.0, (m / (i + 1)) * p[order[i]])
        running = min(running, val)
        adj[order[i]] = running
    return adj


def romano_wolf_stepdown(
    D,
    *,
    B: int = 4999,
    alpha: float = 0.05,
    seed: int | None = None,
    bootstrap_indices: np.ndarray | None = None,
) -> dict:
    """Romano-Wolf (2005) stepdown adjusted p-values for a family of contrasts.

    ``D`` is the ``(J, K)`` matrix of per-replication contrasts (columns =
    primary contrasts, rows = common replications). For each contrast the
    studentized statistic is ``T_k = mean(D_k)/(sd(D_k)/sqrt(J))``. The null
    distribution is approximated by the *ordinary* i.i.d. bootstrap (block length
    1 — the replications are i.i.d., plan §2.3 item 4), and the stepdown uses the
    max-combining statistic over the set of hypotheses still in contention.

    Returns FWER-adjusted p-values plus the marginal t statistics and the number
    of rejections at level ``alpha``.
    """
    D = np.asarray(D, dtype=float)
    if D.ndim != 2:
        raise validate.ValidationError("D must be a (J, K) matrix of contrasts.")
    J, K = D.shape
    if J < 2 or K < 2:
        raise validate.ValidationError("Need >= 2 replications and >= 2 contrasts.")
    if not 0.0 < alpha < 1.0:
        raise validate.ValidationError("alpha must lie in (0, 1).")

    means = D.mean(axis=0)
    sds = D.std(axis=0, ddof=1)
    ses = sds / np.sqrt(J)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_obs = np.where(sds > 0, means / ses, 0.0)
        t_obs = np.where(sds == 0, 0.0, t_obs)

    if bootstrap_indices is not None:
        idx = np.asarray(bootstrap_indices, dtype=int)
        if idx.ndim != 2 or idx.shape[0] != J or idx.shape[1] != B:
            raise validate.ValidationError(
                f"bootstrap_indices must have shape (J={J}, B={B}); got {idx.shape}."
            )
        idx = idx.reshape(J, B)
    else:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, J, size=(J, B))

    t_b = np.empty((B, K))
    for b in range(B):
        Db = D[idx[:, b], :]
        mb = Db.mean(axis=0)
        sb = Db.std(axis=0, ddof=1)
        se_b = sb / np.sqrt(J)
        with np.errstate(divide="ignore", invalid="ignore"):
            tb = np.where(sb > 0, (mb - means) / se_b, 0.0)
            tb = np.where(sb == 0, 0.0, tb)
        t_b[b, :] = tb

    # Stepdown: iterate hypotheses from most to least significant.
    # The bootstrap max-statistic is built from |t_b|, so the observed statistic
    # it is compared against must also be |t_obs|. Comparing |t_b| against the
    # SIGNED t_obs makes every negative contrast unrejectable: for a
    # lower-is-better loss "A beats B" means D_j < 0 and t_obs < 0, and
    # Pr(max|t_b| >= negative) == 1, so the whole case-study family would come
    # back with adjusted p = 1.
    abs_t = np.abs(t_obs)
    order = np.argsort(-abs_t)
    p_adj = np.zeros(K)
    running = 0.0
    for step in range(K):
        rem = order[step:]           # hypotheses still in contention
        M_b = np.max(np.abs(t_b[:, rem]), axis=1)
        p_step = float(np.mean(M_b >= abs_t[order[step]]))
        running = max(running, p_step)
        p_adj[order[step]] = min(1.0, running)

    rejections = np.array(p_adj <= alpha, dtype=bool)
    return {
        "p_values": p_adj,
        "t_stats": t_obs,
        "means": means,
        "B": B,
        "alpha": alpha,
        "rejections": rejections,
        "family_size": K,
        "seed": seed,
    }


def p_adjust(p, method: str):
    """Dispatch to a multiplicity correction by name.

    ``method`` in {'bonferroni', 'holm', 'bh', 'benjamini_hochberg', 'none'}.
    """
    method = method.lower()
    if method in ("none", None):
        return _check_pvals(p).copy()
    if method == "bonferroni":
        return bonferroni(p)
    if method in ("holm",):
        return holm(p)
    if method in ("bh", "benjamini_hochberg", "fdr"):
        return benjamini_hochberg(p)
    raise validate.ValidationError(f"Unknown multiplicity method {method!r}.")


def planning_alpha(correction: str, K: int, alpha: float = 0.05, r: int | None = None) -> tuple[float, str]:
    """Map a correction to the planning alpha and its interpretation (§4.3).

    Returns ``(alpha_plan, note)``:
    * 'bonferroni' / 'holm' -> ``alpha/K`` (Holm is asymptotically no worse, but
      planning at ``alpha/K`` is conservative and valid);
    * 'bh' -> ``alpha*r/K`` with ``r`` the pre-declared conservative number of
      claims expected to be rejected (``r = 1`` is the safe default);
    * 'romano_wolf' -> ``alpha`` is returned as a schedule; power must be planned
      by simulation from ``Sigma_D`` (§4.4), so ``alpha_plan = alpha`` is the
      level the stepdown is run at.
    """
    if K < 1:
        raise validate.ValidationError("family size K must be >= 1.")
    if not 0.0 < alpha < 1.0:
        raise validate.ValidationError("alpha must lie in (0, 1).")
    c = (correction or "none").lower()
    if c in ("none",):
        return alpha, "no multiplicity control (caller declared acceptance)"
    if c in ("bonferroni", "holm"):
        return alpha / K, f"plan every contrast at alpha/K = {alpha:.5f}/{K}"
    if c in ("bh", "benjamini_hochberg", "fdr"):
        r = 1 if r is None else r
        if r < 1 or r > K:
            raise validate.ValidationError("r (expected rejections) must be in [1, K].")
        return alpha * r / K, (
            f"BH: plan at alpha*r/K with r = {r} pre-declared expected rejections "
            f"(r = 1 reduces to Bonferroni)."
        )
    if c in ("romano_wolf", "romano-wolf", "rw"):
        return alpha, "Romano-Wolf power is planned by simulation from Sigma_D (§4.4); use this alpha at the stepdown."
    raise validate.ValidationError(f"Unknown correction {correction!r}.")


def family_power_conjunctive(power: np.ndarray) -> float:
    """Conjunctive ('all K') family power = product of marginal powers.

    Under working independence this is the probability that all primary claims
    hold; with correlated contrasts the joint probability is larger, so the
    product is a lower bound. The E1 outer-MC study measures the true joint rate.
    """
    pw = np.asarray(power, dtype=float)
    if np.any((pw < 0) | (pw > 1)):
        raise validate.ValidationError("marginal powers must lie in [0, 1].")
    return float(np.prod(pw))


def family_power_disjunctive(power: np.ndarray) -> float:
    """Disjunctive ('at least one') family power = 1 - product(1 - power_k)."""
    pw = np.asarray(power, dtype=float)
    if np.any((pw < 0) | (pw > 1)):
        raise validate.ValidationError("marginal powers must lie in [0, 1].")
    return float(1.0 - np.prod(1.0 - pw))
