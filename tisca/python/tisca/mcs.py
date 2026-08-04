"""Model Confidence Set layer — `mcs.py`.

Implements §7.2 / plan §2.3: **Hansen's SPA**, **White's Reality Check**, and the
**Model Confidence Set** (Hansen, Lunde & Nason 2011) with both the ``T_R`` and
``T_max`` statistics.

The MCS algorithm below is translated line-for-line from the CRAN ``MCS`` package
(Bernardi & Catania, ``MCS::MCSprocedure``) so that the Python reference
implementation and the R oracle agree to bootstrap MCSE on test losses. Because
TISCA's replications are i.i.d., the block length collapses to 1 (an ordinary
i.i.d. bootstrap); both implementations are run at ``k = 0`` for parity.

``bootstrap_indices`` (an ``(T, B)`` integer matrix of resample indices) lets a
caller drive *both* this implementation and the CRAN package from the *same*
bootstrap resamples, giving an exact (bit-for-bit within floating point)
parity check rather than merely agreement within Monte Carlo error.
"""

from __future__ import annotations

import numpy as np

from . import validate

__all__ = [
    "loss_differential",
    "spa_pvalue",
    "reality_check_pvalue",
    "mcs",
    "mcs_elimination_path",
    "GetD",
    "GetIndices",
]


def GetD(mL) -> tuple[np.ndarray, np.ndarray]:
    """CRAN ``MCS:::GetD`` — pairwise loss differences and model averages.

    ``mL`` is ``(T, m)``. Returns ``(mD_ij_bar, vD_i_bar)`` where
    ``mD_ij_bar[i, j] = mean(L_i - L_j)`` and
    ``vD_i_bar[i] = sum_j mD_ij_bar[i, j] / (m - 1)``.
    """
    mL = validate.validate_loss_array(mL)
    iM = mL.shape[1]
    mD_ij_bar = np.zeros((iM, iM))
    for i in range(iM):
        for j in range(i, iM):
            if i != j:
                v = np.mean(mL[:, i] - mL[:, j])
                mD_ij_bar[i, j] = v
                mD_ij_bar[j, i] = -v
    if iM > 1:
        vD_i_bar = mD_ij_bar.sum(axis=1) / (iM - 1)
    else:
        vD_i_bar = np.zeros(iM)
    return mD_ij_bar, vD_i_bar


def GetIndices(iT: int, k: int, B: int, seed: int | None = None) -> np.ndarray:
    """CRAN ``MCS:::GetIndices`` — bootstrap resample indices with block length ``k``.

    For ``k == 0`` this is an ordinary i.i.d. multinomial resample of the rows.
    Returns an ``(iT, B)`` matrix. Uses a NumPy MT-compatible RNG when ``seed`` is
    given; if an exact match to R's draws is required, supply ``bootstrap_indices``
    directly to :func:`mcs` instead.
    """
    if k > 0:
        iN = int(np.ceil(iT / k))
    else:
        iN = iT
    rng = np.random.default_rng(seed)
    mIndices = np.zeros((iT, B), dtype=int)
    for b in range(B):
        vIndices = rng.integers(1, iT - k + 1, size=iN)  # sample(1:(iT-k), iN, repl) in 1..(iT-k)
        if k > 0:
            blocks = np.concatenate([v + np.arange(k) for v in vIndices])
            mIndices[:, b] = blocks[:iT]
        else:
            # sample(1:iT, iT, replace=TRUE) -> values in 1..iT
            mIndices[:, b] = rng.integers(1, iT + 1, size=iT)
    # Convert to 0-based for numpy indexing.
    if k == 0:
        return mIndices - 1
    return mIndices - 1


def _centered_means(mL, idx):
    """Given resample indices ``idx`` (T,), return the resampled row-mean vector."""
    return np.mean(mL[idx], axis=0)


def _mcs_stats(mL, idx):
    """Compute the per-variable statistics for one data set and its resamples.

    Mirrors the CRAN loop body for a single block. ``mL`` is ``(T, m)``, ``idx``
    is ``(T, B)`` 0-based. Returns ``(mD_ij_bar, vD_i_bar, mD_ij_bar_var,
    vD_i_bar_var, tij, ti, tij_res, ti_res)``.
    """
    iT, iM = mL.shape
    B = idx.shape[1]

    # Observed GetD on the full data.
    mD_ij_bar, vD_i_bar = GetD(mL)

    # Resampled data: mL_res_all shape (T, B, m).
    mL_res = mL[idx]  # (T, B, m) -- row b = mL[idx[:, b], :]
    # Per-resample column means: shape (B, m)
    mean_res = mL_res.mean(axis=0)
    # aD_ij_bar_res[b, i, j] = mean_res[b, i] - mean_res[b, j]
    aD_ij_bar_res = mean_res[:, :, None] - mean_res[:, None, :]  # (B, m, m)
    # mD_i_bar_res[i, b] = sum_j aD_ij_bar_res[i,j,b] / (m-1)
    mD_i_bar_res = aD_ij_bar_res.sum(axis=2) / (iM - 1)  # (B, m) -> transpose to (m, B)

    # Variances across the B resamples, centered on the observed statistic.
    mD_ij_bar_var = np.mean((aD_ij_bar_res - mD_ij_bar[None, :, :]) ** 2, axis=0)  # (m, m)
    vD_i_bar_var = np.mean((mD_i_bar_res.T - vD_i_bar[:, None]) ** 2, axis=1)  # (m,)

    # Statistics, oriented to match the CRAN indexing:
    #   ti_res[i, b]   -> shape (m, B)
    #   tij_res[i, j, b] -> shape (m, m, B)
    a_perm = aD_ij_bar_res.transpose(1, 2, 0)          # (m, m, B)
    mD_i_bar_res_perm = mD_i_bar_res.T                 # (m, B)
    with np.errstate(divide="ignore", invalid="ignore"):
        tij = np.where(mD_ij_bar_var > 0, mD_ij_bar / np.sqrt(mD_ij_bar_var), 0.0)
        ti = np.where(vD_i_bar_var > 0, vD_i_bar / np.sqrt(vD_i_bar_var), 0.0)
        tij_res = np.where(
            mD_ij_bar_var[:, :, None] > 0,
            (a_perm - mD_ij_bar[:, :, None]) / np.sqrt(mD_ij_bar_var)[:, :, None],
            0.0,
        )
        ti_res = np.where(
            vD_i_bar_var[:, None] > 0,
            (mD_i_bar_res_perm - vD_i_bar[:, None]) / np.sqrt(vD_i_bar_var)[:, None],
            0.0,
        )
    return mD_ij_bar, vD_i_bar, mD_ij_bar_var, vD_i_bar_var, tij, ti, tij_res, ti_res


def mcs(
    Loss,
    alpha: float = 0.15,
    B: int = 4999,
    statistic: str = "Tmax",
    k: int | None = None,
    bootstrap_indices: np.ndarray | None = None,
    seed: int | None = None,
    model_names=None,
) -> dict:
    """Model Confidence Set procedure (Hansen, Lunde & Nason 2011).

    ``Loss`` is a ``(T, m)`` matrix, one scalar loss per model per replication,
    lower is better. ``statistic`` is 'Tmax' or 'TR'. Returns a dict with the
    per-model table (average loss, elimination-order p-values, MCS p-values), the
    included/excluded model sets, and the elimination path.

    If ``bootstrap_indices`` (``(T, B)`` 0-based) is supplied, it is used verbatim
    so the same resamples can drive this implementation and the R oracle.
    Otherwise, if ``k`` is None, it defaults to ``k = 0`` (i.i.d. replications,
    block length 1).
    """
    mL = validate.validate_loss_array(Loss)
    iT, iM = mL.shape
    if model_names is None:
        vModels = [f"model_{i + 1}" for i in range(iM)]
    else:
        vModels = list(model_names)
        if len(vModels) != iM:
            raise validate.ValidationError("model_names length must match Loss columns.")
    vModels = np.array(vModels)

    if k is None:
        k = 0
    if k < 0 or k >= iT:
        raise validate.ValidationError("k must be in [0, T).")
    if B <= 0:
        raise validate.ValidationError("B must be positive.")

    if bootstrap_indices is not None:
        idx = np.asarray(bootstrap_indices, dtype=int)
        if idx.ndim != 2 or idx.shape != (iT, B):
            raise validate.ValidationError(
                f"bootstrap_indices must have shape (T={iT}, B={B}); got {idx.shape}."
            )
    else:
        idx = GetIndices(iT, k, B, seed=seed)

    # Both must start as NaN: the model that survives the elimination path never
    # gets a p assigned in the loop, and the backfill below turns NaN into 1.0.
    # Initialising mTab_p to zeros made that backfill a no-op, so the surviving
    # (best) model carried an MCS p-value of 0 -- harmless only for as long as
    # inclusion was decided on p_H0 instead of p_mcs.
    mTab_p = np.full(iM, np.nan)
    mTab_p_H0 = np.full(iM, np.nan)
    mTab_avg = mL.mean(axis=0)
    elim_order = []
    elim_p = []

    working = mL.copy()
    working_idx = vModels.copy()
    working_boot = idx.copy()
    while working.shape[1] > 1:
        iM_now = working.shape[1]
        mD_ij_bar, vD_i_bar, mD_ij_bar_var, vD_i_bar_var, tij, ti, tij_res, ti_res = _mcs_stats(
            working, working_boot
        )
        if statistic == "Tmax":
            T_Max = float(np.max(ti))
            T_Max_res = ti_res.max(axis=0)  # max over models, per bootstrap b
            p = float(np.mean(T_Max_res > T_Max))
            e = int(np.argmax(ti))
        elif statistic == "TR":
            T_R = float(np.max(np.abs(tij)))
            T_R_res = np.abs(tij_res).max(axis=(0, 1))  # max over (i, j), per bootstrap b
            p = float(np.mean(T_R_res > T_R))
            e = int(np.argmax(np.max(tij, axis=1)))
        else:
            raise validate.ValidationError("statistic must be 'Tmax' or 'TR'.")

        eliminated_model = working_idx[e]
        pos = int(np.where(vModels == eliminated_model)[0][0])
        mTab_p_H0[pos] = p
        # CRAN: MCS p-Value for the eliminated model = max over all non-NA
        # p-Value-for-H0 entries computed so far (the running maximum).
        mTab_p[pos] = float(np.nanmax(mTab_p_H0))
        elim_order.append(str(eliminated_model))
        elim_p.append(p)

        working = np.delete(working, e, axis=1)
        working_idx = np.delete(working_idx, e)

    # The final surviving model has no elimination p; its MCS p is set to 1.
    mTab_p_H0[np.isnan(mTab_p_H0)] = 1.0
    mTab_p[np.isnan(mTab_p)] = 1.0

    # Order models by MCS p-value (ascending), as CRAN does.
    order = np.argsort(mTab_p)
    table = np.column_stack(
        [mTab_avg[order], mTab_p_H0[order], mTab_p[order]]
    )
    table_names = vModels[order]

    # Hansen, Lunde & Nason (2011) define the MCS at level alpha as
    # M*_{1-alpha} = { i : p_MCS_i > alpha }, where p_MCS is the RUNNING MAXIMUM
    # of the elimination p-values along the path -- not the per-step p_H0. The
    # loop above eliminates every model regardless of its p-value, so selecting
    # on p_H0 keeps models whose elimination was already implied by an earlier
    # rejection and drops models the MCS should retain.
    included = vModels[mTab_p > alpha]
    excluded = vModels[mTab_p <= alpha]

    return {
        "avg_loss": mTab_avg,
        "p_H0": mTab_p_H0,
        "p_mcs": mTab_p,
        "table": table,
        "table_names": table_names,
        "included": [str(x) for x in included],
        "excluded": [str(x) for x in excluded],
        "alpha": alpha,
        "B": B,
        "statistic": statistic,
        "k": k,
        "elimination_order": elim_order,
        "elimination_pvalues": elim_p,
        "model_names": [str(x) for x in vModels],
        "seed": seed,
    }


def mcs_elimination_path(Loss, **kwargs):
    """Return only the elimination path (order of models removed and their p-values)."""
    res = mcs(Loss, **kwargs)
    return list(zip(res["elimination_order"], res["elimination_pvalues"]))


def _loss_center(mL, base_series: int = 0) -> np.ndarray:
    """Centre a loss matrix by subtracting the champion series' losses."""
    mL = validate.validate_loss_array(mL)
    return mL - mL[:, base_series, None]


def reality_check_pvalue(
    Loss, champion: int = 0, *, B: int = 4999, seed: int | None = None, k: int = 0,
    direction: str = "champion_beats_some",
) -> dict:
    """White's Reality Check (2000), a max-type test on paired loss differentials.

    **Read the null carefully — it is not "the champion beats all benchmarks".**
    Both directions are max-type tests and neither establishes dominance over
    every competitor by rejection:

    ``direction='champion_beats_some'`` (default, the orientation this module
    has always computed): ``d_k = L_k - L_champ``, ``T = max_k mean(d_k)``.
    ``H0: max_k E[L_k - L_champ] <= 0`` — the champion is (weakly) the worst of
    the set. Rejection licenses only "the champion beats **at least one**
    competitor".

    ``direction='nobody_beats_champion'`` (the classical SPA/RC role
    assignment, champion = benchmark): ``d_k = L_champ - L_k``,
    ``H0: max_k E[L_champ - L_k] <= 0`` — no competitor beats the champion.
    Rejection is evidence **against** the champion. Failing to reject is the
    usual "no competitor was shown to beat it" diagnostic.

    The positive claim "the champion beats every benchmark" is an
    intersection-union statement, ``min_k E[L_k - L_champ] > 0``. It is
    established by rejecting **all** K paired contrasts under FWER control
    (:func:`tisca.multiplicity.romano_wolf_stepdown`), not by either max-test
    here. ``conjunctive_beats_all`` in the return value reports that verdict at
    the marginal level as a convenience flag.
    """
    if direction not in ("champion_beats_some", "nobody_beats_champion"):
        raise validate.ValidationError(
            "direction must be 'champion_beats_some' or 'nobody_beats_champion'."
        )
    mL = validate.validate_loss_array(Loss)
    iT, iM = mL.shape
    others = [c for c in range(iM) if c != champion]
    # _loss_center(mL, champion)[:, k] == L_k - L_champ.
    rck = _loss_center(mL, champion)[:, others]  # (T, m-1)
    if direction == "nobody_beats_champion":
        rck = -rck
    dbar = rck.mean(axis=0)
    T_obs = float(np.max(dbar))

    idx = GetIndices(iT, k, B, seed=seed)
    T_b = np.empty(B)
    for b in range(B):
        rck_b = rck[idx[:, b], :]
        dbar_b = rck_b.mean(axis=0)
        T_b[b] = float(np.max(dbar_b - dbar))
    p = float(np.mean(T_b > T_obs))
    return {
        "p_value": p,
        "T_obs": T_obs,
        "B": B,
        "seed": seed,
        "direction": direction,
        "null": (
            "max_k E[L_k - L_champ] <= 0 (champion is weakly worst)"
            if direction == "champion_beats_some"
            else "max_k E[L_champ - L_k] <= 0 (no competitor beats the champion)"
        ),
        "conjunctive_beats_all": bool(np.all(_loss_center(mL, champion)[:, others].mean(axis=0) > 0)),
    }


def spa_pvalue(
    Loss, champion: int = 0, *, B: int = 4999, seed: int | None = None, k: int = 0,
    direction: str = "champion_beats_some",
) -> dict:
    """Hansen's Superior Predictive Ability (2005): the studentized max-test.

    ``d_k = L_k - L_champ`` (default ``direction``), ``T = max_k sqrt(T)
    dbar_k / sd_k``. The null and the interpretation of a rejection are exactly
    as documented in :func:`reality_check_pvalue` — in particular, rejection
    does **not** establish that the champion beats every benchmark.

    Recentring is **full** (``d - dbar``), which is Hansen's liberal variant
    ``SPA_l`` (equivalently the studentized Reality Check). Hansen's recommended
    consistent variant ``SPA_c`` recentres only the series that fail the
    threshold ``dbar_k >= -sqrt(var_k/T) * sqrt(2 log log T)``; ``SPA_l`` is an
    upper bound on the ``SPA_c`` p-value, so it is conservative for rejection.
    Report which variant was used.
    """
    if direction not in ("champion_beats_some", "nobody_beats_champion"):
        raise validate.ValidationError(
            "direction must be 'champion_beats_some' or 'nobody_beats_champion'."
        )
    mL = validate.validate_loss_array(Loss)
    iT, iM = mL.shape
    others = [c for c in range(iM) if c != champion]
    rck = _loss_center(mL, champion)[:, others]   # column k == L_k - L_champ
    if direction == "nobody_beats_champion":
        rck = -rck
    rbar = rck.mean(axis=0)
    sd = rck.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_obs_all = np.where(sd > 0, np.sqrt(iT) * rbar / sd, 0.0)
    T_obs = float(np.max(t_obs_all))

    idx = GetIndices(iT, k, B, seed=seed)
    rck_c = rck - rbar[None, :]
    T_b = np.empty(B)
    for b in range(B):
        rck_b = rck_c[idx[:, b], :]
        rbar_b = rck_b.mean(axis=0)
        sd_b = rck_b.std(axis=0, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            t_b = np.where(sd_b > 0, np.sqrt(iT) * rbar_b / sd_b, 0.0)
        T_b[b] = float(np.max(t_b))
    p = float(np.mean(T_b > T_obs))
    return {
        "p_value": p,
        "T_obs": T_obs,
        "B": B,
        "seed": seed,
        "direction": direction,
        "variant": "SPA_l (full recentring; conservative relative to SPA_c)",
        "null": (
            "max_k E[L_k - L_champ] <= 0 (champion is weakly worst)"
            if direction == "champion_beats_some"
            else "max_k E[L_champ - L_k] <= 0 (no competitor beats the champion)"
        ),
    }


# Convenience alias.
loss_differential = _loss_center
