"""Genuine joint Module-B outer-Monte-Carlo engine.

Each outer repetition contains a ``(J, K+1)`` loss block with one shared method-A
column and K distinct benchmark columns.  All K paired contrasts are tested in
one family.  This module deliberately leaves the scalar :mod:`engine` untouched,
so the committed 1,983-cell E1 results remain a provenance artifact.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass

import numpy as np
from scipy import stats

from .. import multiplicity
from .designs import _inflate, solve_J
from .joint_families import contrasts_from_losses, draw_joint_losses


DEFAULT_JOINT_CONFIG = {
    "family": "normal",
    "rho": 0.0,
    "sigma_a": 1.0,
    "sigma_b": 1.0,
    "theta": 0.0,
    "R": 2000,
    "J0": 50,
    "Jmax": 1000,
    "alpha": 0.05,
    "alpha_adj": None,
    "delta": 0.5,
    "power_target": 0.80,
    "gamma": 0.20,
    "correction": "none",
    "K": 1,
    "design": "D4",
    "matrix": None,
    "seed": 0,
    "B": 999,
    "fixed_J": None,
    "benchmark_residual_corr": 0.0,
    "outer_seed_block": 32,
    "outer_chunk": 32,
}


def _merge(config: dict | None) -> dict:
    out = dict(DEFAULT_JOINT_CONFIG)
    out.update(config or {})
    K = int(out["K"])
    if out.get("alpha_adj") is None:
        out["alpha_adj"] = multiplicity.planning_alpha(
            out.get("correction", "none"), K, alpha=float(out["alpha"]), r=1
        )[0]
    if out["design"] not in ("D3", "D4"):
        raise ValueError("joint Module B preserves the canonical D3/D4 design factor")
    return out


def _theta_vector(value, K: int) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.ndim == 0:
        return np.full(K, float(value))
    value = value.reshape(-1)
    if value.size != K:
        raise ValueError(f"theta must be scalar or length {K}")
    return value


def _delta_vector(value, K: int) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    if value.ndim == 0:
        return np.full(K, float(value))
    value = value.reshape(-1)
    if value.size != K:
        raise ValueError(f"delta must be scalar or length {K}")
    return value


def _seed_for_outer_block(master_seed: int, stream: int, block_id: int) -> int:
    ss = np.random.SeedSequence(int(master_seed), spawn_key=(int(stream), int(block_id)))
    return int(ss.generate_state(1, dtype=np.uint64)[0])


def draw_outer_range(config: dict, start: int, stop: int, J: int, *, stream: int) -> np.ndarray:
    """Addressable joint draws, invariant to outer-repetition shard offsets.

    Randomness is allocated in fixed outer blocks.  A request for repetitions
    ``17:41`` generates and slices the same fixed blocks as a request for ``0:64``;
    changing worker counts, chunks or shard boundaries therefore cannot change a
    row.
    """
    cfg = _merge(config)
    start, stop, J = int(start), int(stop), int(J)
    if start < 0 or stop < start:
        raise ValueError("invalid outer range")
    if start == stop:
        return np.empty((0, J, int(cfg["K"]) + 1))
    block_size = int(cfg.get("outer_seed_block", 32))
    pieces = []
    first = start // block_size
    last = (stop - 1) // block_size
    for block_id in range(first, last + 1):
        block_start = block_id * block_size
        seed = _seed_for_outer_block(cfg["seed"], stream, block_id)
        full = draw_joint_losses(
            cfg["family"], block_size, J, int(cfg["K"]), rho=cfg.get("rho"),
            sigma_a=float(cfg.get("sigma_a", 1.0)),
            sigma_b=float(cfg.get("sigma_b", 1.0)), theta=cfg.get("theta", 0.0),
            matrix=cfg.get("matrix"),
            benchmark_residual_corr=float(cfg.get("benchmark_residual_corr", 0.0)),
            seed=seed,
        )
        lo = max(start, block_start) - block_start
        hi = min(stop, block_start + block_size) - block_start
        pieces.append(full[lo:hi])
    return np.concatenate(pieces, axis=0)


def _paired_stats(D: np.ndarray, J_used: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Studentized paired statistics for ``D`` of shape ``(R, Jmax, K)``."""
    D = np.asarray(D, dtype=float)
    R, Jmax, _K = D.shape
    J_used = np.asarray(J_used, dtype=int).reshape(R)
    mask = np.arange(Jmax)[None, :] < J_used[:, None]
    selected = np.where(mask[..., None], D, 0.0)
    sums = selected.sum(axis=1)
    means = sums / J_used[:, None]
    squares = (selected * selected).sum(axis=1)
    var = np.clip(
        (squares - J_used[:, None] * means * means) / (J_used[:, None] - 1),
        0.0, None,
    )
    sds = np.sqrt(var)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.sqrt(J_used)[:, None] * means / sds
    tstat = np.where(np.isfinite(tstat), tstat, 0.0)
    return tstat, means, sds


def _welch_stats(losses: np.ndarray, J_used: np.ndarray) -> tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Legacy D3 Welch statistics for K benchmarks sharing one observed control."""
    R, Jmax, _ = losses.shape
    K = losses.shape[2] - 1
    J_used = np.asarray(J_used, dtype=int).reshape(R)
    mask = np.arange(Jmax)[None, :] < J_used[:, None]
    A = losses[..., 0]
    B = losses[..., 1:]
    mean_a = np.where(mask, A, 0.0).sum(axis=1) / J_used
    mean_b = np.where(mask[..., None], B, 0.0).sum(axis=1) / J_used[:, None]
    s2a = np.where(mask, (A - mean_a[:, None]) ** 2, 0.0).sum(axis=1) / (J_used - 1)
    s2b = np.where(mask[..., None], (B - mean_b[:, None, :]) ** 2, 0.0).sum(axis=1) / (
        J_used[:, None] - 1
    )
    va = s2a[:, None] / J_used[:, None]
    vb = s2b / J_used[:, None]
    diff = mean_a[:, None] - mean_b
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = diff / np.sqrt(va + vb)
        df = (va + vb) ** 2 / (
            va ** 2 / (J_used[:, None] - 1) + vb ** 2 / (J_used[:, None] - 1)
        )
    df = np.where(np.isfinite(df) & (df > 0), df, J_used[:, None] - 1)
    tstat = np.where(np.isfinite(tstat), tstat, 0.0)
    # Preserve the scalar D3 reporting scale so K=1 reduces to the old harness.
    report_sd = np.sqrt((s2a[:, None] + s2b) / 2.0)
    return tstat, diff, report_sd, df, np.sqrt(s2a[:, None] + s2b)


def _plan_J(config: dict, pilot_or_conf: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cfg = config
    R = pilot_or_conf.shape[0]
    K = int(cfg["K"])
    fixed = cfg.get("fixed_J")
    if fixed is not None:
        J = np.full(R, int(fixed), dtype=int)
        if np.any(J < 2) or np.any(J > int(cfg["Jmax"])):
            raise ValueError("fixed_J must lie in [2, Jmax]")
        return J, np.ones((R, K), dtype=bool)

    delta = _delta_vector(cfg.get("delta", 0.5), K)
    Jmax = int(cfg["Jmax"])
    if cfg["design"] == "D4":
        Dp = contrasts_from_losses(pilot_or_conf)
        s = Dp.std(axis=1, ddof=1)
        s = np.where(np.isfinite(s) & (s > 0), s, 1e-12)
        sigma_eff = _inflate(s, int(cfg["J0"]), float(cfg.get("gamma", 0.20)))
        flat_J, flat_ok = solve_J(
            float(cfg["power_target"]), np.broadcast_to(delta, s.shape).ravel(),
            sigma_eff.ravel(), float(cfg["alpha_adj"]), mode=1, Jmax=Jmax,
        )
    else:
        # D3 plans from the unpaired marginal variances at nominal alpha, as v1 did.
        lead = pilot_or_conf[:, :int(cfg["J0"]), :]
        s_a = lead[..., 0].std(axis=1, ddof=1)
        s_b = lead[..., 1:].std(axis=1, ddof=1)
        sigma_eff = np.sqrt(s_a[:, None] ** 2 + s_b ** 2)
        sigma_eff = np.where(np.isfinite(sigma_eff) & (sigma_eff > 0), sigma_eff, 1e-12)
        flat_J, flat_ok = solve_J(
            float(cfg["power_target"]), np.broadcast_to(delta, sigma_eff.shape).ravel(),
            sigma_eff.ravel(), float(cfg["alpha"]), mode=1, Jmax=Jmax,
        )
    per_contrast = np.asarray(flat_J).reshape(R, K)
    ok = np.asarray(flat_ok).reshape(R, K)
    return np.clip(per_contrast.max(axis=1), 2, Jmax).astype(int), ok


def _adjust_rows(p_raw: np.ndarray, method: str) -> np.ndarray:
    method = method.lower()
    if method == "none":
        return p_raw.copy()
    return np.vstack([multiplicity.p_adjust(row, method) for row in p_raw])


def _rw_seed(master_seed: int, outer_index: int) -> int:
    ss = np.random.SeedSequence(int(master_seed), spawn_key=(70_001, int(outer_index)))
    return int(ss.generate_state(1, dtype=np.uint64)[0] % np.uint64(2**63 - 1))


@dataclass
class JointRawResult:
    outer_index: np.ndarray
    J: np.ndarray
    theta_hat: np.ndarray
    report_sd: np.ndarray
    p_raw: np.ndarray
    p_adjusted: np.ndarray
    reject: np.ndarray
    ci_cover: np.ndarray
    planning_ok: np.ndarray


def _evaluate_chunk(cfg: dict, start: int, stop: int) -> JointRawResult:
    R = stop - start
    K = int(cfg["K"])
    Jmax = int(cfg["Jmax"])
    if cfg["design"] == "D4" and cfg.get("fixed_J") is None:
        pilot = draw_outer_range(cfg, start, stop, int(cfg["J0"]), stream=10)
        J, planning_ok = _plan_J(cfg, pilot)
    else:
        # D3 reuses the leading confirmatory rows as its internal pilot.  A fixed-J
        # D4 validation cell likewise needs no pilot draw.
        J = None
        planning_ok = None

    conf = draw_outer_range(cfg, start, stop, Jmax, stream=20)
    if J is None:
        J, planning_ok = _plan_J(cfg, conf)
    D = contrasts_from_losses(conf)
    if cfg["design"] == "D4":
        tstat, theta_hat, report_sd = _paired_stats(D, J)
        p_raw = 2.0 * stats.t.sf(np.abs(tstat), J[:, None] - 1)
    else:
        tstat, theta_hat, report_sd, df, _plan_scale = _welch_stats(conf, J)
        p_raw = 2.0 * stats.t.sf(np.abs(tstat), df)
    p_raw = np.where(np.isfinite(p_raw), p_raw, 1.0)

    correction = str(cfg.get("correction", "none")).lower()
    if correction in ("romano_wolf", "romano-wolf", "rw") and K > 1:
        adjusted = np.empty_like(p_raw)
        for r in range(R):
            # The SAME row indices resample every contrast column inside the
            # tested implementation, preserving the K-dimensional dependence.
            rw = multiplicity.romano_wolf_stepdown(
                D[r, :J[r], :], B=int(cfg.get("B", 999)), alpha=float(cfg["alpha"]),
                seed=_rw_seed(cfg["seed"], start + r),
            )
            adjusted[r] = rw["p_values"]
    else:
        # Romano-Wolf at K=1 is exactly the marginal test.  This convention also
        # makes K=1 D3 reduce to the scalar legacy result rather than switching its
        # statistic merely because the correction label says Romano-Wolf.
        method = "none" if correction in ("romano_wolf", "romano-wolf", "rw") else correction
        adjusted = _adjust_rows(p_raw, method)
    reject = adjusted <= float(cfg["alpha"])

    crit = stats.t.ppf(1.0 - float(cfg["alpha_adj"]) / 2.0, J - 1)
    half = crit[:, None] * report_sd / np.sqrt(J)[:, None]
    theta_true = _theta_vector(cfg.get("theta", 0.0), K)
    ci_cover = (theta_true[None, :] >= theta_hat - half) & (
        theta_true[None, :] <= theta_hat + half
    )
    return JointRawResult(
        outer_index=np.arange(start, stop, dtype=int), J=J, theta_hat=theta_hat,
        report_sd=report_sd, p_raw=p_raw, p_adjusted=adjusted, reject=reject,
        ci_cover=ci_cover, planning_ok=planning_ok,
    )


def _concat(parts: list[JointRawResult]) -> JointRawResult:
    fields = JointRawResult.__dataclass_fields__
    return JointRawResult(**{
        name: np.concatenate([getattr(part, name) for part in parts], axis=0)
        for name in fields
    })


def run_joint_results(config: dict, *, outer_start: int = 0,
                      outer_stop: int | None = None) -> JointRawResult:
    """Run raw joint outer repetitions, optionally over an addressable slice."""
    cfg = _merge(config)
    start = int(outer_start)
    stop = int(cfg["R"] if outer_stop is None else outer_stop)
    if stop < start:
        raise ValueError("outer_stop precedes outer_start")
    chunk = max(1, int(cfg.get("outer_chunk", 32)))
    parts = []
    for lo in range(start, stop, chunk):
        parts.append(_evaluate_chunk(cfg, lo, min(stop, lo + chunk)))
    if not parts:
        raise ValueError("no outer repetitions requested")
    return _concat(parts)


def rate_mcse(p: float, R: int) -> float:
    """Bernoulli-rate MCSE required by the Module-B reporting contract."""
    if not np.isfinite(p):
        return float("nan")
    p = float(np.clip(p, 0.0, 1.0))
    return float(np.sqrt(p * (1.0 - p) / int(R)))


def summarize_joint(config: dict, raw: JointRawResult, *, wall_seconds: float) -> dict:
    cfg = _merge(config)
    R, K = raw.reject.shape
    theta = _theta_vector(cfg.get("theta", 0.0), K)
    true_alt = ~np.isclose(theta, 0.0)
    true_null = ~true_alt
    marginal = raw.reject.mean(axis=0)
    marginal_mean = float(marginal.mean())

    fwer = float(np.mean(raw.reject[:, true_null].any(axis=1))) if true_null.any() else np.nan
    conjunctive = (
        float(np.mean(raw.reject[:, true_alt].all(axis=1))) if true_alt.any() else np.nan
    )
    disjunctive = (
        float(np.mean(raw.reject[:, true_alt].any(axis=1))) if true_alt.any() else np.nan
    )
    discoveries = raw.reject.sum(axis=1)
    false_discoveries = raw.reject[:, true_null].sum(axis=1)
    fdp = false_discoveries / np.maximum(discoveries, 1)
    is_bh = str(cfg.get("correction", "none")).lower() in ("bh", "fdr", "benjamini_hochberg")
    fdr = float(fdp.mean()) if is_bh else np.nan

    J = raw.J.astype(float)
    q = np.quantile(J, [0.05, 0.50, 0.95])
    errors = raw.theta_hat - theta[None, :]
    ci_coverage = float(raw.ci_cover.mean())
    row = {
        "module": "B_joint",
        "family": cfg["family"],
        "K": K,
        "correction": cfg["correction"],
        "rho": cfg.get("rho"),
        "theta": float(theta[0]) if np.all(theta == theta[0]) else np.nan,
        "theta_vector": json.dumps(theta.tolist(), separators=(",", ":")),
        "design": cfg["design"],
        "R": R,
        "J0": int(cfg["J0"]),
        "Jmax": int(cfg["Jmax"]),
        "alpha": float(cfg["alpha"]),
        "alpha_adj": float(cfg["alpha_adj"]),
        "marginal_level_or_power": marginal_mean,
        "marginal_mcse": rate_mcse(marginal_mean, R),
        "marginal_estimand": "level" if true_null.all() else "power",
        "marginal_by_contrast": json.dumps(marginal.tolist(), separators=(",", ":")),
        "marginal_mcse_by_contrast": json.dumps(
            [rate_mcse(value, R) for value in marginal], separators=(",", ":")
        ),
        "fwer": fwer,
        "fwer_mcse": rate_mcse(fwer, R),
        "conjunctive_power": conjunctive,
        "conjunctive_power_mcse": rate_mcse(conjunctive, R),
        "disjunctive_power": disjunctive,
        "disjunctive_power_mcse": rate_mcse(disjunctive, R),
        "fdr": fdr,
        "fdr_mcse": rate_mcse(fdr, R),
        "mean_fdp": fdr,
        "E_J": float(J.mean()),
        "sd_J": float(J.std(ddof=0)),
        "q05_J": float(q[0]),
        "q50_J": float(q[1]),
        "q95_J": float(q[2]),
        "P_J_eq_Jmax": float(np.mean(J == int(cfg["Jmax"]))),
        "ci_coverage": ci_coverage,
        "ci_coverage_mcse": rate_mcse(ci_coverage, R),
        "bias": float(errors.mean()),
        "rmse": float(np.sqrt(np.mean(errors ** 2))),
        "mean_wall": float(wall_seconds / max(R, 1)),
        "seed": int(cfg["seed"]),
        "bootstrap_B": int(cfg.get("B", 0)),
        "benchmark_residual_corr": float(cfg.get("benchmark_residual_corr", 0.0)),
        "n_true_alternatives": int(true_alt.sum()),
        "n_true_nulls": int(true_null.sum()),
        "planning_target_met_all": float(np.mean(raw.planning_ok.all(axis=1))),
    }
    # Stable rectangular per-contrast columns make the singular headline field
    # auditable without parsing JSON.  Unused columns are explicit NA.
    for k in range(6):
        value = float(marginal[k]) if k < K else np.nan
        row[f"marginal_contrast_{k + 1}"] = value
        row[f"marginal_contrast_{k + 1}_mcse"] = rate_mcse(value, R)
    return row


def run_joint_cell(config: dict, *, outer_start: int = 0,
                   outer_stop: int | None = None) -> tuple[dict, JointRawResult]:
    """Run and summarize one joint Module-B cell."""
    cfg = _merge(config)
    started = time.perf_counter()
    raw = run_joint_results(cfg, outer_start=outer_start, outer_stop=outer_stop)
    wall = time.perf_counter() - started
    return summarize_joint(cfg, raw, wall_seconds=wall), raw


__all__ = [
    "DEFAULT_JOINT_CONFIG",
    "JointRawResult",
    "draw_outer_range",
    "rate_mcse",
    "run_joint_cell",
    "run_joint_results",
    "summarize_joint",
]
