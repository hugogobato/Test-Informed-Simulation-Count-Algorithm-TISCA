"""Vectorised outer-Monte-Carlo engine and entry point (P2-T3).

Repeats the *entire* simulation-design procedure ``R`` times and returns the
unconditional operating characteristics, exactly the object P3-T2/E1 needs:

  unconditional Type I error (theta = 0), power, CI coverage of theta,
  bias / RMSE of theta_hat at the stopping size, E[J], sd(J), quantiles of J,
  P(J = J_max), and mean wall-clock.

The engine is pure NumPy/SciPy, is reproducible through per-replication
``SeedSequence`` children (seed protocol section 4), supports chunked/resumable
output, and exposes a ``--dry-run`` cost estimator through the CLI.
"""

from __future__ import annotations

import time

import numpy as np
from scipy import stats

from . import families as _fam
from .designs import make_design, paired_pvalue

DEFAULT_CONFIG = {
    "family": "normal",
    "rho": 0.0,
    "sigma_a": 1.0,
    "sigma_b": 1.0,
    "theta": 0.0,
    "sigma_D": None,          # oracle sigma (D6); overrides the derived value
    "R": 2000,                # outer Monte Carlo repetitions
    "J0": 50,                 # pilot size
    "Jmax": 1000,             # budgeted maximum confirmatory rows per rep
    "alpha": 0.05,
    "alpha_adj": None,        # multiplicity-adjusted level (None -> alpha)
    "mode": 1,                # planning/test mode (1 = two-sided), M1..M5
    "delta": 0.5,             # planning alternative (effect size, in sigma_D units)
    "power_target": 0.80,
    "gamma": 0.20,            # variance-uncertainty assurance level
    "correction": "none",
    "K": 1,                   # family size (multiplicity, Module B)
    "design": "D6",
    "matrix": None,           # empirical-family loss matrix (M, 2)
    "seed": 0,
    "fixed_J": None,          # D6: force a fixed J
    "mcse": 0.05,             # D5 precision target
    "half_width": None,
    "nmax_looks": 3,          # D2 adaptive max looks
    "batch": 50,              # rows-per-checkpoint for resumable output
    "outdir": None,
    "parquet": False,
}


def _merge(cfg):
    merged = dict(DEFAULT_CONFIG)
    merged.update(cfg or {})
    return merged


_SIGMA_D_CACHE = {}
_SIGMA_D_N = 1_000_000        # MCSE of the sd estimate ~ sigma/sqrt(2n) = 0.07%
_SIGMA_D_SEED = 20260806      # fixed, so the oracle is a constant of the harness


def sigma_D_true(cfg):
    """True paired-difference sd of the configured family. D6 plans from this.

    Closed form for the bivariate normal, where ``sigma_D = sqrt(sigma_a^2 +
    sigma_b^2 - 2 rho sigma_a sigma_b)`` is exact. **That formula is not valid for
    the other six families** and using it everywhere is what made D6 a non-oracle:

      * ``mix`` adds a 2% catastrophic component to method B, so its true sigma_D
        is 1.6x to 4.6x the formula (the formula is 37-78% too small);
      * for ``lognormal``, ``gamma``, ``beta``, ``t3`` and ``empirical`` the copula
        transform attenuates the Pearson correlation relative to the design's rank
        ``rho``, so the formula is 4-24% too small once ``rho >= 0.6``.

    An oracle that plans from a sigma 20-78% below the truth under-sizes J, so the
    "true sigma known" reference design was reporting optimistic power and E[J]
    exactly where the non-normal families were supposed to stress it. For every
    non-normal family sigma_D is therefore estimated once, to 0.07% MCSE, from a
    dedicated large draw at a fixed seed, and cached: it is a deterministic
    constant of the family/rho/scale cell, independent of ``theta``.
    """
    a, b = float(cfg["sigma_a"]), float(cfg["sigma_b"])
    r = None if cfg["rho"] is None else float(cfg["rho"])
    fam = cfg.get("family", "normal")
    if fam == "normal":
        return float(np.sqrt(a * a + b * b - 2 * r * a * b))
    mat = cfg.get("matrix")
    key = (fam, r, a, b, None if mat is None else hash(np.asarray(mat).tobytes()))
    if key not in _SIGMA_D_CACHE:
        block = _fam.sample_batch(fam, 1, _SIGMA_D_N, rho=r, sigma_a=a, sigma_b=b,
                                  theta=0.0, matrix=mat, master_seed=_SIGMA_D_SEED)
        d = block[0, :, 0] - block[0, :, 1]
        _SIGMA_D_CACHE[key] = float(d.std(ddof=1))
    return _SIGMA_D_CACHE[key]


def _draw_block(family, R, J, rho, sigma_a, sigma_b, theta, matrix, master_seed, key):
    """Vectorised (R, J, 2) draw with a seed separated by block ``key``."""
    ss = np.random.SeedSequence(master_seed, spawn_key=(key,))
    return _fam.draw_rep_losses(family, R, J, rho=rho, sigma_a=sigma_a,
                                sigma_b=sigma_b, theta=theta, matrix=matrix,
                                master_seed=ss.generate_state(1)[0])


def _prepare(cfg):
    """Fill in derived fields: sigma_D (if missing) and alpha_adj."""
    cfg = dict(cfg)
    if cfg["sigma_D"] is None:
        cfg["sigma_D"] = sigma_D_true(cfg)
    if cfg["alpha_adj"] is None:
        cfg["alpha_adj"] = cfg["alpha"] / max(1, cfg.get("K", 1))
    return cfg


def run_results(cfg):
    """Run one design over ``R`` outer repetitions; return raw per-rep arrays.

    `cfg` is mutated for derived fields for the caller; returns ``Result``.
    """
    cfg = _prepare(cfg)
    R = cfg["R"]
    J0 = int(cfg["J0"])
    Jmax = int(cfg["Jmax"])
    seed = cfg["seed"]
    # confirmatory pool must be at least the pilot for reuse designs; ensure Jmax >= J0
    Jmax = max(Jmax, J0 + 1)
    pilot = _draw_block(cfg["family"], R, J0, cfg["rho"], cfg["sigma_a"],
                        cfg["sigma_b"], cfg["theta"], cfg["matrix"], seed, key=0)
    conf = _draw_block(cfg["family"], R, Jmax, cfg["rho"], cfg["sigma_a"],
                       cfg["sigma_b"], cfg["theta"], cfg["matrix"], seed, key=1)
    draws = {"pilot": pilot, "conf": conf}
    fn = make_design(cfg["design"], cfg)
    res, meta = fn(cfg, draws, cfg["alpha"])
    return res, meta


def _ops(res, theta_true, cfg):
    """Reduce a per-rep Result to the operating-characteristic row dict."""
    R = len(res)
    if R == 0:
        raise ValueError("no outer repetitions run")
    T = res.col("T")
    p = res.col("p")
    J = res.col("J")
    th = res.col("theta_hat")
    sD = res.col("s_D")
    capped = res.col("capped")
    alpha = cfg["alpha_adj"]

    reject = np.nan_to_num(p) <= alpha
    # CI coverage of theta: theta_true in Dbar +/- t_{1-a/2,J-1} s_D/sqrt(J)
    Jf = np.asarray(J, dtype=float)
    crit = stats.t.ppf(1.0 - alpha / 2.0, np.maximum(Jf, 2) - 1.0)
    theta_hat = np.asarray(th, dtype=float)
    m = ~np.isnan(theta_hat) & ~np.isnan(sD)
    with np.errstate(divide="ignore", invalid="ignore"):
        half = np.where(m, crit * np.asarray(sD, float) / np.sqrt(Jf), np.nan)
        lo = np.where(m, theta_hat - half, np.nan)
        hi = np.where(m, theta_hat + half, np.nan)
    ci_cover = np.where(m, (theta_true >= lo) & (theta_true <= hi), False)

    nz = Jf
    J_mean = float(np.mean(nz))
    J_sd = float(np.std(nz))
    J_q = np.quantile(nz, [0.05, 0.50, 0.95])
    E_theta = float(np.nanmean(theta_hat))
    bias = E_theta - float(theta_true)
    rmse = float(np.sqrt(np.nanmean((theta_hat - theta_true) ** 2)))
    tilt = cfg["theta"] == 0.0
    return {
        "reject_rate": float(np.mean(reject)),
        "t1e_or_power": float(np.mean(reject)),
        "ci_cover": float(np.mean(ci_cover)),
        "E_theta": E_theta,
        "bias": bias,
        "rmse": rmse,
        "E_J": J_mean,
        "sd_J": J_sd,
        "q05_J": float(J_q[0]),
        "q50_J": float(J_q[1]),
        "q95_J": float(J_q[2]),
        "P_J_eq_Jmax": float(np.mean(capped)),
        "pJmax": float(np.mean(nz >= cfg["Jmax"])),
    }


def run_e1(cfg, checkpoint=None, t0=None):
    """Run one E1 cell and return (summary dict, per-rep Result, meta).

    ``checkpoint`` optionally is a callback ``checkpoint(partial_result)`` invoked
    every ``cas`` replication for resumable progress (see ``chunked_run``).
    """
    cfg = _prepare(dict(cfg))
    start = time.time()
    res, meta = run_results(cfg)
    wall = time.time() - start
    summ = _ops(res, cfg["theta"], cfg)
    summ["mean_wall"] = wall / max(1, cfg["R"])
    summ.update({k: cfg[k] for k in ("family", "rho", "sigma_a", "sigma_b", "theta",
                                     "design", "J0", "Jmax", "alpha", "alpha_adj",
                                     "delta", "power_target", "gamma", "correction",
                                     "K", "seed", "R")})
    return summ, res, meta


def chunked_run(cfg, rows_per_chunk=None, out=None, checkpoint_cb=None):
    """Resumable driver: run in row-chunks and write/checkpoint after each.

    ``cfg`` is a single-cell config. ``rows_per_chunk`` defaults to ``batch``.
    If ``out`` is a writable path, each chunk appends a parquet/CSV shard; the
    caller can later concatenate. Returns the final summary.
    """
    cfg = _prepare(dict(cfg))
    R = cfg["R"]
    chunk = rows_per_chunk or cfg.get("batch", 50)
    if not out:
        return run_e1(cfg)[0]
    import os
    import pandas as pd

    def _plan_cfg(r0, r1):
        c = dict(cfg)
        c.update(R=r1 - r0, seed=cfg["seed"] + r0,
                 outdir=None, parquet=False)
        return c

    cols = None
    for r0 in range(0, R, chunk):
        r1 = min(R, r0 + chunk)
        sub_res, _ = run_e1(_plan_cfg(r0, r1))
        row = {k: (float(v) if isinstance(v, np.floating) else v)
               for k, v in sub_res.items()}
        df = pd.DataFrame([row])
        if out.endswith(".parquet"):
            if os.path.exists(out) and os.path.getsize(out) > 0:
                df.to_parquet(out, engine="pyarrow", append=True)
            else:
                df.to_parquet(out, engine="pyarrow")
        else:
            header = not os.path.exists(out) or os.path.getsize(out) == 0
            df.to_csv(out, mode="a", header=header, index=False)
        if checkpoint_cb:
            checkpoint_cb(r1)
    # a single contiguous equivalence isn't exact for a summarised index; the
    # OCs are per-cell, so the summary is best recomputed over the full R in one
    # call. Provide it directly:
    return run_e1(dict(cfg))[0]


def dry_run(cfg):
    """Estimated cost of a cell: measure one back-of-the-envelope run at R=64 and
    extrapolate. Returns the projected wall-clock and a safe notebook count hint.
    """
    rng_cfg = _prepare(dict(cfg))
    small = dict(rng_cfg, R=64)
    t0 = time.time()
    _, meta = run_results(small)
    dt = time.time() - t0
    per_rep = dt / 64.0
    est = per_rep * rng_cfg["R"]
    core_hours = est / 3600.0 * 2.0  # assume mc.cores=2 speedup ~1.62 -> cost factor
    return {
        "measured_per_rep_s": per_rep,
        "projected_cell_s": est,
        "projected_core_hours": est / 3600.0,
        "hint": f"~{core_hours:.2f} vCPU-hours (mc.cores=2); runs in seconds on Colab",
    }
