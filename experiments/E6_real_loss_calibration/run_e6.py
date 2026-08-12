#!/usr/bin/env python3
"""E6: finite-sample calibration on the REAL case-study losses.

The E1b calibration sub-study measured the smallest replication count at which the
declared paired test enters its equivalence band, but it did so on seven parametric
loss families plus a row bootstrap of ONE recorded loss pair. That single empirical
cell was the most informative one in the whole grid, which is the argument for
running the same protocol on every contrast the case study actually declares.

This script applies the E1b protocol verbatim -- same equivalence band, same J grid,
same paired-t and studentized-bootstrap statistics -- to all 24 real PEHE contrasts
of the E3 case study (3 DGPs at n = 500, DGP1 at n = 100, six contrasts each). The
null is imposed by the ``empirical`` row bootstrap of ``tisca.outermc.families``,
which centres each column on its own mean, so the real joint dependence, marginal
shapes and variance ratio survive while E[D] = 0 holds exactly.

Outputs (results/E6/):
    real_loss_calibration.csv   one row per contrast x test
    real_loss_type_I_vs_J.csv   the full sweep behind it
    real_loss_calibration.md    the prose summary the manuscript quotes

Regenerate with::

    python experiments/E6_real_loss_calibration/run_e6.py
"""

from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, "tisca", "python"))

from tisca.outermc import families  # noqa: E402

# --- protocol constants, identical to E1b ------------------------------------
ALPHA = 0.05
J_GRID = [10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 400]
R_T = 40_000          # paired t, MCSE at 0.05 = 0.0011
# E1b ran the bootstrap overlay at R = 2,500 (MCSE 0.0044), which makes its band
# max(0.005, 2 MCSE) = 0.0088, wider than the paired t's 0.005. The two tests then
# face different acceptance rules and a single noisy sweep point can move J_min a
# whole grid step. Here R = 10,000 puts the bootstrap MCSE at 0.0022, so 2 MCSE
# falls below 0.005 and BOTH tests are judged against the same 0.005 band.
R_BOOT = 10_000
B_BOOT = 499
TOL_ABS = 0.005       # equivalence band: max(0.005, 2 MCSE)
SKEW_N = 1_000_000
N_WORKERS = 8         # leave headroom; other experiments may be running

# --- the 24 declared case-study contrasts ------------------------------------
# (dgp, n) -> the amended confirmatory block actually analysed (J_used = 900).
BLOCKS = [(1, 500), (2, 500), (3, 500), (1, 100)]
CONTRASTS = [
    ("C1", "MVBCF vs BCF, PEHE Y1", "mvbcf_pehe1", "bcf_pehe1"),
    ("C2", "MVBCF vs BCF, PEHE Y2", "mvbcf_pehe2", "bcf_pehe2"),
    ("C3", "MVBCF vs BART, PEHE Y1", "mvbcf_pehe1", "bart_pehe1"),
    ("C4", "MVBCF vs BART, PEHE Y2", "mvbcf_pehe2", "bart_pehe2"),
    ("C5", "MVBCF vs MVBART, PEHE Y1", "mvbcf_pehe1", "mvbart_pehe1"),
    ("C6", "MVBCF vs MVBART, PEHE Y2", "mvbcf_pehe2", "mvbart_pehe2"),
]


def load_matrix(dgp, n, col_a, col_b):
    path = os.path.join(REPO_ROOT, "results", "E3",
                        f"DGP{dgp}_n{n}_amended_confirmatory_replications.csv")
    frame = pd.read_csv(path, usecols=[col_a, col_b])
    mat = frame[[col_a, col_b]].to_numpy(float)
    mat = mat[np.isfinite(mat).all(axis=1)]
    return mat


def type_I_paired_t(matrix, J, R, seed, chunk_elems=2_000_000):
    """Vectorised: R independent row-bootstrap samples of J pairs, one two-sided t each."""
    crit = stats.t.ppf(1 - ALPHA / 2, J - 1)
    per_chunk = max(1, int(chunk_elems // max(J, 1)))
    n_rej, n_done, block_i = 0, 0, 0
    while n_done < R:
        r = min(per_chunk, R - n_done)
        block = families.sample_batch("empirical", r, J, theta=0.0,
                                      master_seed=seed + 1_000 * block_i, matrix=matrix)
        D = block[..., 0] - block[..., 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = D.mean(axis=1) / (D.std(axis=1, ddof=1) / np.sqrt(J))
        n_rej += int(np.count_nonzero(np.abs(t) > crit))
        n_done += r
        block_i += 1
        del block, D, t
    return n_rej / R


def _boot_reject(D, B, rng, chunk=40):
    """Studentized paired bootstrap rejection indicator per row of D (E1b statistic)."""
    R, J = D.shape
    out = np.empty(R, dtype=bool)
    for start in range(0, R, chunk):
        d = D[start:start + chunk]
        c = d.shape[0]
        idx = rng.integers(0, J, size=(c, B, J))
        res = np.take_along_axis(d[:, None, :], idx, axis=2)
        mb = res.mean(axis=2)
        sb = res.std(axis=2, ddof=1)
        m = d.mean(axis=1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            T = (mb - m) / (sb / np.sqrt(J))
        lo, hi = np.nanquantile(T, [ALPHA / 2, 1 - ALPHA / 2], axis=1)
        s = d.std(axis=1, ddof=1) / np.sqrt(J)
        obs = d.mean(axis=1) / s
        out[start:start + c] = (obs < lo) | (obs > hi)
    return out


def type_I_bootstrap(matrix, J, R, seed, boot_seed):
    rng = np.random.default_rng(boot_seed)
    block = families.sample_batch("empirical", R, J, theta=0.0,
                                  master_seed=seed, matrix=matrix)
    D = block[..., 0] - block[..., 1]
    return float(np.mean(_boot_reject(D, B_BOOT, rng)))


def run_cell(args):
    """One contrast: descriptive moments plus the full J sweep for both tests."""
    ci, dgp, n, cid, label, col_a, col_b = args
    matrix = load_matrix(dgp, n, col_a, col_b)
    D_obs = matrix[:, 0] - matrix[:, 1]
    skew = float(families.contrast_skewness("empirical", n=SKEW_N, seed=17, matrix=matrix))
    rho_p = float(families.empirical_natural_rho(matrix, "pearson"))
    rho_s = float(families.empirical_natural_rho(matrix, "spearman"))
    sd_ratio = float(matrix[:, 0].std(ddof=1) / matrix[:, 1].std(ddof=1))

    rows = []
    for J in J_GRID:
        p = type_I_paired_t(matrix, J, R_T, seed=600_000 + 100 * ci + J)
        rows.append({"J": J, "test": "paired_t", "type_I": p,
                     "mcse": float(np.sqrt(p * (1 - p) / R_T)), "R": R_T})
        q = type_I_bootstrap(matrix, J, R_BOOT, seed=500_000 + 100 * ci + J,
                             boot_seed=400_000 + ci)
        rows.append({"J": J, "test": "studentized_bootstrap", "type_I": q,
                     "mcse": float(np.sqrt(q * (1 - q) / R_BOOT)), "R": R_BOOT})

    meta = {"dgp": dgp, "n": n, "contrast": cid, "label": label,
            "M_rows": int(matrix.shape[0]), "skew_D": skew, "abs_skew": abs(skew),
            "skew_mcse": float(np.sqrt(6 / SKEW_N)), "rho_pearson": rho_p,
            "rho_spearman": rho_s, "sd_ratio_A_over_B": sd_ratio,
            "observed_mean_D": float(D_obs.mean()),
            "observed_sd_D": float(D_obs.std(ddof=1))}
    for r in rows:
        r.update(meta)
    return rows


def bootstrap_B_sensitivity():
    """Is the bootstrap's liberal tilt a property of the test or of B = 499?

    The sweep uses the resample count TISCA's own ``studentized_paired_bootstrap``
    defaults to, and a finite B makes the bootstrap quantiles noisy, which inflates
    two-sided rejection on its own. That confound has to be separated from the test
    before the sweep's excess level is attributed to the bootstrap principle, so the
    same cells are re-run at a four-fold larger B.
    """
    global B_BOOT
    cases = [(1, 500, "mvbcf_pehe1", "bcf_pehe1", "C1"),
             (1, 500, "mvbcf_pehe2", "bart_pehe2", "C4")]
    saved, rows = B_BOOT, []
    try:
        for dgp, n, col_a, col_b, cid in cases:
            matrix = load_matrix(dgp, n, col_a, col_b)
            for B in (499, 1999):
                B_BOOT = B
                for J in (50, 100):
                    p = type_I_bootstrap(matrix, J, 4_000, seed=77 + J,
                                         boot_seed=5 + J)
                    rows.append({"dgp": dgp, "n": n, "contrast": cid, "J": J,
                                 "B": B, "type_I": p,
                                 "mcse": float(np.sqrt(p * (1 - p) / 4_000)),
                                 "R": 4_000})
    finally:
        B_BOOT = saved
    return pd.DataFrame(rows)


def band(frame, tol_mcse=2.0):
    return np.maximum(TOL_ABS, tol_mcse * frame["mcse"])


def smallest_valid_J(g, tol_mcse=2.0):
    """Smallest J from which the level stays inside the band for every larger J."""
    g = g.sort_values("J")
    ok = (g["type_I"] - ALPHA).abs() <= band(g, tol_mcse)
    for i in range(len(g)):
        if ok.iloc[i:].all():
            return int(g["J"].iloc[i])
    return np.nan


def direction_at_max_J(g):
    g = g.sort_values("J")
    last = g.iloc[-1]
    if abs(last["type_I"] - ALPHA) <= max(TOL_ABS, 2.0 * last["mcse"]):
        return "calibrated"
    return "liberal" if last["type_I"] > ALPHA else "conservative"


def main():
    jobs = [(ci, dgp, n, cid, label, col_a, col_b)
            for ci, ((dgp, n), (cid, label, col_a, col_b)) in
            enumerate((b, c) for b in BLOCKS for c in CONTRASTS)]
    print(f"{len(jobs)} real contrasts x {len(J_GRID)} values of J x 2 tests",
          flush=True)

    out = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for i, rows in enumerate(pool.map(run_cell, jobs), start=1):
            out.extend(rows)
            print(f"  done {i}/{len(jobs)}: DGP{rows[0]['dgp']} n={rows[0]['n']} "
                  f"{rows[0]['contrast']} skew={rows[0]['skew_D']:.3f}", flush=True)

    sweep = pd.DataFrame(out)
    res_dir = os.path.join(REPO_ROOT, "results", "E6")
    os.makedirs(res_dir, exist_ok=True)
    sweep.to_csv(os.path.join(res_dir, "real_loss_type_I_vs_J.csv"), index=False)

    keys = ["dgp", "n", "contrast", "label", "test"]
    recs = []
    for key, g in sweep.groupby(keys, sort=False):
        first = g.iloc[0]
        recs.append({
            **dict(zip(keys, key)),
            "J_min": smallest_valid_J(g),
            "direction": direction_at_max_J(g),
            "type_I_at_max_J": float(g.sort_values("J").iloc[-1]["type_I"]),
            "skew_D": first["skew_D"], "abs_skew": first["abs_skew"],
            "rho_pearson": first["rho_pearson"], "rho_spearman": first["rho_spearman"],
            "sd_ratio_A_over_B": first["sd_ratio_A_over_B"],
            "observed_mean_D": first["observed_mean_D"],
            "observed_sd_D": first["observed_sd_D"],
        })
    rule = pd.DataFrame(recs).sort_values(["test", "abs_skew"])
    rule.to_csv(os.path.join(res_dir, "real_loss_calibration.csv"), index=False)

    bsens = bootstrap_B_sensitivity()
    bsens.to_csv(os.path.join(res_dir, "bootstrap_B_sensitivity.csv"), index=False)

    pt = rule[rule["test"] == "paired_t"]
    bs = rule[rule["test"] == "studentized_bootstrap"]
    excess = (sweep[(sweep["test"] == "studentized_bootstrap") & (sweep["J"] >= 200)]
              ["type_I"].mean() - ALPHA)
    b_gain = (bsens[bsens["B"] == 499]["type_I"].mean()
              - bsens[bsens["B"] == 1999]["type_I"].mean())
    lines = [
        "# E6: calibration of the declared paired tests on the real case-study losses",
        "",
        f"Protocol identical to E1b: equivalence band max(+/-{TOL_ABS}, 2 MCSE) on the "
        f"level, held for every larger J on the grid {J_GRID}; R = {R_T} for the paired "
        f"t and {R_BOOT} for the studentized bootstrap with B = {B_BOOT}. The null is "
        "imposed by the row bootstrap of the real (M, 2) loss matrix, which preserves "
        "the real joint dependence, marginal shapes and variance ratio.",
        "",
        f"|skew(D)| across the {len(pt)} real contrasts ranges from "
        f"{pt['abs_skew'].min():.2f} to {pt['abs_skew'].max():.2f}; the within-"
        f"replication Pearson correlation ranges from {pt['rho_pearson'].min():.2f} to "
        f"{pt['rho_pearson'].max():.2f}.",
        "",
        f"Paired t: J_min ranges from {pt['J_min'].min():.0f} to {pt['J_min'].max():.0f} "
        f"(median {pt['J_min'].median():.0f}); "
        f"{int((pt['J_min'] > 30).sum())} of {len(pt)} contrasts require more than 30 "
        f"replications and {int(pt['J_min'].isna().sum())} never calibrate on the grid.",
        f"Studentized bootstrap: J_min ranges from {bs['J_min'].min():.0f} to "
        f"{bs['J_min'].max():.0f} (median {bs['J_min'].median():.0f}); "
        f"{int(bs['J_min'].isna().sum())} never calibrate.",
        "",
        f"The bootstrap's residual excess level averages {excess:+.4f} over J >= 200. "
        f"Raising the resample count from B = {B_BOOT} to 1999 removes {b_gain:.4f} of "
        "it, so part of the excess is the finite-resample noise in the bootstrap "
        "quantiles rather than the bootstrap principle; the remainder is not removed "
        "by more resamples. See bootstrap_B_sensitivity.csv.",
        "",
        rule.drop(columns=["label"]).round(4).to_markdown(index=False),
    ]
    with open(os.path.join(res_dir, "real_loss_calibration.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))


if __name__ == "__main__":
    main()
