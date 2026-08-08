#!/usr/bin/env python3
"""E3 confirmatory analysis under ANALYSIS_PLAN.md AMENDMENT 1 (plan P3-T5[b]).

This is the code that performs the amendment's reanalysis. **It re-runs nothing.**
``run_cell.R`` produced 1000 replications per cell, labelled with seeds ``0..999``
and drawn from disjoint L'Ecuyer-CMRG substreams of cell master 1. AMENDMENT 1
partitions those existing rows:

  * seeds ``0..J0-1``   -> the pilot. Used ONLY to estimate ``s_D``, ``Sigma_D``
                           and hence ``J*``. **Discarded from all inference.**
  * seeds ``J0..999``   -> the confirmatory set. Every reported estimate, CI and
                           p-value is computed here and nowhere else.

Disjoint substreams make the two blocks independent samples, which is the only
property Algorithm 1 asks of a pilot; the spec never required a particular seed
label. The reserved master-2 block (``1_000_001..``) is not used, because it
exists for one of four cells and was produced by a superseded ``run_cell.R``
(``CALIBRATION.md`` D1). ``shard_table.csv`` marks those rows ``superseded``.

Run ``collect_shards.py`` first, then::

    python experiments/E3_mvbcf_casestudy/analyse_e3.py
    python experiments/E3_mvbcf_casestudy/analyse_e3.py --pilot-size 50   # sensitivity

Outputs under ``results/E3``:

    planning_table.csv        per contrast: s_D, sigma_UB, J_precision, J_power, J_final
    planning_sensitivity.csv  the same at J0 = 100 / 50 / 25 (declared in the amendment)
    sigma_D_matrix.csv        the 6x6 pilot Sigma_D the family procedure plans from
    paired_contrasts.csv      estimate, Monte Carlo CI, then p-value (IJDA #1e order)
    family_inference.csv      Romano-Wolf adjusted p, with Bonferroni and Holm alongside
    coverage_calibration.csv  the C7 calibration family, |cov - 0.95|
    mcs_pehe.csv              Model Confidence Set over the four models on PEHE
    e3_summary.md             the numbers in report order
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import t as _t

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "tisca", "python"))

from tisca import inference, mcs as _mcs, multiplicity, planning, validate  # noqa: E402

RESULTS = os.path.join(_ROOT, "results", "E3")

# ---------------------------------------------------------------------------
# Frozen analysis constants (ANALYSIS_PLAN.md sections 1, 3, 4; not revisited here)
# ---------------------------------------------------------------------------
ALPHA = 0.05
K_FAMILY = 6                       # C1-C6 form one family
ALPHA_ADJ = ALPHA / K_FAMILY       # Bonferroni-adjusted planning level, plan 3.5
GAMMA = 0.20                       # variance-uncertainty assurance, plan section 4
POWER_TARGET = 0.80
MODE = "M1"                        # two-sided equality, both in planning and in test
J_MAX = 1000                       # declared budget cap
DELTA_FRACTION = 0.25              # delta = 0.25 * sd(tau), plan section 4
RW_B = 4999
RW_SEED = 20260806

#: The six primary contrasts. ``outcome`` selects which ``sd(tau)`` scales delta.
CONTRASTS = [
    ("C1", "MVBCF vs BCF, PEHE Y1", "mvbcf_pehe1", "bcf_pehe1", 1),
    ("C2", "MVBCF vs BCF, PEHE Y2", "mvbcf_pehe2", "bcf_pehe2", 2),
    ("C3", "MVBCF vs BART, PEHE Y1", "mvbcf_pehe1", "bart_pehe1", 1),
    ("C4", "MVBCF vs BART, PEHE Y2", "mvbcf_pehe2", "bart_pehe2", 2),
    ("C5", "MVBCF vs MVBART, PEHE Y1", "mvbcf_pehe1", "mvbart_pehe1", 1),
    ("C6", "MVBCF vs MVBART, PEHE Y2", "mvbcf_pehe2", "mvbart_pehe2", 2),
]

MODELS = ["mvbcf", "bcf", "bart", "mvbart"]

CELLS = [(1, 500), (2, 500), (3, 500), (1, 100)]


def sd_tau(dgp: int, outcome: int) -> float:
    """Population ``sd(tau_k(X))`` of the DGP, in closed form.

    ``delta`` is declared relative to the cross-sectional spread of the unit-level
    effect (plan section 4), so this is the scale the planning alternative is
    quoted in. ``run_cell.R:generate_test`` builds tau as a linear combination of
    independent ``U(0,1)`` covariates in every DGP:

        DGP1/DGP2   tau1 = 10(2 X4 + 2 X5)      tau2 = 10(1 X4 + 3 X5)
        DGP3        tau1 = 10(2 X4 + 2 X2)      tau2 = 10(1 X3 + 3 X5)

    so ``Var(tau) = sum_k (10 c_k)^2 / 12`` and the value is identical across the
    three DGPs. Taking the closed form rather than the realised test-set sd is
    deliberate: the test set is redrawn every replication, so its empirical sd is
    a random quantity and would make ``delta`` -- a *declared* planning constant --
    depend on which replication one happened to look at.
    """
    coefficients = {1: (20.0, 20.0), 2: (10.0, 30.0)}[outcome]
    return float(math.sqrt(sum(c * c for c in coefficients) / 12.0))


def load_cell(dgp: int, n: int, results_dir: str) -> pd.DataFrame:
    path = os.path.join(results_dir, f"DGP{dgp}_n{n}_replications.csv")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}; run collect_shards.py first")
    frame = pd.read_csv(path)
    seeds = frame["seed"].to_numpy()
    if len(seeds) != len(set(seeds)):
        raise SystemExit(f"{path}: duplicate seeds")
    if set(seeds) != set(range(len(seeds))):
        raise SystemExit(f"{path}: seeds are not a complete 0..{len(seeds) - 1} block")
    failures = int((frame["converged_flag"] == 0).sum())
    if failures:
        raise SystemExit(f"{path}: {failures} rows with converged_flag = 0")
    return frame.sort_values("seed").reset_index(drop=True)


def split(frame: pd.DataFrame, J0: int):
    """AMENDMENT 1: seeds 0..J0-1 size J; seeds J0..999 are tested. No overlap."""
    pilot = frame[frame["seed"] < J0]
    conf = frame[frame["seed"] >= J0]
    if len(pilot) != J0:
        raise SystemExit(f"pilot has {len(pilot)} rows, expected {J0}")
    return pilot, conf


def contrast_matrix(frame: pd.DataFrame) -> np.ndarray:
    """``(J, 6)`` paired differences ``D = L_MVBCF - L_benchmark``, lower is better."""
    return np.column_stack([frame[a].to_numpy(float) - frame[b].to_numpy(float)
                            for _, _, a, b, _ in CONTRASTS])


def plan_cell(dgp: int, pilot: pd.DataFrame, J0: int, mcse_fraction: float):
    """Solve ``J*`` per contrast from the pilot alone (plan section 4)."""
    D = contrast_matrix(pilot)
    rows = []
    for k, (cid, label, _, _, outcome) in enumerate(CONTRASTS):
        s_tau = sd_tau(dgp, outcome)
        delta = DELTA_FRACTION * s_tau
        target_mcse = mcse_fraction * s_tau
        sd_pilot = float(np.std(D[:, k], ddof=1))
        sigma_ub = validate.validate_pilot_samples(sd_pilot, J0, gamma=GAMMA)
        j_precision = planning.required_J_mcse(sigma_ub, target_mcse)
        j_power = planning.required_J_power(
            MODE, delta, sigma_ub, POWER_TARGET, ALPHA_ADJ, J_max=J_MAX)
        rows.append(dict(
            dgp=dgp, contrast=cid, label=label, J0=J0,
            sd_tau=s_tau, delta=delta, target_mcse=target_mcse,
            sd_pilot=sd_pilot, sigma_ub=sigma_ub,
            inflation=sigma_ub / sd_pilot if sd_pilot > 0 else float("nan"),
            J_precision=j_precision, J_power=j_power,
            J_contrast=max(j_precision, j_power)))
    j_final = planning.combine_J([r["J_contrast"] for r in rows], J_MAX)
    for r in rows:
        r["J_final"] = j_final
        r["capped_at_J_max"] = int(j_final >= J_MAX)
    return rows, D


def test_cell(dgp: int, n: int, conf: pd.DataFrame):
    """Every reported number comes from the confirmatory block only."""
    D = contrast_matrix(conf)
    J = D.shape[0]
    # CI at the multiplicity-adjusted level, so the interval and the family-level
    # decision are the same statement (plan section 3.5 / C4).
    crit = float(_t.ppf(1.0 - ALPHA_ADJ / 2.0, df=J - 1))

    rows = []
    for k, (cid, label, _, _, _) in enumerate(CONTRASTS):
        res = inference.paired_t(D[:, k], alternative="two-sided")
        rows.append(dict(
            dgp=dgp, n=n, contrast=cid, label=label, J_used=J,
            # Estimate and its Monte Carlo CI FIRST, p-value second (IJDA #1e).
            estimate=res["estimate"], mcse=res["mcse"],
            ci_low=res["estimate"] - crit * res["se"],
            ci_high=res["estimate"] + crit * res["se"],
            ci_level=1.0 - ALPHA_ADJ,
            sd_confirmatory=res["sd"], t=res["t"], df=res["df"],
            p_unadjusted=res["p_value"]))

    p_raw = np.array([r["p_unadjusted"] for r in rows])
    rw = multiplicity.romano_wolf_stepdown(D, B=RW_B, alpha=ALPHA_ADJ, seed=RW_SEED)
    p_bonf = multiplicity.bonferroni(p_raw)
    p_holm = multiplicity.holm(p_raw)
    fam = []
    for k, r in enumerate(rows):
        fam.append(dict(
            dgp=dgp, n=n, contrast=r["contrast"], label=r["label"],
            p_unadjusted=p_raw[k],
            # Romano-Wolf is the pre-registered family procedure; the other two are
            # the declared sensitivity, reported at the unadjusted family level ALPHA
            # because each adjusted p already carries the multiplicity.
            p_romano_wolf=float(rw["p_values"][k]),
            p_bonferroni=float(p_bonf[k]),
            p_holm=float(p_holm[k]),
            reject_romano_wolf=int(rw["p_values"][k] <= ALPHA),
            reject_bonferroni=int(p_bonf[k] <= ALPHA),
            reject_holm=int(p_holm[k] <= ALPHA)))
    return rows, fam, D


def coverage_cell(dgp: int, n: int, conf: pd.DataFrame):
    """C7, secondary family: calibration framing |cov - 0.95|, never 'higher is better'."""
    out = []
    for outcome in (1, 2):
        mv = np.abs(conf[f"mvbcf_cov95{outcome}"].to_numpy(float) - 0.95)
        for bench in ("bcf", "bart", "mvbart"):
            bm = np.abs(conf[f"{bench}_cov95{outcome}"].to_numpy(float) - 0.95)
            res = inference.paired_t(mv - bm, alternative="two-sided")
            out.append(dict(
                dgp=dgp, n=n, outcome=outcome, benchmark=bench,
                mvbcf_mean_abs_dev=float(mv.mean()),
                benchmark_mean_abs_dev=float(bm.mean()),
                estimate=res["estimate"], mcse=res["mcse"], p_value=res["p_value"],
                mvbcf_mean_width=float(conf[f"mvbcf_wid95{outcome}"].mean()),
                benchmark_mean_width=float(conf[f"{bench}_wid95{outcome}"].mean())))
    return out


def mcs_cell(dgp: int, n: int, conf: pd.DataFrame):
    """MCS over the four models on PEHE. An inference layer, never a stopping rule.

    The plan also asks for an interval-score MCS. ``run_cell.R`` records only the
    per-replication mean coverage and mean width, not the per-observation interval
    score, so that layer is NOT reconstructible from the committed shards and is
    reported as unavailable rather than approximated from coverage and width.
    """
    out = []
    for outcome in (1, 2):
        loss = np.column_stack([conf[f"{m}_pehe{outcome}"].to_numpy(float)
                                for m in MODELS])
        res = _mcs.mcs(loss, alpha=0.10, B=RW_B, seed=RW_SEED, model_names=MODELS)
        keep = set(res["included"])
        for i, m in enumerate(MODELS):
            out.append(dict(dgp=dgp, n=n, outcome=outcome, model=m,
                            mean_pehe=float(loss[:, i].mean()),
                            in_mcs_90=int(m in keep)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=RESULTS)
    ap.add_argument("--pilot-size", type=int, default=100,
                    help="J0 declared in AMENDMENT 1 (default 100)")
    ap.add_argument("--mcse-fraction", type=float, default=0.05,
                    help="precision target m as a fraction of sd(tau). The plan fixes "
                         "the SCALE (relative to sd(tau)) but not the number, so this "
                         "is an explicit input and its sensitivity is tabulated.")
    ap.add_argument("--skip-mcs", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.results_dir, exist_ok=True)
    plan_rows, contrast_rows, family_rows, cov_rows, mcs_rows, sens_rows = (
        [], [], [], [], [], [])
    sigma_frames = []

    for dgp, n in CELLS:
        frame = load_cell(dgp, n, args.results_dir)
        pilot, conf = split(frame, args.pilot_size)

        rows, D_pilot = plan_cell(dgp, pilot, args.pilot_size, args.mcse_fraction)
        for r in rows:
            r["n"] = n
        plan_rows += rows

        # The 6x6 Sigma_D the family procedure plans from: the amendment's second
        # justification for J0 = 100 is that this object is estimated from J0 rows.
        sigma = pd.DataFrame(np.cov(D_pilot, rowvar=False),
                             index=[c[0] for c in CONTRASTS],
                             columns=[c[0] for c in CONTRASTS])
        sigma.insert(0, "n", n)
        sigma.insert(0, "dgp", dgp)
        sigma.insert(0, "contrast", sigma.index)
        sigma_frames.append(sigma)

        # Pre-declared sensitivity: J* at J0 = 100 / 50 / 25, reported whichever way
        # it falls (amendment, item 3).
        for j0 in (100, 50, 25):
            if j0 > len(frame):
                continue
            sub, _ = split(frame, j0)
            srows, _ = plan_cell(dgp, sub, j0, args.mcse_fraction)
            for r in srows:
                r["n"] = n
            sens_rows += srows

        c_rows, f_rows, _ = test_cell(dgp, n, conf)
        contrast_rows += c_rows
        family_rows += f_rows
        cov_rows += coverage_cell(dgp, n, conf)
        if not args.skip_mcs:
            mcs_rows += mcs_cell(dgp, n, conf)

    def write(name, rows):
        frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
        path = os.path.join(args.results_dir, name)
        frame.to_csv(path, index=False)
        print(f"wrote {path}  ({len(frame)} rows)")
        return frame

    plan = write("planning_table.csv", plan_rows)
    write("planning_sensitivity.csv", sens_rows)
    write("sigma_D_matrix.csv", pd.concat(sigma_frames, ignore_index=True))
    contrasts = write("paired_contrasts.csv", contrast_rows)
    family = write("family_inference.csv", family_rows)
    write("coverage_calibration.csv", cov_rows)
    if mcs_rows:
        write("mcs_pehe.csv", mcs_rows)

    lines = [
        "# E3 confirmatory analysis (ANALYSIS_PLAN.md AMENDMENT 1)",
        "",
        f"J0 = {args.pilot_size} (seeds 0..{args.pilot_size - 1}), discarded from "
        f"inference. Confirmatory set = seeds {args.pilot_size}..999, "
        f"{1000 - args.pilot_size} replications per cell.",
        "",
        f"Planning: mode {MODE} two-sided, alpha_adj = 0.05/{K_FAMILY} = "
        f"{ALPHA_ADJ:.5f}, gamma = {GAMMA}, power target {POWER_TARGET}, "
        f"delta = {DELTA_FRACTION} x sd(tau), precision target "
        f"{args.mcse_fraction} x sd(tau).",
        "",
        "## Planned J* per cell",
        "",
        plan.groupby(["dgp", "n"])["J_final"].first().to_string(),
        "",
        "## Primary family: estimate (Monte Carlo CI) then p-value",
        "",
        contrasts[["dgp", "n", "contrast", "estimate", "mcse", "ci_low", "ci_high",
                   "p_unadjusted"]].to_string(index=False),
        "",
        "## Family-adjusted p-values (Romano-Wolf primary; Bonferroni/Holm sensitivity)",
        "",
        family[["dgp", "n", "contrast", "p_romano_wolf", "p_bonferroni",
                "p_holm"]].to_string(index=False),
        "",
        "Note: the interval-score MCS layer of plan section 3 is not reconstructible "
        "from the committed shards -- run_cell.R records mean coverage and mean width "
        "per replication, not the per-observation interval score -- so only the PEHE "
        "MCS is reported.",
    ]
    summary = os.path.join(args.results_dir, "e3_summary.md")
    with open(summary, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
