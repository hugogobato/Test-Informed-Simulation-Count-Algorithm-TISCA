#!/usr/bin/env python3
"""What the paired design buys in the case study, like for like.

Re-solves the same planning problem twice from the same pilot block: once on
the paired contrast standard deviation, and once on the two-sample standard
deviation an unpaired analysis of the identical losses would use. Every other
setting is held fixed at the frozen analysis constants, so the ratio isolates
the pairing and nothing else.

Also reports the within-replication correlation of the two loss columns in the
confirmatory block, which is the quantity that drives the ratio.

Outputs ``results/E3/pairing_gain.csv``. No model is refit.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2, nct

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# Frozen analysis constants (ANALYSIS_PLAN.md AMENDMENT 1).
ALPHA_ADJ = 0.05 / 6
GAMMA = 0.20
POWER_TARGET = 0.80
DELTA_FRACTION = 0.25
MCSE_FRACTION = 0.05
PILOT_SEEDS = 100
J_CAP = 20000

CELLS = [("DGP1_n500", 1, 500), ("DGP2_n500", 2, 500),
         ("DGP3_n500", 3, 500), ("DGP1_n100", 1, 100)]
BENCHMARKS = ["bcf", "bart", "mvbart"]


def sd_tau(outcome: int) -> float:
    """Population sd(tau_k(X)), closed form (see analyse_e3.py)."""
    coefficients = {1: (20.0, 20.0), 2: (10.0, 30.0)}[outcome]
    return float(np.sqrt(sum(c * c for c in coefficients) / 12.0))


def _smallest_J(delta: float, sigma: float, paired: bool) -> int:
    """Smallest J reaching POWER_TARGET at ALPHA_ADJ, two-sided.

    ``paired`` uses the one-sample t with J-1 df and se = sigma/sqrt(J);
    otherwise the two-sample t with 2(J-1) df and the same per-arm count.
    """
    for J in range(3, J_CAP + 1):
        df = (J - 1) if paired else 2 * (J - 1)
        crit = float(nct.ppf(1 - ALPHA_ADJ / 2, df=df, nc=0))
        ncp = delta / (sigma / np.sqrt(J))
        power = (1 - float(nct.cdf(crit, df=df, nc=ncp))
                 + float(nct.cdf(-crit, df=df, nc=ncp)))
        if power >= POWER_TARGET:
            return J
    raise RuntimeError("power target not reached within the cap")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replications-dir",
                    default=os.path.join(ROOT, "results", "E3", "v2"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "E3",
                                                  "pairing_gain.csv"))
    args = ap.parse_args(argv)

    inflation = float(np.sqrt((PILOT_SEEDS - 1)
                              / chi2.ppf(GAMMA, PILOT_SEEDS - 1)))
    rows = []
    for cell, dgp, n in CELLS:
        frame = pd.read_csv(os.path.join(args.replications_dir,
                                         f"{cell}_replications.csv"))
        pilot = frame[frame["seed"] < PILOT_SEEDS]
        confirm = frame[frame["seed"] >= PILOT_SEEDS]
        for bench in BENCHMARKS:
            for outcome in (1, 2):
                a_col, b_col = f"mvbcf_pehe{outcome}", f"{bench}_pehe{outcome}"
                a, b = pilot[a_col].to_numpy(), pilot[b_col].to_numpy()
                sd_paired = float(np.std(a - b, ddof=1)) * inflation
                sd_unpaired = float(np.sqrt(np.var(a, ddof=1)
                                            + np.var(b, ddof=1))) * inflation
                delta = DELTA_FRACTION * sd_tau(outcome)
                target = MCSE_FRACTION * sd_tau(outcome)
                rows.append({
                    "dgp": dgp, "n": n, "benchmark": bench, "outcome": outcome,
                    "rho_confirmatory": float(np.corrcoef(
                        confirm[a_col], confirm[b_col])[0, 1]),
                    "sd_paired": sd_paired, "sd_unpaired": sd_unpaired,
                    "J_power_paired": _smallest_J(delta, sd_paired, True),
                    "J_power_unpaired": _smallest_J(delta, sd_unpaired, False),
                    "J_precision_paired": int(np.ceil((sd_paired / target) ** 2)),
                    "J_precision_unpaired": int(
                        np.ceil((sd_unpaired / target) ** 2)),
                })

    out = pd.DataFrame(rows)
    out["J_paired"] = out[["J_power_paired", "J_precision_paired"]].max(axis=1)
    out["J_unpaired"] = out[["J_power_unpaired",
                             "J_precision_unpaired"]].max(axis=1)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)

    summary = out.groupby(["dgp", "n"]).agg(
        J_final_paired=("J_paired", "max"),
        J_final_unpaired=("J_unpaired", "max"),
        rho_min=("rho_confirmatory", "min"),
        rho_max=("rho_confirmatory", "max"))
    summary["ratio"] = (summary["J_final_unpaired"]
                        / summary["J_final_paired"]).round(2)
    print(summary.to_string())
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
