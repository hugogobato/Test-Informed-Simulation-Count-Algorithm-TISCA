#!/usr/bin/env python3
"""What the replications past the planned count actually bought (S5.6).

Section 5.6 asks whether the additional replications delivered precision the
study had declared it needed. Answering that requires comparing the realized
Monte Carlo error against the declared target ``m = 0.05 sd(tau)``, not against
the error at some other count, and comparing the declared verdicts at the count
TISCA planned with the verdicts on the full confirmatory block.

The planned count is ``max(design J, calibration J)``, exactly as S5.2 states:
the design requirement comes from ``planning_table.csv`` (the precision and
power layers solved on the pilot) and the calibration requirement from
``results/E6/real_loss_calibration.csv`` (the smallest count at which the
paired-t on that cell's own realized losses holds its nominal level).

Both counts are read from the same seed-sorted confirmatory block, so the
planned-count rows are a prefix of the full-block rows and no model is refit.
Nothing here is a retrospective power calculation: these are statements about
what one realized dataset shows at two counts.

Outputs ``results/E3/planned_count_check.csv``.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tisca", "python"))

from tisca import mcs as _mcs, multiplicity  # noqa: E402

# Frozen analysis constants (ANALYSIS_PLAN.md AMENDMENT 1).
ALPHA_ADJ = 0.05 / 6
DELTA_FRACTION = 0.25
MCSE_FRACTION = 0.05
PILOT_SEEDS = 100
RW_B = 4999
RW_SEED = 20260806
MCS_ALPHA = 0.05

CELLS = [(1, 500), (2, 500), (3, 500), (1, 100)]
MODELS = ["mvbcf", "bcf", "bart", "mvbart"]
CONTRASTS = [
    ("C1", "mvbcf_pehe1", "bcf_pehe1", 1),
    ("C2", "mvbcf_pehe2", "bcf_pehe2", 2),
    ("C3", "mvbcf_pehe1", "bart_pehe1", 1),
    ("C4", "mvbcf_pehe2", "bart_pehe2", 2),
    ("C5", "mvbcf_pehe1", "mvbart_pehe1", 1),
    ("C6", "mvbcf_pehe2", "mvbart_pehe2", 2),
]


def sd_tau(outcome: int) -> float:
    """Population sd(tau_k(X)), closed form (see analyse_e3.py)."""
    coefficients = {1: (20.0, 20.0), 2: (10.0, 30.0)}[outcome]
    return float(math.sqrt(sum(c * c for c in coefficients) / 12.0))


def planned_counts(root: str) -> dict[tuple[int, int], int]:
    """``max(design, calibration)`` per cell, the count S5.2 plans."""
    design = pd.read_csv(os.path.join(root, "results", "E3", "planning_table.csv"))
    calib = pd.read_csv(os.path.join(root, "results", "E6",
                                     "real_loss_calibration.csv"))
    calib = calib[calib["test"] == "paired_t"]
    out = {}
    for dgp, n in CELLS:
        d = int(design[(design["dgp"] == dgp) & (design["n"] == n)]["J_final"].max())
        c = int(calib[(calib["dgp"] == dgp) & (calib["n"] == n)]["J_min"].max())
        out[(dgp, n)] = max(d, c)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replications-dir",
                    default=os.path.join(ROOT, "results", "E3", "v2"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "E3",
                                                  "planned_count_check.csv"))
    args = ap.parse_args(argv)

    planned = planned_counts(ROOT)
    rows = []
    for dgp, n in CELLS:
        frame = pd.read_csv(os.path.join(args.replications_dir,
                                         f"DGP{dgp}_n{n}_replications.csv"))
        confirm = (frame[frame["seed"] >= PILOT_SEEDS]
                   .sort_values("seed").reset_index(drop=True))
        for J, tag in ((planned[(dgp, n)], "planned"), (len(confirm), "full")):
            block = confirm.iloc[:J]
            D = np.column_stack([block[a].to_numpy(float) - block[b].to_numpy(float)
                                 for _, a, b, _ in CONTRASTS])
            rw = multiplicity.romano_wolf_stepdown(D, B=RW_B, alpha=ALPHA_ADJ,
                                                   seed=RW_SEED)
            retained = {}
            for outcome in (1, 2):
                loss = np.column_stack([block[f"{m}_pehe{outcome}"].to_numpy(float)
                                        for m in MODELS])
                res = _mcs.mcs(loss, alpha=MCS_ALPHA, B=RW_B, seed=RW_SEED,
                               model_names=MODELS)
                retained[outcome] = res["included"]
            for k, (name, _a, _b, outcome) in enumerate(CONTRASTS):
                d = D[:, k]
                mcse = float(d.std(ddof=1) / math.sqrt(J))
                target = MCSE_FRACTION * sd_tau(outcome)
                rows.append({
                    "dgp": dgp, "n": n, "J": J, "block": tag, "contrast": name,
                    "outcome": outcome, "estimate": float(d.mean()), "mcse": mcse,
                    "target_mcse": target, "target_over_mcse": target / mcse,
                    "delta": DELTA_FRACTION * sd_tau(outcome),
                    "rejected": int(np.asarray(rw["rejections"])[k]),
                    "mcs_retains_mvbcf_alone": int(retained[outcome] == ["mvbcf"]),
                })

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)

    for tag, g in out.groupby("block", sort=False):
        combos = g.drop_duplicates(["dgp", "n", "outcome"])
        print(f"{tag:8s} target/MCSE {g['target_over_mcse'].min():.2f}"
              f"--{g['target_over_mcse'].max():.2f}, "
              f"{int(g['rejected'].sum())}/24 rejected, "
              f"{int(combos['mcs_retains_mvbcf_alone'].sum())}/8 MCS singletons")
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
