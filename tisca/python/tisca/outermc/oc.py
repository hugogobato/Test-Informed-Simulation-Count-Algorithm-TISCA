"""Operating-characteristic summarisation (E1/P3-T2).

``summarize_ocs`` turns the per-cell summaries produced by ``run_e1`` into a
rectangular results table with one row per cell, each statistic carrying its
Monte Carlo standard error (MCSE). With ``R`` outer repetitions, the MCSE of a
rate near ``p`` is ``sqrt(p(1-p)/R)``: near 0.05 at R=5000 it is 0.0031, near
0.80 it is 0.0057 (the values quoted in REVISION_PLAN.md P3-T2).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def mcse_rate(p, R):
    """MCSE of a binomial rate estimator."""
    p = np.clip(p, 0.0, 1.0)
    return float(np.sqrt(p * (1.0 - p) / max(1, R)))


def summarize_ocs(summaries):
    """Collapse a list of per-cell summary dicts into a tidy results DataFrame.

    Columns: every factor + the operating characteristics `*_est` and `*_mcse`.
    """
    rows = []
    for s in summaries:
        R = s.get("R", np.nan)
        row = dict(s)
        for key, label in [
            ("reject_rate", "reject_rate"),
            ("ci_cover", "ci_cover"),
            ("P_J_eq_Jmax", "pJmax"),
        ]:
            v = s.get(key)
            if isinstance(v, (int, float, np.number)):
                row[f"{label}_mcse"] = mcse_rate(v, R)
        row[f"E_J"] = s.get("E_J")
        row["bias"] = s.get("bias")
        row["rmse"] = s.get("rmse")
        rows.append(row)
    df = pd.DataFrame(rows)
    return df
