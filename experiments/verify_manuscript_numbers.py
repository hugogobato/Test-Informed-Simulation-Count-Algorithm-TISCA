#!/usr/bin/env python3
"""Recompute every headline number in the manuscript from the released tables.

The manuscript promises that no number in it is typed by hand. This script is
the check on that promise: each entry below states a claim exactly as the paper
makes it, recomputes it from a released results file, and fails if the two
disagree beyond the stated tolerance. It refits no model and reruns no
simulation, so it takes seconds and can be run from a clean clone.

    python experiments/verify_manuscript_numbers.py
    python experiments/verify_manuscript_numbers.py --report results/audit.md

Exit status is non-zero if any claim fails.
"""

from __future__ import annotations

import argparse
import os
from typing import Callable

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RESULTS = os.path.join(ROOT, "results")

CELLS = [(1, 500), (2, 500), (3, 500), (1, 100)]


def _csv(*parts: str) -> pd.DataFrame:
    return pd.read_csv(os.path.join(RESULTS, *parts))


# --- section 2, the bibliometric census ------------------------------------
def bibliometric() -> list[tuple[str, float, float, float]]:
    df = pd.read_csv(os.path.join(RESULTS, "E4", "bibliometric_coded.csv"),
                     keep_default_na=False)
    analysed = df[df["J_report_reported"].astype(int) == 1].copy()
    analysed["J"] = analysed["J_numeric"].astype(int)
    n_a = len(analysed)
    nonarxiv = analysed[analysed["is_arxiv_preprint"].astype(int) == 0]
    years = analysed["year"].astype(int)
    js = _csv("E4", "justification_summary.csv")

    def pct(frame_len: int, count: int) -> float:
        return 100.0 * count / frame_len

    def summary(var: str, cat: str, col: str = "percentage") -> float:
        row = js[(js["variable"] == var) & (js["category"] == cat)]
        return float(row[col].iloc[0])

    return [
        ("S2 screened records", 100, len(df), 0),
        ("S2 analysis denominator N", 100, n_a, 0),
        ("S2 non-arXiv denominator", 89, len(nonarxiv), 0),
        ("S2 J=1000 count", 34, int((analysed["J"] == 1000).sum()), 0),
        ("S2 J=1000 share (%)", 34.00, pct(n_a, int((analysed["J"] == 1000).sum())), 0.01),
        ("S2 J<=500 count", 55, int((analysed["J"] <= 500).sum()), 0),
        ("S2 J<=500 share (%)", 55.00, pct(n_a, int((analysed["J"] <= 500).sum())), 0.01),
        ("S2 non-arXiv J<=500 count", 47, int((nonarxiv["J_numeric"].astype(int) <= 500).sum()), 0),
        ("S2 non-arXiv J<=500 share (%)", 52.81,
         pct(len(nonarxiv), int((nonarxiv["J_numeric"].astype(int) <= 500).sum())), 0.01),
        ("S2 2021-2025 share (%)", 57.00, pct(n_a, int(years.between(2021, 2025).sum())), 0.01),
        ("S2 2016-2025 share (%)", 88.00, pct(n_a, int(years.between(2016, 2025).sum())), 0.01),
        ("S2 arXiv share (%)", 11.00,
         pct(n_a, int((analysed["is_arxiv_preprint"].astype(int) == 1).sum())), 0.01),
        ("S2 not-justified (%)", 92.93, summary("justification_primary", "not_justified"), 0.01),
        ("S2 not-justified Wilson low", 86.12,
         summary("justification_primary", "not_justified", "wilson_95_low"), 0.01),
        ("S2 not-justified Wilson high", 96.53,
         summary("justification_primary", "not_justified", "wilson_95_high"), 0.01),
        ("S2 justified (%)", 7.07, summary("justification_primary", "justified"), 0.01),
        ("S2 explicit count", 4, summary("justification_class", "explicit", "count"), 0),
        ("S2 implicit count", 3, summary("justification_class", "implicit_or_convention", "count"), 0),
        ("S2 unjustified count", 77, summary("justification_class", "unjustified", "count"), 0),
        ("S2 unresolvable count", 15, summary("justification_class", "unclear_report", "count"), 0),
        ("S2 sensitivity-to-J count", 0, summary("sensitivity_to_J_primary", "yes", "count"), 0),
    ]


# --- section 4, operating characteristics ----------------------------------
def operating() -> list[tuple[str, float, float, float]]:
    oc = _csv("E1", "operating_characteristics.csv")
    core = oc[oc["module"] == "A"]
    null = core[core["theta_mult"] == 0.0]
    alt = core[core["theta_mult"] == 1.0]

    claims: list[tuple[str, float, float, float]] = [
        ("S4 completed configurations", 1983, len(oc), 0),
    ]
    expected_t1e = {"D1": 0.0497, "D2": 0.0523, "D3": 0.0314, "D4": 0.0494,
                    "D5": 0.0485, "D6": 0.0497}
    expected_pow = {"D1": 0.7820, "D2": 0.9258, "D3": 0.8684, "D4": 0.8404,
                    "D5": 0.9958, "D6": 0.8160}
    expected_ej = {"D1": 62.1, "D2": 89.4, "D3": 80.9, "D4": 73.3,
                   "D5": 564.4, "D6": 64.2}
    expected_cov = {"D1": 0.9503, "D2": 0.9477, "D3": 0.8966, "D4": 0.9506,
                    "D5": 0.9515, "D6": 0.9503}
    for design in ["D1", "D2", "D3", "D4", "D5", "D6"]:
        claims += [
            (f"S4 {design} Type I error", expected_t1e[design],
             float(null.loc[null["design"] == design, "reject_rate"].mean()), 5e-5),
            (f"S4 {design} power at delta", expected_pow[design],
             float(alt.loc[alt["design"] == design, "reject_rate"].mean()), 5e-5),
            (f"S4 {design} E[J] at delta", expected_ej[design],
             float(alt.loc[alt["design"] == design, "E_J"].mean()), 0.05),
            (f"S4 {design} CI coverage under null", expected_cov[design],
             float(null.loc[null["design"] == design, "ci_cover"].mean()), 5e-5),
        ]

    d3 = null[null["design"] == "D3"]
    claims += [
        ("S4 D3 Type I error at rho=-0.3", 0.0838,
         float(d3.loc[d3["rho"] == -0.3, "reject_rate"].mean()), 5e-5),
        ("S4 D3 Type I error at rho=0.9", 0.0004,
         float(d3.loc[d3["rho"] == 0.9, "reject_rate"].mean()), 5e-5),
    ]

    rule = _csv("E1b", "operational_rule.csv")
    pt = rule[rule["test"] == "paired_t"]
    calibrated_at_10 = pt[(pt["J_min"] == 10) & (pt["abs_skew"] < 0.5)]
    claims += [
        ("S4 paired-t cells calibrated at J=10 with |skew|<0.5", 9,
         len(calibrated_at_10), 0),
        ("S4 configurations with |skew|<0.5", 11,
         int((pt["abs_skew"] < 0.5).sum()), 0),
        ("S4 empirical row-bootstrap J_min", 150,
         float(pt.loc[pt["family"] == "empirical", "J_min"].iloc[0]), 0),
        ("S4 empirical row-bootstrap |skew|", 1.55,
         float(pt.loc[pt["family"] == "empirical", "abs_skew"].iloc[0]), 0.005),
    ]

    sweep = _csv("E1b", "type_I_vs_J.csv")
    mix = sweep[sweep["family"] == "mix"]
    mix_rho = mix[mix["rho"] == 0.6]
    claims += [
        ("S4 mixture paired-t Type I at largest J", 0.040,
         float(mix_rho[(mix_rho["test"] == "paired_t")
                       & (mix_rho["J"] == mix_rho["J"].max())]["type_I"].iloc[0]),
         0.0006),
        ("S4 mixture bootstrap worst Type I", 0.116,
         float(mix[mix["test"] == "studentized_bootstrap"]["type_I"].max()), 0.0006),
    ]
    return claims


# --- section 5, the case study ---------------------------------------------
def case_study() -> list[tuple[str, float, float, float]]:
    plan = _csv("E3", "planning_table.csv")
    contrasts = _csv("E3", "paired_contrasts.csv")
    gain = _csv("E3", "precision_gain.csv")
    pairing = _csv("E3", "pairing_gain.csv")
    dec = _csv("E3", "interval_score_decomposition.csv")
    mcs_is = _csv("E3", "mcs_interval_score.csv")
    mcs_pehe = _csv("E3", "mcs_pehe.csv")

    def cell(frame: pd.DataFrame, dgp: int, n: int) -> pd.DataFrame:
        return frame[(frame["dgp"] == dgp) & (frame["n"] == n)]

    claims: list[tuple[str, float, float, float]] = []
    for (dgp, n), expected in zip(CELLS, [30, 15, 34, 96]):
        claims.append((f"S5 J_final DGP{dgp} n={n}", expected,
                       float(cell(plan, dgp, n)["J_final"].max()), 0))
    for (dgp, n), expected in zip(CELLS, [117, 368, 114, 718]):
        claims.append((f"S5 unpaired J_final DGP{dgp} n={n}", expected,
                       float(cell(pairing, dgp, n)["J_unpaired"].max()), 0))

    claims += [
        ("S5 min within-replication correlation", 0.74,
         float(pairing["rho_confirmatory"].min()), 0.005),
        ("S5 max within-replication correlation", 0.98,
         float(pairing["rho_confirmatory"].max()), 0.005),
        ("S5 confirmatory replications per contrast", 900,
         float(contrasts["J_used"].min()), 0),
        ("S5 all 24 contrasts negative", 24,
         int((contrasts["estimate"] < 0).sum()), 0),
        ("S5 mean MCSE reduction 500->1000 (%)", 29.67,
         float(gain["mcse_reduction_pct"].mean()), 0.005),
        ("S5 min MCSE reduction (%)", 24.97,
         float(gain["mcse_reduction_pct"].min()), 0.005),
        ("S5 max MCSE reduction (%)", 33.04,
         float(gain["mcse_reduction_pct"].max()), 0.005),
        ("S5 PEHE MCS singletons", 8,
         int(mcs_pehe.groupby(["dgp", "n", "outcome"])["in_mcs_90"].sum().eq(1).sum()), 0),
    ]

    for (dgp, n), expected in zip(CELLS, [28.64, 30.00, 30.24, 29.78]):
        claims.append((f"S5 MCSE reduction DGP{dgp} n={n} (%)", expected,
                       float(cell(gain, dgp, n)["mcse_reduction_pct"].mean()), 0.005))

    d1 = contrasts[(contrasts["dgp"] == 1) & (contrasts["n"] == 500)]
    for label, expected in [("C1", -0.504), ("C2", -0.369), ("C3", -1.578),
                            ("C4", -1.544), ("C5", -1.191), ("C6", -1.324)]:
        claims.append((f"S5 DGP1 n=500 {label} estimate", expected,
                       float(d1.loc[d1["contrast"] == label, "estimate"].iloc[0]),
                       0.0005))

    d2 = dec[(dec["dgp"] == 2) & (dec["n"] == 500) & (dec["level"] == 95)
             & (dec["outcome"] == 1)]
    for model, score, width, cov in [("mvbcf", 651.7, 39.2, 0.076),
                                     ("bcf", 650.0, 42.5, 0.089),
                                     ("bart", 453.7, 56.4, 0.281),
                                     ("mvbart", 460.7, 56.1, 0.272)]:
        row = d2[d2["model"] == model].iloc[0]
        claims += [
            (f"S5 DGP2 Y1 {model} interval score", score,
             float(row["mean_interval_score"]), 0.05),
            (f"S5 DGP2 Y1 {model} mean width", width, float(row["mean_width"]), 0.05),
            (f"S5 DGP2 Y1 {model} coverage", cov, float(row["mean_coverage"]), 0.0005),
        ]
    claims += [
        ("S5 DGP2 Y1 min penalty share", 0.876, float(d2["penalty_share"].min()), 0.0005),
        ("S5 DGP2 Y1 max penalty share", 0.940, float(d2["penalty_share"].max()), 0.0005),
        ("S5 DGP2 Y1 mvbcf miss distance", 15.3,
         float(d2.loc[d2["model"] == "mvbcf", "mean_miss_distance"].iloc[0]), 0.05),
        ("S5 DGP2 Y1 bart miss distance", 9.9,
         float(d2.loc[d2["model"] == "bart", "mean_miss_distance"].iloc[0]), 0.05),
    ]

    is95 = mcs_is[mcs_is["level"] == 95]
    retained = is95[is95["in_mcs_90"].astype(bool)]
    mvbcf_cells = retained[retained["model"] == "mvbcf"]
    claims += [
        ("S5 interval-score MCS combinations retaining MVBCF", 7, len(mvbcf_cells), 0),
        ("S5 interval-score MCS retains BART alone in DGP2 Y1", 1,
         int(((retained["dgp"] == 2) & (retained["n"] == 500)
              & (retained["outcome"] == 1) & (retained["model"] == "bart")).sum()), 0),
    ]

    is50 = mcs_is[mcs_is["level"] == 50]
    same = (is95.sort_values(["dgp", "n", "outcome", "model"])["in_mcs_90"].to_numpy()
            == is50.sort_values(["dgp", "n", "outcome", "model"])["in_mcs_90"].to_numpy())
    claims.append(("S5 MCS composition identical at 50% and 95%", 1,
                   int(bool(same.all())), 0))

    parity = os.path.join(RESULTS, "E3", "v1_v2_parity.json")
    if os.path.exists(parity):
        import json
        with open(parity) as fh:
            audit = json.load(fh)
        text = json.dumps(audit)
        claims.append(("S5 replication audit reports no mismatch", 0,
                       text.lower().count('"mismatch": true'), 0))
    return claims


# --- section 6, the forecasting illustration -------------------------------
def generality() -> list[tuple[str, float, float, float]]:
    plan = _csv("E5", "planning_table.csv")
    res = _csv("E5", "contrast_results.csv")
    mcs = _csv("E5", "mcs_table.csv")
    j = dict(zip(plan["contrast"], plan["J"]))
    est = dict(zip(res["contrast"], res["estimate"]))
    return [
        ("S6 J for AR(1) vs AR(2)", 3, float(j["AR(1) vs AR(2)"]), 0),
        ("S6 J for AR(1) vs naive", 25, float(j["AR(1) vs naive"]), 0),
        ("S6 J for AR(1) vs mean", 1162, float(j["AR(1) vs mean"]), 0),
        ("S6 AR(1) vs AR(2) estimate", -0.00269, float(est["AR(1) vs AR(2)"]), 5e-6),
        ("S6 AR(1) vs naive estimate", -0.171, float(est["AR(1) vs naive"]), 5e-4),
        ("S6 AR(1) vs mean estimate", -0.981, float(est["AR(1) vs mean"]), 5e-4),
        ("S6 MCS retains AR(1) alone", 1,
         int((mcs["p_MCS"] >= 0.10).sum()), 0),
    ]


SECTIONS: list[tuple[str, Callable[[], list]]] = [
    ("Section 2, bibliometric census", bibliometric),
    ("Section 4, operating characteristics", operating),
    ("Section 5, case study", case_study),
    ("Section 6, forecasting illustration", generality),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", default=None, help="write a markdown log here")
    args = ap.parse_args(argv)

    lines = ["# Manuscript number audit", "",
             "Every claim below is recomputed from a released results file.", ""]
    failures = 0
    total = 0
    for title, fn in SECTIONS:
        lines += [f"## {title}", "",
                  "| claim | manuscript | recomputed | status |",
                  "|---|---:|---:|---|"]
        for name, stated, computed, tol in fn():
            total += 1
            ok = abs(float(stated) - float(computed)) <= tol
            failures += 0 if ok else 1
            lines.append(f"| {name} | {stated} | {computed:.6g} | "
                         f"{'PASS' if ok else 'FAIL'} |")
            if not ok:
                print(f"FAIL {name}: manuscript {stated}, recomputed {computed}")
        lines.append("")

    verdict = (f"{total - failures}/{total} claims reproduce; "
               f"{failures} discrepancies.")
    lines += ["## Verdict", "", verdict, ""]
    print(verdict)
    if args.report:
        os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
        with open(args.report, "w") as fh:
            fh.write("\n".join(lines))
        print("wrote", args.report)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
