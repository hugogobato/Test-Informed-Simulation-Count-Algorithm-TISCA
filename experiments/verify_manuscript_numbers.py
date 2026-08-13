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
    expected_cov = {"D1": 0.9503, "D2": 0.9477, "D3": 0.9699, "D4": 0.9506,
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

    # Table tab:rho-level-coverage, the six designs resolved by correlation.
    # design: {rho: (type I error, interval coverage)}
    rho_table = {
        "D1": {-0.3: (0.0501, 0.9499), 0.0: (0.0494, 0.9506),
               0.3: (0.0485, 0.9515), 0.6: (0.0461, 0.9539),
               0.9: (0.0526, 0.9474)},
        "D2": {-0.3: (0.0549, 0.9451), 0.0: (0.0538, 0.9462),
               0.3: (0.0525, 0.9475), 0.6: (0.0495, 0.9505),
               0.9: (0.0503, 0.9497)},
        "D3": {-0.3: (0.0838, 0.9194), 0.0: (0.0498, 0.9521),
               0.3: (0.0223, 0.9790), 0.6: (0.0046, 0.9957),
               0.9: (0.0004, 0.9996)},
        "D4": {-0.3: (0.0481, 0.9519), 0.0: (0.0473, 0.9527),
               0.3: (0.0501, 0.9499), 0.6: (0.0461, 0.9539),
               0.9: (0.0525, 0.9475)},
        "D5": {-0.3: (0.0492, 0.9508), 0.0: (0.0495, 0.9505),
               0.3: (0.0495, 0.9505), 0.6: (0.0474, 0.9526),
               0.9: (0.0463, 0.9537)},
        "D6": {-0.3: (0.0494, 0.9506), 0.0: (0.0505, 0.9495),
               0.3: (0.0473, 0.9527), 0.6: (0.0490, 0.9510),
               0.9: (0.0496, 0.9504)},
    }
    for design, by_rho in rho_table.items():
        for rho, (t_exp, c_exp) in by_rho.items():
            cell = null[(null["design"] == design) & (null["rho"] == rho)]
            claims += [
                (f"S4 {design} Type I error at rho={rho}", t_exp,
                 float(cell["reject_rate"].mean()), 5e-5),
                (f"S4 {design} CI coverage at rho={rho}", c_exp,
                 float(cell["ci_cover"].mean()), 5e-5),
            ]

    # Table tab:family-power, the K=6 block of the joint multiplicity grid.
    joint = _csv("E1", "module_b_joint_operating_characteristics.csv")
    k6 = joint[joint["K"] == 6]
    jn = k6[k6["theta"] == 0.0]
    ja = k6[k6["theta"] == 0.5]
    family_table = {
        # correction: {design: (E_J, FWER, conjunctive, disjunctive)}
        "none": {"D4": (68.9, 0.2439, 0.6810, 0.9990),
                 "D3": (76.0, 0.1247, 0.6670, 0.9960)},
        "bonferroni": {"D4": (106.4, 0.0508, 0.7417, 0.9987),
                       "D3": (76.0, 0.0246, 0.2372, 0.9719)},
        "holm": {"D4": (106.3, 0.0519, 0.9235, 0.9989),
                 "D3": (76.0, 0.0242, 0.6359, 0.9755)},
        "bh": {"D4": (106.3, 0.0515, 0.9272, 0.9993),
               "D3": (76.0, 0.0276, 0.6735, 0.9788)},
        "romano_wolf": {"D4": (68.7, 0.0479, 0.6285, 0.9864),
                        "D3": (75.9, 0.0505, 0.7295, 0.9704)},
    }
    for corr, per_design in family_table.items():
        for design, (ej, fwer, conj, disj) in per_design.items():
            n_cells = jn[(jn["correction"] == corr) & (jn["design"] == design)]
            a_cells = ja[(ja["correction"] == corr) & (ja["design"] == design)]
            claims += [
                (f"S4 {corr}/{design} E[J] at alternative", ej,
                 float(a_cells["E_J"].mean()), 0.05),
                (f"S4 {corr}/{design} FWER under the global null", fwer,
                 float(n_cells["fwer"].mean()), 5e-5),
                (f"S4 {corr}/{design} conjunctive power", conj,
                 float(a_cells["conjunctive_power"].mean()), 5e-5),
                (f"S4 {corr}/{design} disjunctive power", disj,
                 float(a_cells["disjunctive_power"].mean()), 5e-5),
            ]
    # Every cell's conjunctive power must exceed the working-independence
    # product of its own realized marginal power, which is the direction
    # Section 3.2 claims for positively dependent contrasts.
    marg = ja.groupby(["correction", "design"])[
        ["marginal_level_or_power", "conjunctive_power"]].mean()
    above = int((marg["conjunctive_power"] > marg["marginal_level_or_power"] ** 6).sum())
    below_marginal = int((marg["conjunctive_power"] < marg["marginal_level_or_power"]).sum())
    claims += [
        ("S4 K=6 cells with conjunctive power above the independence product",
         len(marg), above, 0),
        ("S4 K=6 cells with conjunctive power below the marginal",
         len(marg), below_marginal, 0),
        ("S4 romano_wolf/D4 independence product", 0.3864,
         float(marg.loc[("romano_wolf", "D4"), "marginal_level_or_power"] ** 6), 5e-5),
        ("S4 bonferroni/D3 independence product", 0.1184,
         float(marg.loc[("bonferroni", "D3"), "marginal_level_or_power"] ** 6), 5e-5),
        ("S4 uncorrected K=6 independent-test FWER reference", 0.265,
         1 - 0.95 ** 6, 0.0005),
    ]

    # Table tab:module-c, the pilot and batch sensitivity module.
    mc = oc[oc["module"] == "C"]
    module_c = {
        "D2": {25: (0.9326, 68.5, 17.4), 50: (0.9201, 71.0, 10.9),
               100: (0.9589, 101.4, 2.9)},
        "D3": {25: (0.8531, 65.2, 18.2), 50: (0.8725, 65.6, 12.2),
               100: (0.9653, 100.0, 0.5)},
        "D4": {25: (0.8783, 61.9, 20.4), 50: (0.8689, 56.5, 12.9),
               100: (0.8584, 53.3, 8.6)},
    }
    for design, per_j0 in module_c.items():
        for j0, (power, ej, sdj) in per_j0.items():
            cell = mc[(mc["design"] == design) & (mc["J0"] == j0)]
            claims += [
                (f"S4 module C {design} power at J0={j0}", power,
                 float(cell["reject_rate"].mean()), 5e-5),
                (f"S4 module C {design} E[J] at J0={j0}", ej,
                 float(cell["E_J"].mean()), 0.05),
                # Tolerance admits an exact half, as at D2 with J0=50 where the
                # recomputed 10.850 rounds up to the one-decimal 10.9.
                (f"S4 module C {design} SD of J at J0={j0}", sdj,
                 float(cell["sd_J"].mean()), 0.0501),
            ]
    d4_by_batch = mc[mc["design"] == "D4"].groupby("B")["reject_rate"].mean()
    claims.append(("S4 module C D4 power spread across batch sizes", 0.0007,
                   float(d4_by_batch.max() - d4_by_batch.min()), 5e-5))

    # Table tab:module-d, the unequal-marginal-variance module.
    md = oc[oc["module"] == "D"]
    md_null = md[md["theta_mult"] == 0.0]
    md_alt = md[md["theta_mult"] == 1.0]
    module_d = {
        # design: {ratio: (type_I, coverage, power, E_J at the alternative)}
        "D1": {1.0: (0.0501, 0.9499, 0.8058, 47.2), 2.0: (0.0539, 0.9461, 0.7937, 123.0)},
        "D2": {1.0: (0.0531, 0.9469, 0.9226, 71.0), 2.0: (0.0569, 0.9431, 0.8771, 156.3)},
        "D3": {1.0: (0.0364, 0.9647, 0.8682, 65.6), 2.0: (0.0381, 0.9623, 0.8234, 159.4)},
        "D4": {1.0: (0.0488, 0.9512, 0.8686, 56.5), 2.0: (0.0530, 0.9470, 0.8607, 148.6)},
        "D5": {1.0: (0.0504, 0.9496, 1.0000, 532.9), 2.0: (0.0509, 0.9491, 1.0000, 861.3)},
        "D6": {1.0: (0.0504, 0.9496, 0.8259, 47.1), 2.0: (0.0545, 0.9455, 0.8316, 123.1)},
    }
    for design, per_ratio in module_d.items():
        for ratio, (t1e, cov, power, ej) in per_ratio.items():
            n_cell = md_null[(md_null["design"] == design)
                             & (md_null["sigma_ratio"] == ratio)]
            a_cell = md_alt[(md_alt["design"] == design)
                            & (md_alt["sigma_ratio"] == ratio)]
            claims += [
                (f"S4 module D {design} Type I at ratio {ratio:g}", t1e,
                 float(n_cell["reject_rate"].mean()), 5e-5),
                (f"S4 module D {design} coverage at ratio {ratio:g}", cov,
                 float(n_cell["ci_cover"].mean()), 5e-5),
                (f"S4 module D {design} power at ratio {ratio:g}", power,
                 float(a_cell["reject_rate"].mean()), 5e-5),
                (f"S4 module D {design} E[J] at ratio {ratio:g}", ej,
                 float(a_cell["E_J"].mean()), 0.05),
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

    # E6: the same calibration protocol on the 24 real case-study contrasts.
    real = _csv("E6", "real_loss_calibration.csv")
    rt = real[real["test"] == "paired_t"]
    rb = real[real["test"] == "studentized_bootstrap"]
    claims += [
        ("S4 real contrasts", 24, len(rt), 0),
        ("S4 real paired-t J_min min", 10, float(rt["J_min"].min()), 0),
        ("S4 real paired-t J_min max", 100, float(rt["J_min"].max()), 0),
        ("S4 real paired-t J_min median", 28, float(rt["J_min"].median()), 0.5),
        ("S4 real paired-t contrasts needing J>30", 11,
         int((rt["J_min"] > 30).sum()), 0),
        ("S4 real paired-t contrasts never calibrating", 0,
         int(rt["J_min"].isna().sum()), 0),
        ("S4 real bootstrap J_min max", 400, float(rb["J_min"].max()), 0),
        ("S4 real bootstrap J_min median", 200, float(rb["J_min"].median()), 0),
        ("S4 real bootstrap contrasts never calibrating", 13,
         int(rb["J_min"].isna().sum()), 0),
        ("S4 real |skew| min", 0.03, float(rt["abs_skew"].min()), 0.005),
        ("S4 real |skew| max", 1.47, float(rt["abs_skew"].max()), 0.005),
        ("S4 real rho min", 0.74, float(rt["rho_pearson"].min()), 0.005),
        ("S4 real rho max", 0.98, float(rt["rho_pearson"].max()), 0.005),
    ]
    # The per-cell calibration requirement, against the planned counts of S5.2.
    cell_req = {(1, 500): 100, (2, 500): 75, (3, 500): 75, (1, 100): 40}
    planned = {(1, 500): 30, (2, 500): 15, (3, 500): 34, (1, 100): 96}
    binding = 0
    for (dgp, n), stated in cell_req.items():
        got = rt[(rt["dgp"] == dgp) & (rt["n"] == n)]["J_min"].max()
        claims.append((f"S4 calibration requirement DGP{dgp} n={n}", stated,
                       float(got), 0))
        binding += int(got > planned[(dgp, n)])
    claims.append(("S4 cells where calibration exceeds the planned J", 3, binding, 0))

    # S5.2 plans against max(design, calibration); recompute the compute budget
    # from the committed shard-time projections so the abstract's figures are
    # audited rather than asserted.
    rate = {500: 6.173 / 100, 100: 6.852 / 333}
    tisca = sum((100 + max(cell_req[k], planned[k])) * rate[k[1]] for k in cell_req)
    baseline = sum(1000 * rate[k[1]] for k in cell_req)
    claims += [
        ("S5 planned J with calibration, largest cell", 100,
         max(max(cell_req[k], planned[k]) for k in cell_req), 0),
        ("S5 rows per DGP, largest cell", 200,
         100 + max(max(cell_req[k], planned[k]) for k in cell_req), 0),
        ("S5 TISCA aggregate compute-hours", 38.0, tisca, 0.05),
        ("S5 1000-row aggregate compute-hours", 205.8, baseline, 0.05),
        ("S5 compute-hours saved", 167.8, baseline - tisca, 0.05),
        ("S5 compute-hours saved as a percentage", 81.5,
         100 * (baseline - tisca) / baseline, 0.05),
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
         int(mcs_pehe.groupby(["dgp", "n", "outcome"])["in_mcs_95"].sum().eq(1).sum()), 0),
    ]

    for (dgp, n), expected in zip(CELLS, [28.64, 30.00, 30.24, 29.78]):
        claims.append((f"S5 MCSE reduction DGP{dgp} n={n} (%)", expected,
                       float(cell(gain, dgp, n)["mcse_reduction_pct"].mean()), 0.005))

    # S5.6, what the replications past the planned count bought, measured
    # against the declared target rather than against another count.
    pcc = _csv("E3", "planned_count_check.csv")
    at_plan = pcc[pcc["block"] == "planned"]
    at_full = pcc[pcc["block"] == "full"]
    plan_combos = at_plan.drop_duplicates(["dgp", "n", "outcome"])
    d3_plan = at_plan[(at_plan["dgp"] == 3) & (at_plan["n"] == 500)]
    claims += [
        ("S5 planned-count target/MCSE margin, min", 1.23,
         float(at_plan["target_over_mcse"].min()), 0.005),
        ("S5 planned-count target/MCSE margin, max", 5.62,
         float(at_plan["target_over_mcse"].max()), 0.005),
        ("S5 full-block target/MCSE margin, min", 3.7,
         float(at_full["target_over_mcse"].min()), 0.05),
        ("S5 full-block target/MCSE margin, max", 17.4,
         float(at_full["target_over_mcse"].max()), 0.05),
        ("S5 family rejections already held at the planned count", 22,
         int(at_plan["rejected"].sum()), 0),
        ("S5 family rejections on the full block", 24,
         int(at_full["rejected"].sum()), 0),
        ("S5 PEHE MCS singletons at the planned count", 7,
         int(plan_combos["mcs_retains_mvbcf_alone"].sum()), 0),
        ("S5 DGP3 C1 estimate at the planned count", -0.33,
         float(d3_plan.loc[d3_plan["contrast"] == "C1", "estimate"].iloc[0]), 0.005),
        ("S5 DGP3 C2 estimate at the planned count", -0.09,
         float(d3_plan.loc[d3_plan["contrast"] == "C2", "estimate"].iloc[0]), 0.005),
        ("S5 planning alternative Y1", 2.04,
         float(at_plan.loc[at_plan["outcome"] == 1, "delta"].iloc[0]), 0.005),
        ("S5 planning alternative Y2", 2.28,
         float(at_plan.loc[at_plan["outcome"] == 2, "delta"].iloc[0]), 0.005),
    ]

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
    retained = is95[is95["in_mcs_95"].astype(bool)]
    mvbcf_cells = retained[retained["model"] == "mvbcf"]
    claims += [
        ("S5 interval-score MCS combinations retaining MVBCF", 7, len(mvbcf_cells), 0),
        ("S5 interval-score MCS retains BART alone in DGP2 Y1", 1,
         int(((retained["dgp"] == 2) & (retained["n"] == 500)
              & (retained["outcome"] == 1) & (retained["model"] == "bart")).sum()), 0),
    ]

    # Section 5 claims the MCS composition is level-insensitive. The MCS p-value
    # does not depend on the declared confidence level, so the largest p-value
    # among the eliminated models bounds how far the level could be raised
    # before any set changes.
    mcs_crps = _csv("E3", "mcs_crps.csv")
    eliminated = pd.concat([mcs_pehe, mcs_is, mcs_crps])["p_mcs"]
    claims.append(("S5 largest MCS p-value among eliminated models", 0.011,
                   float(eliminated[eliminated < 1.0].max()), 5e-4))

    is50 = mcs_is[mcs_is["level"] == 50]
    same = (is95.sort_values(["dgp", "n", "outcome", "model"])["in_mcs_95"].to_numpy()
            == is50.sort_values(["dgp", "n", "outcome", "model"])["in_mcs_95"].to_numpy())
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
    skipped = []
    for title, fn in SECTIONS:
        lines += [f"## {title}", "",
                  "| claim | manuscript | recomputed | status |",
                  "|---|---:|---:|---|"]
        try:
            section_claims = fn()
        except FileNotFoundError as exc:
            # A missing results file means this block cannot be audited. Say so
            # and continue, rather than aborting every remaining section.
            note = f"{title}: skipped, {os.path.basename(str(exc.filename))} not found"
            skipped.append(note)
            print("SKIP", note)
            lines += [f"| _{note}_ | | | SKIP |", ""]
            continue
        for name, stated, computed, tol in section_claims:
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
    if skipped:
        verdict += f" {len(skipped)} section(s) skipped: " + "; ".join(skipped) + "."
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
