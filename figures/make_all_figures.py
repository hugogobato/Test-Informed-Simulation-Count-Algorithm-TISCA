#!/usr/bin/env python3
"""Regenerate every manuscript figure from the committed results CSVs (P4-T2).

Single entry point; the acceptance criterion for P4-T2 is that every figure in
the paper is regenerable by ONE script from the committed results tables:

    results/E4/bibliometric_coded.csv        -> Fig 1, 2, 3a, 3b
    results/E3/planning_table.csv            -> Fig 4  (planning curves)
    results/E3/paired_contrasts.csv          -> Fig 5  (forest plots)
    results/E1/operating_characteristics.csv -> Fig 6, 7
    results/E1b/type_I_vs_J.csv              -> Fig 8 (summary) and 8a-8c
    results/E1b/operational_rule.csv         -> Fig 8 (left panel)
    results/E6/real_loss_calibration.csv     -> Fig E6 (real-loss calibration)
    results/E3/v2/DGP*_n*_replications.csv   -> Fig 9  (MCS paths)
    results/E3/mcs_interval_score.csv        -> Fig 10 (interval-score MCS)

Every curve drawn here is recomputed from the committed tables; nothing is
hand-typed into a caption.

Usage::

    python figures/make_all_figures.py                 # write to ./figures
    python figures/make_all_figures.py --paper-dir ..  # also publish Fig*.png

No model is refit and no new simulation is run by this script.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd
from scipy.stats import nct as _nct

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "tisca", "python"))
from tisca import mcs as _mcs  # noqa: E402

FIGURES = HERE
RESULTS = os.path.join(ROOT, "results")

# Print-size legibility (SoutoNeto section 4) and the colour-blind-safe palette
# used by every other figure in the revision package.
matplotlib.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
    "axes.spines.right": False, "figure.autolayout": False,
})
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]

# ---- frozen analysis constants (ANALYSIS_PLAN.md AMENDMENT 1) -------------
ALPHA = 0.05
K_FAMILY = 6
ALPHA_ADJ = ALPHA / K_FAMILY
GAMMA = 0.20
POWER_TARGET = 0.80
MODE = "M1"
DELTA_FRACTION = 0.25
MCSE_FRACTION = 0.05
J_MAX = 1000
RW_B = 4999
RW_SEED = 20260806
MCS_ALPHA = ALPHA  # 95% MCS, the level every other interval and test in the paper uses

CELLS = [(1, 500), (2, 500), (3, 500), (1, 100)]
CONTRASTS = [
    ("C1", "MVBCF vs BCF, PEHE Y1", "mvbcf_pehe1", "bcf_pehe1", 1),
    ("C2", "MVBCF vs BCF, PEHE Y2", "mvbcf_pehe2", "bcf_pehe2", 2),
    ("C3", "MVBCF vs BART, PEHE Y1", "mvbcf_pehe1", "bart_pehe1", 1),
    ("C4", "MVBCF vs BART, PEHE Y2", "mvbcf_pehe2", "bart_pehe2", 2),
    ("C5", "MVBCF vs MVBART, PEHE Y1", "mvbcf_pehe1", "mvbart_pehe1", 1),
    ("C6", "MVBCF vs MVBART, PEHE Y2", "mvbcf_pehe2", "mvbart_pehe2", 2),
]
MODELS = ["mvbcf", "bcf", "bart", "mvbart"]
DESIGNS = ["D1", "D2", "D3", "D4", "D5", "D6"]
PAPER_FIGURES = ["Fig1.png", "Fig2.png", "Fig3a.png", "Fig3b.png", "Fig4.png",
                 "Fig5.png", "Fig6.png", "Fig7.png", "Fig8.png", "Fig8a.png",
                 "Fig8b.png", "Fig8c.png", "Fig9.png", "Fig10.png",
                 "Fig_E3_precision_gain.png", "Fig_E4_justification.png",
                 "Fig_E6_real_loss_calibration.png"]

RETURN = object()


def download(fn: str) -> None:
    try:
        from google.colab import files  # type: ignore
        files.download(fn)
        print("Downloaded:", fn)
    except Exception:
        print("(not on Colab / download skipped):", fn)


def _sd_tau(outcome: int) -> float:
    """Population sd(tau_k(X)), closed form (see analyse_e3.py)."""
    coefficients = {1: (20.0, 20.0), 2: (10.0, 30.0)}[outcome]
    return float(np.sqrt(sum(c * c for c in coefficients) / 12.0))


def _power_2sided(J: int, delta: float, sigma: float, alpha: float) -> float:
    crit = float(_nct.ppf(1 - alpha / 2, df=J - 1, nc=0))
    ncp = np.sqrt(J) * delta / sigma
    return (1 - float(_nct.cdf(crit, df=J - 1, nc=ncp))
            + float(_nct.cdf(-crit, df=J - 1, nc=ncp)))


# ---------------------------------------------------------------------------
# Fig 1-3: bibliometric study (P3-T1 / results/E4)
# ---------------------------------------------------------------------------
def fig_bibliometrics() -> None:
    bib = os.path.join(RESULTS, "E4", "bibliometric_coded.csv")
    df = pd.read_csv(bib, keep_default_na=False)
    analysed = df[df["J_report_reported"].astype(int) == 1].copy()
    analysed["J"] = analysed["J_numeric"].astype(int)
    N_A = len(analysed)
    nonarxiv = analysed[analysed["is_arxiv_preprint"].astype(int) == 0]

    # Fig 1 -- J distribution, largest share first.
    vc = analysed["J"].value_counts(normalize=True) * 100
    vc = vc.sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.barh([str(int(x)) for x in vc.index], vc.values, color="#4C72B0")
    ax.invert_yaxis()
    ax.set_xlabel("Percentage of analysed studies (%)")
    ax.set_ylabel("Number of simulations $J$")
    ax.set_title(f"Distribution of $J$ across analysed studies ($N = {N_A}$)")
    for i, v in enumerate(vc.values):
        ax.text(v + 0.4, i, f"{v:.1f}%", va="center")
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig1.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)

    # Fig 2 -- publisher / venue distribution (top 12).
    vp = analysed["publisher"].mask(analysed["publisher"].eq(""),
                                    analysed["venue"]).replace("", "unlisted")
    vc2 = vp.value_counts(normalize=True) * 100
    vc2 = vc2.head(12).sort_values()
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.barh(vc2.index, vc2.values, color="#C44E52")
    ax.set_xlabel("Percentage of analysed studies (%)")
    ax.set_ylabel("Publisher / venue")
    ax.set_title(f"Publisher distribution ($N = {N_A}$)")
    for i, v in enumerate(vc2.values):
        ax.text(v + 0.3, i, f"{v:.1f}%", va="center")
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig2.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)

    # Fig 3a -- year distribution (all analysed).
    years = analysed["year"].dropna().astype(int)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    years.value_counts().sort_index().plot(kind="bar", ax=ax, color="#55A868",
                                           width=0.75)
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Number of analysed studies")
    ax.set_title(f"Year distribution ($N = {len(years)}$)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig3a.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)

    # Fig 3b -- non-arXiv subset (v1 caption wrongly said 89; it is 88).
    years_na = nonarxiv["year"].dropna().astype(int)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    years_na.value_counts().sort_index().plot(kind="bar", ax=ax,
                                              color="#8172B3", width=0.75)
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Number of studies (non-arXiv)")
    ax.set_title(f"Year distribution, non-arXiv ($N = {len(years_na)}$)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig3b.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)

    # Fig E4 -- justification audit.  This belongs in the appendix because it
    # documents the coding denominator and the sample composition rather than
    # reporting a methodological result.
    summary_path = os.path.join(RESULTS, "E4", "justification_summary.csv")
    summary = pd.read_csv(summary_path, keep_default_na=False)
    # The coding itself is three-way: unresolvable count-to-source mappings are
    # coded unjustified at the row level, not pooled here.
    classes = ["explicit", "implicit_or_convention", "unjustified"]
    labels = ["Explicit criterion", "Implicit or convention", "No stated reason"]
    rows = summary[(summary["variable"] == "justification_class") &
                   (summary["category"].isin(classes))].copy()
    counts = {row["category"]: int(row["count"]) for _, row in rows.iterrows()}
    values = [counts.get(c, 0) for c in classes]
    # Read the denominator from the summary file so the figure cannot drift
    # away from the coded corpus.
    denom = int(rows["denominator"].iloc[0])
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    bars = ax.barh(labels[::-1], values[::-1],
                   color=["#C44E52", "#DD8452", "#4C72B0"])
    ax.set_xlabel(f"Number of records (denominator = {denom})")
    ax.set_title("Reason given for the reported outer replication count")
    for bar, value in zip(bars, values[::-1]):
        ax.text(value + 0.8, bar.get_y() + bar.get_height() / 2,
                f"{value} ({100 * value / denom:.1f}%)", va="center")
    ax.set_xlim(0, max(values) * 1.22)
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig_E4_justification.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig E3: marginal precision gain from additional case-study replications
# ---------------------------------------------------------------------------
def fig_precision_gain() -> None:
    j_grid = np.arange(100, 1001, 100)
    colours = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]
    labels = [f"DGP{dgp}, $n={n}$" for dgp, n in CELLS]
    normalized = []
    for dgp, n in CELLS:
        path = os.path.join(RESULTS, "E3",
                            f"DGP{dgp}_n{n}_confirmatory_replications.csv")
        frame = pd.read_csv(path).sort_values("seed").reset_index(drop=True)
        mean_mcses = []
        for j in j_grid:
            block = frame.iloc[:j]
            mcses = []
            for _, _, mv_col, benchmark_col, _ in CONTRASTS:
                diff = block[mv_col].to_numpy(float) - block[benchmark_col].to_numpy(float)
                mcses.append(np.std(diff, ddof=1) / np.sqrt(j))
            mean_mcses.append(float(np.mean(mcses)))
        normalized.append(np.asarray(mean_mcses) / mean_mcses[0])

    shard_path = os.path.join(ROOT, "experiments", "E3_mvbcf_casestudy",
                              "shard_table_v2.csv")
    shard = pd.read_csv(shard_path)
    shard["hours_per_replication"] = shard["projected_hours"] / shard["replication_count"]
    rates = shard.groupby("n")["hours_per_replication"].median().to_dict()

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0))
    ax = axes[0]
    for colour, label, curve in zip(colours, labels, normalized):
        ax.plot(j_grid, curve, marker="o", ms=4, lw=1.8, color=colour, label=label)
    ax.plot(j_grid, np.sqrt(100 / j_grid), color="0.25", ls="--", lw=1.2,
            label="square-root benchmark")
    ax.set_xlabel("cumulative replications $J$")
    ax.set_ylabel("mean MCSE / mean MCSE at $J=100$")
    ax.set_title("Marginal precision")
    ax.set_xticks(j_grid)
    ax.set_ylim(0, 1.08)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for colour, label, (_, n) in zip(colours, labels, CELLS):
        ax.plot(j_grid, j_grid * rates[n], marker="o", ms=4, lw=1.8,
                color=colour, label=label)
    ax.set_xlabel("cumulative replications $J$")
    ax.set_ylabel("projected aggregate compute-hours")
    ax.set_title("Cost of additional replications")
    ax.set_xticks(j_grid)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig_E3_precision_gain.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 4: planned power and precision curves (headline cell, DGP1 n=500)
# ---------------------------------------------------------------------------
def fig_case_study_planning() -> None:
    plan = pd.read_csv(os.path.join(RESULTS, "E3", "planning_table.csv"))
    cell = plan[(plan["dgp"] == 1) & (plan["n"] == 500)]
    J_grid = np.arange(5, 401)

    # The third panel is deliberately a pilot-sensitivity diagnostic.  It keeps
    # the figure tied to the committed planning tables and makes clear that the
    # plotted power crossing is not a claim that one particular grid point is
    # a universal minimum.
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0),
                             gridspec_kw={"width_ratios": [1.05, 1.05, 0.9]})

    ax = axes[0]
    for k, (cid, _, _, _, outcome) in enumerate(CONTRASTS):
        r = cell[cell["contrast"] == cid]
        delta = DELTA_FRACTION * _sd_tau(outcome)
        sig_lo = float(r["sd_pilot"].iloc[0])
        sig_hi = float(r["sigma_ub"].iloc[0])
        p_lo = np.array([_power_2sided(int(j), delta, sig_lo, ALPHA_ADJ)
                         for j in J_grid])
        p_hi = np.array([_power_2sided(int(j), delta, sig_hi, ALPHA_ADJ)
                         for j in J_grid])
        ax.plot(J_grid, p_lo, color=PALETTE[k % 6], lw=1.4,
                label=f"{cid}  $J_{{power}}$={int(r['J_power'].iloc[0])}")
        ax.fill_between(J_grid, p_lo, p_hi, color=PALETTE[k % 6], alpha=0.15,
                        linewidth=0)
    ax.axhline(POWER_TARGET, color="crimson", lw=1.4, ls="--",
               label="target 0.80")
    ax.set_xlabel("replications $J$")
    ax.set_ylabel("planned power (analytic noncentral $t$)")
    ax.set_title("Power layer at $\\alpha = 0.05/6$")
    ax.legend(frameon=False, fontsize=9, ncol=2)

    ax = axes[1]
    for k, (cid, _, _, _, outcome) in enumerate(CONTRASTS):
        r = cell[cell["contrast"] == cid]
        sig = float(r["sigma_ub"].iloc[0])
        target = MCSE_FRACTION * _sd_tau(outcome)
        mcse = sig / np.sqrt(J_grid)
        ax.plot(J_grid, mcse, color=PALETTE[k % 6], lw=1.4,
                label=f"{cid} ($\\sigma_{{\\mathrm{{UB}}}}$)")
        ax.axhline(target, color=PALETTE[k % 6], lw=0.9, ls=":")
    ax.axvline(float(cell["J_final"].iloc[0]), color="black", lw=1.5, ls="--",
               label=f"$J_{{final}}$ = {int(cell['J_final'].iloc[0])}")
    ax.set_xlabel("replications $J$")
    ax.set_ylabel("planned MCSE  $\\sigma_\\mathrm{UB}/\\sqrt{J}$")
    ax.set_title("Precision layer")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=9, ncol=2)

    # Pilot-size sensitivity is the available stability diagnostic for the
    # case-study planning tables.  The confirmatory E3 analysis uses J0=100;
    # J0=25 and J0=50 are declared sensitivity analyses, not new fits.
    sens = pd.read_csv(os.path.join(RESULTS, "E3", "planning_sensitivity.csv"))
    sens = sens[(sens["dgp"] == 1) & (sens["n"] == 500)]
    summary = (sens.groupby("J0", as_index=False)["J_final"]
               .agg(["min", "median", "max"]).reset_index())
    ax = axes[2]
    ax.plot(summary["J0"], summary["median"], marker="o", color="#4C72B0",
            label="median $J_{final}$")
    ax.fill_between(summary["J0"], summary["min"], summary["max"],
                    color="#4C72B0", alpha=0.18, label="min--max across contrasts")
    ax.axhline(float(cell["J_final"].iloc[0]), color="black", lw=1.2, ls="--",
               label="$J_{final}$ at $J_0=100$")
    ax.set_xlabel("pilot size $J_0$")
    ax.set_ylabel("planned $J_{final}$")
    ax.set_title("Pilot-size sensitivity")
    # Anchor at zero. The panel's message is that the planned count barely moves
    # with the pilot size; an autoscaled axis magnifies a spread of a few
    # replications into a dramatic-looking curve and says the opposite.
    ax.set_ylim(0, 60)
    ax.set_xticks(sorted(summary["J0"].unique()))
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig4.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 5: paired-contrast forest plots (estimate + MC CI first, p second)
# ---------------------------------------------------------------------------
def fig_contrast_forest() -> None:
    cc = pd.read_csv(os.path.join(RESULTS, "E3", "paired_contrasts.csv"))
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2))
    for ax, (dgp, n) in zip(axes.ravel(), CELLS):
        g = cc[(cc["dgp"] == dgp) & (cc["n"] == n)].iloc[::-1]
        y = np.arange(len(g))
        ax.errorbar(g["estimate"], y, fmt="o", color="#4C72B0", ms=5,
                    xerr=[g["estimate"] - g["ci_low"],
                          g["ci_high"] - g["estimate"]],
                    capsize=3, lw=1.4)
        ax.axvline(0, color="0.35", lw=1.0)
        ax.set_yticks(y)
        ax.set_yticklabels([lab for _, lab, _, _, _ in
                            CONTRASTS][::-1], fontsize=11)
        ax.set_title(f"DGP{dgp}, $n={n}$: PEHE paired differences")
        ax.set_xlabel("MVBCF $-$ benchmark (lower is better)")
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig5.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 6: operating characteristics of the TISCA v2 procedure (E1 module A)
# ---------------------------------------------------------------------------
def fig_operating_characteristics() -> None:
    oc = pd.read_csv(os.path.join(RESULTS, "E1", "operating_characteristics.csv"))
    oc = oc[oc["module"] == "A"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.0))

    ax = axes[0]
    null = oc[oc["theta_mult"] == 0.0]
    for i, d in enumerate(DESIGNS):
        v = null.loc[null["design"] == d, "reject_rate"]
        ax.scatter(np.full(len(v), i) + np.random.default_rng(i).normal(0, 0.07, len(v)),
                   v, s=9, alpha=0.35, color=PALETTE[i], linewidths=0)
        ax.scatter([i], [v.mean()], marker="_", s=700, color="black", zorder=3)
    ax.axhline(ALPHA, color="crimson", lw=1.2, ls="--", label="nominal 0.05")
    ax.set_xticks(range(len(DESIGNS)))
    ax.set_xticklabels(DESIGNS)
    ax.set_ylabel("unconditional Type I error")
    ax.set_title("Achieved level at $\\theta = 0$")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    alt = oc[oc["theta_mult"] == 1.0]
    for i, d in enumerate(DESIGNS):
        s = alt[alt["design"] == d]
        ax.scatter(s["E_J"], s["reject_rate"], s=16, alpha=0.55,
                   color=PALETTE[i], label=d, linewidths=0)
    ax.axhline(POWER_TARGET, color="crimson", lw=1.2, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("$E[J]$ (log scale)")
    ax.set_ylabel("achieved power")
    ax.set_title("Power at $\\theta = \\delta$ vs replications spent")
    ax.legend(frameon=False, ncol=2, fontsize=10)

    ax = axes[2]
    piv = (alt[alt["rho"].notna()]
           .pivot_table(index="rho", columns="design", values="E_J",
                        aggfunc="mean"))
    for i, d in enumerate(DESIGNS):
        if d in piv:
            ax.plot(piv.index, piv[d], marker="o", color=PALETTE[i], label=d)
    ax.set_xlabel("design correlation $\\rho$")
    ax.set_ylabel("$E_J$ at the planning alternative")
    ax.set_title("Cost of ignoring the pairing")
    ax.legend(frameon=False, ncol=1, fontsize=9)

    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig6.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 7: distribution of the chosen J over repetitions of the whole procedure
# ---------------------------------------------------------------------------
def fig_J_distribution() -> None:
    oc = pd.read_csv(os.path.join(RESULTS, "E1", "operating_characteristics.csv"))
    oc = oc[oc["module"] == "A"].copy()
    alt = oc[oc["theta_mult"] == 1.0].copy()

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    for i, d in enumerate(DESIGNS):
        s = alt[alt["design"] == d]
        xs = np.full(len(s), i) + np.random.default_rng(i).uniform(-0.16, 0.16, len(s))
        ax.errorbar(xs, s["q50_J"], yerr=[s["q50_J"] - s["q05_J"],
                                          s["q95_J"] - s["q50_J"]],
                    fmt="none", ecolor=PALETTE[i], alpha=0.6, elinewidth=1.3,
                    capsize=0)
        ax.scatter(xs, s["q50_J"], s=12, color=PALETTE[i], alpha=0.8,
                   linewidths=0)
        ax.scatter([i], [s["q50_J"].mean()], marker="_", s=500, color="black",
                   zorder=3)
    ax.set_xticks(range(len(DESIGNS)))
    ax.set_xticklabels(DESIGNS)
    ax.set_ylabel("$J$ chosen by the procedure")
    ax.set_title("Distribution of $J$ over $R$ repetitions of the whole "
                 "procedure (median $\\pm$ [0.05, 0.95] quantiles)")
    ax.set_yscale("log")
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig7.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 8: Type I error vs J, panelled by |skew(D)| (E1b)
# ---------------------------------------------------------------------------
_TESTS = [("paired_t", "paired $t$"),
          ("studentized_bootstrap", "stud. bootstrap")]


def _e1b_cells(sweep: pd.DataFrame) -> pd.DataFrame:
    """The 15 observed (family, rho) calibration cells, ordered by skewness."""
    cells = sweep[["family", "rho"]].drop_duplicates()
    sk = sweep[["family", "rho", "abs_skew"]].drop_duplicates()
    return cells.merge(sk, on=["family", "rho"], how="left").sort_values("abs_skew")


def _e1b_cell_rows(sweep: pd.DataFrame, fam: str, rho) -> pd.DataFrame:
    mask = sweep["family"].eq(fam)
    if pd.isna(rho):
        mask &= sweep["rho"].isna()
    else:
        mask &= sweep["rho"].eq(float(rho))
    return sweep[mask]


def _draw_calibration_cell(ax, sub: pd.DataFrame, title: str) -> None:
    for i, (test, lab) in enumerate(_TESTS):
        g = sub[sub["test"] == test].sort_values("J")
        ax.plot(g["J"], g["type_I"], marker="o", ms=4.0, color=PALETTE[i],
                label=lab)
        ax.fill_between(g["J"], g["type_I"] - 2 * g["mcse"],
                        g["type_I"] + 2 * g["mcse"], color=PALETTE[i],
                        alpha=0.18, linewidth=0)
    ax.axhline(ALPHA, color="crimson", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 0.125)
    ax.set_title(title, fontsize=11)


def fig_type_I_vs_skewness() -> None:
    """Fig 8 (main text) and Fig 8a-8c (appendix).

    The 15-cell grid is unreadable at single-column width, so the main text
    carries the two panels that state the finding -- the replication count each
    test needs to become calibrated, against the skewness of the paired
    contrast, and the one family that never calibrates -- and the full grid is
    split across three appendix figures of five cells each.
    """
    sweep = pd.read_csv(os.path.join(RESULTS, "E1b", "type_I_vs_J.csv"))
    rule = pd.read_csv(os.path.join(RESULTS, "E1b", "operational_rule.csv"))
    cells_sorted = _e1b_cells(sweep)

    # ---- Fig 8: the two-panel summary -------------------------------------
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.4, 4.6))

    # One marker per configuration; consecutive markers are different loss
    # families, so they are deliberately not joined by a line.
    for i, (test, lab) in enumerate(_TESTS):
        g = rule[rule["test"] == test].sort_values("abs_skew")
        ax.scatter(g["abs_skew"], g["J_min"], s=70, color=PALETTE[i],
                   marker="o" if i == 0 else "^", label=lab, zorder=3,
                   linewidths=0)
    ax.axhline(30, color="0.45", ls=":", lw=1.6, label="conventional $J>30$")
    ax.set_yscale("log")
    ax.set_xlabel("$|\\gamma_D|$, absolute skewness of $D_j$")
    ax.set_ylabel("smallest calibrated $J$")
    ax.set_title("What the requirement depends on")
    ax.legend(frameon=False, fontsize=11, loc="upper left")

    mix = _e1b_cell_rows(sweep, "mix", 0.6)
    _draw_calibration_cell(bx, mix, "Catastrophic-failure mixture,\n"
                                    "$\\rho=0.6$: neither test calibrates")
    bx.set_xlabel("replications $J$ (log scale)")
    bx.set_ylabel("unconditional Type I error")
    bx.legend(frameon=False, fontsize=11, loc="upper left")

    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig8.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)

    # ---- Fig 8a-8c: the full grid, five cells per figure -------------------
    groups = [cells_sorted.iloc[0:5], cells_sorted.iloc[5:10],
              cells_sorted.iloc[10:15]]
    for suffix, group in zip("abc", groups):
        fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.8), sharex=True,
                                 sharey=True)
        axes = np.atleast_1d(axes).ravel()
        for ax, item in zip(axes, group.itertuples(index=False)):
            sub = _e1b_cell_rows(sweep, item.family, item.rho)
            rho_lab = "NA" if pd.isna(item.rho) else f"{item.rho:g}"
            _draw_calibration_cell(
                ax, sub,
                f"{item.family}, $\\rho$={rho_lab}, "
                f"$|\\gamma_D|={item.abs_skew:.2f}$")
        legend_ax = axes[len(group)]
        legend_ax.axis("off")
        legend_ax.grid(False)
        handles, labels = axes[0].get_legend_handles_labels()
        handles.append(plt.Line2D([0], [0], color="crimson", ls="--", lw=1.0))
        labels.append("nominal 0.05")
        legend_ax.legend(handles, labels, frameon=False, fontsize=12,
                         loc="center")
        for ax in axes[len(group) + 1:]:
            ax.set_visible(False)
        fig.supxlabel("replications $J$ (log scale)")
        fig.supylabel("unconditional Type I error")
        fig.tight_layout()
        out = os.path.join(FIGURES, f"Fig8{suffix}.png")
        fig.savefig(out); plt.close()
        print("wrote", out)
        download(out)


# ---------------------------------------------------------------------------
# Fig E6: the same calibration protocol on the 24 real case-study contrasts
# ---------------------------------------------------------------------------
def fig_real_loss_calibration() -> None:
    """Calibration requirement measured on the losses the case study recorded.

    Left, the smallest calibrated J against the skewness of the real paired
    contrast, one marker per declared contrast; right, the one real contrast on
    which the studentized bootstrap never enters the band while the paired t
    does. Both panels use the E1b protocol unchanged, so the only difference
    from Fig 8 is that the loss law is the recorded one rather than a
    parametric family.
    """
    rule = pd.read_csv(os.path.join(RESULTS, "E6", "real_loss_calibration.csv"))
    sweep = pd.read_csv(os.path.join(RESULTS, "E6", "real_loss_type_I_vs_J.csv"))
    j_max = float(sweep["J"].max())

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.4, 4.6))

    for i, (test, lab) in enumerate(_TESTS):
        g = rule[rule["test"] == test].sort_values("abs_skew")
        # A cell that never calibrates has no finite requirement; it is drawn
        # as an open marker above the grid rather than silently dropped.
        finite = g[g["J_min"].notna()]
        never = g[g["J_min"].isna()]
        ax.scatter(finite["abs_skew"], finite["J_min"], s=70, color=PALETTE[i],
                   marker="o" if i == 0 else "^", label=lab, zorder=3,
                   linewidths=0)
        if len(never):
            ax.scatter(never["abs_skew"], np.full(len(never), j_max * 1.6),
                       s=90, facecolors="none", edgecolors=PALETTE[i],
                       marker="o" if i == 0 else "^", zorder=3, linewidths=1.6,
                       label=f"{lab}, never calibrates")
    ax.axhline(30, color="0.45", ls=":", lw=1.6, label="conventional $J>30$")
    ax.set_yscale("log")
    ax.set_xlabel("$|\\gamma_D|$, absolute skewness of the real $D_j$")
    ax.set_ylabel("smallest calibrated $J$")
    ax.set_title("24 declared case-study contrasts")
    ax.legend(frameon=False, fontsize=10, loc="lower right")

    # The right panel isolates the bootstrap's own excess from skewness, so it
    # shows the LEAST skewed contrast the bootstrap still fails to calibrate on:
    # there the paired t is essentially exact at every J, and whatever separates
    # the two curves cannot be attributed to the shape of D.
    worst = rule[(rule["test"] == "studentized_bootstrap")
                 & rule["J_min"].isna()].sort_values("abs_skew")
    if len(worst):
        key = worst.iloc[0]
    else:                       # every cell calibrated: show the hardest one
        key = rule[rule["test"] == "studentized_bootstrap"] \
            .sort_values("J_min").iloc[-1]
    sub = sweep[(sweep["dgp"] == key["dgp"]) & (sweep["n"] == key["n"])
                & (sweep["contrast"] == key["contrast"])]
    _draw_calibration_cell(
        bx, sub, f"DGP{int(key['dgp'])}, $n$={int(key['n'])}, "
                 f"{key['label']}\n$|\\gamma_D|={key['abs_skew']:.2f}$")
    bx.set_xlabel("replications $J$ (log scale)")
    bx.set_ylabel("unconditional Type I error")
    bx.legend(frameon=False, fontsize=11, loc="upper left")

    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig_E6_real_loss_calibration.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


# ---------------------------------------------------------------------------
# Fig 9: MCS composition / elimination path for the case study (E3)
# ---------------------------------------------------------------------------
def fig_mcs_paths() -> None:
    rows = []
    for dgp, n in CELLS:
        frame = pd.read_csv(os.path.join(RESULTS, "E3", "v2",
                                         f"DGP{dgp}_n{n}_replications.csv"))
        frame = frame[frame["seed"] >= 100]  # confirmatory block (AMENDMENT 1)
        for outcome in (1, 2):
            loss = np.column_stack([frame[f"{m}_pehe{outcome}"].to_numpy(float)
                                    for m in MODELS])
            res = _mcs.mcs(loss, alpha=MCS_ALPHA, B=RW_B, seed=RW_SEED,
                           model_names=MODELS)
            included = set(res["included"])
            order = {m: i + 1 for i, m in enumerate(res["elimination_order"])}
            for i, m in enumerate(MODELS):
                rows.append(dict(
                    dgp=dgp, n=n, outcome=outcome, model=m,
                    mean=float(loss[:, i].mean()),
                    elim_step=order.get(m, 0), included=int(m in included)))
    md = pd.DataFrame(rows)

    grid = [(dgp, n, o) for dgp, n in CELLS for o in (1, 2)]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for r_i, (dgp, n, out) in enumerate(grid):
        g = md[(md["dgp"] == dgp) & (md["n"] == n) & (md["outcome"] == out)]
        for m_i, m in enumerate(MODELS):
            r = g[g["model"] == m].iloc[0]
            if r["included"]:
                ax.add_patch(plt.Rectangle((m_i - 0.4, r_i - 0.4), 0.8, 0.8,
                                           fc="#55A868", ec="0.4"))
                ax.text(m_i, r_i, f"in ({r['mean']:.2f})", ha="center",
                        va="center", fontsize=8)
            else:
                ax.add_patch(plt.Rectangle((m_i - 0.4, r_i - 0.4), 0.8, 0.8,
                                           fc="#FFFFFF", ec="0.5"))
                ax.text(m_i, r_i, f"step {r['elim_step']} ({r['mean']:.2f})",
                        ha="center", va="center", fontsize=8)
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels(MODELS)
    ax.set_yticks(range(len(grid)))
    ax.set_yticklabels([f"DGP{dgp} n={n}, Y{o}" for dgp, n, o in grid])
    ax.set_xlim(-0.6, len(MODELS) - 0.4)
    ax.set_ylim(-0.5, len(grid) - 0.5)
    ax.invert_yaxis()
    ax.set_title("MCS composition on PEHE (95% confidence)\n"
                 "in = retained, step = bootstrap elimination order "
                 "(mean PEHE in brackets)")
    # Keep the tick labels -- they say which cell and which model each box is --
    # but drop the frame and grid, which carry no information here.
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig9.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


def fig_interval_score() -> None:
    """Fig 10: the exact 95% CATE interval score (P5.5-T1).

    Left, the MCS composition on the score itself. Right, the same scores split
    into their two components, ``score = width + penalty``, as shares. The split
    is the point of the figure: a model can lose on the score by being wide, by
    missing the truth, or by both, and only the split says which. In DGP2 the
    penalty share reaches 0.94 for every method, which is what a confounded DGP
    does to interval statements that no coverage-only summary makes comparable.
    """
    md = pd.read_csv(os.path.join(RESULTS, "E3", "mcs_interval_score.csv"))
    md = md[md["level"] == 95]
    dc = pd.read_csv(os.path.join(RESULTS, "E3", "interval_score_decomposition.csv"))
    dc = dc[dc["level"] == 95]

    grid = [(dgp, n, o) for dgp, n in CELLS for o in (1, 2)]
    labels = [f"DGP{dgp} $n$={n}, Y{o}" for dgp, n, o in grid]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.6, 5.8),
                                 gridspec_kw={"width_ratios": [1.0, 1.25]})

    for r_i, (dgp, n, out) in enumerate(grid):
        g = md[(md["dgp"] == dgp) & (md["n"] == n) & (md["outcome"] == out)]
        for m_i, m in enumerate(MODELS):
            r = g[g["model"] == m].iloc[0]
            kept = bool(r["in_mcs_95"])
            ax.add_patch(plt.Rectangle((m_i - 0.45, r_i - 0.45), 0.9, 0.9,
                                       fc="#55A868" if kept else "#FFFFFF",
                                       ec="0.4"))
            ax.text(m_i, r_i, f"{'in' if kept else 'out'}\n{r['mean_interval_score']:.1f}",
                    ha="center", va="center", fontsize=8.5)
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels(MODELS)
    ax.set_yticks(range(len(grid)))
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.6, len(MODELS) - 0.4)
    ax.set_ylim(-0.6, len(grid) - 0.4)
    ax.invert_yaxis()
    ax.set_title("MCS on the exact 95% interval score\n(95% confidence; mean score shown)")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    height = 0.19
    for m_i, m in enumerate(MODELS):
        offs = (m_i - 1.5) * height
        widths, penalties, ys = [], [], []
        for r_i, (dgp, n, out) in enumerate(grid):
            r = dc[(dc["dgp"] == dgp) & (dc["n"] == n) & (dc["outcome"] == out)
                   & (dc["model"] == m)].iloc[0]
            total = r["mean_width"] + r["mean_penalty"]
            widths.append(r["mean_width"] / total)
            penalties.append(r["mean_penalty"] / total)
            ys.append(r_i + offs)
        bx.barh(ys, widths, height=height, color=PALETTE[m_i], alpha=0.95,
                edgecolor="none", label=m if m_i == 0 else m)
        bx.barh(ys, penalties, height=height, left=widths, color=PALETTE[m_i],
                alpha=0.28, edgecolor="none")
    bx.set_yticks(range(len(grid)))
    bx.set_yticklabels([])
    bx.set_ylim(-0.6, len(grid) - 0.4)
    bx.invert_yaxis()
    bx.set_xlim(0, 1)
    bx.set_xlabel("share of the mean interval score")
    bx.set_title("Interval width (solid) and non-coverage penalty (pale)\n"
                 "bars within a row are mvbcf, bcf, bart, mvbart top to bottom")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETTE[i], alpha=0.95)
               for i in range(len(MODELS))]
    bx.legend(handles, MODELS, loc="lower center", ncol=4, frameon=False,
              bbox_to_anchor=(0.5, -0.24))
    fig.tight_layout()
    out = os.path.join(FIGURES, "Fig10.png")
    fig.savefig(out); plt.close()
    print("wrote", out)
    download(out)


def main(argv=None) -> int:
    global FIGURES
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures-dir", default=FIGURES)
    ap.add_argument("--paper-dir", default=None,
                    help="also copy the manuscript-named PNGs here")
    ap.add_argument("--skip", nargs="*", default=[],
                    help="figure function names to skip")
    args = ap.parse_args(argv)

    FIGURES = os.path.abspath(args.figures_dir)
    os.makedirs(FIGURES, exist_ok=True)

    for name in ["fig_bibliometrics", "fig_precision_gain", "fig_case_study_planning",
                 "fig_contrast_forest", "fig_operating_characteristics",
                 "fig_J_distribution", "fig_type_I_vs_skewness",
                 "fig_real_loss_calibration", "fig_mcs_paths",
                 "fig_interval_score"]:
        if name in args.skip:
            continue
        fn = globals()[name]
        print("--", name)
        fn()

    if args.paper_dir:
        os.makedirs(args.paper_dir, exist_ok=True)
        for name in PAPER_FIGURES:
            src = os.path.join(FIGURES, name)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(args.paper_dir, name))
                print("copied", name, "->", args.paper_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
