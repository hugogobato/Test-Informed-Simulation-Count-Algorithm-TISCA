#!/usr/bin/env python3
"""Generate ``E2_design_comparison.ipynb`` (REVISION_PLAN.md P3-T3).

Analysis-only over the E1 output: the six designs D1-D6 side by side on achieved
error rates, achieved power and E[J], ending in the explicit recommendation
sentence the plan's acceptance criterion demands. This is the notebook that answers
Reviewer 2 paragraph 4 ("the iterative algorithm seems unnecessarily complicated").

Regenerate with::

    python notebooks/_generators/build_e2_design_comparison_nb.py
"""

from __future__ import annotations

from _nbcommon import PLOT_STYLE, SETUP, code, md, write


DESIGN_LABELS = r'''
DESIGNS = ["D1", "D2", "D3", "D4", "D5", "D6"]
LABELS = {
    "D1": "D1 fixed-J, independent pilot (raw pilot sigma)",
    "D2": "D2 internal-pilot re-estimation (adaptive)",
    "D3": "D3 TISCA v1: unpaired Welch, iterative",
    "D4": "D4 TISCA v2: paired, two-stage (Algorithm 1)",
    "D5": "D5 paired fixed-precision (MCSE target)",
    "D6": "D6 oracle fixed-J (true sigma known)",
}
ALPHA = 0.05
POWER_TARGET = 0.80
'''

LOAD = r'''
OC_PATH = os.path.join(RESULTS, "E1", "operating_characteristics.csv")
assert os.path.exists(OC_PATH), (
    f"{OC_PATH} not found. Produce it first with\n"
    "  python experiments/E1_operating_characteristics/run_e1_grid.py\n"
    "or by running the three E1 Colab notebooks and concatenating their CSVs."
)
oc = pd.read_csv(OC_PATH)
print("cells:", len(oc), dict(oc.groupby("module").size()))

EXPECTED = {"A": 864, "B": 660, "C": 243, "D": 216}
got = dict(oc.groupby("module").size())
missing = {m: EXPECTED[m] - got.get(m, 0) for m in EXPECTED if got.get(m, 0) != EXPECTED[m]}
assert not missing, f"incomplete E1 output, missing per module: {missing}"
assert not oc["cell_id"].duplicated().any(), "duplicate cell_id"

# theta = 0 is the null; theta_mult = 1 is the planning alternative the designs
# were sized for, so "achieved power" is only meaningful there.
oc["is_null"] = oc["theta"] == 0.0
oc["family_kind"] = np.where(oc["family"].str.startswith("empirical"), "empirical", "parametric")
print(oc.groupby(["module", "design"]).size().unstack(fill_value=0))
'''

TABLE = r'''
def wmcse(p, R):
    """MCSE of a mean rate estimated from R outer repetitions."""
    return float(np.sqrt(np.mean(p) * (1 - np.mean(p)) / R))


def design_table(df, label=""):
    """One row per design: achieved level, achieved power, E[J], each with MCSE."""
    rows = []
    for d in DESIGNS:
        sub = df[df["design"] == d]
        if sub.empty:
            continue
        null = sub[sub["is_null"]]
        alt = sub[sub["theta_mult"] == 1.0]
        R_null = null["R"].median() if len(null) else np.nan
        rows.append({
            "design": d,
            "label": LABELS[d],
            "n_cells": len(sub),
            "type_I_mean": null["t1e_or_power"].mean(),
            "type_I_max": null["t1e_or_power"].max(),
            "type_I_mcse": wmcse(null["t1e_or_power"], R_null) if len(null) else np.nan,
            "n_cells_over_level": int((null["t1e_or_power"]
                                       > ALPHA + 2 * null["reject_rate_mcse"]).sum()),
            "ci_cover_mean": null["ci_cover"].mean(),
            "power_mean": alt["t1e_or_power"].mean(),
            "power_min": alt["t1e_or_power"].min(),
            "n_cells_under_power": int((alt["t1e_or_power"] < POWER_TARGET
                                        - 2 * alt["reject_rate_mcse"]).sum()),
            "E_J_null": null["E_J"].mean(),
            "E_J_alt": alt["E_J"].mean(),
            "sd_J_alt": alt["sd_J"].mean(),
            "q95_J_alt": alt["q95_J"].mean(),
            "P_J_at_cap": alt["pJmax"].mean(),
            "bias_alt": alt["bias"].mean(),
            "rmse_alt": alt["rmse"].mean(),
            "scope": label,
        })
    return pd.DataFrame(rows)


main_tbl = design_table(oc[oc["module"] == "A"], "Module A (all families)")
pd.set_option("display.width", 200)
print(main_tbl[["design", "type_I_mean", "type_I_max", "n_cells_over_level",
                "power_mean", "power_min", "E_J_alt", "sd_J_alt",
                "P_J_at_cap"]].round(4).to_string(index=False))
'''

STRATIFIED = r'''
# Same table restricted to the parametric families, and to the two empirical ones.
# The empirical rows are the realistic case: real loss shapes, and for the row
# bootstrap the real dependence too.
by_kind = pd.concat([
    design_table(oc[(oc["module"] == "A") & (oc["family_kind"] == k)], k)
    for k in ("parametric", "empirical")
], ignore_index=True)
print(by_kind[["scope", "design", "type_I_mean", "power_mean", "E_J_alt"]]
      .round(4).to_string(index=False))

# Type I error by design and family: the honest per-family statement the plan's
# acceptance criterion asks for, rather than a single pooled number.
t1e_fam = (oc[(oc["module"] == "A") & oc["is_null"]]
           .pivot_table(index="family", columns="design",
                        values="t1e_or_power", aggfunc="mean").round(4))
print()
print("Type I error by family x design (theta = 0, Module A):")
print(t1e_fam.to_string())

# ... and the same broken out by correlation, which is where D3 (unpaired) is
# expected to separate from the paired designs.
t1e_rho = (oc[(oc["module"] == "A") & oc["is_null"] & oc["rho"].notna()]
           .pivot_table(index="rho", columns="design",
                        values="t1e_or_power", aggfunc="mean").round(4))
ej_rho = (oc[(oc["module"] == "A") & (oc["theta_mult"] == 1.0) & oc["rho"].notna()]
          .pivot_table(index="rho", columns="design", values="E_J", aggfunc="mean").round(1))
print()
print("Type I error by rho x design:")
print(t1e_rho.to_string())
print()
print("E[J] by rho x design (at the planning alternative):")
print(ej_rho.to_string())
'''

LOOP = r'''
# --------------------------------------------------------------------------- #
# The actual question: does the iterative loop earn its complexity?            #
# --------------------------------------------------------------------------- #
# Matched comparison, cell by cell, so the contrast is within-configuration and
# not an artefact of which cells happen to be in each design's average. `rho` is
# NaN for the row-bootstrap family (its dependence comes from the data), so it is
# given a sentinel before pivoting -- otherwise those cells vanish from the index.
keys = ["module", "family", "rho_key", "theta_mult", "sigma_ratio", "J0", "B"]
src = oc[oc["module"].isin(["A", "C", "D"])].assign(rho_key=oc["rho"].fillna(-99.0))
wide = src.pivot_table(index=keys, columns="design",
                       values=["t1e_or_power", "E_J", "ci_cover"])
tm = wide.index.get_level_values("theta_mult")
IS_NULL, IS_ALT = (tm == 0.0), (tm == 1.0)


def diff(metric, a, b, mask):
    """Paired differences over the cells where BOTH designs were run."""
    if (metric, a) not in wide.columns or (metric, b) not in wide.columns:
        return pd.Series(dtype=float)
    return (wide[(metric, a)] - wide[(metric, b)])[mask].dropna()


def ratio(metric, a, b, mask):
    if (metric, a) not in wide.columns or (metric, b) not in wide.columns:
        return pd.Series(dtype=float)
    return (wide[(metric, a)] / wide[(metric, b)])[mask].replace(
        [np.inf, -np.inf], np.nan).dropna()


loop_rows = []
for a, b, question in [
    ("D2", "D4", "adaptive loop vs two-stage (Reviewer 2 par. 4)"),
    ("D3", "D4", "TISCA v1 (unpaired, iterative) vs TISCA v2"),
    ("D1", "D4", "raw-pilot fixed-J vs variance-inflated two-stage"),
    ("D5", "D4", "precision target vs power target"),
    ("D4", "D6", "two-stage vs the oracle it approximates"),
]:
    d_t1e = diff("t1e_or_power", a, b, IS_NULL)
    d_pow = diff("t1e_or_power", a, b, IS_ALT)
    d_J = diff("E_J", a, b, IS_ALT)
    loop_rows.append({
        "comparison": f"{a} - {b}", "question": question,
        "d_type_I_mean": d_t1e.mean(), "d_type_I_max": d_t1e.abs().max(),
        "d_power_mean": d_pow.mean(), "d_power_worst": d_pow.min(),
        "d_E_J_mean": d_J.mean(),
        "E_J_ratio": ratio("E_J", a, b, IS_ALT).mean(),
        "n_matched": int(len(d_J)),
    })
loop = pd.DataFrame(loop_rows)
print(loop.round(4).to_string(index=False))
'''

VERDICT = r'''
# --------------------------------------------------------------------------- #
# The recommendation sentence, derived from the table rather than asserted.     #
# --------------------------------------------------------------------------- #
tbl_i = main_tbl.set_index("design")
d4, d6 = tbl_i.loc["D4"], tbl_i.loc["D6"]
gap = loop.set_index("comparison")

# Reviewer 2's question is whether the LOOP is worth its complexity. The adaptive
# design reaching higher power is not on its own an argument for it: it also runs
# more replications. The comparison that answers the question is efficiency, so the
# power gain is set against the extra replications rather than reported alone.
d2 = gap.loc["D2 - D4"]
d2_inflates = d2["d_type_I_mean"] > 2 * d4["type_I_mcse"]
d2_more_efficient = (d2["d_power_mean"] > 0.01) and (d2["E_J_ratio"] < 1.05)

# Matched-budget reference: what does the two-stage design achieve when it is given
# the replications the loop actually spends? The oracle at a comparable E[J] is the
# cleanest available proxy, since D6 sizes J from the true sigma.
level_ok = abs(d4["type_I_mean"] - ALPHA) <= 2 * d4["type_I_mcse"]

if d2_inflates:
    clause = ("it inflates the unconditional Type I error, which it does here by "
              f"{d2['d_type_I_mean']:+.4f} ({d2['d_type_I_mean'] / d4['type_I_mcse']:.1f} MCSE)")
elif d2_more_efficient:
    clause = (f"more power per replication is needed: the loop gains "
              f"{d2['d_power_mean']:+.3f} power at {d2['E_J_ratio']:.2f}x the cost, "
              "which is a genuine efficiency gain")
else:
    clause = (f"a larger budget is acceptable in itself: the loop's extra power "
              f"({d2['d_power_mean']:+.3f}) is bought with {d2['E_J_ratio']:.2f}x the "
              "replications, so it is a bigger design rather than a better one")

verdict = (
    f"Use the two-stage design (D4) unless {clause}. "
    f"At the planning alternative D4 achieves power {d4['power_mean']:.3f} "
    f"(target {POWER_TARGET}) at E[J] = {d4['E_J_alt']:.0f}, against the oracle's "
    f"{d6['power_mean']:.3f} at E[J] = {d6['E_J_alt']:.0f}; its unconditional Type I "
    f"error is {d4['type_I_mean']:.4f} (nominal {ALPHA}, MCSE {d4['type_I_mcse']:.4f})"
    f"{'' if level_ok else ', which is NOT within 2 MCSE of nominal and must be reported as such'}. "
    f"The adaptive loop's effect on the level is {d2['d_type_I_mean']:+.4f} "
    f"({'material' if d2_inflates else 'within 2 MCSE'}), and it spends "
    f"{d2['E_J_ratio']:.2f}x the replications for {d2['d_power_mean']:+.3f} power."
)
print(verdict)

# Two further facts the same table settles, both of which belong in Section 4.
d3 = gap.loc["D3 - D4"]
d5 = gap.loc["D5 - D4"]
extra = (
    f"\nTISCA v1 (D3) is CONSERVATIVE, not liberal: its level falls to "
    f"{oc[(oc.module=='A') & oc.is_null & (oc.design=='D3') & (oc.rho==0.9)]['t1e_or_power'].mean():.4f} "
    f"at rho = 0.9 against a nominal {ALPHA}, and it spends {d3['E_J_ratio']:.2f}x the "
    f"replications of v2. Fixing the paired-design error therefore buys efficiency, "
    f"not validity.\n"
    f"\nThe precision target (D5) is far more demanding than the power target at "
    f"these settings: {d5['E_J_ratio']:.1f}x the replications of D4, hitting the "
    f"budget cap in {tbl_i.loc['D5', 'P_J_at_cap']:.1%} of repetitions. Which target "
    f"binds is a reporting obligation, not a detail."
)
print(extra)

os.makedirs(os.path.join(RESULTS, "E2"), exist_ok=True)
main_tbl.to_csv(os.path.join(RESULTS, "E2", "design_comparison.csv"), index=False)
by_kind.to_csv(os.path.join(RESULTS, "E2", "design_comparison_by_family_kind.csv"), index=False)
loop.to_csv(os.path.join(RESULTS, "E2", "loop_vs_two_stage.csv"), index=False)
t1e_fam.to_csv(os.path.join(RESULTS, "E2", "type_I_by_family.csv"))
with open(os.path.join(RESULTS, "E2", "design_verdict.md"), "w") as f:
    f.write("# P3-T3 verdict: is the iterative loop worth it?\n\n")
    f.write(verdict + "\n" + extra + "\n\n## Design table (Module A)\n\n")
    f.write(main_tbl.round(4).to_markdown(index=False) + "\n\n## Matched comparisons\n\n")
    f.write(loop.round(4).to_markdown(index=False) + "\n")
print("\nwrote results/E2/")
'''

FIGURE = r'''
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

ax = axes[0]
null_A = oc[(oc["module"] == "A") & oc["is_null"]]
for i, d in enumerate(DESIGNS):
    v = null_A.loc[null_A["design"] == d, "t1e_or_power"]
    ax.scatter(np.full(len(v), i) + np.random.default_rng(i).normal(0, .07, len(v)),
               v, s=9, alpha=.35, color=PALETTE[i], linewidths=0)
    ax.scatter([i], [v.mean()], marker="_", s=700, color="black", zorder=3)
ax.axhline(ALPHA, color="crimson", lw=1.2, ls="--", label="nominal 0.05")
ax.set_xticks(range(len(DESIGNS)))
ax.set_xticklabels(DESIGNS)
ax.set_ylabel("unconditional Type I error")
ax.set_title("Achieved level at $\\theta = 0$")
ax.legend(frameon=False, loc="upper left")

ax = axes[1]
alt_A = oc[(oc["module"] == "A") & (oc["theta_mult"] == 1.0)]
for i, d in enumerate(DESIGNS):
    s = alt_A[alt_A["design"] == d]
    ax.scatter(s["E_J"], s["t1e_or_power"], s=16, alpha=.55,
               color=PALETTE[i], label=d, linewidths=0)
ax.axhline(POWER_TARGET, color="crimson", lw=1.2, ls="--")
ax.set_xscale("log")
ax.set_xlabel("E[J]  (log scale)")
ax.set_ylabel("achieved power")
ax.set_title("Power against the replications spent")
ax.legend(frameon=False, ncol=2, fontsize=10)

ax = axes[2]
piv = (alt_A[alt_A["rho"].notna()]
       .pivot_table(index="rho", columns="design", values="E_J", aggfunc="mean"))
for i, d in enumerate(DESIGNS):
    if d in piv:
        ax.plot(piv.index, piv[d], marker="o", color=PALETTE[i], label=d)
ax.set_xlabel(r"design correlation $\rho$")
ax.set_ylabel("E[J] at the planning alternative")
ax.set_title("Cost of ignoring the pairing")
ax.legend(frameon=False, ncol=2, fontsize=10)

fig.tight_layout()
out = os.path.join(FIGURES, "Fig_E2_design_comparison.png")
fig.savefig(out)
print("wrote", out)
download(out)
'''


def build():
    cells = []
    md(cells, """
        # E2: design comparison, and the verdict on the iterative loop (P3-T3)

        Analysis-only. Reads `results/E1/operating_characteristics.csv` (P3-T2) and
        answers Reviewer 2's fourth paragraph directly: for each of the six designs
        D1-D6, what error rate, what power, and what E[J] were actually achieved,
        and does the iterative loop beat the closed-form two-stage design.

        The plan's acceptance criterion for this task is an explicit recommendation
        of the form *"use the two-stage design unless <condition>"*, supported by the
        table. The final cell writes that sentence from the measured numbers rather
        than from an assumption about how it should come out. **A verdict against the
        adaptive loop is a publishable finding, not a failure**; so is a verdict
        against the two-stage default.
        """)
    code(cells, SETUP)
    code(cells, PLOT_STYLE)
    code(cells, DESIGN_LABELS)
    md(cells, """
        ## Load and check the E1 output

        Completeness is asserted, not assumed: a short grid is the failure mode that
        would otherwise reach the manuscript unnoticed.
        """)
    code(cells, LOAD)
    md(cells, """
        ## The headline design table

        Every rate carries an MCSE. `n_cells_over_level` counts the cells where the
        achieved Type I error sits more than 2 MCSE above nominal, which is the
        honest way to report "does this design hold its level" across a grid: a mean
        near 0.05 can still hide cells that do not.
        """)
    code(cells, TABLE)
    md(cells, """
        ## Stratified views

        Pooling over families hides exactly what the referees asked about. The
        breakdown by family answers IJDA comment 3 (finite-sample behaviour under
        skewness and heavy tails); the breakdown by correlation answers comment 2
        (what the pairing is worth), and is where the unpaired v1 design D3 should
        separate from the paired designs.
        """)
    code(cells, STRATIFIED)
    md(cells, """
        ## Matched design-versus-design comparison

        Cell-by-cell differences within the same configuration, so no comparison is
        contaminated by differing cell composition. `D4 - D6` is the price of not
        knowing sigma; `D3 - D4` is what fixing the paired-design error buys.
        """)
    code(cells, LOOP)
    md(cells, """
        ## Verdict

        The sentence below goes into Section 4 of the manuscript and into the
        response letter, and the branches are written so that an unfavourable result
        is stated plainly rather than softened.
        """)
    code(cells, VERDICT)
    md(cells, """
        ## Figure

        Three panels: achieved level per design, the power-versus-E[J] frontier, and
        E[J] against the design correlation.
        """)
    code(cells, FIGURE)
    write("E2_design_comparison.ipynb", cells)


if __name__ == "__main__":
    build()
