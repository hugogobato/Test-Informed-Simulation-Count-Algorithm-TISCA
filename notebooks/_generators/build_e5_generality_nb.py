#!/usr/bin/env python3
"""Generate ``E5_generality_demo.ipynb`` (REVISION_PLAN.md P3-T6).

Answers SoutoNeto section 3: "Why is TISCA given in a context of estimating the
treatment effect? The approach seems flexible." The demonstration is one-step-ahead
forecasting of a stationary AR(1) process, chosen because **the ranking of the
competing forecasters is known in closed form**, so every number the v2 workflow
produces can be checked against the truth rather than merely reported.

With y_t = phi y_{t-1} + eps_t, eps ~ N(0, sigma^2), |phi| < 1, the population
one-step-ahead mean squared errors are

    AR(1) at the true parameter   sigma^2
    naive (random walk)           sigma^2 [1 + (1 - phi) / (1 + phi)]
    sample-mean forecast          sigma^2 / (1 - phi^2)

and an AR(2) fitted to AR(1) data is asymptotically tied with AR(1) but pays one
extra estimation penalty of sigma^2 / n_train at any finite training length.

That last contrast is the interesting one, and it turned out to be a better example
than the tie it was meant to be. Its true value is about -0.0025, roughly 5% of the
declared practical margin -- and with the J the design selects, it is detected at
p < 1e-30. The demonstration therefore reproduces IJDA comment 1 exactly: with
enough replications an arbitrarily small and practically irrelevant difference
becomes statistically significant, which is why the manuscript reports the estimate
and its interval first and tests against a declared margin (mode M3) before calling
anything an improvement.

Regenerate with::

    python notebooks/_generators/build_e5_generality_nb.py
"""

from __future__ import annotations

from _nbcommon import PLOT_STYLE, SETUP, code, md, write


DGP = r'''
from tisca import estimands, inference, mcs, multiplicity, planning

PHI, SIGMA = 0.7, 1.0
N_TRAIN, N_TEST = 400, 200
METHODS = ["AR(1)", "AR(2)", "naive", "mean"]
CONTROL = 0                      # AR(1) is the proposed method
ALPHA = 0.05

# Two reference values per forecaster, and the distinction between them matters.
#
# ASYMPTOTIC: the one-step-ahead MSE at the true parameters. Here AR(1) and AR(2)
# are exactly tied, because the second lag has a true coefficient of zero.
ASYMPTOTIC_MSE = {
    "AR(1)": SIGMA ** 2,
    "AR(2)": SIGMA ** 2,
    "naive": SIGMA ** 2 * (1 + (1 - PHI) / (1 + PHI)),
    "mean":  SIGMA ** 2 / (1 - PHI ** 2),
}

# FINITE-SAMPLE: what a forecaster fitted on N_TRAIN observations actually achieves.
# A linear predictor with p estimated parameters pays an estimation penalty of
# approximately sigma^2 * p / n, so AR(2) is genuinely -- if barely -- worse than
# AR(1) at any finite training length: one extra estimated coefficient costs
# sigma^2 / N_TRAIN. The sample-mean forecast pays the variance of ybar, which for
# an AR(1) is Var(y) * (1 + phi) / ((1 - phi) n). These are O(1/n) approximations
# and are labelled as such; the notebook checks the estimates against them with a
# tolerance rather than claiming they are exact.
N_PARAMS = {"AR(1)": 2, "AR(2)": 3, "naive": 0, "mean": 1}
VAR_Y = SIGMA ** 2 / (1 - PHI ** 2)
FINITE_MSE = {
    "AR(1)": SIGMA ** 2 * (1 + N_PARAMS["AR(1)"] / N_TRAIN),
    "AR(2)": SIGMA ** 2 * (1 + N_PARAMS["AR(2)"] / N_TRAIN),
    "naive": ASYMPTOTIC_MSE["naive"],
    "mean":  VAR_Y + VAR_Y * (1 + PHI) / ((1 - PHI) * N_TRAIN),
}
TRUE_CONTRAST = {m: FINITE_MSE["AR(1)"] - FINITE_MSE[m] for m in METHODS[1:]}

# The practical margin: a difference in one-step-ahead MSE smaller than 5% of the
# innovation variance is not worth acting on. It is declared HERE, before any data
# is seen, and it is the same number used for the planning alternative and for the
# minimum-effect test below (IJDA comment 4: planning, testing and claiming must be
# aligned, and the MDE is a planning device, not a licence to claim importance).
MARGIN = 0.05 * SIGMA ** 2

print("asymptotic MSE  :", {k: round(v, 5) for k, v in ASYMPTOTIC_MSE.items()})
print("finite-sample   :", {k: round(v, 5) for k, v in FINITE_MSE.items()})
print("true contrasts  :", {k: round(v, 5) for k, v in TRUE_CONTRAST.items()})
print(f"practical margin: {MARGIN} (5% of sigma^2)")
print("\nNote the AR(1)-vs-AR(2) contrast: true value "
      f"{TRUE_CONTRAST['AR(2)']:+.5f}, which is real but is "
      f"{abs(TRUE_CONTRAST['AR(2)']) / MARGIN:.1%} of the practical margin.")
'''

SIM = r'''
def simulate(seed):
    """One replication: simulate a series, fit four forecasters, return four losses.

    The loss is the mean squared one-step-ahead forecast error over the test block,
    i.e. the replicate-level quantity L_j of the estimand table: one scalar per
    method per replication, lower is better, evaluated on the SAME simulated series
    for every method. That common series is what makes the contrasts paired.
    """
    rng = np.random.default_rng(seed)
    n = N_TRAIN + N_TEST + 50                       # burn-in for stationarity
    e = rng.normal(0.0, SIGMA, size=n)
    y = np.empty(n)
    y[0] = e[0] / np.sqrt(1 - PHI ** 2)             # stationary start
    for t in range(1, n):
        y[t] = PHI * y[t - 1] + e[t]
    y = y[50:]
    train, test = y[:N_TRAIN], y[N_TRAIN - 2:]      # keep 2 lags for the AR(2) start

    def ols(lags):
        X = np.column_stack([train[lags - k - 1:len(train) - k - 1] for k in range(lags)])
        X = np.column_stack([np.ones(len(X)), X])
        b, *_ = np.linalg.lstsq(X, train[lags:], rcond=None)
        return b

    b1, b2 = ols(1), ols(2)
    ybar = train.mean()

    losses = []
    idx = np.arange(2, len(test))                   # forecast origins inside `test`
    actual = test[idx]
    losses.append(np.mean((actual - (b1[0] + b1[1] * test[idx - 1])) ** 2))
    losses.append(np.mean((actual - (b2[0] + b2[1] * test[idx - 1]
                                     + b2[2] * test[idx - 2])) ** 2))
    losses.append(np.mean((actual - test[idx - 1]) ** 2))
    losses.append(np.mean((actual - ybar) ** 2))
    return losses


probe = np.array([simulate(1_000_001 + j) for j in range(30)])
print("30-replication means:", dict(zip(METHODS, probe.mean(0).round(5))))
print("finite-sample ref   :", {k: round(v, 5) for k, v in FINITE_MSE.items()})
'''

PLAN = r'''
# --------------------------------------------------------------------------- #
# Stage 1: the independent pilot, on the reserved seed block.                   #
# --------------------------------------------------------------------------- #
J0 = 50
pilot = np.array([simulate(1_000_001 + j) for j in range(J0)])
D_pilot = pilot[:, [1, 2, 3]] - pilot[:, [CONTROL]]     # AR(2)/naive/mean minus AR(1)
# sign convention: D = L_proposed - L_benchmark, negative favours the proposal
D_pilot = -D_pilot

# The planning alternative is the practical margin declared in the first cell, and
# the precision target is one fifth of it, so that an interval is never wider than
# the difference it has to resolve.
DELTA = -MARGIN
TARGET_MCSE = MARGIN / 5

rows = []
for k, name in enumerate(METHODS[1:]):
    sd = float(D_pilot[:, k].std(ddof=1))
    sigma_ub = planning.inflate_std(sd, J0, gamma=0.20)
    alpha_adj, note = multiplicity.planning_alpha("bonferroni", K=3, alpha=ALPHA)
    J_prec = planning.required_J_mcse(sigma_ub, target_mcse=TARGET_MCSE)
    J_pow = planning.required_J_power(mode="M1", delta=DELTA, sigma=sigma_ub,
                                      alpha=alpha_adj, target_power=0.80)
    rows.append({"contrast": f"AR(1) vs {name}", "sd_pilot": sd, "sigma_UB": sigma_ub,
                 "J_precision": J_prec, "J_power": J_pow,
                 "J": max(J_prec, J_pow), "alpha_adj": alpha_adj,
                 "true_contrast": TRUE_CONTRAST[name]})
plan_tbl = pd.DataFrame(rows)
J_FINAL = int(plan_tbl["J"].max())
print(plan_tbl.round(5).to_string(index=False))
print(f"\nJ_final = max over contrasts and targets = {J_FINAL}   ({note})")
'''

CONFIRM = r'''
# --------------------------------------------------------------------------- #
# Stage 2: the confirmatory run on a disjoint seed block; the pilot is discarded. #
# --------------------------------------------------------------------------- #
L = np.array([simulate(1 + j) for j in range(J_FINAL)])
D = -(L[:, [1, 2, 3]] - L[:, [CONTROL]])
print("confirmatory replications:", L.shape)

alpha_adj = float(plan_tbl["alpha_adj"].iloc[0])
res = []
for k, name in enumerate(METHODS[1:]):
    t = inference.paired_t(D[:, k], alternative="two-sided")
    b = inference.studentized_paired_bootstrap(D[:, k], B=4999, seed=7 + k)
    # The interval is reported at the MULTIPLICITY-ADJUSTED level, so that it agrees
    # with the test actually being conducted (IJDA comment 4b: same sidedness, same
    # critical value in planning, testing and reporting).
    lo, hi = estimands.nc_ci(t["estimate"], t["se"], level=1 - alpha_adj, df=t["df"])
    res.append({
        "contrast": f"AR(1) vs {name}",
        "estimate": t["estimate"], "mcse": t["mcse"],
        "ci_lo": lo, "ci_hi": hi, "half_width": (hi - lo) / 2,
        "p_paired_t": t["p_value"], "p_bootstrap": b["p_value"],
        "true": TRUE_CONTRAST[name],
        "covers_truth": bool(lo <= TRUE_CONTRAST[name] <= hi),
        "mcse_target_met": bool(t["mcse"] <= TARGET_MCSE),
    })
conf = pd.DataFrame(res)
print(conf.round(5).to_string(index=False))

# Family-level inference: Romano-Wolf exploits the (strong) dependence between the
# three contrasts, all of which share the AR(1) column.
rw = multiplicity.romano_wolf_stepdown(D, B=4999, alpha=ALPHA, seed=11)
conf["p_romano_wolf"] = rw["p_values"]
conf["p_bonferroni"] = multiplicity.bonferroni(conf["p_paired_t"].to_numpy())
conf["p_holm"] = multiplicity.holm(conf["p_paired_t"].to_numpy())
print()
print(conf[["contrast", "p_paired_t", "p_bonferroni", "p_holm",
            "p_romano_wolf"]].round(6).to_string(index=False))
'''

MCS_CELL = r'''
# --------------------------------------------------------------------------- #
# The Model Confidence Set over all four forecasters.                           #
# --------------------------------------------------------------------------- #
# This is the layer that answers "which methods are indistinguishable from the
# best" without selecting a benchmark post hoc. The truth here is that AR(1) and
# AR(2) are tied and both beat naive and mean, so a correct MCS retains exactly
# {AR(1), AR(2)} -- a checkable claim, which is the whole point of choosing a
# problem with a known answer.
res_mcs = mcs.mcs(L, alpha=0.15, B=4999, statistic="Tmax", seed=21, model_names=METHODS)
# `table` is (m, 3) in elimination order; `table_names` labels its ROWS.
tbl = pd.DataFrame(np.asarray(res_mcs["table"]),
                   columns=["avg_loss", "p_elimination", "p_MCS"])
tbl.insert(0, "model", list(res_mcs["table_names"]))
print(tbl.round(5).to_string(index=False))
print()
print("MCS (85%) retains :", res_mcs["included"])
print("eliminated        :", res_mcs["excluded"])
# What SHOULD the MCS contain? Not the asymptotic tie: at this training length
# AR(1) is genuinely best and AR(2) is genuinely second, so with J large enough to
# resolve a gap of sigma^2/n the set can legitimately narrow to {AR(1)}. The
# checkable claim is the weaker and more honest one: the two forecasters that are
# worse by a practically meaningful margin must be eliminated.
print("naive and mean both eliminated:",
      not ({"naive", "mean"} & set(res_mcs["included"])))
print("elimination order (worst first):", list(res_mcs["elimination_order"]))
'''

VERIFY = r'''
# --------------------------------------------------------------------------- #
# Did the design deliver what it promised?                                      #
# --------------------------------------------------------------------------- #
# The reference values are O(1/n) approximations, so agreement is checked with a
# tolerance rather than by asserting interval coverage of an exact constant.
conf["abs_err_vs_reference"] = (conf["estimate"] - conf["true"]).abs()
order_est = list(np.argsort(L.mean(axis=0)))
order_true = list(np.argsort([FINITE_MSE[m] for m in METHODS]))

# The minimum-effect test (mode M3): is the contrast worse than the practical
# margin, not merely different from zero? This is the claim the manuscript should
# make, and it separates the three contrasts in exactly the way a plain two-sided
# test of equality does not.
m3 = []
for k, name in enumerate(METHODS[1:]):
    d = D[:, k]
    t_m3 = inference.paired_t(d + MARGIN, alternative="less")   # H0: E[D] >= -margin
    m3.append({"contrast": f"AR(1) vs {name}", "estimate": float(d.mean()),
               "margin": -MARGIN, "p_minimum_effect": t_m3["p_value"],
               "beats_margin": bool(t_m3["p_value"] < ALPHA / 3)})
m3 = pd.DataFrame(m3)
print("Minimum-effect tests at the declared practical margin:")
print(m3.round(6).to_string(index=False))

checks = pd.DataFrame([
    {"claim": "every contrast met its MCSE target",
     "value": f"max MCSE = {conf['mcse'].max():.5f} vs target {TARGET_MCSE:.5f}",
     "pass": bool(conf["mcse_target_met"].all())},
    {"claim": "estimates match the finite-sample reference within 1% of sigma^2",
     "value": f"max |error| = {conf['abs_err_vs_reference'].max():.5f}",
     "pass": bool(conf["abs_err_vs_reference"].max() < 0.01 * SIGMA ** 2)},
    {"claim": "the ranking of all four forecasters is recovered exactly",
     "value": f"estimated {[METHODS[i] for i in order_est]}",
     "pass": order_est == order_true},
    {"claim": "both practically inferior methods are rejected at the margin",
     "value": f"naive p = {m3.loc[1, 'p_minimum_effect']:.2e}, "
              f"mean p = {m3.loc[2, 'p_minimum_effect']:.2e}",
     "pass": bool(m3.loc[1:2, "beats_margin"].all())},
    {"claim": "the practically negligible difference is NOT claimed as important",
     "value": f"AR(1) vs AR(2): two-sided p = {conf.loc[0, 'p_paired_t']:.2e} "
              f"(significant) but minimum-effect p = "
              f"{m3.loc[0, 'p_minimum_effect']:.3f} (not important)",
     "pass": bool(not m3.loc[0, "beats_margin"])},
    {"claim": "the MCS retains only models close to the best",
     "value": f"{sorted(res_mcs['included'])}",
     "pass": "mean" not in res_mcs["included"] and "naive" not in res_mcs["included"]},
])
print()
print(checks.to_string(index=False))

os.makedirs(os.path.join(RESULTS, "E5"), exist_ok=True)
pd.DataFrame(L, columns=METHODS).assign(seed=np.arange(J_FINAL)).to_csv(
    os.path.join(RESULTS, "E5", "generality_demo_losses.csv"), index=False)
plan_tbl.to_csv(os.path.join(RESULTS, "E5", "planning_table.csv"), index=False)
conf.to_csv(os.path.join(RESULTS, "E5", "contrast_results.csv"), index=False)
m3.to_csv(os.path.join(RESULTS, "E5", "minimum_effect_tests.csv"), index=False)
tbl.to_csv(os.path.join(RESULTS, "E5", "mcs_table.csv"), index=False)
checks.to_csv(os.path.join(RESULTS, "E5", "verification.csv"), index=False)
print("\nwrote results/E5/")
'''

DRAFT = r'''
draft = f"""
### Draft for Section 6: a demonstration outside causal inference

Nothing in the design layer of Section 3 refers to treatment effects. To make that
concrete, and to make it checkable, we apply the identical workflow to one-step-ahead
forecasting of a stationary AR(1) process, where the competing forecasters can be
ranked in closed form. With y_t = {PHI} y_(t-1) + e_t and e_t ~ N(0, {SIGMA}), the
population one-step-ahead mean squared errors are {ASYMPTOTIC_MSE['AR(1)']:.3f} for an
AR(1) forecast at the true parameter, {ASYMPTOTIC_MSE['naive']:.3f} for the random walk
and {ASYMPTOTIC_MSE['mean']:.3f} for the sample-mean forecast; an AR(2) fitted to AR(1)
data is asymptotically tied with the AR(1) but pays one further estimation penalty of
approximately sigma^2 / n at any finite training length, here
{abs(TRUE_CONTRAST['AR(2)']):.4f}. Before seeing any data we declared a practical
margin of {MARGIN} on the loss scale, five per cent of the innovation variance, and used
it both as the planning alternative and as the margin of the final claim.

The pilot of {J0} independent replications, after variance-uncertainty inflation and
Bonferroni planning at alpha/3, required J = {J_FINAL} replications, set by the
precision target rather than the power target and driven by the widest contrast
(AR(1) against the sample mean). On the confirmatory block the ranking of all four
forecasters was recovered exactly, every estimate agreed with its finite-sample
reference to within {conf['abs_err_vs_reference'].max():.4f}, and every Monte Carlo
standard error met its target.

The AR(1)-versus-AR(2) contrast is the instructive one. Its estimate is
{conf.loc[0, 'estimate']:+.4f} with a Monte Carlo interval of
[{conf.loc[0, 'ci_lo']:+.4f}, {conf.loc[0, 'ci_hi']:+.4f}], and the two-sided test of
equality rejects overwhelmingly (p = {conf.loc[0, 'p_paired_t']:.1e}). It would be
wrong to report that as evidence that AR(1) is the better forecaster in any useful
sense: the difference is {abs(conf.loc[0, 'estimate']) / MARGIN:.0%} of the margin we
declared to matter, and the minimum-effect test against that margin does not reject
(p = {m3.loc[0, 'p_minimum_effect']:.3f}). The two genuinely inferior forecasters are
rejected at the margin (p < {max(m3.loc[1, 'p_minimum_effect'], m3.loc[2, 'p_minimum_effect']):.0e}),
and the Model Confidence Set eliminates both. This is the distinction of Section 3.4
made visible on a problem where the truth is known: a sufficiently large J makes any
non-zero difference significant, so the design must be sized for precision, and the
claim must be tested against a margin declared in advance.

The point is not that forecasting is a novel application. It is that the estimand
table, the paired contrasts, the two-stage sizing and the multiple-comparison layer
transfer without modification, and that on a problem whose answer is known they return
it -- including the part of the answer that says a statistically certain difference is
not worth acting on.
"""
with open(os.path.join(RESULTS, "E5", "section_6_draft.md"), "w") as f:
    f.write(draft)
print(draft)
'''

FIGURE = r'''
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))

ax = axes[0]
pos = np.arange(len(METHODS))
ax.boxplot([L[:, k] for k in range(len(METHODS))], positions=pos, widths=.55,
           showfliers=False, medianprops=dict(color="black"))
ax.scatter(pos, [FINITE_MSE[m] for m in METHODS], color="crimson", zorder=4,
           marker="D", s=45, label="finite-sample reference")
ax.scatter(pos, [ASYMPTOTIC_MSE[m] for m in METHODS], color="0.35", zorder=4,
           marker="x", s=40, label="asymptotic MSE")
ax.set_xticks(pos)
ax.set_xticklabels(METHODS)
ax.set_ylabel("mean squared one-step-ahead error")
ax.set_title(f"Replicate-level losses, J = {J_FINAL}")
ax.legend(frameon=False)

ax = axes[1]
y = np.arange(len(conf))
ax.errorbar(conf["estimate"], y,
            xerr=[conf["estimate"] - conf["ci_lo"], conf["ci_hi"] - conf["estimate"]],
            fmt="o", color=PALETTE[0], capsize=4, label="estimate and Monte Carlo CI")
ax.scatter(conf["true"], y, marker="D", color="crimson", zorder=4,
           label="finite-sample reference")
ax.axvline(-MARGIN, color="0.45", ls=":", lw=1.2, label="practical margin")
ax.axvline(0, color="0.4", lw=1)
ax.set_yticks(y)
ax.set_yticklabels(conf["contrast"])
ax.invert_yaxis()
ax.set_xlabel(r"paired contrast $E[D]$, negative favours AR(1)")
ax.set_title("Estimate first, p-value second")
ax.legend(frameon=False, fontsize=10)

ax = axes[2]
ax.plot(plan_tbl["contrast"], plan_tbl["J_precision"], marker="o",
        color=PALETTE[1], label="precision target")
ax.plot(plan_tbl["contrast"], plan_tbl["J_power"], marker="s",
        color=PALETTE[2], label="power target")
ax.axhline(J_FINAL, color="0.35", ls="--", label=f"$J_{{final}}$ = {J_FINAL}")
ax.set_ylabel("required $J$")
ax.set_title("Which target binds")
ax.tick_params(axis="x", rotation=20)
ax.legend(frameon=False, fontsize=10)

fig.tight_layout()
out = os.path.join(FIGURES, "Fig_E5_generality.png")
fig.savefig(out)
print("wrote", out)
download(out)
'''


def build():
    cells = []
    md(cells, """
        # E5: the same workflow outside causal inference (P3-T6)

        SoutoNeto's third point asks why the method is presented only for treatment
        effects when it appears to be general. This notebook runs the complete TISCA
        v2 workflow (estimand, independent pilot, variance-inflated two-stage sizing,
        paired contrasts, Romano-Wolf family inference, Model Confidence Set) on a
        forecasting problem, and checks every output against a closed-form answer.

        The competitors are an AR(1) fit, an AR(2) fit, a random walk and a
        sample-mean forecast, on data generated by a stationary AR(1). AR(1) and
        AR(2) are asymptotically tied and the other two are strictly worse, so the
        correct Model Confidence Set is known in advance. Runtime: about a minute.
        """)
    code(cells, SETUP)
    code(cells, PLOT_STYLE)
    code(cells, DGP)
    md(cells, """
        ## The replication function

        One simulated series per replication, all four forecasters evaluated on that
        same series. This is what makes the contrasts paired, and it is the reason
        the design layer transfers unchanged from the causal setting.
        """)
    code(cells, SIM)
    md(cells, """
        ## Stage 1: pilot and planning
        """)
    code(cells, PLAN)
    md(cells, """
        ## Stage 2: confirmatory run and paired inference
        """)
    code(cells, CONFIRM)
    md(cells, """
        ## Model Confidence Set
        """)
    code(cells, MCS_CELL)
    md(cells, """
        ## Verification against the analytic answer

        These are the checks that make the section worth including: a workflow that
        cannot recover a known ranking has not been demonstrated, only exercised.
        """)
    code(cells, VERIFY)
    md(cells, """
        ## Draft text for Section 6
        """)
    code(cells, DRAFT)
    md(cells, """
        ## Figure
        """)
    code(cells, FIGURE)
    write("E5_generality_demo.ipynb", cells)


if __name__ == "__main__":
    build()
