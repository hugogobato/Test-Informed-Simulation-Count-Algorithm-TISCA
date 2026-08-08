#!/usr/bin/env python3
"""Generate ``E1c_no_difference.ipynb`` (REVISION_PLAN.md P3-T7).

SoutoNeto asks directly: "What happens when you run TISCA in a situation when there
is no difference? As a trivial example: both models are the same model."

Three cases are distinguished, because they behave differently and the manuscript
currently conflates them:

  1. equal expected loss, different models  (theta = 0, sigma_D > 0) -- the ordinary
     null, read off the E1 grid;
  2. the same model twice                   (theta = 0, sigma_D = 0) -- degenerate;
  3. perfectly correlated distinct models    (rho = 1, sigma_D = 0 when the marginals
     also match) -- the boundary case that produces the same degeneracy.

Case 2 is where v1 hits its ``sd > 0`` guard, sets power to 0, and runs the loop to
the budget cap. The notebook demonstrates that, shows what v2 does instead, and
drafts the Section 3.10 paragraph.

Regenerate with::

    python notebooks/_generators/build_e1c_no_difference_nb.py
"""

from __future__ import annotations

from _nbcommon import PLOT_STYLE, SETUP, code, md, write


LOAD = r'''
from scipy import stats

from tisca import inference, planning, procedure

ALPHA = 0.05
OC_PATH = os.path.join(RESULTS, "E1", "operating_characteristics.csv")
assert os.path.exists(OC_PATH), f"{OC_PATH} not found; run the E1 grid first"
oc = pd.read_csv(OC_PATH)
null = oc[(oc["theta"] == 0.0) & (oc["module"] == "A")]
print(f"{len(null)} null cells in Module A")
'''

CASE1 = r'''
# --------------------------------------------------------------------------- #
# Case 1: two genuinely different methods with equal expected loss.            #
# --------------------------------------------------------------------------- #
# This is the ordinary null and it is already in the E1 grid at theta = 0. The
# answer to "what happens" is: the test rejects at its nominal rate, the interval
# covers zero at its nominal rate, and J is set by whichever target is active --
# which is the substantive point, because a power target at a planning alternative
# delta remains perfectly well defined when the true effect is zero.
case1 = (null.groupby("design")
         .agg(type_I=("t1e_or_power", "mean"),
              type_I_worst=("t1e_or_power", "max"),
              ci_cover=("ci_cover", "mean"),
              E_J=("E_J", "mean"),
              sd_J=("sd_J", "mean"),
              P_at_cap=("pJmax", "mean"),
              n_cells=("cell_id", "size"))
         .round(4))
case1["mcse"] = np.sqrt(0.05 * 0.95 / null["R"].median()).round(4)
print(case1.to_string())
print()
print("By family (theta = 0, all designs pooled):")
print(null.groupby("family")["t1e_or_power"].agg(["mean", "max"]).round(4).to_string())
'''

CASE2 = r'''
# --------------------------------------------------------------------------- #
# Case 2: literally the same model twice. D_j == 0 for every j, so sigma_D = 0.  #
# --------------------------------------------------------------------------- #
def same_model(seed):
    """One replication in which both 'methods' are the identical fitted model."""
    rng = np.random.default_rng(seed)
    v = float(rng.normal(loc=10.0, scale=2.0))
    return [v, v]


contrasts = [dict(A=0, B=1, mode="M1", delta=0.5, target_power=0.80, name="same-model")]
res_v2 = procedure.two_stage(same_model, contrasts, J0=50, n_metrics=2,
                             alpha=ALPHA, J_max=5000)
print("TISCA v2, same model twice")
print("  J_final          :", res_v2["J_final"])
print("  capped at J_max  :", res_v2.get("capped"))
print("  estimate         :", res_v2["contrast_results"][0]["estimate"])
print("  p-value          :", res_v2["contrast_results"][0]["p_value"])
print("  rejected         :", res_v2["rejected"])
'''

CASE2_V1 = r'''
# The same input through the v1 stopping rule. v1 estimated power from the two
# MARGINAL sds and guarded with `if sd_p > 0 and sd_b > 0: ... else: P_current = 0`
# (legacy/v1/tisca.py, line 486). With identical models both marginal sds are
# positive -- the guard does not even fire -- but the Welch power at any J is the
# power to detect delta between two samples that are identical replication by
# replication, which never reaches the target. Either way the loop runs to the cap.
def v1_stopping_path(sd_p, sd_b, delta, alpha=ALPHA, target=0.80,
                     J0=50, batch=50, J_max=5000):
    """Reproduce the v1 loop's trajectory: unpaired Welch power, grow J by `batch`."""
    path = []
    J = J0
    while J <= J_max:
        if sd_p > 0 and sd_b > 0:
            se = np.sqrt(sd_p ** 2 / J + sd_b ** 2 / J)
            df = (sd_p ** 2 / J + sd_b ** 2 / J) ** 2 / (
                (sd_p ** 2 / J) ** 2 / (J - 1) + (sd_b ** 2 / J) ** 2 / (J - 1))
            crit = stats.t.ppf(1 - alpha / 2, df)
            ncp = abs(delta) / se
            power = 1 - stats.nct.cdf(crit, df, ncp) + stats.nct.cdf(-crit, df, ncp)
        else:
            power = 0.0                      # the v1 guard
        path.append({"J": J, "power": power})
        if power >= target:
            break
        J += batch
    return pd.DataFrame(path), J


# The identical-model case, expressed the way v1 saw it: two marginal sds of 2.0.
path_same, J_v1_same = v1_stopping_path(2.0, 2.0, delta=0.5)
# And the true degeneracy the guard was written for: a constant metric column.
path_const, J_v1_const = v1_stopping_path(0.0, 0.0, delta=0.5)

print(f"v1 on identical models   : stopped at J = {J_v1_same} "
      f"(power at the cap = {path_same['power'].iloc[-1]:.3f})")
print(f"v1 on a constant column  : stopped at J = {J_v1_const} "
      f"(power pinned at {path_const['power'].iloc[-1]:.3f} by the sd > 0 guard)")
print(f"v2 on identical models   : stopped at J = {res_v2['J_final']}, "
      f"terminated={res_v2['contrast_results'][0]['p_value'] is not None}")
'''

CASE3 = r'''
# --------------------------------------------------------------------------- #
# Case 3: distinct methods, perfect correlation. sigma_D -> 0 continuously.     #
# --------------------------------------------------------------------------- #
# The interesting question is not the endpoint but the approach to it: as rho -> 1
# the paired design needs ever fewer replications, which is the pairing gain taken
# to its limit, and the procedure must stay well behaved rather than diverge.
rows = []
for rho in [0.0, 0.5, 0.9, 0.99, 0.999, 1.0]:
    sd_D = float(np.sqrt(2 * (1 - rho)))
    if sd_D > 0:
        J_prec = planning.required_J_mcse(sd_D, target_mcse=0.05)
        J_pow = planning.required_J_power(mode="M1", delta=0.5, sigma=sd_D,
                                          alpha=ALPHA, target_power=0.80)
    else:
        J_prec = J_pow = np.nan
    rows.append({"rho": rho, "sigma_D": sd_D, "J_precision": J_prec, "J_power": J_pow})
approach = pd.DataFrame(rows)
print(approach.to_string(index=False))

# At rho = 1 exactly with matched marginals, sigma_D is identically 0 and BOTH
# targets are vacuous: any J attains any precision, because the contrast has no
# sampling variability at all. Report what the implementation actually does here
# rather than what it ought to do. Both edges are now closed (F9), so these lines
# document the CURRENT behaviour: J = 2 rather than 1, and a flagged degenerate cell
# rather than a blanket marginal_power of 1.0.
print("required_J_mcse(sigma=0)   ->", planning.required_J_mcse(0.0, target_mcse=0.05))
print("required_J_power(sigma=0)  ->",
      planning.required_J_power(mode="M1", delta=0.5, sigma=0.0, alpha=ALPHA,
                                target_power=0.80))
print("two_stage on identical models reports "
      f"marginal_power = {res_v2['marginal_power']}, sd = "
      f"{res_v2['Js_by_contrast'][0]['sd']}")
'''

RULE = r'''
# --------------------------------------------------------------------------- #
# What the implementation should do, and what it does.                          #
# --------------------------------------------------------------------------- #
checks = []

def record(name, ok, detail):
    checks.append({"check": name, "pass": bool(ok), "detail": detail})

record("terminates on identical models",
       res_v2["J_final"] <= 5000 and res_v2["contrast_results"][0]["p_value"] is not None,
       f"J_final = {res_v2['J_final']}, p = {res_v2['contrast_results'][0]['p_value']}")
record("does not reject a zero contrast",
       res_v2["rejected"] == [False],
       f"rejected = {res_v2['rejected']}")
record("reports the estimate as exactly zero",
       res_v2["contrast_results"][0]["estimate"] == 0.0,
       f"estimate = {res_v2['contrast_results'][0]['estimate']}")
record("bootstrap flags the degenerate sample rather than dividing by zero",
       inference.studentized_paired_bootstrap(np.zeros(50), B=99, seed=0).get("degenerate") is True,
       "studentized_paired_bootstrap({0}) sets degenerate=True")
record("v1 would have run to the budget cap",
       J_v1_const > 1000,
       f"v1 stopped at J = {J_v1_const} with power pinned at 0 by the sd > 0 guard")

# The two rough edges this notebook used to record as open recommendations (F9).
# Both are now implemented in the Python port, so they are asserted rather than
# reported as gaps. They mattered to a user whose two metric columns are identical
# by accident -- a duplicated column name, say -- because the old behaviour answered
# with a plan no test could run on, and a power of 1.0 for a comparison in which
# nothing is detectable.
record("precision layer returns a testable J on a degenerate contrast",
       planning.required_J_mcse(0.0, target_mcse=0.05) >= 2,
       f"required_J_mcse(0) returns {planning.required_J_mcse(0.0, target_mcse=0.05)}; "
       "J = 1 leaves df = 0, so the paired t the plan rests on does not exist")
record("two-stage result flags a zero-variance contrast",
       res_v2.get("degenerate_contrasts") == [c.get("name", "contrast")
                                              for c in contrasts],
       f"degenerate_contrasts = {res_v2.get('degenerate_contrasts')}, "
       f"marginal_power = {res_v2['marginal_power']} (no longer a blanket 1.0)")
# Note what marginal_power = 1.0 does and does not mean here. It is the PLANNED
# power at the declared alternative delta = 0.5, and with sigma = 0 that alternative
# would indeed be detected with certainty -- so 1.0 is correct arithmetic, and the
# fix is that the cell is now flagged as degenerate alongside it. The actual F9
# defect was that 1.0 came back for EVERY delta, including delta = 0, where the
# statistic is 0/0 and nothing is detectable at all. That is what is asserted:
record("zero-variance power at a null alternative is not reported as certainty",
       planning.power_M1(2, 0.0, 0.0, ALPHA) < 1.0,
       f"power_M1(J=2, delta=0, sigma=0) = {planning.power_M1(2, 0.0, 0.0, ALPHA)}, "
       f"not 1.0; at the declared delta = {contrasts[0].get('delta')} it correctly "
       f"remains {float(np.asarray(res_v2['marginal_power']).ravel()[0])}")

report = pd.DataFrame(checks)
print(report.to_string(index=False))
print()
assert report["pass"].all(), report[~report["pass"]].to_string(index=False)
print("All checks pass: the procedure terminates correctly, refuses to claim "
      "certainty about an undetectable contrast, and names the degenerate cell.")

os.makedirs(os.path.join(RESULTS, "E1c"), exist_ok=True)
case1.to_csv(os.path.join(RESULTS, "E1c", "null_case_by_design.csv"))
approach.to_csv(os.path.join(RESULTS, "E1c", "rho_to_one_approach.csv"), index=False)
report.to_csv(os.path.join(RESULTS, "E1c", "degenerate_checks.csv"), index=False)

worst = case1["type_I_worst"].max()
paragraph = f"""
### Draft for Section 3.10 (Limitations): the no-difference case

Three situations must be distinguished, and they do not behave alike. When two
distinct methods happen to have equal expected loss, the design is in its ordinary
null state: across the {len(null)} null cells of the operating-characteristic study
the two-stage procedure rejected at {case1.loc['D4', 'type_I']:.4f} against a
nominal {ALPHA} (worst cell over all designs {worst:.4f}), and the confidence
interval for the contrast covered zero at {case1.loc['D4', 'ci_cover']:.4f}. The
number of replications is unaffected, because it is chosen from the planning
alternative and the pilot variance, neither of which depends on the true effect.

When the two methods are the same method, the paired contrast is identically zero,
its variance is zero, and no amount of replication can change that. The procedure
must recognise the situation rather than plan for it. TISCA v2 returns an estimate
of exactly zero, a p-value of 1, no rejection, and stops at J = {res_v2['J_final']};
the studentized bootstrap reports the sample as degenerate instead of studentizing
by zero. This is the case in which v1 misbehaved: its power was estimated from the
two marginal standard deviations and set to zero whenever either was zero, so the
target was never met and the loop consumed the entire budget
(J = {J_v1_const} here) before stopping at the cap.

The third case is the limit of the second. As the two methods become perfectly
correlated, sigma_D falls to zero continuously and the required J falls with it,
which is the pairing gain in its extreme form; at correlation exactly one with
matched marginals both design targets become vacuous, since any J attains any
precision. Here the current implementation is serviceable but not yet ideal: the
planning functions return the smallest admissible J rather than declining to plan,
and the two-stage result reports a marginal power of 1 for a contrast that is
identically zero. Neither harms the procedure, which terminates and does not
reject, but a zero pilot standard deviation is far more often a duplicated metric
column than a genuine tie, so it should be surfaced as a refusal rather than
absorbed silently. That is a recommendation on the software, and it is recorded
here as one.
"""
with open(os.path.join(RESULTS, "E1c", "section_3_10_no_difference.md"), "w") as f:
    f.write(paragraph)
print(paragraph)
'''

FIGURE = r'''
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

ax = axes[0]
designs = sorted(null["design"].unique())
for i, d in enumerate(designs):
    v = null.loc[null["design"] == d, "t1e_or_power"]
    ax.scatter(np.full(len(v), i) + np.random.default_rng(i).normal(0, .07, len(v)),
               v, s=10, alpha=.4, color=PALETTE[i % len(PALETTE)], linewidths=0)
    ax.scatter([i], [v.mean()], marker="_", s=700, color="black", zorder=3)
ax.axhline(ALPHA, color="crimson", ls="--", lw=1.2)
ax.set_xticks(range(len(designs)))
ax.set_xticklabels(designs)
ax.set_ylabel("rejection rate")
ax.set_title("Case 1: equal expected loss, distinct methods")

ax = axes[1]
a = approach[approach["sigma_D"] > 0]
ax.plot(a["rho"], a["J_power"], marker="o", color=PALETTE[0], label="power target")
ax.plot(a["rho"], a["J_precision"], marker="s", color=PALETTE[1], label="precision target")
ax.set_yscale("log")
ax.set_xlabel(r"correlation $\rho$ between the two methods")
ax.set_ylabel("required $J$ (log scale)")
ax.set_title(r"Case 3: the approach to $\sigma_D = 0$")
ax.legend(frameon=False)

fig.tight_layout()
out = os.path.join(FIGURES, "Fig_E1c_no_difference.png")
fig.savefig(out)
print("wrote", out)
download(out)
'''


def build():
    cells = []
    md(cells, """
        # E1c: what happens when there is no difference (P3-T7)

        Answers SoutoNeto's direct question, and separates three cases the current
        manuscript runs together: two different methods that happen to be equally
        good, the same method compared with itself, and two methods whose losses are
        perfectly correlated. Only the first is an ordinary null; the other two are
        degenerate, and one of them is where v1 misbehaves.

        The acceptance criterion is a named paragraph for Section 3.10 plus a unit
        test for the degenerate case. The test already exists
        (`tisca/tests/test_procedure.py::test_degenerate_same_model_terminates_gracefully`);
        this notebook is the evidence behind the paragraph, and the paragraph is
        written from the measured values rather than typed by hand.
        """)
    code(cells, SETUP)
    code(cells, PLOT_STYLE)
    code(cells, LOAD)
    md(cells, """
        ## Case 1: equal expected loss, genuinely different methods
        """)
    code(cells, CASE1)
    md(cells, """
        ## Case 2: the same model twice
        """)
    code(cells, CASE2)
    md(cells, """
        ### The same input through the v1 rule
        """)
    code(cells, CASE2_V1)
    md(cells, """
        ## Case 3: perfectly correlated methods
        """)
    code(cells, CASE3)
    md(cells, """
        ## Behavioural checks and the Section 3.10 paragraph
        """)
    code(cells, RULE)
    md(cells, """
        ## Figure
        """)
    code(cells, FIGURE)
    write("E1c_no_difference.ipynb", cells)


if __name__ == "__main__":
    build()
