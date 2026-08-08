#!/usr/bin/env python3
"""Generate ``E1b_skewness_calibration.ipynb`` (REVISION_PLAN.md P3-T4).

Produces the evidence that replaces the deleted "J > 30" claim (IJDA comment 3):
the actual Type I error of the paired t as a function of J and of the standardised
third moment of D_j, with the studentized paired bootstrap overlaid, and the
Berry-Esseen bound shown as the explanatory frame rather than as a guarantee.

The E1 grid does not vary J directly (every design chooses its own J), so this
notebook runs its own small fixed-J sweep on top of the same loss families, and
reads E1 only for the design-level context.

Regenerate with::

    python notebooks/_generators/build_e1b_skewness_nb.py
"""

from __future__ import annotations

from _nbcommon import PLOT_STYLE, SETUP, code, md, write


CONFIG = r'''
from scipy import stats

from tisca import inference
from tisca.outermc import families

ALPHA = 0.05
J_GRID = [10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 400]
R_T = 40_000          # repetitions for the paired t (MCSE at 0.05 = 0.0011)
R_BOOT = 2_500        # repetitions for the bootstrap overlay (MCSE = 0.0044)
B_BOOT = 499          # resamples per repetition
RHOS = [0.0, 0.6]

MATRIX = pd.read_csv(os.path.join(REPO_ROOT, "legacy", "Paper_Experiments",
                                  "DGP1_500_results.csv"))[["mvbcf_pehe1", "bcf_pehe1"]].to_numpy(float)

# (family, rho) cells. The row-bootstrap family carries the data's own dependence,
# so it appears once with rho = None rather than at each design rho.
CELLS = [(f, r) for f in ["normal", "lognormal", "gamma", "beta", "t3", "mix",
                          "empirical_copula"] for r in RHOS]
CELLS += [("empirical", None)]
print(len(CELLS), "family x rho cells x", len(J_GRID), "values of J")
'''

SKEW = r'''
# Standardised third moment of D per cell, estimated once from a large independent
# draw. This is the x-axis of the figure, so it is measured rather than assumed.
skew_rows = []
for fam, rho in CELLS:
    g = families.contrast_skewness(fam, rho=rho, n=1_000_000, seed=17,
                                   matrix=MATRIX if fam.startswith("empirical") else None)
    skew_rows.append({"family": fam, "rho": rho, "skew_D": g,
                      "abs_skew": abs(g), "skew_mcse": np.sqrt(6 / 1_000_000)})
skew = pd.DataFrame(skew_rows).sort_values("abs_skew")
print(skew.round(4).to_string(index=False))
'''

SWEEP_T = r'''
def type_I_paired_t(family, rho, J, R, seed, chunk_elems=2_000_000):
    """Vectorised: R independent samples of J paired draws, one two-sided t each.

    Drawn in chunks of at most ``chunk_elems`` pairs. Only the rejection rate is
    needed, so materialising all R x J pairs at once buys nothing and costs a great
    deal: at R = 40,000 and J = 400 the single-shot version peaked at 2.4 GB, which
    is enough to push a 16 GB workstation into swap. Chunking holds it near 0.3 GB.
    """
    crit = stats.t.ppf(1 - ALPHA / 2, J - 1)
    per_chunk = max(1, int(chunk_elems // max(J, 1)))
    mat = MATRIX if family.startswith("empirical") else None
    n_rej, n_done, block_i = 0, 0, 0
    while n_done < R:
        r = min(per_chunk, R - n_done)
        block = families.sample_batch(family, r, J, rho=rho, theta=0.0,
                                      master_seed=seed + 1_000 * block_i, matrix=mat)
        D = block[..., 0] - block[..., 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = D.mean(axis=1) / (D.std(axis=1, ddof=1) / np.sqrt(J))
        n_rej += int(np.count_nonzero(np.abs(t) > crit))
        n_done += r
        block_i += 1
        del block, D, t
    return n_rej / R


rows = []
for ci, (fam, rho) in enumerate(CELLS):
    for J in J_GRID:
        p = type_I_paired_t(fam, rho, J, R_T, seed=900_000 + 100 * ci + J)
        rows.append({"family": fam, "rho": rho, "J": J, "test": "paired_t",
                     "type_I": p, "mcse": np.sqrt(p * (1 - p) / R_T), "R": R_T})
tt = pd.DataFrame(rows)
print(tt.pivot_table(index="J", columns=["family", "rho"], values="type_I", dropna=False)
      .round(4).to_string())
'''

SWEEP_BOOT = r'''
def _boot_type_I(D, B, rng, chunk=40):
    """Studentized paired bootstrap rejection indicator for each row of D.

    Vectorised over repetitions in chunks so the (chunk, B, J) index array stays
    small. The statistic is the same one ``inference.studentized_paired_bootstrap``
    computes: T_b = (mean_b - mean) / se_b, with the two-sided interval read off
    the T_b quantiles. Agreement with that reference implementation is checked
    below rather than assumed.
    """
    R, J = D.shape
    out = np.empty(R, dtype=bool)
    for start in range(0, R, chunk):
        d = D[start:start + chunk]                       # (c, J)
        c = d.shape[0]
        idx = rng.integers(0, J, size=(c, B, J))
        res = np.take_along_axis(d[:, None, :], idx, axis=2)     # (c, B, J)
        mb = res.mean(axis=2)
        sb = res.std(axis=2, ddof=1)
        m = d.mean(axis=1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            T = (mb - m) / (sb / np.sqrt(J))
        lo, hi = np.nanquantile(T, [ALPHA / 2, 1 - ALPHA / 2], axis=1)
        s = d.std(axis=1, ddof=1) / np.sqrt(J)
        obs = d.mean(axis=1) / s
        out[start:start + c] = (obs < lo) | (obs > hi)
    return out


rows_b = []
for ci, (fam, rho) in enumerate(CELLS):
    rng = np.random.default_rng(700_000 + ci)
    for J in J_GRID:
        block = families.sample_batch(fam, R_BOOT, J, rho=rho, theta=0.0,
                                      master_seed=800_000 + 100 * ci + J,
                                      matrix=MATRIX if fam.startswith("empirical") else None)
        D = block[..., 0] - block[..., 1]
        p = float(np.mean(_boot_type_I(D, B_BOOT, rng)))
        rows_b.append({"family": fam, "rho": rho, "J": J, "test": "studentized_bootstrap",
                       "type_I": p, "mcse": np.sqrt(p * (1 - p) / R_BOOT), "R": R_BOOT})
bt = pd.DataFrame(rows_b)
print(bt.pivot_table(index="J", columns=["family", "rho"], values="type_I", dropna=False)
      .round(4).to_string())
'''

CHECK = r'''
# The vectorised bootstrap above is a reimplementation, so check it against the
# package function on a handful of samples before any of it is believed.
rng = np.random.default_rng(1)
block = families.sample_batch("lognormal", 24, 60, rho=0.3, theta=0.0, master_seed=5)
D = block[..., 0] - block[..., 1]
mine = _boot_type_I(D, 999, np.random.default_rng(2))
theirs = np.array([
    (lambda r: r["ci"][0] > 0 or r["ci"][1] < 0)(
        inference.studentized_paired_bootstrap(d, B=999, seed=3))
    for d in D
])
agree = float(np.mean(mine == theirs))
print(f"vectorised vs tisca.inference agreement on {len(D)} samples: {agree:.3f}")
assert agree >= 0.85, "vectorised bootstrap disagrees with the reference implementation"
'''

RULE = r'''
sweep = pd.concat([tt, bt], ignore_index=True).merge(skew, on=["family", "rho"], how="left")


def cell_mask(frame, fam, rho):
    """Select one (family, rho) cell, treating the row-bootstrap's NaN rho correctly."""
    m = frame["family"] == fam
    return m & (frame["rho"].isna() if pd.isna(rho) else (frame["rho"] == rho))


# An EQUIVALENCE BAND, not an MCSE band. With R = 40,000 the MCSE is 0.0011, so a
# +/- 2 MCSE rule tolerates only +/- 0.0022 and flags ordinary Monte Carlo wobble as
# miscalibration: it declared that the paired t on a bivariate NORMAL contrast --
# where the test is exact at every J -- "needs J = 150". The band is therefore the
# larger of a practically meaningful deviation (0.005, i.e. a level in 0.045-0.055)
# and 2 MCSE. This is the same equivalence-band logic the manuscript adopts for
# coverage in change C6, applied to its own operating characteristics.
TOL_ABS = 0.005


def band(g, tol_mcse=2.0):
    return np.maximum(TOL_ABS, tol_mcse * g["mcse"])


def smallest_valid_J(g, tol_mcse=2.0):
    """Smallest J from which the level stays inside the equivalence band for all larger J."""
    g = g.sort_values("J")
    ok = (g["type_I"] - ALPHA).abs() <= band(g, tol_mcse)
    for i in range(len(g)):
        if ok.iloc[i:].all():
            return int(g["J"].iloc[i])
    return np.nan          # never calibrated anywhere in the grid -- reported as such


def direction_at_max_J(g):
    """Is the cell conservative, liberal, or calibrated at the largest J tested?

    The direction matters far more than the magnitude: a conservative test wastes
    replications, a liberal one invalidates the claim.
    """
    g = g.sort_values("J")
    p = g["type_I"].iloc[-1]
    tol = max(TOL_ABS, 2 * g["mcse"].iloc[-1])
    if p < ALPHA - tol:
        return "conservative"
    if p > ALPHA + tol:
        return "liberal"
    return "calibrated"


agg = (sweep.groupby(["test", "family", "rho"], dropna=False)
       .apply(lambda g: pd.Series({"J_min": smallest_valid_J(g),
                                   "direction": direction_at_max_J(g),
                                   "type_I_at_max_J": g.sort_values("J")["type_I"].iloc[-1]}),
              include_groups=False)
       .reset_index())
rule = agg.merge(skew, on=["family", "rho"], how="left")

# Berry-Esseen frame: the bound on the normal approximation error of the
# STANDARDISED MEAN is C * E|D - mu|^3 / (sigma^3 sqrt(J)) with C = 0.4748. It
# explains the sqrt(J) rate and the dependence on the third moment; it is NOT a
# guarantee for the t statistic (which studentizes by an estimated sd) and it is
# far from tight here. Shown for that reason only.
C_BE = 0.4748
rule["J_berry_esseen_0.005"] = np.ceil((C_BE * rule["abs_skew"] / 0.005) ** 2)
print(rule.sort_values(["test", "abs_skew"]).round(4).to_string(index=False))

os.makedirs(os.path.join(RESULTS, "E1b"), exist_ok=True)
sweep.to_csv(os.path.join(RESULTS, "E1b", "type_I_vs_J.csv"), index=False)
rule.to_csv(os.path.join(RESULTS, "E1b", "operational_rule.csv"), index=False)

# The sentence has to separate three things the old J > 30 rule conflated: cells
# that calibrate quickly, cells that need a large J, and cells that never calibrate
# within the grid at all -- and among the last, CONSERVATIVE and LIBERAL failures
# are not interchangeable. A conservative test costs replications; a liberal one
# costs the claim.
pt = rule[rule["test"] == "paired_t"].copy()
bs = rule[rule["test"] == "studentized_bootstrap"].copy()


def describe(df):
    reached = df[df["J_min"].notna()].sort_values("J_min")
    lib = df[df["direction"] == "liberal"]
    cons = df[df["direction"] == "conservative"]
    return reached, lib, cons


pt_ok, pt_lib, pt_cons = describe(pt)
bs_ok, bs_lib, bs_cons = describe(bs)
mild = pt_ok[pt_ok["abs_skew"] < 0.5]
worst = pt_ok.iloc[-1]


def cells(df):
    return ", ".join(f"{r.family}"
                     + ("" if pd.isna(r.rho) else f" (rho={r.rho:g})")
                     for r in df.itertuples()) or "none"


sentence = (
    f"Replacing the J > 30 heuristic. Calibration is judged against an equivalence "
    f"band of +/-{TOL_ABS} on the level (or 2 MCSE, whichever is wider), held for "
    f"every larger J in the grid. Of the {len(pt)} loss cells, the {len(mild)} with "
    f"|skew(D)| < 0.5 are calibrated for the paired t by J = {int(mild['J_min'].max())}, "
    f"and the most skewed cell of all -- the row bootstrap of the real MVBCF losses, "
    f"skew(D) = {float(pt[pt.family == 'empirical']['skew_D'].iloc[0]):+.2f} -- needs "
    f"J = {int(pt[pt.family == 'empirical']['J_min'].iloc[0])}. "
    f"Two cells never calibrate anywhere in the grid, and the direction differs by "
    f"test in a way that matters: the paired t is CONSERVATIVE on the catastrophic-"
    f"failure mixture ({cells(pt_cons)}), reaching only "
    f"{pt_cons['type_I_at_max_J'].min():.4f} at J = {max(J_GRID)}, whereas the "
    f"studentized bootstrap is LIBERAL on the same cells ({cells(bs_lib)}), reaching "
    f"{bs_lib['type_I_at_max_J'].max():.4f} against a nominal {ALPHA}. "
    f"The bootstrap is therefore not a general remedy for non-normality: it repairs "
    f"skewness but over-rejects under rare catastrophic failures, exactly where the "
    f"t errs in the safe direction. The operating requirement is a function of the "
    f"shape of D, which a pilot can estimate, and not a universal constant."
)
print()
print(sentence)
with open(os.path.join(RESULTS, "E1b", "operational_rule.md"), "w") as f:
    f.write("# P3-T4: finite-sample calibration and the rule replacing J > 30\n\n")
    f.write(sentence + "\n\n" + rule.round(4).to_markdown(index=False) + "\n")
'''

FIGURE = r'''
# Panelled by skewness: one panel per cell, ordered by |skew(D)|, both tests overlaid.
cells_sorted = (skew.sort_values("abs_skew")[["family", "rho"]]
                .itertuples(index=False, name=None))
cells_sorted = list(cells_sorted)
ncol = 5
nrow = int(np.ceil(len(cells_sorted) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.7 * nrow),
                         sharex=True, sharey=True)
axes = np.atleast_1d(axes).ravel()

for ax, (fam, rho) in zip(axes, cells_sorted):
    g = sweep[cell_mask(sweep, fam, rho)]
    for i, (test, lab) in enumerate([("paired_t", "paired $t$"),
                                     ("studentized_bootstrap", "stud. bootstrap")]):
        h = g[g["test"] == test].sort_values("J")
        ax.plot(h["J"], h["type_I"], marker="o", ms=3.5, color=PALETTE[i], label=lab)
        ax.fill_between(h["J"], h["type_I"] - 2 * h["mcse"], h["type_I"] + 2 * h["mcse"],
                        color=PALETTE[i], alpha=.18, linewidth=0)
    ax.axhline(ALPHA, color="crimson", lw=1.0, ls="--")
    ax.axvline(30, color="0.55", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 0.115)
    g_sk = float(skew[cell_mask(skew, fam, rho)]["skew_D"].iloc[0])
    lab = "data" if pd.isna(rho) else f"{rho:g}"
    ax.set_title(f"{fam}, $\\rho$={lab}\n$\\gamma_D$={g_sk:+.2f}", fontsize=10)

for ax in axes[len(cells_sorted):]:
    ax.set_visible(False)
axes[0].legend(frameon=False, fontsize=9, loc="upper right")
fig.supxlabel("replications $J$ (log scale); dotted line = the deleted $J>30$ heuristic")
fig.supylabel("unconditional Type I error")
fig.tight_layout()
out = os.path.join(FIGURES, "Fig_E1b_type_I_vs_J.png")
fig.savefig(out)
print("wrote", out)

# Companion panel: worst-case level over the J grid against |skew(D)|, which is the
# relationship the operational rule is stated in terms of.
fig2, ax = plt.subplots(figsize=(6.2, 4.2))
for i, (test, lab) in enumerate([("paired_t", "paired $t$"),
                                 ("studentized_bootstrap", "stud. bootstrap")]):
    r = rule[rule["test"] == test].sort_values("abs_skew")
    ax.plot(r["abs_skew"], r["J_min"], marker="o", color=PALETTE[i], label=lab)
ax.set_xlabel(r"$|\gamma_D|$, standardised third moment of the paired contrast")
ax.set_ylabel("smallest $J$ with a calibrated level")
ax.set_yscale("log")
ax.axhline(30, color="0.55", ls=":", label="the $J>30$ heuristic")
ax.legend(frameon=False)
ax.set_title("What the requirement actually depends on")
fig2.tight_layout()
out2 = os.path.join(FIGURES, "Fig_E1b_J_min_vs_skew.png")
fig2.savefig(out2)
print("wrote", out2)
download(out)
download(out2)
'''


def build():
    cells = []
    md(cells, """
        # E1b: skewness and finite-sample calibration (P3-T4)

        IJDA comment 3 removes the manuscript's "J > 30" justification and asks for
        the replacement to be *measured*: the actual Type I error of the paired
        contrast as a function of J and of the standardised third moment of D_j,
        with a studentized bootstrap alternative, and with Berry-Esseen used to
        explain the rate rather than to claim a guarantee.

        Two tests are swept over a J grid for every loss family: the paired t and the
        studentized paired bootstrap. The output is a figure panelled by skewness and
        an operational rule stated in terms of a quantity the pilot can estimate.

        Runtime is a few minutes; the bootstrap sweep dominates it.
        """)
    code(cells, SETUP)
    code(cells, PLOT_STYLE)
    code(cells, CONFIG)
    md(cells, """
        ## How skewed is each cell's contrast?

        `contrast_skewness` draws a million pairs per cell, so its own MCSE is
        `sqrt(6/n)` = 0.0024. Note the two empirical families differ sharply here:
        the row bootstrap keeps the real data's `skew(D)` near -1.55, while the
        copula variant is much closer to symmetric, because differencing two
        similarly shaped marginals under an exchangeable copula cancels most of the
        asymmetry.
        """)
    code(cells, SKEW)
    md(cells, """
        ## Sweep 1: the paired t
        """)
    code(cells, SWEEP_T)
    md(cells, """
        ## Sweep 2: the studentized paired bootstrap
        """)
    code(cells, SWEEP_BOOT)
    md(cells, """
        ## Cross-check the vectorised bootstrap against the package implementation
        """)
    code(cells, CHECK)
    md(cells, """
        ## The rule that replaces "J > 30"

        `J_min` is the smallest grid value from which the level stays within 2 MCSE
        of nominal for **every** larger J, so a single lucky crossing does not
        qualify. The Berry-Esseen column is the J at which the classical bound on
        the normal approximation of the standardised mean falls below 0.005. It is
        an upper bound for a different statistic (it does not cover studentization)
        and it is loose here; it is tabulated to show that the requirement scales
        with the third moment, which is the point the manuscript needs.
        """)
    code(cells, RULE)
    md(cells, """
        ## Figures
        """)
    code(cells, FIGURE)
    write("E1b_skewness_calibration.ipynb", cells)


if __name__ == "__main__":
    build()
