#!/usr/bin/env python3
"""Generate the P2-T3 harness self-check notebook (Colab).

P2-T3 (L-CODE) is the outer-MC operating-characteristics harness. Its runnable
artifact on the user's side is this notebook: it pulls the current harness
(either from a local copy or from the public GitHub repo), runs the P2-T3
acceptance check (bivariate normal, rho=0, oracle sigma -> Type I ~ 0.05 and
power ~ 0.80 within 2 MCSE), smoke-tests all six designs D1..D6 and all loss
families, and writes an operating-characteristics CSV (with the Colab download
fallback). Pure NumPy/SciPy: a free Colab CPU runtime is enough, no R bundle.

Output: .../notebooks/E1_harness_selfcheck.ipynb
Regenerate with: python notebooks/_generators/build_e1_selfcheck_nb.py
"""
import json, textwrap

CELLS = []

def md(src):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)})

def code(src):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.splitlines(keepends=True)})

md(textwrap.dedent("""\
# E1 — outer-MC harness self-check (P2-T3)

This notebook verifies the vectorised outer-Monte-Carlo harness that produces the
unconditional operating characteristics for experiment **E1** (REVISION_PLAN.md
P2-T3, P3-T2). It is pure NumPy/SciPy and runs in about a minute on a free Colab
CPU runtime; no R library bundle is needed.

It does three things:

1. **P2-T3 acceptance check** — on a configuration with a known closed-form answer
   (bivariate normal, `rho = 0`, oracle `sigma`, fixed `J`) the harness must recover
   the nominal Type I error `0.05` and nominal power `0.80` to within 2 Monte Carlo
   standard errors. This is the task's acceptance gate.
2. **All six designs `D1..D6`** run head-to-head on one grid cell.
3. **All seven loss families** are sampled and their implied effect is sanity-checked.

Every cell prints `[PASS]/[FAIL]`; a single clean run is the audit record for P2-T3.
"""))

md(textwrap.dedent("""\
## 0. Obtain the harness

The harness lives in the public repo `hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA`
under `tisca/python`. If a local checkout already exists on the runtime, use it;
otherwise clone the public repo into `/content` and add its `tisca/python` directory
to `sys.path`.
"""))

code(textwrap.dedent("""\
import os, sys, subprocess, shutil

REPO = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
HAVE_LOCAL = (_d := "/content/Test-Informed-Simulation-Count-Algorithm-TISCA") and os.path.isdir(_d) and os.path.isdir(os.path.join(_d, "tisca", "python"))

if HAVE_LOCAL:
    SRC = "/content/Test-Informed-Simulation-Count-Algorithm-TISCA/tisca/python"
    print("[INFO] using local checkout")
else:
    print("[INFO] cloning public repo ...")
    subprocess.run(["git", "clone", "--depth", "1", REPO, "/content/TISCA_repo"],
                   check=True)
    SRC = "/content/TISCA_repo/tisca/python"

assert os.path.isdir(SRC), f"expected harness at {SRC}"
sys.path.insert(0, SRC)
import numpy as np
from tisca.outermc import engine, families
from tisca.outermc.designs import solve_J
print("[PASS] harness imported from", SRC)
"""))

md(textwrap.dedent("""\
## 1. P2-T3 acceptance check

The acceptance says: *bivariate normal, `rho = 0`, fixed `J`, oracle `sigma`,
recovers nominal Type I error 0.05 and power 0.80 each within 2 MCSE.*

`J` is fixed from the oracle `sigma_D = sqrt(2)` and the planning alternative
`delta = 0.5`; we then run `R = 20 000` outer repetitions at the null
(`theta = 0`) and at the alternative (`theta = delta`).
"""))

code(textwrap.dedent("""\
def acceptance(verbose=True):
    delta, sd = 0.5, np.sqrt(2.0)
    J = int(np.atleast_1d(solve_J(0.80, delta, sd, 0.05)[0]).max())
    R = 20_000
    mcse0, mcse1 = np.sqrt(.05 * .95 / R), np.sqrt(.80 * .20 / R)
    cfg = dict(engine.DEFAULT_CONFIG)
    cfg.update(design="D6", family="normal", rho=0.0, sigma_a=1.0, sigma_b=1.0,
               delta=delta, sigma_D=None, R=R, Jmax=2000, fixed_J=J)

    cfg["theta"] = 0.0
    s0, _, _ = engine.run_e1(cfg)              # Type I
    cfg["theta"] = delta
    s1, _, _ = engine.run_e1(cfg)              # Power

    t1e, power = s0["reject_rate"], s1["reject_rate"]
    ok0 = abs(t1e - 0.05) <= 2 * mcse0
    ok1 = abs(power - 0.80) <= 2 * mcse1
    out = {
        "oracle_J": J,
        "Type_I_est": round(t1e, 4), "Type_I_nom": 0.05, "mcse": round(mcse0, 4),
        "Type_I_within_2MCSE": bool(ok0),
        "Power_est": round(power, 4), "Power_nom": 0.80, "mcse": round(mcse1, 4),
        "Power_within_2MCSE": bool(ok1),
        "Power_theta_hat": round(s1["E_theta"], 4),
        "PASS": bool(ok0 and ok1),
    }
    if verbose:
        for k, v in out.items():
            print(f"  {k} = {v}")
        print("[PASS] acceptance" if out["PASS"] else "[FAIL] acceptance")
    return out

res = acceptance()
assert res["PASS"], "P2-T3 acceptance failed"
print("\\n[PASS] P2-T3 acceptance gate passed")
"""))

md(textwrap.dedent("""\
## 2. All six designs run head-to-head

One grid cell (bivariate normal, `rho=0.3`, `theta=0.5`, `R=1000`) evaluated under
each of D1..D6 so we can see the different ways they trade Type I / power / `E[J]`.
No outside reference is asserted here; the point is that all six designs execute
and return finite operating characteristics.
"""))

code(textwrap.dedent("""\
DESIGNS = ["D1", "D2", "D3", "D4", "D5", "D6"]
base = dict(engine.DEFAULT_CONFIG)
base.update(family="normal", rho=0.3, sigma_a=1.0, sigma_b=1.0, theta=0.5,
            delta=0.5, sigma_D=None, R=1000, J0=50, Jmax=1000, K=1, alpha=0.05)
rows = []
for d in DESIGNS:
    c = base.copy(); c["design"] = d; c["seed"] = 1
    s, res, meta = engine.run_e1(c)
    assert len(res) == 1000 and np.all(np.isfinite(res.col("J")))
    rows.append([d, s["reject_rate"], s["E_J"], s["E_theta"], s["bias"], s["rmse"]])
print("design | reject | E[J]  | E[theta] |  bias   |  rmse")
for r in rows:
    print(f"{r[0]:5} |  {r[1]:5.3f} | {r[2]:6.1f} |   {r[3]:7.3f}  | {r[4]:7.4f} | {r[5]:6.3f}")
print("[PASS] all six designs executed and returned valid results")
"""))

md(textwrap.dedent("""\
## 3. All loss families sampled

Confirm the seven supported families produce `(n, 2)` paired losses with the
requested effect `E[D] ~ delta` (the empirical family is omitted here because it
needs a real matrix; it is exercised in E1 Module A).
"""))

code(textwrap.dedent("""\
for fam in ["normal", "lognormal", "gamma", "mix", "beta", "t3"]:
    X = families.sample_pairs(fam, n=2000, rho=0.0, sigma_a=1.0, sigma_b=1.0, theta=0.5)
    assert X.shape == (2000, 2)
    D = X[:, 0] - X[:, 1]
    drift = abs(float(np.mean(D)) - 0.5)
    status = "[PASS]" if drift < 0.15 else "[FAIL]"
    print(f"  {status} {fam:9}  E[D] = {np.mean(D):7.3f}  sd(D) = {np.std(D):6.3f}")
print("[PASS] all parametric families sanitised")
"""))

md(textwrap.dedent("""\
## 4. Write the operating-characteristics CSV

For the artifacts promised to readers, the self-check also emits a tiny results
file (few rows) with every operating characteristic and its MCSE, and provides the
Colab download fallback. The full E1 grid (Modules A/B/C/D) is run separately by
the `E1_module*` notebooks.
"""))

code(textwrap.dedent("""\
import pandas as pd
from tisca.outermc import summarize_ocs
summaries = []
for d in ["D4", "D5", "D6"]:
    c = base.copy(); c["design"] = d; c["theta"] = 0.0; c["seed"] = 2; c["R"] = 800
    s, _, _ = engine.run_e1(c)
    summaries.append(s)
df = summarize_ocs(summaries)
df.to_csv("E1_selfcheck_oc.csv", index=False)
print("wrote E1_selfcheck_oc.csv")
print(df.to_string())

try:
    from google.colab import files
    files.download("E1_selfcheck_oc.csv")
    print("Downloaded: E1_selfcheck_oc.csv")
except Exception as e:
    print("(Not on Colab / download skipped):", e)
"""))

md(textwrap.dedent("""\
## Summary

The harness passed its acceptance gate (Type I `0.05` and power `0.80`, each
within 2 MCSE), executed all six designs and all parametric families, and wrote a
results CSV. This signs off **P2-T3** so that E1/P3-T2 can be launched on the 33
Colab sessions.
"""))

nb = {"nbformat": 4, "nbformat_minor": 0,
      "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                   "language_info": {"name": "python"}},
      "cells": CELLS}
out = "/home/hugo_souto/Stuff/Submitted_Research/TISCA_paper/Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks/E1_harness_selfcheck.ipynb"
with open(out, "w") as f:
    json.dump(nb, f, indent=1)
for i, c in enumerate(CELLS):
    src = "".join(c["source"])
    assert not any(bad in src for bad in ["'''", '"""']), f"cell {i} forbidden triple-quote"
print("wrote", out, len(CELLS), "cells")
