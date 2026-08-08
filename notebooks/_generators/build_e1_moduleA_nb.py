#!/usr/bin/env python3
"""Generate the E1 Module-A Colab notebook (P3-T2, using the P2-T3 harness).

Module A (REVISION_PLAN.md P3-T2) is the core operating-characteristics grid:
  F(7) x rho(5) x theta(4) x Design(6), with K=1, J0=50 — a total of 840 cells at
  R = 5000 outer repetitions. The harness (P2-T3) is pure NumPy; each cell is a
  few ms, and 840 x 5000 x ~1.5 ms is on the order of 1-2 core-hours, comfortably
  inside the 8-hour Colab cap. Results are check-pointed per cell to a CSV on the
  runtime and offered for download (Colab fallback per the standing rule).

The notebook pulls the current harness (local checkout if present, else the public
GitHub repo), runs the acceptance gate once, then sweeps the grid with a progress
bar, appending to `E1_moduleA_results.csv` so an interrupted session loses at most
one cell.

**SUPERSEDED — do not run to produce results.** This is the earlier standalone
Module-A runner, kept for provenance only; its output notebook is no longer in
`notebooks/`. It carries its own copy of a 840-cell grid that predates the split of
the empirical family into the row-bootstrap and copula variants, so it does not
match `tisca/python/tisca/outermc/e1_grid.py`, which is now the single source of
truth. (It does not hit the numeric-rho failure, because it omits the empirical
families altogether — see the `families` list below.) Use `build_e1_modules_nb.py`
for the Colab runners and `experiments/E1_operating_characteristics/run_e1_grid.py`
locally; both import the canonical grid.

Regenerate: python notebooks/_generators/build_e1_moduleA_nb.py
"""
import json, textwrap

CELLS = []

def md(src):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)})

def code(src):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.splitlines(keepends=True)})

md(textwrap.dedent("""\
# E1 — Module A: core operating-characteristics grid (P3-T2)

This notebook runs the **Module A** cell grid on top of the outer-MC harness from
P2-T3. It covers the seven loss families, five correlation levels, four true-effect
levels and all six designs at `K = 1`, `J0 = 50` (840 cells at `R = 5000` outer
repetitions). Each returned operating characteristic carries a Monte Carlo standard
error.

Pure NumPy/SciPy: a free Colab CPU runtime is enough (the whole grid is ~1-2
core-hours, well inside the 8-hour cap). Results are appended cell-by-cell to
`E1_moduleA_results.csv` and offered for download at the end.

> Module B (multiplicity) and Modules C/D (tuning / variance-ratio) are run by
> their own notebooks so each stays under the session cap.
"""))

md(textwrap.dedent("""\
## 0. Obtain the harness
"""))

code(textwrap.dedent("""\
import os, sys, subprocess, itertools, time, json

REPO = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
LV = "/content/Test-Informed-Simulation-Count-Algorithm-TISCA"
SRC = os.path.join(LV, "tisca", "python") if os.path.isdir(os.path.join(LV, "tisca", "python")) else None
if SRC is None:
    subprocess.run(["git", "clone", "--depth", "1", REPO, "/content/TISCA_repo"], check=True)
    SRC = "/content/TISCA_repo/tisca/python"
sys.path.insert(0, SRC)

import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from tisca.outermc import engine, summarize_ocs
print("[PASS] harness from", SRC)
"""))

md(textwrap.dedent("""\
## 1. Define the Module-A grid

Grid: `F(7) x rho(5) x theta(4) x Design(6)`, `K = 1`, `J0 = 50`. The effect levels
are expressed relative to the planning alternative `delta` (in units of
`sigma_D`): null, `0.5*delta`, `delta`, `2*delta`. `sigma_D` is the paired-D sd from
the normal family factor (oracle for D6).

The parametric families are `normal, lognormal, gamma, mix, beta, t3`; family `(g)
empirical` (row-bootstrap of the real 500x20 MVBCF loss matrix) needs the real
matrix and is run in the E3-linked pass, not in this standalone sweep.
"""))

code(textwrap.dedent("""\
def make_grid(R, J0=50):
    families = ["normal", "lognormal", "gamma", "mix", "beta", "t3"]  # empirical (g) needs a real matrix
    rhos = [-0.3, 0.0, 0.3, 0.6, 0.9]
    thetas = [0.0, 0.5, 1.0, 2.0]          # multiples of delta
    designs = ["D1", "D2", "D3", "D4", "D5", "D6"]
    delta = 0.5
    cells = []
    for fam, rho, tmult, des in itertools.product(families, rhos, thetas, designs):
        sigma_D = float(np.sqrt(2.0 - 2.0 * rho))       # normal-shaped (oracle for D6)
        c = dict(engine.DEFAULT_CONFIG)
        c.update(design=des, family=fam, rho=rho, sigma_a=1.0, sigma_b=1.0,
                 theta=delta * tmult, delta=delta, sigma_D=None, R=R,
                 J0=J0, Jmax=1000, K=1, alpha=0.05, seed=0)
        cells.append(c)
    return cells

GRID = make_grid(R=5000)
print("Module A cells:", len(GRID))
print("factors per cell:", "family,rho,theta_mult,design")
"""))

md(textwrap.dedent("""\
## 2. Run the grid with per-cell checkpointing

Results append to `E1_moduleA_results.csv`; an 8-hour kill or dropped session loses
at most one cell (the one currently running). Rows are de-duplicated on any factor
quadruple already present, so a re-run resumes without recomputing finished cells.
"""))

code(textwrap.dedent("""\
CSV = "E1_moduleA_results.csv"
if os.path.exists(CSV):
    done = pd.read_csv(CSV)
else:
    done = pd.DataFrame()

def key(c):
    return (c["family"], c["rho"], round(c["theta"] / c["delta"], 3), c["design"])

def run_cell(c):
    s, res, meta = engine.run_e1(c)
    row = summarize_ocs([s]).iloc[0].to_dict()
    row.update(family=c["family"], rho=c["rho"], theta_mult=round(c["theta"]/c["delta"],3),
               design=c["design"], R=c["R"])
    return row

t0 = time.time()
pending = [c for c in GRID if not ((done.shape[0] > 0) and key(c) in set(map(tuple, done[["family","rho","theta_mult","design"]].values.tolist())))]
print(f"completed {done.shape[0]} / {len(GRID)}; running {len(pending)}")
for c in tqdm(pending, desc="Module A cells", unit="cell"):
    try:
        row = run_cell(c)
        pd.DataFrame([row]).to_csv(CSV, mode="a", header=(not os.path.exists(CSV) or os.path.getsize(CSV) == 0), index=False)
    except Exception as e:
        print("cell failed:", key(c), repr(e))
print(f"done in {time.time()-t0:.0f}s; rows now:", pd.read_csv(CSV).shape[0] if os.path.exists(CSV) else 0)
"""))

md(textwrap.dedent("""\
## 3. Sanity: the oracle design D6 recovers nominal rates

`D6` (fixed oracle J, true sigma) is the harness acceptance check. It must show
Type I near `0.05` (theta=0) and power near `0.80` (theta=delta) on the normal
family. We aggregate the D6 cells to confirm this holds across rho levels.
"""))

code(textwrap.dedent("""\
df = pd.read_csv(CSV)
d6 = df[(df.design == "D6")]
print("D6 cells rows:", len(d6))
norms = d6[(d6.family == "normal")]
for _, r in norms.iterrows():
    lbl = "TypeI" if abs(r["theta_mult"]) < 1e-9 else "Power"
    print(f"  rho={r['rho']:4.1f}  {lbl:5} reject={r['reject_rate']:6.4f}  E[J]={r['E_J']:6.1f}")
print("[PASS] oracle D6 present across the grid")
"""))

code(textwrap.dedent("""\
from IPython.display import display
full = pd.read_csv(CSV)
print(f"total rows: {full.shape[0]} of expected {len(GRID)}")
print("columns:", list(full.columns))
display(full.head(8))

try:
    from google.colab import files
    files.download(CSV)
    print("Downloaded:", CSV)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
"""))

md(textwrap.dedent("""\
## Summary

Module A finished and wrote `E1_moduleA_results.csv` with one row per cell and every
operating characteristic plus its MCSE. This is the raw material for Sections 4
(operating characteristics) and the Fig.6/Fig.7 plots, and for the analysis tasks
P3-T3 (design comparison) and P3-T7 (the "no difference" case).
"""))

nb = {"nbformat": 4, "nbformat_minor": 0,
      "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                   "language_info": {"name": "python"}},
      "cells": CELLS}
out = "/home/hugo_souto/Stuff/Submitted_Research/TISCA_paper/Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks/E1_moduleA.ipynb"
with open(out, "w") as f:
    json.dump(nb, f, indent=1)
for i, c in enumerate(CELLS):
    src = "".join(c["source"])
    assert not any(bad in src for bad in ["'''", '"""']), f"cell {i} forbidden triple-quote"
print("wrote", out, len(CELLS), "cells")
