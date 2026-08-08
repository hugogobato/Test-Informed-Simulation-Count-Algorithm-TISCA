#!/usr/bin/env python3
"""Generate ``E3_seed_verification.ipynb`` (P3-T5 acceptance; P0-T2 acceptance).

Closes the last two open boxes on the E3 case study:

* **P3-T5 acceptance** — "a re-run of 3 randomly chosen seeds reproducing the stored
  rows exactly". The confirmatory shards ran across 10 different hostnames, so this
  is a cross-machine reproducibility check as well as a within-machine one.
* **P0-T2 acceptance / CALIBRATION.md §3.5** — "running replication 37 twice, on
  different machines and different worker counts, yields bit-identical metrics": the
  four-way identity test over ``mc.cores`` in {1, 2} crossed with shard-aligned
  versus shard-offset seed ranges.

**CALIBRATION.md §3.6 — which of `stochtree::bcf`, `dbarts`, `bartCause` and
`fast_bart()` actually honour an *explicit* seed — is NOT answered here, and the
argument that it was is wrong.** The claim used to be that per-model column
prefixes settle it: if a re-run reproduces `bcf_*` exactly then `bcf` honoured its
seed. It does not follow. `run_cell.R` positions a fixed L'Ecuyer substream in the
*global* RNG state immediately before each fit, so a model that ignored its seed
argument entirely and drew from the global stream would reproduce its columns just
as exactly as one that honoured it. Both hypotheses predict identical output, so
the test has no power to separate them.

What the prefix comparison does establish is real and worth having: **run-level
reproducibility per model** — re-running a seed reproduces that model's columns,
which rules out wall-clock or PID seeding and uncontrolled nondeterminism inside
each fit. The columns are reported under that name.

The test that would isolate the seed argument is the one CALIBRATION.md actually
specifies: fit the same model twice on one replication with the **same explicit
seed** but from **two different global RNG states**. Only agreement there
attributes the reproducibility to the seed argument. That run does not exist, so
G2 stays open and the §3.6 box stays unticked.

Regenerate with::

    python notebooks/_generators/build_e3_seed_verification_nb.py \\
        --bundle-folder-url '<drive folder>' --bundle-sha256 '<sha>'
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

from _nbcommon import SETUP, md, code, write

_SPEC = importlib.util.spec_from_file_location(
    "_e3gen", Path(__file__).with_name("build_e3_notebooks.py"))
_E3 = importlib.util.module_from_spec(_SPEC)
sys.modules["_e3gen"] = _E3
_SPEC.loader.exec_module(_E3)


CONFIG = r'''
import glob
import os
import subprocess

import numpy as np
import pandas as pd

# --- what to verify -------------------------------------------------------- #
# Test A: exact reproduction of stored rows. n = 500 costs ~9.6 min/replication at
# one core, so three seeds is ~29 minutes -- which is exactly what the acceptance
# criterion asks for and no more.
TEST_A_DGP, TEST_A_N, TEST_A_SEEDS = 1, 500, None      # None -> drawn below
N_TEST_A_SEEDS = 3
TEST_A_RNG_SEED = 20260806                             # fixed: the choice is auditable

# Test B: the four-way identity test. Run at n = 100 (~2.7 min/replication) because
# it needs 4 configurations x 2 seeds = 8 replications, and the property under test
# (stream construction) does not depend on the training size.
TEST_B_DGP, TEST_B_N = 1, 100
TEST_B_SEEDS = [37, 38]        # contiguous, so a shard-offset run is well defined

WORK = "/content/e3_verify"
os.makedirs(WORK, exist_ok=True)


def find_stored(dgp, n, mode="confirmatory"):
    """Locate the committed shard CSVs for one cell."""
    pats = [
        os.path.join(REPO_ROOT, "results", "E3", f"DGP{dgp}_n{n}_{mode}_replications.csv"),
        os.path.join(REPO_ROOT, "notebooks", "E3_shards",
                     f"E3_DGP{dgp}_n{n}_{mode}_shard*.csv"),
        f"/content/E3_DGP{dgp}_n{n}_{mode}_shard*.csv",
        f"/content/*DGP{dgp}_n{n}*.csv",
    ]
    files = []
    for p in pats:
        files = sorted(glob.glob(p))
        if files:
            break
    if not files:
        raise FileNotFoundError(
            f"No stored rows for DGP{dgp} n={n}. Either commit the shard CSVs to "
            "notebooks/E3_shards/, or upload them to /content before running this "
            "notebook.")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    print(f"stored rows for DGP{dgp} n={n}: {len(df)} from {len(files)} file(s)")
    return df


STORED_A = find_stored(TEST_A_DGP, TEST_A_N)
rng = np.random.default_rng(TEST_A_RNG_SEED)
TEST_A_SEEDS = sorted(rng.choice(np.sort(STORED_A["seed"].unique()),
                                 size=N_TEST_A_SEEDS, replace=False).tolist())
print("seeds selected for exact re-run:", TEST_A_SEEDS)
print("hostnames that produced the stored rows:",
      sorted(STORED_A["hostname"].astype(str).unique()))
'''

RUNNER = r'''
def run_seeds(dgp, n, s_start, s_end, out_name, cores=1, mode="confirmatory"):
    """Invoke run_cell.R over one contiguous seed range and return the rows."""
    out = os.path.join(WORK, out_name)
    if os.path.exists(out):
        os.remove(out)
    cmd = ["Rscript", "/content/e3/run_cell.R", str(dgp), str(n),
           str(s_start), str(s_end), "--out", out,
           "--cores", str(cores), "--mode", mode]
    env = dict(os.environ, TISCA_MVBCF_CPP="/content/e3/MVBCF_Code.cpp")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=env)
    print(f"  {' '.join(cmd[2:])} -> {time.time() - t0:.0f}s")
    if p.returncode != 0:
        print(p.stdout[-3000:])
        print(p.stderr[-3000:])
        raise RuntimeError("run_cell.R failed")
    return pd.read_csv(out)


# Columns that MUST differ or are allowed to differ between runs: they record the
# execution environment, not the replication's result.
ENV_COLS = {"hostname", "git_sha", "session_hash", "replication_seconds",
            "fit_seconds_mvbcf", "fit_seconds_bcf1", "fit_seconds_bcf2",
            "fit_seconds_bart1", "fit_seconds_bart2", "fit_seconds_mvbart",
            "error_message"}


def compare_rows(a, b, label_a, label_b):
    """Exact comparison of every non-environment column, keyed on `seed`."""
    a = a.set_index("seed").sort_index()
    b = b.set_index("seed").sort_index()
    common = sorted(set(a.index) & set(b.index))
    assert common, f"no shared seeds between {label_a} and {label_b}"
    cols = [c for c in a.columns if c in b.columns and c not in ENV_COLS]
    rows = []
    for c in cols:
        x, y = a.loc[common, c], b.loc[common, c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            both_nan = x.isna() & y.isna()
            d = (x - y).abs().where(~both_nan, 0.0)
            identical = bool(((d == 0) | both_nan).all())
            maxdiff = float(d.max())
        else:
            identical = bool((x.astype(str) == y.astype(str)).all())
            maxdiff = np.nan
        rows.append({"column": c, "identical": identical, "max_abs_diff": maxdiff,
                     "model": c.split("_")[0] if "_" in c else "meta"})
    out = pd.DataFrame(rows)
    n_bad = int((~out["identical"]).sum())
    print(f"{label_a} vs {label_b}: {len(cols)} columns compared over "
          f"{len(common)} seed(s); {n_bad} differ")
    if n_bad:
        print(out[~out["identical"]].sort_values("max_abs_diff", ascending=False)
              .head(20).to_string(index=False))
    return out
'''

TEST_A = r'''
# --------------------------------------------------------------------------- #
# Test A: do three randomly chosen seeds reproduce their stored rows exactly?   #
# --------------------------------------------------------------------------- #
frames = []
for s in TEST_A_SEEDS:
    frames.append(run_seeds(TEST_A_DGP, TEST_A_N, s, s, f"A_seed{s}.csv", cores=1))
rerun_A = pd.concat(frames, ignore_index=True)
cmp_A = compare_rows(STORED_A[STORED_A["seed"].isin(TEST_A_SEEDS)], rerun_A,
                     "stored", "re-run")
'''

TEST_A_MODELS = r'''
# Per-model RUN-LEVEL REPRODUCIBILITY. Read this for what it is: a model whose
# columns come back identical was reproduced by a re-run of the same seed. That is
# NOT the same as CALIBRATION.md section 3.6's question, which asks whether the model
# honours its *explicit seed argument*. run_cell.R sets a fixed L'Ecuyer substream in
# the global RNG state right before each fit, so a model that ignored its seed
# argument and drew from the global stream would reproduce these columns just as
# exactly. The two explanations are observationally identical here, so section 3.6
# is NOT answered by this table -- see the closing note.
per_model = (cmp_A[cmp_A["model"].isin(["mvbcf", "bcf", "bart", "mvbart"])]
             .groupby("model")
             .agg(columns=("column", "size"),
                  n_differing=("identical", lambda v: int((~v).sum())),
                  max_abs_diff=("max_abs_diff", "max")))
per_model["reproduces_under_rerun"] = per_model["n_differing"] == 0
print(per_model.to_string())
print()
meta = cmp_A[~cmp_A["model"].isin(["mvbcf", "bcf", "bart", "mvbart"])]
print("non-model columns differing:",
      meta.loc[~meta["identical"], "column"].tolist() or "none")
'''

TEST_B = r'''
# --------------------------------------------------------------------------- #
# Test B: the four-way identity test (P0-T2 acceptance, CALIBRATION.md 3.5).     #
# --------------------------------------------------------------------------- #
# mc.cores in {1, 2} x {shard-aligned, shard-offset}. "Shard-offset" means the same
# seeds are requested as part of a range that starts earlier, which is exactly what
# happens when a cell is split into shards; the stream for index j must not depend
# on where its shard began.
lo, hi = min(TEST_B_SEEDS), max(TEST_B_SEEDS)
variants = {}
variants["cores1_aligned"] = run_seeds(TEST_B_DGP, TEST_B_N, lo, hi,
                                       "B_c1_aligned.csv", cores=1)
variants["cores2_aligned"] = run_seeds(TEST_B_DGP, TEST_B_N, lo, hi,
                                       "B_c2_aligned.csv", cores=2)
off = max(0, lo - 2)
variants["cores1_offset"] = run_seeds(TEST_B_DGP, TEST_B_N, off, hi,
                                      "B_c1_offset.csv", cores=1)
variants["cores2_offset"] = run_seeds(TEST_B_DGP, TEST_B_N, off, hi,
                                      "B_c2_offset.csv", cores=2)

ref = "cores1_aligned"
cmp_B = {}
for k, v in variants.items():
    if k == ref:
        continue
    cmp_B[k] = compare_rows(variants[ref], v, ref, k)
'''

VERDICT = r'''
# --------------------------------------------------------------------------- #
# Verdict                                                                       #
# --------------------------------------------------------------------------- #
checks = []
checks.append({
    "check": "P3-T5: 3 random seeds reproduce their stored rows exactly",
    "detail": f"seeds {TEST_A_SEEDS}, {len(cmp_A)} columns each; "
              f"{int((~cmp_A['identical']).sum())} differing",
    "pass": bool(cmp_A["identical"].all()),
})
checks.append({
    "check": "cross-machine: the stored rows came from a different host",
    "detail": f"stored hosts {sorted(STORED_A['hostname'].astype(str).unique())} "
              f"vs this session",
    "pass": True,       # informational; the exactness check above is the real test
})
for k, c in cmp_B.items():
    checks.append({
        "check": f"P0-T2 four-way identity: {ref} == {k}",
        "detail": f"{int((~c['identical']).sum())} of {len(c)} columns differ",
        "pass": bool(c["identical"].all()),
    })
for m, r in per_model.iterrows():
    checks.append({
        "check": f"run-level reproducibility of {m}",
        "detail": f"{int(r['n_differing'])} of {int(r['columns'])} columns differ, "
                  f"max |diff| = {r['max_abs_diff']}",
        "pass": bool(r["reproduces_under_rerun"]),
    })
report = pd.DataFrame(checks)
print(report.to_string(index=False))
print()
print("ALL PASS" if report["pass"].all() else "*** FAILURES ABOVE -- do not tick "
      "the CALIBRATION.md boxes ***")

os.makedirs(os.path.join(RESULTS, "E3"), exist_ok=True)
report.to_csv(os.path.join(RESULTS, "E3", "seed_verification.csv"), index=False)
cmp_A.to_csv(os.path.join(RESULTS, "E3", "seed_verification_columns.csv"), index=False)
download(os.path.join(RESULTS, "E3", "seed_verification.csv"))

print()
print("--- paste into experiments/E3_mvbcf_casestudy/CALIBRATION.md ---")
print(f"""
### Seed verification, {pd.Timestamp.today().date()}

- [{'x' if report['pass'].all() else ' '}] P3-T5 acceptance: seeds {TEST_A_SEEDS} of
      DGP{TEST_A_DGP} n={TEST_A_N} re-run and compared against the stored rows on
      {len(cmp_A)} non-environment columns. Differing columns:
      {int((~cmp_A['identical']).sum())}.
- [{'x' if all(c['identical'].all() for c in cmp_B.values()) else ' '}] P0-T2 §3.5
      four-way identity test (mc.cores 1 vs 2, shard-aligned vs shard-offset) on
      DGP{TEST_B_DGP} n={TEST_B_N}, seeds {TEST_B_SEEDS}.
- [ ] P0-T2 §3.6 per-model seed HONOURING: **still open**. What this notebook
      establishes is run-level reproducibility per model:
      {dict(per_model['reproduces_under_rerun'])}.
      That does not isolate the explicit seed argument, because run_cell.R
      deterministically positions the global RNG stream before every fit. The
      discriminating run -- same model, same explicit seed, two different global
      RNG states -- has not been done.
""")
'''


def build(bundle_source, bundle_sha):
    cells = []
    md(cells, """
        # E3 seed verification (P3-T5 and P0-T2 acceptance)

        Two acceptance criteria are still open on the case study, and both are
        reproducibility claims that were asserted but never tested:

        1. **P3-T5** requires that three randomly chosen seeds be re-run and
           reproduce their stored rows *exactly*. The confirmatory shards ran across
           ten different Colab hosts, so this doubles as a cross-machine check.
        2. **P0-T2** requires that a replication give bit-identical metrics under
           different worker counts and different shard boundaries.

        `CALIBRATION.md` §3.6 also asks which model implementations honour an
        explicit seed. **This notebook does not answer that**, and an earlier
        version wrongly claimed it did. Re-running a seed reproduces a model's
        columns whether the model honoured its seed argument or ignored it and drew
        from the global stream, which `run_cell.R` positions deterministically
        before every fit. What is reported per model is therefore **run-level
        reproducibility**, and the §3.6 box stays unticked.

        Runtime is about 50 minutes: 3 replications at n = 500 for Test A and 8 at
        n = 100 for Test B. One Colab session, comfortably inside the cap.

        **A failure here is important and must not be worked around.** If a model
        does not honour its seed, the case-study data is not reproducible and the
        Data Availability statement cannot claim that it is.
        """)
    # SETUP first: it defines REPO_ROOT, RESULTS, FIGURES and download(), which
    # every later cell uses. The E3 shared cells (hardware report, R install, bundle
    # restore, driver download, fast_bart compile) follow; its leading markdown cell
    # is dropped because it describes a shard run.
    code(cells, SETUP)
    row = {"round": "verification", "dgp": 1, "n": 500, "mode": "confirmatory",
           "seed_start": 0, "seed_end": 0}
    for c in _E3._shared_cells(bundle_source, bundle_sha, row, None)[1:]:
        cells.append(c)
    code(cells, "import time\n" + CONFIG)
    md(cells, "## Runner and comparison helpers")
    code(cells, RUNNER)
    md(cells, """
        ## Test A: exact reproduction of stored rows

        Environment columns (`hostname`, `git_sha`, `session_hash`, and every timing
        column) are excluded from the comparison because they record where the run
        happened, not what it produced. Everything else must match exactly, not
        approximately: a tolerance here would hide precisely the drift the check
        exists to detect.
        """)
    code(cells, TEST_A)
    md(cells, "### Per-model verdict (CALIBRATION.md §3.6)")
    code(cells, TEST_A_MODELS)
    md(cells, "## Test B: the four-way identity test")
    code(cells, TEST_B)
    md(cells, "## Verdict, and the block to paste into CALIBRATION.md")
    code(cells, VERDICT)
    write("E3_seed_verification.ipynb", cells)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    _E3._bundle_args(ap)
    args = ap.parse_args()
    source, sha = _E3._validated_bundle(args)
    build(source, sha)


if __name__ == "__main__":
    main()
