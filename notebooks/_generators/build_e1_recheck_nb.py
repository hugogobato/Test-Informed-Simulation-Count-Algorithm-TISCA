#!/usr/bin/env python3
"""Generate the sharded Colab runners that re-run the E1 grid after the D3 interval fix.

Why this exists
---------------
``designs.design_v1_welch`` (D3) and ``joint_engine._welch_stats`` emitted their
reporting scale as ``sqrt((s2A + s2B) / 2)``, the average marginal sd, while both
engines build an interval as ``theta_hat +/- t * s / sqrt(J)``. D3's own Welch test
uses ``SE = sqrt(s2A/J + s2B/J) = sqrt(s2A + s2B)/sqrt(J)``, so the interval was
narrower than the design's own test by exactly ``sqrt(2)``. The published grid shows
the contradiction directly: at ``rho = 0``, where an unpaired analysis of paired rows
is exactly valid, D3's measured level is a correct 0.0498 while its coverage is 0.837.

Only D3's ``ci_cover`` is affected. Levels, power and ``E[J]`` are computed from
p-values and ``J``, which the fix does not touch, so a re-run doubles as a regression
check: every non-coverage column must come back identical.

What the notebooks do
---------------------
Each shard clones the public repo, applies the two patches with assertions (idempotent
if the fix is already upstream), runs a numerical self-check that reproduces both the
artefact and the repair before spending any compute, then executes its slice of the
canonical 1,983-cell grid with a process pool sized by available RAM.

Cells are assigned to shards by ``index % n_shards``, so the expensive Module-A cells
are spread evenly rather than landing in one session.

Merge the downloaded CSVs with::

    python experiments/E1_operating_characteristics/merge_e1_shards.py ~/Downloads/E1_recheck_shard*.csv

Regenerate the notebooks with::

    python notebooks/_generators/build_e1_recheck_nb.py --shards 4
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "tisca" / "python"))

from tisca.outermc import e1_grid  # noqa: E402  (needs the path above)

TOTAL_CELLS = sum(e1_grid.EXPECTED[m] for m in "ABCD")
DEFAULT_TIMINGS = ROOT.parent / "results" / "E1" / "operating_characteristics.csv"


def _assign_shards(cell_ids, n_shards, timings_path):
    """Split the grid into shards of roughly equal RUNNING TIME, not equal cell count.

    Cell cost is very uneven: Module A averages ~13 s and Module B ~1 s, and the
    slowest single cell (a beta-family D3 cell at R=5,000, Jmax=1,000) takes 94 s
    against a median of 3.5 s. Plain round-robin therefore leaves the slowest shard
    1.7x the fastest at 32 shards, and every session waits for the slowest one.

    The previous run recorded ``cell_seconds`` per cell, so the split uses longest-
    processing-time-first: sort by measured cost descending and give each cell to the
    least-loaded shard. Without a timings file this degrades to round-robin over a
    cost-sorted order, which is still better than round-robin over grid order.
    """
    costs = {}
    if timings_path and os.path.exists(timings_path):
        frame = pd.read_csv(timings_path, usecols=["cell_id", "cell_seconds"])
        # 95 of the recorded cell_seconds are negative: the original parallel run was
        # timed with a wall clock that skewed under load on WSL. They are provenance
        # metadata only and affect no published quantity, but as a cost signal they are
        # unusable, so non-positive entries are treated as unknown and get the median.
        frame = frame[frame["cell_seconds"] > 0]
        costs = dict(zip(frame["cell_id"].astype(str), frame["cell_seconds"]))
    default = float(np.median(list(costs.values()))) if costs else 1.0
    ordered = sorted(cell_ids, key=lambda c: -costs.get(c, default))

    shards = [[] for _ in range(n_shards)]
    loads = [0.0] * n_shards
    for cell_id in ordered:
        k = min(range(n_shards), key=lambda i: loads[i])
        shards[k].append(cell_id)
        loads[k] += costs.get(cell_id, default)
    return shards, loads, bool(costs)


def _md(cells: list[dict], source: str) -> None:
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": textwrap.dedent(source).splitlines(keepends=True),
    })


def _code(cells: list[dict], source: str) -> None:
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": textwrap.dedent(source).splitlines(keepends=True),
    })


# The patches are expressed as exact old -> new source substitutions rather than as a
# diff, so that a silent upstream edit makes the notebook fail loudly instead of
# applying a hunk to a file that has moved underneath it.
SETUP = r'''
import multiprocessing as mp
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

REPO_URL = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
CLONED_REPO = "/content/TISCA_repo"

if not os.path.isdir(os.path.join(CLONED_REPO, "tisca", "python")):
    subprocess.run(["git", "clone", "--depth", "1", REPO_URL, CLONED_REPO], check=True)
REPO_ROOT = CLONED_REPO
SOURCE_ROOT = os.path.join(REPO_ROOT, "tisca", "python")
assert os.path.isdir(SOURCE_ROOT), f"TISCA package not found at {SOURCE_ROOT}"

# ---- the D3 reporting-scale fix -------------------------------------------------
# Applied to the checkout before the package is imported. Each patch asserts that it
# found either the old text (and rewrites it) or the new text (already upstream); a
# file that matches neither stops the run rather than producing numbers of unknown
# provenance.
PATCHES = [
    (os.path.join(SOURCE_ROOT, "tisca", "outermc", "designs.py"),
     "                 np.sqrt((s2A + s2B) / 2), capped), dict(",
     "                 np.sqrt(s2A + s2B), capped), dict("),
    (os.path.join(SOURCE_ROOT, "tisca", "outermc", "joint_engine.py"),
     "    report_sd = np.sqrt((s2a[:, None] + s2b) / 2.0)",
     "    report_sd = np.sqrt(s2a[:, None] + s2b)"),
]

for path, old, new in PATCHES:
    with open(path) as fh:
        text = fh.read()
    if old in text:
        with open(path, "w") as fh:
            fh.write(text.replace(old, new))
        print("[PATCHED]", os.path.basename(path))
    elif new in text:
        print("[ALREADY FIXED]", os.path.basename(path))
    else:
        raise RuntimeError(
            f"{path}: neither the old nor the fixed reporting scale was found. "
            "The upstream file has changed; regenerate this notebook."
        )

sys.path.insert(0, SOURCE_ROOT)

OUTPUT_ROOT = "/content/TISCA_E1_recheck"
os.makedirs(OUTPUT_ROOT, exist_ok=True)

# ---- the worker, written to a real module -----------------------------------------
# multiprocessing.Pool sends tasks by pickling the callable by qualified name, so a
# function defined in a notebook cell cannot be a worker even under fork: the child
# fails with "attribute lookup _run_one on __main__ failed". Writing the worker to an
# importable module is what makes the pool usable from a notebook at all.
WORKER_SOURCE = r"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.environ["TISCA_SOURCE_ROOT"])
from tisca.outermc import engine, summarize_ocs

_MATRIX = None


# The pre-specified real loss pair driving the two empirical families.
def load_empirical_matrix():
    path = os.path.join(os.environ["TISCA_REPO_ROOT"], "legacy", "Paper_Experiments",
                        "DGP1_500_results.csv")
    raw = pd.read_csv(path)
    matrix = raw[["mvbcf_pehe1", "bcf_pehe1"]].to_numpy(dtype=float)
    matrix = matrix[np.all(np.isfinite(matrix), axis=1)]
    assert matrix.ndim == 2 and matrix.shape[1] == 2 and matrix.shape[0] >= 2
    return matrix


# Execute one grid cell and return its operating-characteristic row.
def run_one(cell):
    global _MATRIX
    cfg = dict(cell["config"])
    if cfg.get("matrix") == "__MATRIX__":
        if _MATRIX is None:
            _MATRIX = load_empirical_matrix()
        cfg["matrix"] = _MATRIX
    t0 = time.time()
    summary, _, _ = engine.run_e1(cfg)
    row = summarize_ocs([summary]).iloc[0].to_dict()
    row.update(cell["factors"])
    row.update(cell_id=cell["cell_id"], module=cell["module"],
               projected_R=cfg["R"], bootstrap_B=cfg.get("B", np.nan),
               cell_seconds=round(time.time() - t0, 3))
    return row
"""

os.environ["TISCA_SOURCE_ROOT"] = SOURCE_ROOT
os.environ["TISCA_REPO_ROOT"] = REPO_ROOT
WORKER_PATH = os.path.join(OUTPUT_ROOT, "e1_worker.py")
with open(WORKER_PATH, "w") as fh:
    fh.write(WORKER_SOURCE)
sys.path.insert(0, OUTPUT_ROOT)

import e1_worker

from tisca.outermc import e1_grid, engine, summarize_ocs
from tisca import multiplicity

print("cpu_count:", os.cpu_count())
print("outputs:", OUTPUT_ROOT)
'''


SELFCHECK = r'''
# Reproduce the defect and its repair before spending any compute. With equal group
# sizes the two interval scales differ by exactly sqrt(2), so at rho = 0 -- where an
# unpaired analysis of paired rows is valid and the Welch test is correctly sized --
# the as-coded interval must land near 2*Phi(1.96/sqrt(2)) - 1 = 0.834 while the
# consistent interval must land near 0.95. If that separation does not appear, the
# patch above did not do what this notebook claims it does.
from scipy import stats

rng = np.random.default_rng(0)
R_CHK, J_CHK = 200_000, 60
crit = stats.t.ppf(0.975, J_CHK - 1)
print(f"{'rho':>5} {'as-coded':>10} {'consistent':>12}")
checked = {}
for rho in (-0.3, 0.0, 0.6):
    cov = np.array([[1.0, rho], [rho, 1.0]])
    x = rng.multivariate_normal([0.0, 0.0], cov, size=(R_CHK, J_CHK))
    A, B = x[..., 0], x[..., 1]
    s2A, s2B = A.var(1, ddof=1), B.var(1, ddof=1)
    dbar = A.mean(1) - B.mean(1)
    as_coded = float(np.mean(np.abs(dbar) <= crit * np.sqrt((s2A + s2B) / 2) / np.sqrt(J_CHK)))
    consistent = float(np.mean(np.abs(dbar) <= crit * np.sqrt((s2A + s2B) / J_CHK)))
    checked[rho] = (as_coded, consistent)
    print(f"{rho:>5.1f} {as_coded:>10.4f} {consistent:>12.4f}")

assert abs(checked[0.0][0] - 0.834) < 0.01, checked[0.0]
assert abs(checked[0.0][1] - 0.950) < 0.01, checked[0.0]
assert checked[0.6][1] > 0.99, checked[0.6]     # over-covers at positive rho
assert checked[-0.3][1] < 0.93, checked[-0.3]   # under-covers at negative rho
print("[PASS] the sqrt(2) artefact and its repair both reproduce")
'''


GRID = r'''
MATRIX = e1_worker.load_empirical_matrix()

# The grid is imported, never retyped: e1_grid.py is the single source of truth for
# all 1,983 cells and the local runner reads the same object.
GRID_ALL = e1_grid.make_grid("ABCD", matrix="__MATRIX__",
                             planning_alpha=multiplicity.planning_alpha)
for _m, _n in e1_grid.EXPECTED.items():
    _got = sum(c["module"] == _m for c in GRID_ALL)
    assert _got == _n, f"module {_m}: {_got} cells, expected {_n}"


def resolve_oracle_sigma(grid, matrix):
    """Fill each cell's true sigma_D once, in the parent.

    engine.sigma_D_true estimates the non-normal families' sigma from a 1,000,000-draw
    sample and caches inside one process, so leaving the call to the workers repeats
    every estimate in every worker.
    """
    table = {}
    for cell in grid:
        cfg = cell["config"]
        key = (cfg["family"], cfg["rho"], cfg["sigma_a"], cfg["sigma_b"])
        if key not in table:
            probe = dict(cfg)
            probe["matrix"] = matrix if cfg["matrix"] == "__MATRIX__" else cfg["matrix"]
            table[key] = engine.sigma_D_true(probe)
        cfg["sigma_D"] = table[key]
    return table


_table = resolve_oracle_sigma(GRID_ALL, MATRIX)
print(f"resolved {len(_table)} oracle sigma_D values")

# SHARD_CELL_IDS is the explicit, time-balanced slice this session owns (see
# _assign_shards in the generator). The grid itself is still imported rather than
# retyped; only the partition is listed here, so a missing or duplicated cell is
# caught by the assertion below and again by the merge script across all shards.
_by_id = {c["cell_id"]: c for c in GRID_ALL}
missing = sorted(set(SHARD_CELL_IDS) - set(_by_id))
assert not missing, f"shard names {len(missing)} cells absent from the grid: {missing[:3]}"
SHARD = [_by_id[cid] for cid in SHARD_CELL_IDS]
print(f"shard {SHARD_INDEX + 1}/{N_SHARDS}: {len(SHARD)} of {len(GRID_ALL)} cells")
print("by module:", {m: sum(c["module"] == m for c in SHARD) for m in "ABCD"})
print(f"budgeted work: {EXPECTED_CPU_SECONDS:.0f} CPU-s "
      f"(<= ~{EXPECTED_CPU_SECONDS / 60:.0f} min at 1 core; an upper bound, the "
      f"source timings were recorded under load)")
'''


RUN = r'''
def available_gb():
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return float("inf")


# Each worker holds an (R, Jmax, 2) block plus working copies: measured peak RSS is
# ~0.52 GB on the heaviest cells, so the pool is capped by memory as well as cores.
# A 40-core VM that budgets by cores alone will swap.
MEM_PER_WORKER = 0.6
by_mem = max(1, int((available_gb() * 0.7) // MEM_PER_WORKER))
JOBS = max(1, min(os.cpu_count() or 2, by_mem))
print(f"{available_gb():.1f} GB available -> {JOBS} workers "
      f"(cpu {os.cpu_count()}, memory cap {by_mem})")

OUTPUT_FILE = os.path.join(OUTPUT_ROOT, OUTPUT_NAME)

rows = []
done = set()
if os.path.exists(OUTPUT_FILE) and os.path.getsize(OUTPUT_FILE) > 0:
    previous = pd.read_csv(OUTPUT_FILE)
    rows = previous.to_dict("records")
    done = set(previous["cell_id"].astype(str))
    print(f"resuming: {len(done)} cells already checkpointed")
pending = [c for c in SHARD if c["cell_id"] not in done]
print(f"{len(pending)} cells pending")

started = time.time()
# The checkpoint rewrites the whole frame rather than appending one row at a time.
# Module-B rows carry factor columns (K, correction) that Module-A rows do not, so
# appending each row in its own column order silently misaligns the CSV: values land
# under the wrong headers and cell_id stops matching. Building one DataFrame unifies
# the columns, and a few hundred rows are cheap enough to rewrite per cell.
if pending:
    with mp.get_context("fork").Pool(JOBS) as pool:
        for row in tqdm(pool.imap_unordered(e1_worker.run_one, pending, chunksize=1),
                        total=len(pending), unit="cell"):
            rows.append(row)
            pd.DataFrame(rows).to_csv(OUTPUT_FILE, index=False)
print(f"elapsed {time.time() - started:.0f}s")
'''


VERIFY = r'''
RESULTS = pd.read_csv(OUTPUT_FILE)
expected = {c["cell_id"] for c in SHARD}
got = set(RESULTS["cell_id"].astype(str))
missing = sorted(expected - got)
assert not missing, f"{len(missing)} cells missing, first={missing[:3]}"
assert not RESULTS["cell_id"].duplicated().any(), "duplicate cell_id in checkpoint"
print(f"[PASS] {len(RESULTS)} rows, all {len(expected)} shard cells present")

# The fix moves D3 coverage and nothing else, so D3 null cells are the one place the
# change should be visible in this shard.
d3_null = RESULTS[(RESULTS["design"] == "D3") & (RESULTS["theta"] == 0)]
if len(d3_null):
    print("\nD3 null cells, coverage by rho (was 0.79-0.99 rising in rho, "
          "artefactually centred near 0.90):")
    print(d3_null.groupby("rho")[["ci_cover", "reject_rate"]].mean().round(4))

try:
    from google.colab import files
    files.download(OUTPUT_FILE)
    print("Downloaded:", OUTPUT_FILE)
except Exception as e:
    print("(Not on Colab / download skipped):", e)
'''


def _nb(cells: list[dict]) -> dict:
    return {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {
            "kernelspec": {"name": "python3", "display_name": "Python 3"},
            "language_info": {"name": "python"},
        },
        "cells": cells,
    }


def _notebook(shard_index: int, n_shards: int, cell_ids: list[str],
              cpu_seconds: float) -> dict:
    name = f"E1_recheck_shard{shard_index + 1}_of{n_shards}_results.csv"
    approx = len(cell_ids)
    cells: list[dict] = []
    _md(cells, f"""
        # E1 re-run after the D3 interval fix, shard {shard_index + 1} of {n_shards}

        Re-runs roughly {approx:,} of the {TOTAL_CELLS:,} canonical E1 cells with the
        D3 reporting scale corrected.

        **What was wrong.** Both engines build an interval as
        `theta_hat +/- t * s / sqrt(J)`, but the unpaired design D3 emitted `s` as the
        average marginal sd `sqrt((sA^2 + sB^2)/2)` while its own Welch test uses
        `SE = sqrt(sA^2 + sB^2) / sqrt(J)`. The interval was therefore narrower than
        the design's own test by exactly `sqrt(2)`. In the published grid D3 has a
        correct level of 0.0498 at `rho = 0` alongside a coverage of 0.837, which is
        the signature of that error rather than a property of unpaired analysis.

        **What changes.** Only D3's `ci_cover`. Levels, power and `E[J]` come from
        p-values and `J`, which the fix does not touch, so those columns must come back
        identical to the previous run; the merge script checks exactly that.

        Run all cells, then download `{name}` when prompted. Shards are independent and
        can run in parallel sessions; the run cell also parallelises across the cores of
        this VM, sized by available RAM.

        This shard carries about {cpu_seconds:.0f} CPU-seconds of work by the previous
        run's timings, so expect at most ~{cpu_seconds / 60:.0f} min on one core and
        ~{cpu_seconds / 120:.0f} min on a typical 2-core Colab VM. Both are upper bounds:
        those timings were recorded under a fully loaded local machine, and a measured
        re-run of one shard came in at roughly half the predicted cost.
    """)
    _code(cells, SETUP)
    _md(cells, """
        ## Self-check: reproduce the defect and its repair

        Before any grid compute, confirm numerically that the two interval scales differ
        the way the diagnosis says they do. This is a closed-form prediction, so a
        failure here means the diagnosis is wrong, not that the run was unlucky.
    """)
    _code(cells, SELFCHECK)
    _md(cells, """
        ## The canonical grid and this shard's slice

        Cells are imported from `tisca/python/tisca/outermc/e1_grid.py`, the single
        source of truth for all 1,983 cells. Only the partition is listed below, and it
        is balanced by measured running time rather than by cell count: Module A averages
        ~13 s per cell against Module B's ~1 s, so equal-sized shards would not be
        equal-length shards.
    """)
    ids = ",\n    ".join(repr(c) for c in cell_ids)
    _code(cells, f"SHARD_INDEX = {shard_index}\nN_SHARDS = {n_shards}\n"
                 f'OUTPUT_NAME = "{name}"\n'
                 f"EXPECTED_CPU_SECONDS = {cpu_seconds:.1f}\n"
                 f"SHARD_CELL_IDS = [\n    {ids},\n]\n"
                 f"assert len(SHARD_CELL_IDS) == len(set(SHARD_CELL_IDS))\n")
    _code(cells, GRID)
    _md(cells, """
        ## Run the shard

        One checkpointed row per completed cell, so an interrupted session loses at most
        the cells in flight; re-running this cell resumes.
    """)
    _code(cells, RUN)
    _md(cells, "## Completeness check and download")
    _code(cells, VERIFY)
    return _nb(cells)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shards", type=int, default=4,
                    help="number of Colab sessions to split the grid across")
    ap.add_argument("--timings", default=str(DEFAULT_TIMINGS),
                    help="CSV with cell_id,cell_seconds used to balance the shards")
    args = ap.parse_args(argv)
    if args.shards < 1:
        raise SystemExit("--shards must be at least 1")

    grid = e1_grid.make_grid("ABCD", matrix="__MATRIX__")
    cell_ids = [c["cell_id"] for c in grid]
    assert len(cell_ids) == TOTAL_CELLS == len(set(cell_ids))

    shards, loads, measured = _assign_shards(cell_ids, args.shards, args.timings)
    covered = [cid for shard in shards for cid in shard]
    assert sorted(covered) == sorted(cell_ids), "shards do not partition the grid"

    source = "measured cell_seconds" if measured else "uniform cost (no timings file)"
    print(f"{args.shards} shards balanced by {source}: "
          f"{min(len(s) for s in shards)}-{max(len(s) for s in shards)} cells each, "
          f"{min(loads):.0f}-{max(loads):.0f} CPU-s "
          f"(imbalance {max(loads) / max(min(loads), 1e-9):.2f}x)")

    for index in range(args.shards):
        notebook = _notebook(index, args.shards, shards[index], loads[index])
        name = f"E1_recheck_shard{index + 1}_of{args.shards}.ipynb"
        out = ROOT / name
        for cell_index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] == "code":
                compile("".join(cell["source"]), f"{name}:cell-{cell_index}", "exec")
        out.write_text(json.dumps(notebook, indent=1) + "\n")
        print(f"wrote {out.name}: {len(shards[index])} cells, {loads[index]:.0f} CPU-s")


if __name__ == "__main__":
    main()
