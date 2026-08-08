#!/usr/bin/env python3
"""Run the full E1 grid locally, in parallel, and write the P3-T2 results CSV.

The Colab runners in ``notebooks/`` execute the same grid one cell at a time in a
single session. On a multi-core workstation the whole study is a few minutes, so
this is the path to prefer: it removes the session-loss risk entirely and produces
``results/E1/operating_characteristics.csv`` in one shot.

Usage::

    python experiments/E1_operating_characteristics/run_e1_grid.py            # all modules
    python experiments/E1_operating_characteristics/run_e1_grid.py -m AD -j 8
    python experiments/E1_operating_characteristics/run_e1_grid.py --dry-run

Resumable: completed ``cell_id``s in the output file are skipped, so an interrupted
run continues where it stopped. Each worker is pinned to one BLAS thread, because
the engine is already vectorised and nested threading only causes contention.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                                  # noqa: E402
import pandas as pd                                                 # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "tisca", "python"))

from tisca import multiplicity                                      # noqa: E402
from tisca.outermc import e1_grid, engine, summarize_ocs            # noqa: E402

MATRIX_CSV = os.path.join(_ROOT, "legacy", "Paper_Experiments", "DGP1_500_results.csv")
MATRIX_COLS = ["mvbcf_pehe1", "bcf_pehe1"]
DEFAULT_OUT = os.path.join(_ROOT, "results", "E1", "operating_characteristics.csv")

_MATRIX = None          # per-worker cache; the (500, 2) matrix is cheap to reload


def load_matrix(path=MATRIX_CSV):
    """The pre-specified real loss pair driving the two empirical families."""
    raw = pd.read_csv(path)
    missing = [c for c in MATRIX_COLS if c not in raw.columns]
    if missing:
        raise ValueError(f"{path} lacks {missing}")
    mat = raw[MATRIX_COLS].to_numpy(dtype=float)
    mat = mat[np.all(np.isfinite(mat), axis=1)]
    if mat.shape[0] < 2:
        raise ValueError("empirical loss matrix has fewer than 2 usable rows")
    return mat


def _run_one(cell):
    """Execute one cell and return its operating-characteristic row."""
    global _MATRIX
    cfg = dict(cell["config"])
    if cfg.get("matrix") == "__MATRIX__":
        if _MATRIX is None:
            _MATRIX = load_matrix()
        cfg["matrix"] = _MATRIX
    t0 = time.time()
    summary, _, _ = engine.run_e1(cfg)
    row = summarize_ocs([summary]).iloc[0].to_dict()
    row.update(cell["factors"])
    row.update(cell_id=cell["cell_id"], module=cell["module"],
               projected_R=cfg["R"], bootstrap_B=cfg.get("B", np.nan),
               cell_seconds=round(time.time() - t0, 3))
    return row


def build_grid(modules):
    """Grid with the matrix replaced by a token, so cells stay picklable and small."""
    grid = e1_grid.make_grid(modules, matrix="__MATRIX__",
                             planning_alpha=multiplicity.planning_alpha)
    return grid


def resolve_oracle_sigma(grid, matrix, out_csv=None):
    """Fill each cell's ``sigma_D`` once, in the parent, and return the lookup table.

    ``engine.sigma_D_true`` estimates the non-normal families' sigma from a
    1,000,000-draw sample. Its cache lives inside one process, so leaving the call
    to the workers repeats every estimate in every worker. Resolving the ~70 unique
    (family, rho, sigma_a, sigma_b) keys up front removes that entirely, and makes
    the oracle an explicit, auditable input to the grid instead of a side effect:
    D6's whole claim is that it plans from the true sigma, so the value it used
    belongs in the results directory.
    """
    table, rows = {}, []
    for cell in grid:
        cfg = cell["config"]
        key = (cfg["family"], cfg["rho"], cfg["sigma_a"], cfg["sigma_b"])
        if key not in table:
            probe = dict(cfg)
            probe["matrix"] = matrix if cfg["matrix"] == "__MATRIX__" else cfg["matrix"]
            table[key] = engine.sigma_D_true(probe)
            rows.append({"family": key[0], "rho": key[1], "sigma_a": key[2],
                         "sigma_b": key[3], "sigma_D_true": table[key]})
        cfg["sigma_D"] = table[key]
    if out_csv:
        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        pd.DataFrame(rows).to_csv(out_csv, index=False)
    return table


def available_gb():
    """Currently available RAM in GB, from /proc/meminfo (MemAvailable)."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return float("inf")


def _choose_jobs(args):
    """Cap the worker count by memory, not by core count.

    Each worker holds an ``(R, Jmax, 2)`` block plus working copies: measured peak
    RSS is ~0.52 GB on the heaviest cells. Sizing the pool from cores alone is what
    put this machine into swap, and the user may have other experiments running, so
    the default takes at most ``--mem-fraction`` of *currently available* memory.
    """
    by_mem = max(1, int((available_gb() * args.mem_fraction) // args.mem_per_worker))
    by_cpu = max(1, (os.cpu_count() or 2) - 4)
    if args.jobs is not None:
        if args.jobs > by_mem:
            print(f"[E1] WARNING: -j {args.jobs} needs about "
                  f"{args.jobs * args.mem_per_worker:.1f} GB; only "
                  f"{available_gb() * args.mem_fraction:.1f} GB is budgeted. "
                  f"Consider -j {by_mem}.")
        return args.jobs
    jobs = min(by_cpu, by_mem)
    print(f"[E1] {available_gb():.1f} GB available -> {jobs} workers "
          f"(cpu cap {by_cpu}, memory cap {by_mem} at "
          f"{args.mem_per_worker} GB/worker)")
    return jobs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-m", "--modules", default="ABCD", help="subset of ABCD (default all)")
    p.add_argument("-j", "--jobs", type=int, default=None,
                   help="worker processes (default: memory-aware, see --mem-per-worker)")
    p.add_argument("--mem-per-worker", type=float, default=0.6,
                   help="GB budgeted per worker (measured peak RSS is ~0.52 GB on the "
                        "heaviest cells: R=5000 x Jmax=1000 blocks)")
    p.add_argument("--mem-fraction", type=float, default=0.5,
                   help="fraction of currently AVAILABLE RAM this run may use (default 0.5, "
                        "leaving room for other experiments on the same machine)")
    p.add_argument("-o", "--out", default=DEFAULT_OUT)
    p.add_argument("--dry-run", action="store_true", help="report the plan, run nothing")
    p.add_argument("--fresh", action="store_true", help="ignore an existing output file")
    args = p.parse_args(argv)
    args.jobs = _choose_jobs(args)

    grid = build_grid(args.modules.upper())
    counts = {m: sum(c["module"] == m for c in grid) for m in "ABCD" if m in args.modules.upper()}
    print(f"[E1] {len(grid)} cells {counts}, {args.jobs} workers")
    for m, n in counts.items():
        assert n == e1_grid.EXPECTED[m], f"module {m}: {n} cells, expected {e1_grid.EXPECTED[m]}"

    sig_csv = os.path.join(os.path.dirname(args.out), "oracle_sigma_D.csv")
    table = resolve_oracle_sigma(grid, load_matrix(), out_csv=sig_csv)
    print(f"[E1] resolved {len(table)} oracle sigma_D values -> {sig_csv}")

    done = set()
    if os.path.exists(args.out) and os.path.getsize(args.out) > 0 and not args.fresh:
        done = set(pd.read_csv(args.out)["cell_id"].astype(str))
        print(f"[E1] resuming: {len(done)} cells already in {args.out}")
    pending = [c for c in grid if c["cell_id"] not in done]

    if args.dry_run:
        print(f"[E1] would run {len(pending)} cells")
        return 0
    if not pending:
        print("[E1] nothing to do")
        return 0

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    started = time.time()
    rows, failures = [], []

    if args.jobs <= 1:
        for i, cell in enumerate(pending, 1):
            rows.append(_run_one(cell))
            if i % 100 == 0:
                print(f"  {i}/{len(pending)} cells, {time.time() - started:.0f}s", flush=True)
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(args.jobs) as pool:
            for i, row in enumerate(pool.imap_unordered(_run_one, pending, chunksize=4), 1):
                rows.append(row)
                if i % 100 == 0:
                    print(f"  {i}/{len(pending)} cells, {time.time() - started:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    if done:
        df = pd.concat([pd.read_csv(args.out), df], ignore_index=True)
    df = df.sort_values(["module", "cell_id"]).reset_index(drop=True)
    df.to_csv(args.out, index=False)

    # Completeness is asserted here rather than trusted: a silently short grid is
    # the failure mode that would survive into the manuscript.
    ids = set(df["cell_id"].astype(str))
    expected = {c["cell_id"] for c in grid}
    if ids != expected:
        raise RuntimeError(f"incomplete: {len(expected - ids)} missing, "
                           f"{len(ids - expected)} unexpected")
    if df["cell_id"].duplicated().any():
        raise RuntimeError("duplicate cell_id in output")
    if failures:
        raise RuntimeError(f"{len(failures)} cells failed")

    print(f"[E1] wrote {len(df)} rows to {args.out} in {time.time() - started:.0f}s")
    print(df.groupby("module").size().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
