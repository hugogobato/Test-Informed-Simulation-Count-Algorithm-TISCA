#!/usr/bin/env python3
"""Merge the Colab E1 re-check shards and verify the D3 interval fix changed only D3 coverage.

The re-run exists to correct one quantity: D3's confidence-interval coverage, which the
old reporting scale understated by a factor of sqrt(2) (see
``notebooks/_generators/build_e1_recheck_nb.py``). Every other operating characteristic
is computed from p-values and the selected count ``J``, which the fix does not touch, so
a correct re-run must reproduce them exactly. This script asserts that rather than
assuming it: if a level, a power, or an ``E[J]`` moved, something other than the
intended fix changed and the merge is refused.

Usage::

    python experiments/E1_operating_characteristics/merge_e1_shards.py ~/Downloads/E1_recheck_shard*.csv
    python experiments/E1_operating_characteristics/merge_e1_shards.py shards/*.csv --write

Without ``--write`` the script reports and writes nothing.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "tisca", "python"))

from tisca.outermc import e1_grid  # noqa: E402

DEFAULT_OUT = os.path.join(_ROOT, "results", "E1", "operating_characteristics.csv")

#: Columns the fix must not move. ci_cover is deliberately absent: it is the target.
INVARIANT = ["reject_rate", "t1e_or_power", "E_theta", "bias", "rmse",
             "E_J", "sd_J", "q05_J", "q50_J", "q95_J", "pJmax"]
TOL = 1e-9


def load_shards(paths):
    frames = []
    for path in paths:
        frame = pd.read_csv(path)
        if "cell_id" not in frame.columns:
            raise SystemExit(f"{path}: no cell_id column")
        frame["_source"] = os.path.basename(path)
        frames.append(frame)
        print(f"  {os.path.basename(path):<45} {len(frame):>5} rows")
    return pd.concat(frames, ignore_index=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("shards", nargs="+", help="shard CSVs downloaded from Colab")
    ap.add_argument("-o", "--out", default=DEFAULT_OUT)
    ap.add_argument("--baseline", default=None,
                    help="pre-fix CSV to check invariants against (default: --out if it exists)")
    ap.add_argument("--write", action="store_true", help="write the merged CSV")
    args = ap.parse_args(argv)

    print(f"reading {len(args.shards)} shard files")
    merged = load_shards(args.shards)

    dupes = merged["cell_id"][merged["cell_id"].duplicated()].unique()
    if len(dupes):
        raise SystemExit(f"{len(dupes)} duplicate cell_id across shards, "
                         f"first={list(dupes[:3])}. Overlapping shards?")

    expected = {c["cell_id"] for c in e1_grid.make_grid("ABCD", matrix="__MATRIX__")}
    got = set(merged["cell_id"].astype(str))
    missing = sorted(expected - got)
    extra = sorted(got - expected)
    if missing:
        raise SystemExit(f"{len(missing)} cells missing, first={missing[:3]}")
    if extra:
        raise SystemExit(f"{len(extra)} unexpected cells, first={extra[:3]}")
    print(f"[PASS] all {len(expected)} cells present, no duplicates")

    merged = merged.drop(columns=["_source"]).sort_values(
        ["module", "cell_id"]).reset_index(drop=True)

    baseline_path = args.baseline or (args.out if os.path.exists(args.out) else None)
    if baseline_path and os.path.exists(baseline_path):
        old = pd.read_csv(baseline_path).set_index("cell_id")
        new = merged.set_index("cell_id")
        shared = old.index.intersection(new.index)
        print(f"\ncomparing {len(shared)} cells against {os.path.basename(baseline_path)}")

        moved = []
        for col in INVARIANT:
            if col not in old.columns or col not in new.columns:
                continue
            a = old.loc[shared, col].to_numpy(float)
            b = new.loc[shared, col].to_numpy(float)
            both_nan = np.isnan(a) & np.isnan(b)
            delta = np.where(both_nan, 0.0, np.abs(a - b))
            n_moved = int(np.count_nonzero(delta > TOL))
            if n_moved:
                moved.append((col, n_moved, float(np.nanmax(delta))))
        if moved:
            print("\n[FAIL] columns the fix should not touch have moved:")
            for col, n, worst in moved:
                print(f"  {col:<14} {n:>5} cells changed, max |delta| = {worst:.3e}")
            raise SystemExit(
                "Refusing to merge. The re-run differs from the baseline in quantities "
                "the D3 reporting-scale fix cannot affect, so the shards were produced "
                "by different code or a different seed than the baseline."
            )
        print("[PASS] every invariant column reproduces the baseline exactly")

        cov_delta = (new.loc[shared, "ci_cover"].to_numpy(float)
                     - old.loc[shared, "ci_cover"].to_numpy(float))
        changed = np.abs(cov_delta) > TOL
        design = old.loc[shared, "design"].to_numpy()
        non_d3 = changed & (design != "D3")
        if non_d3.any():
            raise SystemExit(f"[FAIL] coverage moved on {int(non_d3.sum())} non-D3 cells; "
                             "the fix should touch D3 only")
        print(f"[PASS] coverage changed on {int(changed.sum())} cells, all D3")

        d3 = new[(new["design"] == "D3") & (new["theta"] == 0)]
        d3_old = old[(old["design"] == "D3") & (old["theta"] == 0)]
        print("\nD3 null-cell coverage by rho (old -> new):")
        a = d3_old.groupby("rho")["ci_cover"].mean()
        b = d3.groupby("rho")["ci_cover"].mean()
        for rho in sorted(set(a.index) & set(b.index)):
            print(f"  rho={rho:>5}: {a[rho]:.4f} -> {b[rho]:.4f}")
        print(f"  mean:      {a.mean():.4f} -> {b.mean():.4f}")

    if args.write:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        merged.to_csv(args.out, index=False)
        print(f"\nwrote {len(merged)} rows to {args.out}")
    else:
        print("\n(dry run: pass --write to update "
              f"{os.path.relpath(args.out, _ROOT)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
