#!/usr/bin/env python3
"""E1 outer-MC operating-characteristics runner (P2-T3 harness driver).

Runs one or more E1 cells and writes ``results/E1/operating_characteristics.csv``
with one row per cell, each operating characteristic carrying an MCSE.

Usage (cell or JSON grid):
  python run_e1.py --cell D4 normal 0.0 1 0.5 0 5000 --out results/E1/oc.csv
  python run_e1.py --grid e1_moduleA.json --dry-run
  python run_e1.py --cell D6 normal 0.0 1 0.5 0 20000 --acceptance

Positional order for ``--cell``: DESIGN FAMILY RHO THETA SIGMA_RATIO SEED R

Options:
  --J0 J0              pilot size (default 50)
  --Jmax J             budgeted max confirmatory rows per rep (default 1000)
  --alpha A            test level (default 0.05)
  --delta D            planning alternative (default 0.5)
  --power P            power target (default 0.80)
  --gamma G            variance-uncertainty assurance (default 0.20)
  --alpha-adj A        multiplicity-adjusted level (default: alpha/K)
  --K K                family size (default 1)
  --mcse M             D5 MCSE precision target (default 0.05)
  --fixed-j J          D6: force a fixed J
  --out PATH           output CSV (chunked/resumable)
  --chunk C            rows per checkpoint (default 50)
  --dry-run            print projected cost only, run nothing
  --acceptance         run the P2-T3 acceptance check (Type I .05, power .80)
  --matrix PATH        empirical family loss matrix (.npy, Mx2)
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                  "..", "..", "tisca", "python")))

from tisca.outermc import engine  # noqa: E402


def run_acceptance():
    """P2-T3 acceptance: bivariate normal, rho=0, oracle sigma -> recovers the
    nominal Type I (0.05) and power (0.80) within 2 Monte Carlo standard errors."""
    from tisca.outermc.designs import solve_J
    delta, sd = 0.5, np.sqrt(2.0)
    J = int(np.atleast_1d(solve_J(0.80, delta, sd, 0.05)[0]).max())
    R = 20_000
    mcse0, mcse1 = np.sqrt(.05 * .95 / R), np.sqrt(.80 * .20 / R)
    out = {"note": f"P2-T3 acceptance (fixed oracle J={J}, bivariate normal, rho=0)"}
    for theta, nom, mcse, lab in [(0.0, .05, mcse0, "Type_I"), (delta, .80, mcse1, "Power")]:
        c = dict(engine.DEFAULT_CONFIG)
        c.update(design="D6", family="normal", rho=0.0, sigma_a=1.0, sigma_b=1.0,
                 theta=theta, delta=delta, sigma_D=None, R=R, Jmax=2000, fixed_J=J)
        s, _, _ = engine.run_e1(c)
        est = s["reject_rate"]
        ok = abs(est - nom) <= 2 * mcse
        out[f"{lab}_est"] = round(est, 4)
        out[f"{lab}_nom"] = nom
        out[f"{lab}_mcse"] = round(mcse, 4)
        out[f"{lab}_ok"] = bool(ok)
    ok_all = out["Type_I_ok"] and out["Power_ok"]
    out["PASS"] = ok_all
    print(json.dumps(out, indent=2))
    return 0 if ok_all else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cell", nargs="+",
                   metavar="DESIGN_FAMILY_RHO_THETA_SIGMARATIO_SEED_R")
    p.add_argument("--grid", help="path to a JSON list of per-cell dicts")
    p.add_argument("--J0", type=int, default=50)
    p.add_argument("--Jmax", type=int, default=1000)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--alpha-adj", type=float, default=None)
    p.add_argument("--delta", type=float, default=0.5)
    p.add_argument("--power", type=float, default=0.80)
    p.add_argument("--gamma", type=float, default=0.20)
    p.add_argument("--K", type=int, default=1)
    p.add_argument("--mcse", type=float, default=0.05)
    p.add_argument("--matrix")
    p.add_argument("--out")
    p.add_argument("--chunk", type=int, default=50)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--fixed-j", type=int, default=None)
    p.add_argument("--acceptance", action="store_true")
    args = p.parse_args(argv)

    if args.acceptance:
        return run_acceptance()

    cells = []
    if args.cell:
        ds, fam, rho, th, sr, seed, R = args.cell
        cells.append(dict(design=ds, family=fam, rho=float(rho), theta=float(th),
                          sigma_ratio=float(sr), seed=int(seed), R=int(R)))
    elif args.grid:
        cells = json.load(open(args.grid))
    else:
        p.error("provide --cell or --grid")

    def _full(c0):
        c2 = dict(engine.DEFAULT_CONFIG)
        c2.update(design=c0.get("design", "D4"), family=c0["family"],
                  rho=c0["rho"], theta=c0["theta"], sigma_a=c0["sigma_ratio"],
                  sigma_b=1.0, delta=args.delta, power_target=args.power,
                  sigma_D=None, J0=args.J0, Jmax=args.Jmax, alpha=args.alpha,
                  alpha_adj=args.alpha_adj, K=args.K, mcse=args.mcse, gamma=args.gamma,
                  seed=c0["seed"], R=c0["R"], fixed_J=args.fixed_j)
        if args.alpha_adj is None:
            c2["alpha_adj"] = args.alpha / max(1, args.K)
        if args.matrix:
            c2["matrix"] = np.load(args.matrix)
        return c2

    if args.dry_run:
        print(json.dumps(engine.dry_run(_full(cells[-1])), indent=2))
        return 0

    out = args.out or os.path.join(os.path.dirname(__file__), "..", "..", "results", "E1", "operating_characteristics.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    import pandas as pd

    rows = []
    for c in cells:
        c2 = _full(c)
        summarize, res, _ = engine.run_e1(c2)
        from tisca.outermc import summarize_ocs
        row = summarize_ocs([summarize]).iloc[0].to_dict()
        rows.append(row)
    df = pd.DataFrame(rows)
    if out.endswith(".parquet"):
        df.to_parquet(out, engine="pyarrow")
    else:
        df.to_csv(out, index=False)
    print(f"[E1] wrote {len(df)} rows to {out}")
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
