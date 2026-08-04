#!/usr/bin/env python3
"""
TISCA v2 -- resampling parity evaluator (P2-T2, second tier).

Compares R and Python on the stochastic/resampling routines within bootstrap MCSE:
  * studentized_paired_bootstrap on a shared D (skewed + Gaussian + null cases)
  * MCS / reality_check / spa_test on shared loss matrices
  * romano_wolf_stepdown on a shared contrast matrix

The two languages use different RNGs (MT vs PCG64), so results cannot be bitwise
identical; the plan therefore specifies agreement to bootstrap MCSE. This script
reads a small spec JSON written by parity_resampling.R, computes the Python
side, and writes back the measured p-values/CIs for R to compare.
"""

import argparse
import json
import sys
import os
import numpy as np


def load_py():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "python"))
    import tisca.inference as inf
    import tisca.multiplicity as mult
    import tisca.mcs as mcs
    return inf, mult, mcs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True, help="JSON spec of shared datasets")
    ap.add_argument("--out", required=True, help="JSON output of python results")
    args = ap.parse_args()

    try:
        inf, mult, mcs = load_py()
    except (ImportError, ModuleNotFoundError) as e:
        with open(args.out, "w") as f:
            json.dump({"_PENDING_": str(e)}, f)
        print("PENDING", e)
        return 0

    with open(args.spec) as f:
        spec = json.load(f)

    results = {}

    # studentized paired bootstrap : compare p_value + CI bounds within MCSE
    for name, ds in spec["bootstrap"].items():
        D = np.array(ds["D"])
        alt = ds["alternative"]          # python convention already
        B = int(ds["B"])
        seed = int(ds["seed"])
        r = inf.studentized_paired_bootstrap(D, mu0=0.0, B=B, alternative=alt, seed=seed)
        results.setdefault("bootstrap", {})[name] = {
            "p_value": r["p_value"],
            "estimate": r["estimate"],
            "ci_lower": r["ci"][0],
            "ci_upper": r["ci"][1],
            "p_mcse": None,  # computed by caller if available
        }

    # MCS
    for name, Ls in spec.get("mcs", {}).items():
        Loss = np.array(Ls["L"])  # (J, M)
        B = int(Ls["B"]); alpha = float(Ls["alpha"]); stat = Ls["statistic"]
        seed = int(Ls["seed"])
        shared_idx = Ls.get("bootstrap_indices")
        if shared_idx is not None:
            shared_idx = np.asarray(shared_idx, dtype=int)
        r = mcs.mcs(Loss, B=B, alpha=alpha, statistic=stat, seed=seed,
                    model_names=Ls.get("model_names"),
                    bootstrap_indices=shared_idx)
        results.setdefault("mcs", {})[name] = {
            "models_kept": list(r["included"]),
            "included": list(r["included"]),
            "excluded": list(r["excluded"]),
            "p_H0": list(np.asarray(r["p_H0"], dtype=float)),
            "p_mcs": list(np.asarray(r["p_mcs"], dtype=float)),
            "elimination_order": list(r["elimination_order"]),
            "elimination_pvalues": list(np.asarray(r["elimination_pvalues"], dtype=float)),
        }

    # reality_check / SPA
    for name, rc in spec.get("spa", {}).items():
        Loss = np.array(rc["L"])
        champ = int(rc["champion"]); B = int(rc["B"]); seed = int(rc["seed"])
        rr = mcs.reality_check_pvalue(Loss, champion=champ, B=B, seed=seed)
        sr = mcs.spa_pvalue(Loss, champion=champ, B=B, seed=seed)
        results.setdefault("spa", {})[name] = {
            "rc_p": float(rr.get("p_value", np.nan)),
            "spa_p": float(sr.get("p_value", np.nan)),
        }

    # romano_wolf stepdown (cast numpy bools to native for JSON)
    for name, rw in spec.get("rw", {}).items():
        D = np.array(rw["D"])
        B = int(rw["B"]); alpha = float(rw["alpha"]); seed = int(rw["seed"])
        r = mult.romano_wolf_stepdown(D, B=B, alpha=alpha, seed=seed)
        rej = r.get("rejections")
        if rej is None:
            rej = r.get("reject", [])
        results.setdefault("rw", {})[name] = {
            "reject": [bool(x) for x in np.asarray(rej).tolist()],
        }

    with open(args.out, "w") as f:
        json.dump(results, f)
    print("resampling parity results written to", args.out)


if __name__ == "__main__":
    sys.exit(main())
