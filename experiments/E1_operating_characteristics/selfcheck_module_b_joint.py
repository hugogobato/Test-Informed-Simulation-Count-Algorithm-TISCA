#!/usr/bin/env python3
"""Deterministic pre-flight gate for the expensive joint Module-B grid.

This is intentionally the same test file used by the normal project test suite,
not a second set of notebook-only assertions.  It covers the eleven validation
requirements attached to G1 and exits non-zero on the first failed gate.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / "tisca" / "tests" / "test_outer_mc_joint.py"


def main() -> int:
    checks = [
        "K=1 reduction to the scalar Module-B estimand",
        "exact K+1 loss-to-contrast mapping",
        "positive-semidefinite implied Sigma_D on all 33 covariance configurations",
        "increasing unadjusted FWER as K increases",
        "Bonferroni, Holm and Romano-Wolf global-null FWER control",
        "BH FDR stored separately from FWER",
        "joint power sensitivity to cross-contrast dependence",
        "Romano-Wolf preservation of K-dimensional dependence",
        "fixed-seed bit reproducibility",
        "outer-shard offset invariance",
        "no repeated two-column block for K greater than one",
    ]
    print("[joint Module B self-check] gates:")
    for index, check in enumerate(checks, 1):
        print(f"  {index:2d}. {check}")
    command = [sys.executable, "-m", "pytest", str(TEST), "-q"]
    result = subprocess.run(command, cwd=ROOT, env=os.environ.copy(), check=False)
    if result.returncode:
        print("[FAIL] joint Module B self-check")
        return result.returncode
    print(f"[PASS] all {len(checks)} joint Module B gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
