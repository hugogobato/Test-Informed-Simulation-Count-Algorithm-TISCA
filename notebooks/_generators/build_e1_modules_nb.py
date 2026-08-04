#!/usr/bin/env python3
"""Generate the three planned E1 Phase-3 Colab runners.

The generated notebooks are:

* ``E1_modA_C_D.ipynb``: Modules A, C, and D, 840 + 243 + 216 cells;
* ``E1_modB_shard1.ipynb``: the first 300 Module-B cells;
* ``E1_modB_shard2.ipynb``: the remaining 300 Module-B cells.

All three notebooks use the pure-Python P2-T3 harness, checkpoint one completed
cell at a time under ``/content``, and offer the CSV through
``google.colab.files.download`` at the end. They assume the Colab session runs
to completion; a runtime failure loses the local checkpoint.

Regenerate all three with::

    python notebooks/_generators/build_e1_modules_nb.py
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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


COMMON_SETUP = r'''
import itertools
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

REPO_URL = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
LOCAL_REPO = "/content/Test-Informed-Simulation-Count-Algorithm-TISCA"
CLONED_REPO = "/content/TISCA_repo"

if os.path.isdir(os.path.join(LOCAL_REPO, "tisca", "python")):
    REPO_ROOT = LOCAL_REPO
    SOURCE_ROOT = os.path.join(LOCAL_REPO, "tisca", "python")
elif os.path.isdir(os.path.join(CLONED_REPO, "tisca", "python")):
    REPO_ROOT = CLONED_REPO
    SOURCE_ROOT = os.path.join(CLONED_REPO, "tisca", "python")
else:
    subprocess.run(["git", "clone", "--depth", "1", REPO_URL, CLONED_REPO], check=True)
    REPO_ROOT = CLONED_REPO
    SOURCE_ROOT = os.path.join(CLONED_REPO, "tisca", "python")

assert os.path.isdir(SOURCE_ROOT), f"TISCA Python package not found at {SOURCE_ROOT}"
sys.path.insert(0, SOURCE_ROOT)

OUTPUT_ROOT = "/content/TISCA_E1"
os.makedirs(OUTPUT_ROOT, exist_ok=True)
print("Outputs/checkpoints:", OUTPUT_ROOT)

from tisca.outermc import engine, summarize_ocs
from tisca import multiplicity

ALPHA = 0.05
DELTA = 0.5
JMAX = 1000
MATRIX_CANDIDATES = [
    os.path.join(REPO_ROOT, "legacy", "Paper_Experiments", "DGP1_500_results.csv"),
    os.path.join(OUTPUT_ROOT, "E1_empirical_loss_matrix.npy"),
    os.path.join(OUTPUT_ROOT, "E1_empirical_loss_matrix.csv"),
    "/content/E1_empirical_loss_matrix.npy",
    "/content/E1_empirical_loss_matrix.csv",
]


def load_empirical_matrix():
    """Derive the empirical M x 2 pair from the committed 500 x 20 matrix."""
    path = next((p for p in MATRIX_CANDIDATES if os.path.exists(p)), None)
    if path is None or not os.path.exists(path):
        raise FileNotFoundError(
            "The committed legacy/Paper_Experiments/DGP1_500_results.csv was not "
            "found in the repository checkout."
        )
    if path.lower().endswith(".csv"):
        raw = pd.read_csv(path)
        required = ["mvbcf_pehe1", "bcf_pehe1"]
        if all(column in raw.columns for column in required):
            matrix = raw[required].to_numpy(dtype=float)
            print("[PASS] empirical pair:", required, "from", path)
            return matrix
        matrix = pd.read_csv(path, header=None).to_numpy(dtype=float)
    else:
        matrix = np.load(path)
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != 2 or matrix.shape[0] < 2:
        raise ValueError(f"empirical matrix must have shape (M, 2), got {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("empirical matrix contains non-finite values")
    print("[PASS] empirical matrix:", matrix.shape, "from", path)
    return matrix


EMPIRICAL_MATRIX = load_empirical_matrix()


def _append(row, path):
    frame = pd.DataFrame([row])
    header = not os.path.exists(path) or os.path.getsize(path) == 0
    frame.to_csv(path, mode="a", header=header, index=False)


def run_grid(grid, output_name):
    """Run and checkpoint a deterministic grid, resuming completed cell IDs."""
    output_file = os.path.join(OUTPUT_ROOT, output_name)
    error_file = output_file.replace("_results.csv", "_errors.csv")
    done = set()
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        old = pd.read_csv(output_file)
        if "cell_id" not in old.columns:
            raise ValueError(f"existing checkpoint lacks cell_id: {output_file}")
        done = set(old["cell_id"].astype(str))

    pending = [c for c in grid if c["cell_id"] not in done]
    print(f"{output_name}: {len(done)} completed, {len(pending)} pending, {len(grid)} expected")
    failures = []
    started = time.time()
    for cell in tqdm(pending, desc=output_name, unit="cell"):
        try:
            cfg = dict(cell["config"])
            summary, _, _ = engine.run_e1(cfg)
            row = summarize_ocs([summary]).iloc[0].to_dict()
            row.update(cell["factors"])
            row.update(
                cell_id=cell["cell_id"],
                module=cell["module"],
                projected_R=cell["config"]["R"],
                bootstrap_B=cell["config"].get("B", np.nan),
            )
            _append(row, output_file)
        except Exception as exc:
            print("[FAIL]", cell["cell_id"], repr(exc))
            failures.append({"cell_id": cell["cell_id"], "error": repr(exc)})

    if failures:
        pd.DataFrame(failures).to_csv(error_file, index=False)
        raise RuntimeError(f"{len(failures)} cells failed; see {error_file}")

    result = pd.read_csv(output_file)
    if set(result["cell_id"].astype(str)) != {c["cell_id"] for c in grid}:
        missing = sorted({c["cell_id"] for c in grid} - set(result["cell_id"].astype(str)))
        raise RuntimeError(f"checkpoint incomplete; missing {len(missing)} cells, first={missing[:3]}")
    if result["cell_id"].duplicated().any():
        raise RuntimeError("duplicate cell_id detected in checkpoint")
    print(f"[PASS] {len(result)} rows in {output_file}; elapsed {time.time() - started:.0f}s")
    return output_file, result


def download_fallback(output_file):
    try:
        from google.colab import files
        files.download(output_file)
        print("Downloaded:", output_file)
    except Exception as e:
        print("(Not on Colab / download skipped):", e)
'''


def _cfg(**kwargs: object) -> dict:
    cfg = {
        "design": "D4",
        "family": "normal",
        "rho": 0.0,
        "sigma_a": 1.0,
        "sigma_b": 1.0,
        "theta": 0.0,
        "sigma_D": None,
        "R": 5000,
        "J0": 50,
        "Jmax": 1000,
        "alpha": 0.05,
        "alpha_adj": 0.05,
        "mode": 1,
        "delta": 0.5,
        "power_target": 0.80,
        "gamma": 0.20,
        "correction": "none",
        "K": 1,
        "matrix": "EMPIRICAL_MATRIX" if kwargs.get("family") == "empirical" else None,
        "seed": 0,
        "fixed_J": None,
        "mcse": 0.05,
        "batch": 50,
        "B": 50,
    }
    cfg.update(kwargs)
    return cfg


def _cell(module: str, index: int, factors: dict, config: dict) -> dict:
    """Build a JSON-safe cell descriptor; replace the matrix symbol in notebook code."""
    parts = [module] + [f"{k}={v}" for k, v in factors.items()]
    return {
        "cell_id": f"{module}_{index:04d}",
        "module": module,
        "factors": {"cell_index": index, **factors},
        "config": config,
        "label": "|".join(parts),
    }


def _acd_grid_code() -> str:
    return r'''
def make_acd_grid():
    cells = []
    index = 0

    # Module A: F(7) x rho(5) x theta(4) x Design(6), K=1, J0=50, B=50.
    for family, rho, theta_mult, design in itertools.product(
        ["normal", "lognormal", "gamma", "mix", "beta", "t3", "empirical"],
        [-0.3, 0.0, 0.3, 0.6, 0.9],
        [0.0, 0.5, 1.0, 2.0],
        ["D1", "D2", "D3", "D4", "D5", "D6"],
    ):
        factors = dict(module_cell="A", family=family, rho=rho,
                       theta_mult=theta_mult, design=design, J0=50, B=50,
                       sigma_ratio=1.0)
        config = {
            "design": design, "family": family, "rho": rho,
            "sigma_a": 1.0, "sigma_b": 1.0, "theta": DELTA * theta_mult,
            "sigma_D": None, "R": 5000, "J0": 50, "Jmax": JMAX,
            "alpha": ALPHA, "alpha_adj": ALPHA, "mode": 1, "delta": DELTA,
            "power_target": 0.80, "gamma": 0.20, "correction": "none",
            "K": 1, "matrix": EMPIRICAL_MATRIX if family == "empirical" else None,
            "seed": 100000 + index, "fixed_J": None, "mcse": 0.05,
            "batch": 50, "B": 50,
        }
        cells.append(dict(cell_id=f"A_{index:04d}", module="A",
                          factors=factors, config=config))
        index += 1

    # Module C: J0(3) x B(3) x F(3) x rho(3) x Design(3), theta=delta.
    for J0, B, family, rho, design in itertools.product(
        [25, 50, 100], [25, 50, 100], ["normal", "lognormal", "gamma"],
        [-0.3, 0.3, 0.9], ["D2", "D3", "D4"],
    ):
        factors = dict(module_cell="C", family=family, rho=rho,
                       theta_mult=1.0, design=design, J0=J0, B=B,
                       sigma_ratio=1.0)
        config = {
            "design": design, "family": family, "rho": rho,
            "sigma_a": 1.0, "sigma_b": 1.0, "theta": DELTA,
            "sigma_D": None, "R": 5000, "J0": J0, "Jmax": JMAX,
            "alpha": ALPHA, "alpha_adj": ALPHA, "mode": 1, "delta": DELTA,
            "power_target": 0.80, "gamma": 0.20, "correction": "none",
            "K": 1, "matrix": None, "seed": 200000 + index,
            "fixed_J": None, "mcse": 0.05, "batch": B, "B": B,
        }
        cells.append(dict(cell_id=f"C_{index - 840:04d}", module="C",
                          factors=factors, config=config))
        index += 1

    # Module D: sigma-ratio(2) x F(3) x rho(3) x theta(2) x Design(6).
    for sigma_ratio, family, rho, theta_mult, design in itertools.product(
        [1.0, 2.0], ["normal", "lognormal", "gamma"], [-0.3, 0.3, 0.9],
        [0.0, 1.0], ["D1", "D2", "D3", "D4", "D5", "D6"],
    ):
        factors = dict(module_cell="D", family=family, rho=rho,
                       theta_mult=theta_mult, design=design, J0=50, B=50,
                       sigma_ratio=sigma_ratio)
        config = {
            "design": design, "family": family, "rho": rho,
            "sigma_a": sigma_ratio, "sigma_b": 1.0,
            "theta": DELTA * theta_mult, "sigma_D": None, "R": 5000,
            "J0": 50, "Jmax": JMAX, "alpha": ALPHA, "alpha_adj": ALPHA,
            "mode": 1, "delta": DELTA, "power_target": 0.80, "gamma": 0.20,
            "correction": "none", "K": 1, "matrix": None,
            "seed": 300000 + index, "fixed_J": None, "mcse": 0.05,
            "batch": 50, "B": 50,
        }
        cells.append(dict(cell_id=f"D_{index - 1083:04d}", module="D",
                          factors=factors, config=config))
        index += 1

    assert len(cells) == 1299, len(cells)
    assert sum(c["module"] == "A" for c in cells) == 840
    assert sum(c["module"] == "C" for c in cells) == 243
    assert sum(c["module"] == "D" for c in cells) == 216
    return cells


ACD_GRID = make_acd_grid()
print("Module A/C/D cells:", {m: sum(c["module"] == m for c in ACD_GRID) for m in "ACD"})
'''


def _module_b_grid_code() -> str:
    return r'''
def make_module_b_grid():
    cells = []
    index = 0
    for K, correction, rho, family, theta, design in itertools.product(
        [1, 3, 6],
        ["none", "bonferroni", "holm", "bh", "romano_wolf"],
        [-0.3, 0.0, 0.3, 0.6, 0.9],
        ["normal", "empirical"],
        [0.0, DELTA],
        ["D3", "D4"],
    ):
        alpha_plan, alpha_note = multiplicity.planning_alpha(
            correction, K, alpha=ALPHA, r=1
        )
        factors = dict(module_cell="B", family=family, rho=rho,
                       theta_mult=round(theta / DELTA, 3), design=design,
                       J0=50, B=999, K=K, correction=correction,
                       alpha_plan=alpha_plan, alpha_note=alpha_note)
        config = {
            "design": design, "family": family, "rho": rho,
            "sigma_a": 1.0, "sigma_b": 1.0, "theta": theta,
            "sigma_D": None, "R": 2000, "J0": 50, "Jmax": JMAX,
            "alpha": ALPHA, "alpha_adj": alpha_plan, "mode": 1,
            "delta": DELTA, "power_target": 0.80, "gamma": 0.20,
            "correction": correction, "K": K,
            "matrix": EMPIRICAL_MATRIX if family == "empirical" else None,
            "seed": 400000 + index, "fixed_J": None, "mcse": 0.05,
            "batch": 50, "B": 999,
        }
        cells.append(dict(cell_id=f"B_{index:04d}", module="B",
                          factors=factors, config=config))
        index += 1
    assert len(cells) == 600, len(cells)
    return cells


MODULE_B_GRID = make_module_b_grid()
print("Module B total cells:", len(MODULE_B_GRID))
'''


def _notebook_acd() -> dict:
    cells: list[dict] = []
    _md(cells, """
        # E1 Modules A/C/D: operating characteristics, tuning, and variance ratio

        This is the combined Phase-3 runner specified in `REVISION_PLAN.md` P3-T2.
        It executes Module A (840 cells), Module C (243 cells), and Module D
        (216 cells), for 1,299 cells total. Each cell uses 5,000 outer repetitions,
        is checkpointed under `/content`, and is intended to run to completion in
        one Colab session. The empirical family is derived automatically from the
        committed 500 x 20 `legacy/Paper_Experiments/DGP1_500_results.csv`, using
        the pre-specified `mvbcf_pehe1` versus `bcf_pehe1` pair.
    """)
    _code(cells, COMMON_SETUP)
    _md(cells, """
        ## Define the three grids

        Module A spans all seven families, five correlations, four effect levels,
        and six designs. Module C varies pilot size, checkpoint batch, family,
        correlation, and the three adaptive/two-stage designs at the planning
        alternative. Module D varies the marginal variance ratio, family,
        correlation, effect, and design.
    """)
    _code(cells, _acd_grid_code())
    _code(cells, """
        OUTPUT_FILE, RESULTS = run_grid(ACD_GRID, "E1_modA_C_D_results.csv")
        print(RESULTS.groupby("module").size())
        print(RESULTS.head())
    """)
    _md(cells, """
        ## Basic completion checks

        The checkpoint must contain exactly one row per declared cell. The D6
        normal cells in Module A provide the oracle sanity check. Module B is
        deliberately separate because its bootstrap-heavy correction grid is
        split into two notebooks.
    """)
    _code(cells, """
        assert len(RESULTS) == 1299
        assert RESULTS["cell_id"].is_unique
        d6 = RESULTS[(RESULTS["module"] == "A") & (RESULTS["design"] == "D6") &
                     (RESULTS["family"] == "normal")]
        assert len(d6) == 20, len(d6)
        print("[PASS] A/C/D complete; D6 normal oracle rows:", len(d6))
        download_fallback(OUTPUT_FILE)
    """)
    return _nb(cells)


def _notebook_b(shard: int) -> dict:
    cells: list[dict] = []
    _md(cells, f"""
        # E1 Module B, shard {shard} of 2: multiplicity

        This notebook runs Module B from P3-T2. The full factorial contains 600
        cells: K(3) x correction(5) x rho(5) x family(2) x theta(2) x design(2).
        Shard {shard} runs exactly 300 cells, uses 2,000 outer repetitions, and
        uses B=999 as the bootstrap budget recorded with every row. The scalar
        P2-T3 harness receives the correction-specific planning level from
        `tisca.multiplicity.planning_alpha`, including the Romano-Wolf schedule.
        Results are checkpointed under `/content` one cell at a time. This notebook
        assumes the session completes; a runtime failure loses the local checkpoint.

        The empirical family is derived automatically from the committed 500 x 20
        `legacy/Paper_Experiments/DGP1_500_results.csv`, using the pre-specified
        `mvbcf_pehe1` versus `bcf_pehe1` pair.
    """)
    _code(cells, COMMON_SETUP)
    _code(cells, _module_b_grid_code())
    start = 0 if shard == 1 else 300
    end = 300 if shard == 1 else 600
    _code(cells, f"""
        SHARD = {shard}
        START, END = {start}, {end}
        SHARD_GRID = MODULE_B_GRID[START:END]
        assert len(SHARD_GRID) == 300
        assert SHARD_GRID[0]["cell_id"] == f"B_{{START:04d}}"
        assert SHARD_GRID[-1]["cell_id"] == f"B_{{END - 1:04d}}"
        print("Module B shard", SHARD, "cells:", len(SHARD_GRID))
    """)
    _code(cells, f"""
        OUTPUT_FILE, RESULTS = run_grid(SHARD_GRID, "E1_modB_shard{shard}_results.csv")
        assert len(RESULTS) == 300
        assert RESULTS["cell_id"].is_unique
        print("[PASS] Module B shard {shard} complete")
        download_fallback(OUTPUT_FILE)
    """)
    return _nb(cells)


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


def main() -> None:
    outputs = {
        "E1_modA_C_D.ipynb": _notebook_acd(),
        "E1_modB_shard1.ipynb": _notebook_b(1),
        "E1_modB_shard2.ipynb": _notebook_b(2),
    }
    for name, notebook in outputs.items():
        out = ROOT / name
        out.write_text(json.dumps(notebook, indent=1) + "\n")
        for index, cell in enumerate(notebook["cells"]):
            source = "".join(cell["source"])
            if cell["cell_type"] == "code":
                compile(source, f"{name}:cell-{index}", "exec")
        print("wrote", out, len(notebook["cells"]), "cells")


if __name__ == "__main__":
    main()
