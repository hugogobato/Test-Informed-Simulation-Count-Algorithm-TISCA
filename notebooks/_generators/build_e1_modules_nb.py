#!/usr/bin/env python3
"""Generate the three planned E1 Phase-3 Colab runners.

The generated notebooks are:

* ``E1_modA_C_D.ipynb``: Modules A, C, and D, 864 + 243 + 216 = 1,323 cells;
* ``E1_modB_shard1.ipynb``: the first half of the 660 Module-B cells;
* ``E1_modB_shard2.ipynb``: the remainder.

**The grid is not defined here.** It lives in ``tisca/python/tisca/outermc/
e1_grid.py`` and these notebooks import it, exactly as the local runner
``experiments/E1_operating_characteristics/run_e1_grid.py`` does. That module is
the single source of truth for the 1,983 cells.

The generator used to carry its own copy of the grid as a string literal, and the
two drifted: this file was still emitting the pre-repair 1,299-cell version, in
which the empirical family was one entry rather than the two it was split into
(row bootstrap plus copula variant), and in which every empirical cell was handed
a numeric ``rho``. Since the row bootstrap reproduces its matrix's own dependence
and refuses an imposed ``rho``, every one of those cells would now fail at run
time. Importing the canonical grid is what stops that recurring.

All three notebooks use the pure-Python P2-T3 harness, checkpoint one completed
cell at a time under ``/content``, and offer the CSV through
``google.colab.files.download`` at the end. They assume the Colab session runs
to completion; a runtime failure loses the local checkpoint.

Regenerate all three with::

    python notebooks/_generators/build_e1_modules_nb.py
"""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "tisca" / "python"))

from tisca.outermc import e1_grid  # noqa: E402  (needs the path above)

#: Module B is split across two Colab sessions; the boundary is derived, not typed.
_B_TOTAL = e1_grid.EXPECTED["B"]
B_SHARD_BOUNDS = {1: (0, _B_TOTAL // 2), 2: (_B_TOTAL // 2, _B_TOTAL)}
ACD_TOTAL = sum(e1_grid.EXPECTED[m] for m in "ACD")


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

from tisca.outermc import e1_grid, engine, summarize_ocs
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


GRID_CODE = r'''
# The grid comes from the repository, never from a copy pasted into this notebook:
# tisca/python/tisca/outermc/e1_grid.py is the single source of truth for all
# 1,983 cells, and experiments/E1_operating_characteristics/run_e1_grid.py reads
# exactly the same object. `matrix` is threaded in here because the two empirical
# families need the real loss pair; the row-bootstrap variant carries rho=None,
# which is why the grid must not be retyped by hand.
FULL_GRID = e1_grid.make_grid("ABCD", matrix=EMPIRICAL_MATRIX,
                              planning_alpha=multiplicity.planning_alpha)
ACD_GRID = [c for c in FULL_GRID if c["module"] in ("A", "C", "D")]
MODULE_B_GRID = [c for c in FULL_GRID if c["module"] == "B"]

for _m, _n in e1_grid.EXPECTED.items():
    _got = sum(c["module"] == _m for c in FULL_GRID)
    assert _got == _n, f"module {_m}: {_got} cells, expected {_n}"
print("cells per module:", {m: sum(c["module"] == m for c in FULL_GRID)
                            for m in "ABCD"})
print("A/C/D:", len(ACD_GRID), " B:", len(MODULE_B_GRID))
'''


def _notebook_acd() -> dict:
    a, c, d = (e1_grid.EXPECTED[m] for m in "ACD")
    cells: list[dict] = []
    _md(cells, f"""
        # E1 Modules A/C/D: operating characteristics, tuning, and variance ratio

        This is the combined Phase-3 runner specified in `REVISION_PLAN.md` P3-T2.
        It executes Module A ({a} cells), Module C ({c} cells), and Module D
        ({d} cells), for {ACD_TOTAL:,} cells total. Each cell uses 5,000 outer
        repetitions, is checkpointed under `/content`, and is intended to run to
        completion in one Colab session. The empirical families are derived
        automatically from the committed 500 x 20
        `legacy/Paper_Experiments/DGP1_500_results.csv`, using the pre-specified
        `mvbcf_pehe1` versus `bcf_pehe1` pair.
    """)
    _code(cells, COMMON_SETUP)
    _md(cells, """
        ## Load the canonical grid

        The cells are **imported** from `tisca/python/tisca/outermc/e1_grid.py`,
        not restated here, so this notebook and the local runner
        `experiments/E1_operating_characteristics/run_e1_grid.py` cannot drift
        apart. Module A spans the six parametric families plus the two empirical
        variants, five correlations, four effect levels, and six designs. Module C
        varies pilot size, checkpoint batch, family, correlation, and the three
        adaptive/two-stage designs at the planning alternative. Module D varies the
        marginal variance ratio, family, correlation, effect, and design.
    """)
    _code(cells, GRID_CODE)
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
    _code(cells, f"""
        assert len(RESULTS) == {ACD_TOTAL}, len(RESULTS)
        assert RESULTS["cell_id"].is_unique
        d6 = RESULTS[(RESULTS["module"] == "A") & (RESULTS["design"] == "D6") &
                     (RESULTS["family"] == "normal")]
        assert len(d6) == 20, len(d6)
        print("[PASS] A/C/D complete; D6 normal oracle rows:", len(d6))
        download_fallback(OUTPUT_FILE)
    """)
    return _nb(cells)


def _notebook_b(shard: int) -> dict:
    start, end = B_SHARD_BOUNDS[shard]
    size = end - start
    cells: list[dict] = []
    _md(cells, f"""
        # E1 Module B, shard {shard} of 2: multiplicity

        This notebook runs Module B from P3-T2. The full factorial contains
        {_B_TOTAL} cells: K(3) x correction(5) x rho(5) x family(2) x theta(2) x
        design(2), plus the row-bootstrap empirical family, whose dependence is
        fixed by its own data and which therefore carries no `rho` factor.
        Shard {shard} runs exactly {size} cells, uses 2,000 outer repetitions, and
        uses B=999 as the bootstrap budget recorded with every row. The scalar
        P2-T3 harness receives the correction-specific planning level from
        `tisca.multiplicity.planning_alpha`, including the Romano-Wolf schedule.
        Results are checkpointed under `/content` one cell at a time. This notebook
        assumes the session completes; a runtime failure loses the local checkpoint.

        The empirical families are derived automatically from the committed
        500 x 20 `legacy/Paper_Experiments/DGP1_500_results.csv`, using the
        pre-specified `mvbcf_pehe1` versus `bcf_pehe1` pair.
    """)
    _code(cells, COMMON_SETUP)
    _code(cells, GRID_CODE)
    _code(cells, f"""
        SHARD = {shard}
        START, END = {start}, {end}
        SHARD_GRID = MODULE_B_GRID[START:END]
        assert len(SHARD_GRID) == {size}, len(SHARD_GRID)
        assert SHARD_GRID[0]["cell_id"] == f"B_{{START:04d}}"
        assert SHARD_GRID[-1]["cell_id"] == f"B_{{END - 1:04d}}"
        print("Module B shard", SHARD, "cells:", len(SHARD_GRID))
    """)
    _code(cells, f"""
        OUTPUT_FILE, RESULTS = run_grid(SHARD_GRID, "E1_modB_shard{shard}_results.csv")
        assert len(RESULTS) == {size}, len(RESULTS)
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
