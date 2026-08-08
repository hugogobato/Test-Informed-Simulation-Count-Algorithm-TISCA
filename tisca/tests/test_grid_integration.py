"""Integration cover for the defects the unit tests structurally could not catch.

The audit's closing observation was that the suite passed while several *integration*
faults were live. Each of those faults lived in a seam no unit test crossed:

* the notebook generator carried a private copy of the E1 grid, so it kept emitting
  the pre-repair 1,299-cell version while the canonical grid had moved to 1,983 --
  nothing compared the two;
* that stale copy handed a numeric ``rho`` to the row-bootstrap empirical family,
  which raises, so every empirical cell in the Colab notebooks would have died at run
  time -- but only when actually run on Colab;
* ``sample_batch("empirical", matrix=...)`` raised on its own default ``rho = 0.0``,
  which no test exercised because every test spelled out ``rho=None``.

These tests cross those seams: they build the real grid, push representative cells
(including both empirical variants) through the real engine, and check the generated
notebooks against the canonical counts.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from tisca import multiplicity
from tisca.outermc import e1_grid, engine

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_MATRIX_CSV = os.path.join(_ROOT, "legacy", "Paper_Experiments", "DGP1_500_results.csv")
_NOTEBOOKS = os.path.join(_ROOT, "notebooks")


def _matrix():
    import pandas as pd

    if not os.path.exists(_MATRIX_CSV):
        pytest.skip("legacy DGP1_500_results.csv not present")
    raw = pd.read_csv(_MATRIX_CSV)
    return raw[["mvbcf_pehe1", "bcf_pehe1"]].to_numpy(dtype=float)


def _grid():
    return e1_grid.make_grid("ABCD", matrix=_matrix(),
                             planning_alpha=multiplicity.planning_alpha)


def test_canonical_grid_shape_and_empirical_rho_convention():
    """The grid is 1,983 cells and the row bootstrap never carries a numeric rho."""
    grid = _grid()
    assert len(grid) == sum(e1_grid.EXPECTED.values()) == 1983
    for module, expected in e1_grid.EXPECTED.items():
        assert sum(c["module"] == module for c in grid) == expected, module
    assert len({c["cell_id"] for c in grid}) == len(grid), "cell_id must be unique"

    row_bootstrap = [c for c in grid if c["config"]["family"] == "empirical"]
    assert row_bootstrap, "the row-bootstrap family must be present"
    assert all(c["config"]["rho"] is None for c in row_bootstrap), (
        "family 'empirical' reproduces its matrix's own dependence and refuses an "
        "imposed rho; a numeric rho here fails at run time, not at build time")
    copula = [c for c in grid if c["config"]["family"] == "empirical_copula"]
    assert copula and all(isinstance(c["config"]["rho"], float) for c in copula)

    # Every empirical cell must actually receive the loss matrix.
    for c in grid:
        if c["config"]["family"].startswith("empirical"):
            assert c["config"]["matrix"] is not None, c["cell_id"]


@pytest.mark.parametrize("family", ["empirical", "empirical_copula", "t3", "mix"])
def test_representative_cells_execute_end_to_end(family):
    """Push real grid cells through the real engine. This is the seam that broke."""
    grid = _grid()
    cells = [c for c in grid if c["config"]["family"] == family][:2]
    assert cells, family
    for cell in cells:
        cfg = dict(cell["config"])
        cfg["R"] = 200                      # keep it a test, not a study
        summary, _, _ = engine.run_e1(cfg)
        assert np.isfinite(summary["reject_rate"])
        assert np.isfinite(summary["E_J"])
        assert 0.0 <= summary["reject_rate"] <= 1.0


def test_oracle_sigma_resolves_for_every_grid_key():
    """D6 plans from sigma_D, so every distinct cell key must resolve to a finite one."""
    grid = _grid()
    keys = {(c["config"]["family"], c["config"]["rho"],
             c["config"]["sigma_a"], c["config"]["sigma_b"]) for c in grid}
    assert len(keys) > 40
    for family, rho, sa, sb in keys:
        cfg = {"family": family, "rho": rho, "sigma_a": sa, "sigma_b": sb,
               "matrix": _matrix() if family.startswith("empirical") else None}
        sigma = engine.sigma_D_true(cfg)
        assert np.isfinite(sigma) and sigma > 0, (family, rho, sa, sb, sigma)


@pytest.mark.parametrize("name,module,expected", [
    ("E1_modA_C_D.ipynb", "ACD", sum(e1_grid.EXPECTED[m] for m in "ACD")),
    ("E1_modB_shard1.ipynb", "B", e1_grid.EXPECTED["B"] // 2),
    ("E1_modB_shard2.ipynb", "B", e1_grid.EXPECTED["B"] - e1_grid.EXPECTED["B"] // 2),
])
def test_generated_notebooks_agree_with_the_canonical_grid(name, module, expected):
    """The Colab notebooks must import the grid and assert the canonical counts.

    A notebook that restates the grid inline is the failure mode being guarded
    against, so both facts are checked: that ``e1_grid`` is imported, and that the
    count the notebook asserts matches the count the module declares.
    """
    path = os.path.join(_NOTEBOOKS, name)
    if not os.path.exists(path):
        pytest.skip(f"{name} not generated")
    source = "".join("".join(c["source"]) for c in json.load(open(path))["cells"]
                     if c["cell_type"] == "code")
    assert "e1_grid" in source, f"{name} must import the canonical grid module"
    assert "e1_grid.make_grid" in source, f"{name} must build the grid from e1_grid"
    assert str(expected) in source, (
        f"{name} does not assert the canonical count {expected}; the notebook and "
        f"e1_grid.EXPECTED have drifted apart")
    compile(source, name, "exec")
