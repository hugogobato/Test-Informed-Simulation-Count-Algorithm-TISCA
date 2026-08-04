# `notebooks/` — Google Colab notebooks

Notebooks that run on Google Colab (Python runtime, R via `Rscript` calls) to
build and execute the heavy parts of the project across the 33-session execution
model described in `REVISION_PLAN.md` §5.

| Notebook | Task | Purpose |
|---|---|---|
| `P0T4_build_rlib_bundle.ipynb` | **P0-T4** | compile the heavy R packages **once** into `tisca_rlib.tar.gz`, record `sessionInfo()`, publish, and restore-and-verify in a fresh session |
| `P1_math_verification.ipynb` | **P1-T2/T1 (light)** | numerically verify the v2 planning formulas (power modes M1-M5, the plan's §1.2 J table, variance-inflation, rooted-vs-unrooted estimand, the `scipy.stats.t` bug magnitude). Pure NumPy/SciPy, seconds to minutes; no R bundle needed. |

`P0T4_build_rlib_bundle.ipynb` is the main notebook to run once during Phase 0.
Once it has produced and published the tarball, paste the direct-download URL
into Cell 9; every later worker notebook will restore the bundle with a short
`wget <URL> && tar xzf && .libPaths(<dest>/tisca_rlib/rlib)`.

`P1_math_verification.ipynb` is the Phase-1 notebook you run yourself (single
run, on the laptop or a free Colab CPU runtime, no R bundle). Every cell prints
`[PASS]`/`[FAIL]`; a clean run is the audit record for the Phase-1 specs in
`docs/`.

The notebooks are generated from `notebooks/_generators/` (`build_p0t4_bundle_nb.py`,
`build_p1_verification_nb.py`). Regenerate with
`python notebooks/_generators/<name>.py`. Keep the generator and the `.ipynb` in
sync when you edit either.
