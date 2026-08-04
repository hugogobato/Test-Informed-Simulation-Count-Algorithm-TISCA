# `notebooks/` — Google Colab notebooks

Notebooks that run on Google Colab (Python runtime, R via `Rscript` calls) to
build and execute the heavy parts of the project across the 33-session execution
model described in `REVISION_PLAN.md` §5.

| Notebook | Task | Purpose |
|---|---|---|
| `P0T4_build_rlib_bundle.ipynb` | **P0-T4** | compile the heavy R packages **once** into `tisca_rlib.tar.gz`, record `sessionInfo()`, publish, and restore-and-verify in a fresh session |

`P0T4_build_rlib_bundle.ipynb` is the only notebook you need to run yourself
during Phase 0. Once it has produced and published the tarball, paste the
direct-download URL into Cell 9; every later worker notebook will restore the
bundle with a short `wget <URL> && tar xzf && .libPaths(<dest>/tisca_rlib/rlib)`.

The notebooks are generated from `notebooks/_generators/` (only
`build_p0t4_bundle_nb.py` exists for now). Regenerate with:
`python notebooks/_generators/build_p0t4_bundle_nb.py`. Keep the generator and
the `.ipynb` in sync when you edit either.
