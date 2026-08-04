# `legacy/` — Original submission code (v1)

This directory preserves the code that accompanied the **original submission** of:

> Souto, H. G. and Louzada Neto, F. "Beyond Arbitrary Replications: A Principled
> Approach to Simulation Design in Causal Inference." International Journal of
> Data Science and Analytics (major revision, 2026).

**It is retained for audit and historical reference only, and is entirely superseded.**

Everything in this directory corresponds to **TISCA v1** as submitted. It is
superseded by the restructured `tisca/`, `experiments/`, `results/`, `figures/`
and `docs/` elsewhere in this repository, which implement **TISCA v2** as
described in `REVISION_PLAN.md`. Do not use anything under `legacy/` as the
basis for new analyses. In particular:

- `legacy/Paper_Experiments/` contains the original DGP1 simulation notebooks,
  `.R` scripts and result CSVs. Note that `DGP1_100_results.csv` has **99 rows,
  not 100** and is an orphan run (not nested in the others); several rows were
  generated with an unrecoverable `as.numeric(Sys.time())` seed. These defects
  are documented in §1.4 of the revision plan and are fixed in v2.
- `legacy/Usage_Examples/` contains the v1 usage examples (R and Python).
- `legacy/Bibliometric_Study/` contains the v1 bibliometric spreadsheet and its
  analysis notebook, including the hand-summed-percentage error documented in
  §1.5 of the revision plan.
- `legacy/v1/` contains the monolithic v1 source files `TISCA.R` and `tisca.py`.

Known issues in the v1 code (all addressed in v2) are enumerated in §1.6 of
`REVISION_PLAN.md`: the `welch_t_test()` with `paired = FALSE` (the defect at
the heart of review issue #2), per-column NA dropping, the
`scipy.stats.t.cdf(..., loc, ncp)` misread of the noncentral t, and the two
different power estimators (analytic vs 500-draw inner Monte Carlo).

The original authors' study that the case study reproduces is
**not** stored here — it is linked from the paper and referenced from the
methods, not copied.
