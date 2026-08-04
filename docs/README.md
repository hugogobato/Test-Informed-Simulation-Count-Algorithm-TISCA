# `docs/` — TISCA v2 specifications

Stable design and reproducibility documents for the revised project.

| Document | Task | Contents |
|---|---|---|
| `seed_rng_protocol.md` | P0-T2 (specified; **acceptance test pending**) | Seed allocation, stream separation (data vs fit), model-owned RNG, pilot/confirmatory disjoint blocks, sharding/completeness |
| `estimand_table.md` | P1-T1 (final for review, rev. 2) | Performance-measure estimand rows; full table + main.tex equation-fix text |
| `tisca_v2_spec.md` | P1-T2 (draft for review, rev. 2) | Formal spec of the v2 procedures (modes M1-M5, targets, multiplicity, variance inflation, two-stage/adaptive, inference, validation) |
| `power_target_guidance.md` | P1-T3 (draft for review) | When to use the power layer vs precision; the "do not use TISCA for X" list |

The plan references `docs/estimand_table.md` and `docs/tisca_v2_spec.md` as
the gatekeeping specs for Phases 1 and 2.

`seed_rng_protocol.md` is **specified but not yet accepted**: acceptance is the
four-way identity test in its §3.5, which has to run against the real E3 driver.
Its §3.3 stream construction has been verified in R 4.3.3 in isolation.

The Phase-1 numerical check backing the formulas in `estimand_table.md` and
`tisca_v2_spec.md` is `notebooks/P1_math_verification.ipynb` (generated from
`notebooks/_generators/build_p1_verification_nb.py`). It recomputes the
Section 1.2 `J` table from the raw replication CSV, computes the
**spec-compliant** `J_final` under the v2 defaults (multiplicity + variance
inflation, spec §4.5), verifies M1-M5 against Monte Carlo (M5 against the real
studentized TOST rule and against exact quadrature, not against its own
approximation), demonstrates what the variance inflation actually buys
(assurance, not a smaller bias in `E[J]`), checks the rooted-vs-unrooted and
delta-method estimand rows, exhibits the counterexample behind the calibration
row, and quantifies the `scipy.stats.t` bug. It is pure NumPy/SciPy, so no R
bundle is needed. Every section prints `[PASS]`/`[FAIL]`; all sections pass as
of rev. 2.
