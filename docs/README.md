# `docs/` — TISCA v2 specifications

Stable design and reproducibility documents for the revised project.

| Document | Task | Contents |
|---|---|---|
| `seed_rng_protocol.md` | P0-T2 | Seed allocation, stream separation (data vs fit), pilot/confirmatory disjoint blocks, sharding/completeness |
| `estimand_table.md` | P1-T1 (final for review) | Performance-measure estimand rows; full table + main.tex equation-fix text |
| `tisca_v2_spec.md` | P1-T2 (draft for review) | Formal spec of the v2 procedures (modes M1-M5, targets, multiplicity, variance inflation, two-stage/adaptive, inference, validation) |
| `power_target_guidance.md` | P1-T3 (draft for review) | When to use the power layer vs precision; the "do not use TISCA for X" list |

The plan references `docs/estimand_table.md` and `docs/tisca_v2_spec.md` as
the gatekeeping specs for Phases 1 and 2. `seed_rng_protocol.md` is finalised.

The Phase-1 numerical check backing the formulas in `estimand_table.md` and
`tisca_v2_spec.md` is `notebooks/P1_math_verification.ipynb` (generated from
`notebooks/_generators/build_p1_verification_nb.py`). It verifies the M1-M5
power functions against Monte Carlo, reproduces the Section 1.2 J values to
within one replication, checks the variance-inflation factor and the
rooted-vs-unrooted estimand distinction, and quantifies the `scipy.stats.t` bug.
It is pure NumPy/SciPy, so no R bundle is needed.
