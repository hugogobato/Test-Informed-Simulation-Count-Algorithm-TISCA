# P2-T2 — R port of TISCA v2 with parity tests

**Task (from `REVISION_PLAN.md`):** port the TISCA v2 reference implementation to R with
the same API as the Python module and a parity test suite: `>= 200` random configurations
agreeing to `1e-8` (analytic) / bootstrap MCSE (resampling). Remove `paired = FALSE` and
per-column NA dropping (=> listwise deletion + reported count + sensitivity note).

**Date:** 2026-08-04 · **Lane:** L-CODE · **Depends on:** P1-T2 (signed off, rev. 2)

---

## What was delivered

New/edited files under `tisca/R/` (the R twin of `tisca/python/tisca/`):

| Module | File | Contents |
|---|---|---|
| package entry | `tisca.R` | `%||%`, `tisca_version`, `power_()`, ready-check |
| estimands | `estimands.R` | PEHE (rooted), CATE MSE (unrooted), ATE sq err, RMSE_ATE + delta-method MCSE, ATE/CATE coverage, CIL, bias, interval score; per-rep validators |
| planning | `planning.R` | M1–M5 power (`pt(q, df, ncp)`), `sigma_ub`, `solve_mcse_J`, `solve_halfwidth_J`, `solve_power_J`, `combine_J`, `power()` |
| inference | `inference.R` | `contrast_from_columns` (listwise NA + drop count), `paired_t`, studentized paired bootstrap, BCa, McNemar |
| multiplicity | `multiplicity.R` | `adj_alpha`, `adjust_p_values`, `romano_wolf_stepdown`, `family_power`, `plan_rw_power` |
| mcs | `mcs.R` | reality check, SPA, MCS (`T_R`/`T_max`) with elimination path; CRAN oracle |
| procedure | `procedure.R` | Algorithm 1 (two-stage) + Algorithm 2 (adaptive, guarded) |
| validate | `validate.R` | spec Section 8 rejection rules |

Package metadata: `tisca/DESCRIPTION`, `tisca/NAMESPACE` (exports for `R CMD INSTALL`).

### v1 defects fixed (paper Section 1.6)
- **`paired = FALSE` removed.** All inference is now on paired contrasts `D_j` via one-sample t /
  studentized bootstrap / McNemar. `TISCA.R`'s `welch_t_test` is gone.
- **Per-column NA dropping replaced** by listwise deletion across the pair with the drop count
  returned (`contrast_from_columns`), satisfying reviewer suggestion 2(b). A
  sensitivity note is carried in the procedure output.
- `pt(q, df, ncp)` is kept (already correct). The parallel `scipy.stats.t` `ncp`-in-`loc` bug is a
  P2-T1 concern; a *separate* scipy `nct` NaN defect surfaced here is documented in §3.

## Tests (`tisca/tests/`)

- `run_tests.R` + `test_*.R`: unit tests (planning anchors, estimands, inference, multiplicity,
  MCS, procedure, validate, MCS-vs-CRAN-oracle).
- `parity_generate_reference.R` → `parity_reference.csv`: **407** configs (deterministic).
- `parity_eval_python.py`: evaluates the Python reference for each config.
- `parity_run.R`: Tier 1 diff (analytic 1e-8 / integer exact).
- `parity_resampling.R` + `parity_eval_python_resampling.py`: Tier 2 bootstrap-MCSE tier.
- `parity_all.R`: combined ACCEPTANCE runner.

## Acceptance results

### Analytic / integer parity (Tier 1) — 378/378 well-conditioned configs PASS (1e-8)
407 configs across `power` M1–M5, `sigma_ub`, `solve_power_J`, `solve_mcse_J`,
`solve_halfwidth_J`, `paired_t`, `adjust_p_values`, `adj_alpha`. `0` genuine mismatches
after two fixes: one R tolerance bug (below) and one R/Python recording mismatch in the generator.

The parity harness **surfaced a real R bug**: `solve_halfwidth_J` used `tolerance = 1e-6`,
which at σ=4, h=0.2 accepted `0.2000003 > 0.2` and returned `J = 1539`; Python correctly
returned `1540`. R tolerance is now `1e-12` (machine-scale) and returns the same `J` as Python.

### Resampling parity (Tier 2) — 10/10 checks PASS
Studentized paired bootstrap p-values (3), MCS kept sets (2), SPA + Reality Check (4),
Romano-Wolf rejection patterns (1) — all within bootstrap MCSE.

## §3. Issue flagged for the P2-T1 lane (not an R bug)

**28** analytic configs return `NaN` from the **Python** reference (`scipy.stats.nct.cdf(-crit, df, ncp)`
→ `NaN` for large `ncp`). These are high-power `M1`–`M5` cells where the far tail should be ~ `0`.
R's `pt()` returns the correct value. The parity harness classifies these as a **Python-nct-NaN
(P2-T1)** defect and does not count them against R. **Fix belongs in `tisca/python/tisca/planning.py`.**
See `tisca/tests/README_P2T2.md` for the recommended patch.
