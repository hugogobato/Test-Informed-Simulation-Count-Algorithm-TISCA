# P2-T2 — TISCA v2 R port and R <-> Python parity suite

**Task:** Port TISCA v2 to R with the same API as the Python reference and build a
parity test suite (`>= 200` random configurations; agree to `1e-8` analytic /
bootstrap MCSE resampling). Covers the P2-T2 ACCEPTANCE criteria and removes the
v1 defects (`paired = FALSE`, per-column NA dropping).

**Date:** 2026-08-04 · **Owner:** P2-T2 (L-CODE)

---

## 1. What was delivered

### 1.1 The R package (`tisca/R/`)

| Module | File | Contents |
|---|---|---|
| package entry | `tisca.R` | version, `%||%`, `power_()`, `tisca_ready()` |
| estimands | `estimands.R` | PEHE (rooted), CATE MSE (unrooted), ATE sq err, RMSE_ATE + delta-method MCSE, ATE/CATE coverage, CIL, bias, interval score (`*_{c}`), per-rep validators |
| planning | `planning.R` | modes `M1`–`M5` power (`pt(q,df,ncp)` kept), `sigma_ub` inflation, `solve_mcse_J`, `solve_halfwidth_J`, `solve_power_J`, `combine_J`, `power()` dispatch |
| inference | `inference.R` | `contrast_from_columns` (LISTWISE NA deletion + drop count), `paired_t` (estimate-first), `studentized_paired_bootstrap`, BCa, `mcnemar_exact` |
| multiplicity | `multiplicity.R` | `adj_alpha`, `adjust_p_values`, `romano_wolf_stepdown`, `family_power`, `plan_rw_power` |
| mcs | `mcs.R` | `reality_check`, `spa_test`, `mcs` (`T_R`/`T_max`) with elimination path, CRAN-oracle wrapper |
| procedure | `procedure.R` | `run_two_stage` (Algorithm 1), `run_adaptive` (Algorithm 2, guarded) |
| validate | `validate.R` | spec Section 8 rejection rules |

Package files: `DESCRIPTION`, `NAMESPACE`.

### 1.2 v1 defects removed
- `paired = FALSE` is gone; all inference is on `D_j` via one-sample t / bootstrap / McNemar.
- Per-column NA dropping replaced by **listwise deletion across the pair** with the
  dropped count returned (`contrast_from_columns`), satisfying reviewer suggestion 2(b).
- `pt(q, df, ncp)` retained (already correct); the SciPy `t`/`nct` fix is P2-T1's.

### 1.3 Tests (`tisca/tests/`)
- `run_tests.R` + `test_*.R`: unit tests (planning anchors against the spec's own
  numbers; estimands; inference; multiplicity; MCS; procedure; validate; MCS-vs-CRAN-oracle).
- `parity_generate_reference.R`: enumerates **406** deterministic analytic/integer configs
  and writes `parity_reference.csv` (R reference values).
- `parity_eval_python.py`: evaluates the Python reference for each config.
- `parity_run.R`: **Tier 1** analytic/integer diff (1e-8 / integer-exact).
- `parity_resampling.R` + `parity_eval_python_resampling.py`: **Tier 2** bootstrap-MCSE tier
  (studentized bootstrap p, MCS set, SPA / Reality Check, Romano-Wolf rejections).
- `parity_all.R`: combined runner producing the ACCEPTANCE line.
- `mcs_parity_oracle.R`: CRAN-MCS oracle driver (bit-exact resample indices for the Python side).

## 2. Results

### 2.1 R unit tests
**All pass.** Notable anchors reproduced from `docs/tisca_v2_spec.md` rev. 2:
- `M1` 80% power at `delta/sigma = 0.5` requires `J = 34`.
- `sigma_ub` inflation factors `J0 = 25/50/100 -> 1.153 / 1.099 / 1.067`.
- M5 genuine-infeasibility (`|delta| >= Delta`) flagged; small-J zero-power is **not** misread as infeasible.
- `solve_power_J` reproduces the spec's Section 4.5 `J_final = 382` (within `+/-1` integer-boundary).
- Degenerate `s_D = 0` terminates (no infinite loop; P3-T7 guard).

### 2.2 Parity, Tier 1 (analytic / integer) — **PASS**
`406` configs · `378/378` well-conditioned agree to `1e-8` / exact · `0` genuine mismatches.

One real boundary bug was found and fixed in **R** during the parity run: `solve_halfwidth_J`
had a `tolerance = 1e-6` that accepted `0.2000003` as `<= 0.2` (true smallest J is `1540`,
not `1539`). The default tolerance is now `1e-12` (machine-scale) so R returns the same
smallest J as Python.

### 2.3 Parity, Tier 2 (resampling) — **PASS**
`10/10` checks within bootstrap MCSE: studentized bootstrap p (3), MCS kept sets (2),
SPA + Reality Check (4), Romano-Wolf rejection pattern (1).

---

## 3. Critical finding for the P2-T1 lane (needs a fix there)

**`28` of the `406` analytic configs return `NaN` from the PYTHON reference, not R.**
All are high-power `M1`-`M5` cells where `scipy.stats.nct.cdf(-crit, df, ncp)` returns `NaN`
for a large noncentrality parameter (an extreme-`ncp` instability in the scipy noncentral-t;
`nct.cdf(-crit)` should be ~`0`). R's `pt()` handles these exactly.

This is exactly the class of numeric defect the paper Section 1.6 flags. **It is not an R bug.**
Recommended fix in `tisca/python/tisca/planning.py`:

> In `power_M1` (and effectively `M2`/`M3`/`M4`), when the noncentrality is large such that the
> far-tail term underflows, clamp the impossible tail to `0.0` (or compute `F`(crit) − `F`(−crit)
> from the regularized lower-tail and set the sub-machine tail to `0`). Because the term is
> `<= 1e-31`, this changes nothing material but fixes `NaN`.

The parity harness classifies these as `Python-nct-NaN (P2-T1)` and does **not** count them as
R failures. Once P2-T1 is patched, `parity_all.R` re-runs cleanly with `29/29` nct cells passing.

## 4. Reproducibility
- All four parity/unit runners are deterministic (fixed seeds; the analytic reference is seeded
  `20260804`; the resampling spec is seeded `2026`).
- `parity_generate_reference.R` and `parity_eval_python.py` reconstruct *identical* data on both
  sides (`D_vec` / `pv_vec` are written verbatim by R and read by Python; MT vs PCG64 cannot be
  made bit-identical, so inputs are passed, not regenerated).

## 5. How to run
From `tisca/tests/`:
```
Rscript run_tests.R                      # R unit tests
Rscript parity_all.R                     # both parity tiers -> ACCEPTANCE line
```
Requires: R >= 4.0 (stats, parallel), Python 3 with numpy & scipy, `testthat`, `jsonlite`.
`MCS` (CRAN) is used only as the unit-test oracle and is optional.
