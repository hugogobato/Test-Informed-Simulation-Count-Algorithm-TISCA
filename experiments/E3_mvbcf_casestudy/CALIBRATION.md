# E3 — MVBCF case study: calibration and execution log

**Status:** artifacts staged; **execution in progress on the user's Colab sessions**.

## Prepared artifacts (this task, 2026-08-04)

| Artifact | Notes |
|---|---|
| `run_cell.R` | the P3-T5(a) driver, derived from `GitHub_DGP1.R` (see its header for the deviation list). Its stochtree BCF branch now uses the corrected paper-equivalent specification: `num_gfr = 0`, 50/20 trees, 500 burn-in and 500 posterior draws, fixed standardized forest variances 1 and 0.140625, and the paper-equivalent residual prior. Validated locally: parses; the DGP formula blocks, the L'Ecuyer-CMRG stream construction, and the 171-column schema pass `ALL_VALIDATION_PASSED`. **Note what that validation does and does not cover** (see the audit below): it exercises the schema and the seed machinery, never the model calls. |

## Stochtree BCF calibration decision, 2026-08-05

The paired Round 0 benchmark compares this corrected `stochtree::bcf` branch
with the original paper `bcf` implementation on the same generated data,
propensity scores, and replication seeds. The published Table 2 bands remain
reported as a diagnostic reference. The operational E3 choice is to retain
the corrected stochtree implementation when it converges and tracks the paper
implementation, even if its metrics do not fall inside those historical
bands. The shard notebooks therefore report the calibration verdict without
blocking Round 1; incomplete checkpoints and non-converged rows still stop
execution.

## Pre-execution audit, 2026-08-04 (before any Colab session was spent)

The driver was line-by-line diffed against the upstream `GitHub_DGP1.R` and
`MVBCF_Code.cpp`. Six defects were found and fixed. Recording them here because
four of them would have produced a full campaign of unusable shards, and none
would have been visible from the shard CSVs alone.

| # | Defect | Effect if run |
|---|---|---|
| 1 | `fast_bart()` called with 17 arguments; the signature (`MVBCF_Code.cpp:659`) takes 18. `X_mod_test` was omitted, shifting every later positional argument. | Every MVBCF fit aborts; `converged_flag = 0` on 100% of rows. The case study has no treatment model. |
| 2 | `predict(bartc_mod, type = "icate")` returns (draws × observations); the driver averaged it as (observations × draws) in all four places (`rowMeans`/`colMeans` swapped, `apply` margin 1 instead of 2, `rowSds` instead of `colSds`). | Every BART metric (PEHE, bias, CRPS, both coverages, both widths, ATE interval) computed on the wrong axis. Silent: 1000 test points and 500 draws both recycle without error. |
| 3 | Pilot shards: `seeds_emit <- s_start + 1000000L + 1L` collapses the index range to a **length-1 vector**. | Every Round 0 pilot shard runs 1 replication instead of 50. The P3-T5(e) calibration gate would have been decided on a single draw. |
| 4 | `colnames(X)` never set before `bartc()`; upstream sets `V1..Vp` at exactly that point so `predict(..., X_test)` can match columns by name. | `bartCause` matches positionally and warns, or errors, per replication. |
| 5 | All seven model fits reset the RNG to the identical stream state, so `sample.int()` returned the same seed for all of them and `model_seed_*` was constant. | Not a bias, but the recorded per-model seeds were meaningless and the two BCF outcome fits shared one stream. Replaced with fixed per-model L'Ecuyer substreams. |
| 6 | Propensity `dbarts::bart()` ran at 500 burn-in / 500 draws; upstream uses the package defaults (100/1000). | An undocumented deviation in the model that feeds every downstream fit, inside the very comparison the P3-T5(e) gate makes against the published Table 2. |

Two further changes, both deliberate deviations now documented in the driver:

* **The upstream `mvbart` CRPS is wrong and was being copied faithfully.**
  `GitHub_DGP1.R:262` uses `colSds(z_1_0_preds[,k,])` on a (test obs × draws)
  slice, which returns one sd *per draw* where `crps_norm` needs one *per test
  observation*; with `n_test = 1000` and `n_mcmc = 500` it silently recycles.
  Every other model in the upstream script supplies a per-observation sd. Since
  CRPS is promoted to the MCS loss for the uncertainty comparison (plan §1.7b,
  §2.3), the driver now uses `rowSds` and says so.
* **Dimension guards added** on all four model outputs (`stopifnot` on
  `n_test`/`n_mcmc`). Defect 2 was invisible precisely because both dimensions
  are plausible; the guard turns any future orientation change into a loud
  failure instead of a shard of quiet garbage.

Verified locally after the fixes: the script parses; a pilot shard now builds 50
seeds against 50 streams with labels `1000001..1000050`; pilot (master 2) and
confirmatory (master 1) streams have zero overlap; streams are shard-offset
invariant (index 30 is identical whether the shard starts at 0 or at 30); the
seven per-model substreams are distinct; and `stochtree::bcf`'s `tau_hat_test`
was confirmed by direct run to be (n_test × n_mcmc), so the BCF branch's
`rowMeans`/`colMeans` were already correct. `bartCause` and `skewBART` are not
installed locally, so defects 2 and 4 are corrected from the upstream script's
own usage rather than by execution — **the Round 0 pilot must confirm them**,
and the new `stopifnot` guards make that automatic.
| `ANALYSIS_PLAN.md` | the P3-T5(b) **pre-registration**; committed before confirmatory analysis. Freezes family, two-sided α, δ (relative to `sd(τ(X))`), the Romano–Wolf family procedure, and the `stochtree::bcf` calibration band. |
| `notebooks/E3_mvbcf_shard.ipynb` | P3-T5(d) parameterised worker (generator `build_e3_notebooks.py`); edit Cell 2 per session. |
| `notebooks/E3_round0_pilot_calibration.ipynb` | P3-T5(e) Round 0 pilot (50 independent seeds) + the calibration gate vs McJames et al. Table 2. |

Per the revision plan, the single number that governs the schedule is the Colab
2-vCPU speedup. It was measured once on a short `dbarts` proxy; the plan asks to
**re-confirm on the real workload during the Round 0 pilots** (they already are
the real workload — just log per-replication wall-clock). Keep the following
records here as they are produced:

- [x] P0-T3 result (recorded in `REVISION_PLAN.md` §P0-T3): `SPEEDUP = 1.62` on
      a `dbarts` proxy (4 fits: 4.5 s serial vs 2.8 s at `mc.cores = 2`) →
      `11.0` usable core-hours per 8-hour notebook. Decision: `mc.cores = 2`.
- [ ] Round 0 pilots: **re-measured** throughput and peak RAM on the real
      workload, per cell. Record CPU model, `nproc`, R version, and wall-clock
      per replication here (IJDA #14a requires this in the paper).
- [ ] P3-T5(e) `stochtree::bcf` calibration gate against the published Table 2
      (BCF PEHE Y1/Y2, BCF tau-95 coverage Y1/Y2): 50-replication means, SEs,
      and verdict.
- [ ] P3-T5(f) `fast_bart()` vs `mvbcf::run_mvbcf()` equivalence check.
- [ ] P0-T4 library bundle: exact `sessionInfo()`, the `DEPENDENCIES.csv`
      including a non-`NA` `remote_sha` for `skewBART` and `mvbcf`, the
      snapshotted `renv.lock`, the published tarball URL, and the timing of the
      fresh-account restore (the ACCEPTANCE test must run on an account that did
      **not** build the bundle).
- [ ] P0-T2 acceptance: the four-way identity test of `docs/seed_rng_protocol.md`
      §3.5 (`mc.cores` 1 vs 2, shard-aligned vs shard-offset) passing on the real
      driver. **P0-T2 is not accepted until this is green.**
- [ ] P0-T2 §3.6: which of `stochtree::bcf`, `dbarts`, `bartCause` and
      `fast_bart()` actually honour an explicit seed (fit each twice on one
      replication and compare). Any model that does not is a reproducibility
      hole that has to be closed before Round 1.

Relevant plan sections: §P0-T3, §P0-T4, §P3-T5(e,f). Acceptances: the
`stochtree` gate passes before any confirmatory replication runs.
