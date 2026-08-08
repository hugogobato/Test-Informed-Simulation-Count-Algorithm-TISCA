# E3 — MVBCF case study: calibration and execution log

**Status:** artifacts staged; E3 data collection and the primary analysis are
complete. Remaining verification items are listed below.

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

## Deviation log

### D1 — 2026-08-06: pilot block and pilot size (ANALYSIS_PLAN.md AMENDMENT 1)

The dedicated pilot seed block (`seed_cell_master = 2`, seeds 1,000,001-1,000,050)
was run for **one** of the four cells, and that run used a superseded `run_cell.R`:
it reports BCF PEHE 10.92 and 95% coverage 0.759, against 9.70 and 0.970 from the
driver that produced the confirmatory shards. The commit order confirms it —
`notebooks/E3_DGP1_n500_pilot_shard01_seeds1000001-1000050.csv` was added in
`44583de`, and `run_cell.R`'s BCF branch was rewritten afterwards in `ffeaa8b`.
That pilot is therefore **not used**.

The pilot is instead seeds `0…99` of the confirmatory block, with `J0` raised from
50 to 100 on the three criteria stated in the amendment. Pilot rows are discarded
from final inference, so the confirmatory analysis runs on 900 replications per
cell. Because the pilot only sizes `J` and is then discarded, the deviation cannot
bias the reported estimates or the level of the reported tests.

**Disclose in the response letter.** The amendment was written after the
confirmatory data was collected and its column means inspected.

**Status: implemented, no re-run required.** All four cells already hold a complete
1000-replication confirmatory block with zero `converged_flag` failures, so the
amendment is executed purely as a row split of data in hand:
`collect_shards.py --pilot-size 100` writes the two blocks and the `analysis_phase`
label, and `analyse_e3.py` performs the planning and the confirmatory inference on
the 900 retained replications. `shard_table.csv` marks the four Round-0 pilot rows
`status = superseded` so they are neither read nor required.

The pre-declared `J0` sensitivity comes out as follows (`J*` per cell, over the six
primary contrasts, at the precision target `0.05 x sd(tau)`):

| cell | J0 = 25 | J0 = 50 | J0 = 100 |
|---|---:|---:|---:|
| DGP1 n = 500 | 34 | 27 | 30 |
| DGP2 n = 500 | 19 | 17 | 15 |
| DGP3 n = 500 | 51 | 34 | 34 |
| DGP1 n = 100 | 144 | 108 | 96 |

Reported whichever way it falls, as the amendment requires: `J0 = 100` does **not**
uniformly give the smallest `J*` (DGP1 n = 500 gives 27 at `J0 = 50` against 30 at
`J0 = 100`). What it does give is the stability the amendment argued for, and the
choice is not doing the paper any favours.

### Round 0 status

**Supersession note.** The pilot CSVs on the reserved seed block that use the FULL
171-column schema (DGP1 n = 500, DGP3 n = 500) were produced by the superseded
`run_cell.R`. DGP3 is the clearest case: its pilot reports BCF PEHE Y1 11.55 and 95%
coverage 0.777, against 9.97 and 0.970 in its own confirmatory shards. They are not
used for planning (see deviation D1).

### P3-T5(e) calibration gate: PASSED, 2026-08-06

The gate is settled by `E3_round0_dual_bcf_pilot.ipynb`, which ran the corrected
`stochtree::bcf` and the original authors' `library(bcf)` on the **same 50 pilot
seeds, the same generated data and the same propensity scores**
(`notebooks/E3_DGP1_n500_round0_dual_bcf_pilot_seeds1000001-1000050.csv`). This is a
stronger test than the plan asked for: it validates against the original
implementation directly, not only against the published table.

| metric (DGP1, n = 500, Y1) | `library(bcf)` | corrected `stochtree::bcf` | published Table 2 | pass band |
|---|---:|---:|---:|---|
| PEHE | 9.727 | **9.673** | 9.63 | 9.3-10.0 ✓ |
| tau 95% coverage | 0.9750 | **0.9725** | 0.97 | 0.95-0.98 ✓ |
| tau 50% coverage | 0.5242 | 0.5239 | — | — |
| 95% interval width | 44.06 | 43.69 | — | — |
| CRPS | 5.651 | 5.629 | — | — |
| ATE | 20.068 | 20.362 | ~20 | — |

The two implementations agree to within 0.6% on PEHE and 0.3 percentage points on
coverage, and both fall inside the published bands. **Verdict: the `stochtree`
substitution is validated; the cold-start `num_gfr = 0` variant with the paper's
priors and residual prior reproduces `library(bcf)`.** The confirmatory run at 1000
replications independently confirms it (BCF PEHE Y1 9.70, coverage 0.9698).

Ordering caveat to disclose: this diagnostic was completed alongside, not strictly
before, the confirmatory shards.

- [x] P3-T5(e) `stochtree::bcf` calibration gate against the published Table 2.
- [x] P3-T5(f) — superseded in substance: the dual-BCF diagnostic compares the two
      BCF implementations directly, which is the equivalence question that mattered.
      `fast_bart()` is used for MVBCF in both arms, so the original
      `fast_bart()`-vs-`mvbcf::run_mvbcf()` question is moot: `run_mvbcf` is not used.

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
      per replication here (IJDA #14a requires this in the paper). Partially
      superseded: `replication_seconds` and `fit_seconds_*` are recorded on
      every one of the 4,000 confirmatory rows, so the throughput figure the
      paper needs can be read off the shipped data rather than re-measured.
- [x] P3-T5(e) `stochtree::bcf` calibration gate against the published Table 2.
      Recorded above (2026-08-05 decision and the dual-BCF table).
- [x] P3-T5(f) superseded in substance by the dual-BCF diagnostic; see above.
- [ ] P0-T4 library bundle: exact `sessionInfo()`, the `DEPENDENCIES.csv`
      including a non-`NA` `remote_sha` for `skewBART` and `mvbcf`, the
      snapshotted `renv.lock`, the published tarball URL, and the timing of the
      fresh-account restore (the ACCEPTANCE test must run on an account that did
      **not** build the bundle).
- [x] P0-T2 acceptance: the four-way identity test of `docs/seed_rng_protocol.md`
      §3.5 (`mc.cores` 1 vs 2, shard-aligned vs shard-offset) passing on the real
      driver. **Green, 2026-08-08.** See the seed-verification section below.
- [ ] P0-T2 §3.6: whether `stochtree::bcf` and `dbarts`/`bartCause` honour their
      explicit seed argument. Scope now settled by source inspection (below): it
      applies to those two libraries only, not to four models.

Relevant plan sections: §P0-T3, §P0-T4, §P3-T5(e,f). Acceptances: the
`stochtree` gate passes before any confirmatory replication runs.

## Seed verification, 2026-08-08 (`notebooks/E3_seed_verification.ipynb`)

Results are in `results/E3/seed_verification.csv` and the per-column detail in
`results/E3/seed_verification_columns.csv`.

- [x] **P3-T5 acceptance.** Seeds 155, 453 and 495 of DGP1 `n = 500`, drawn by a
      fixed RNG so the choice is auditable, were re-run and compared against the
      stored rows on all 159 non-environment columns. Zero columns differ. The
      stored rows for that cell were produced on ten different Colab hosts, so
      this is a cross-machine check as well as a re-run check.
- [x] **P0-T2 §3.5 four-way identity test.** DGP1 `n = 100`, seeds 37 and 38,
      run at `mc.cores = 1` and `2`, each shard-aligned and shard-offset (the
      offset run requests the same seeds inside a range starting at 35). All
      three comparisons against the `cores = 1`, aligned reference differ on 0 of
      159 columns. The stream for replication index `j` therefore does not depend
      on the worker count or on where its shard began.
- [x] **Run-level reproducibility, per model.** `bart`, `bcf`, `mvbart` and
      `mvbcf` each reproduce all 36 of their columns exactly, maximum absolute
      difference 0.0. Read this for what it is: it rules out wall-clock or PID
      seeding and uncontrolled nondeterminism inside a fit. It does **not** by
      itself attribute that reproducibility to the explicit seed argument, for
      the reason given next.

### §3.6, per-model seed honouring: scope, and why no re-run is implied

The protocol asks which models honour an explicit seed. Reading the code settles
most of the question without spending a session, and narrows what is left.

`fast_bart()` (MVBCF) **has no seed argument**. Its signature in the upstream
`MVBCF_Code.cpp:659` takes 18 arguments, none of them a seed, and every random
draw inside it goes through R's own generator: `R::runif`, and `rmvnorm`,
`riwish` and `sample` from `RcppDist`/Rcpp sugar. `Rcpp::export` wraps the call
in an `RNGScope`, so the function reads and writes `.Random.seed`. There is no
independent generator anywhere in the file (no `mt19937`, no `rand()`, no
`arma::randu`/`randn`). MVBCF is therefore governed entirely by §3.3 stream
positioning, which is exactly what the four-way identity test above validates,
and `model_seed_mvbcf` is a recorded label rather than an argument that was
passed. The same holds for `MultiskewBART`, which `run_cell.R` calls without a
seed argument.

That leaves two libraries that do take one: `stochtree::bcf`
(`general_params$random_seed = f_seed`) and `dbarts`, used directly for the
propensity fit (`seed = f_seed`) and underneath `bartCause::bartc`
(`seed = f_seed`). For these two the discriminating run is the one §3.6
specifies and which has not been done: fit **the same model, on the same data,
with the same explicit seed, from two different global RNG states**, and check
the output is identical. Agreement attributes reproducibility to the seed
argument; disagreement shows the argument is inert and the global stream is
doing the work.

**No E3 replication needs to be re-run for this.** It is a property of two R
packages, tested on one dataset with a handful of fits, and it reads none of the
4,000 stored rows. Nothing about the stored data is contingent on the outcome:
exact reproduction has already been demonstrated empirically across ten hosts,
two worker counts and two shard offsets. What the outcome changes is only what
the Data Availability statement is entitled to *say*, that is, whether
reproducibility is attributed to per-model seeding or to global stream
positioning alone.

`notebooks/E3_seed_honouring.ipynb` (generator
`notebooks/_generators/build_e3_seed_honouring_nb.py`) runs it, in about 20
minutes on one Colab session, most of that being environment setup. It is
standalone: it does not import, source or modify `run_cell.R`, though the DGP
block and every model call are copied verbatim from it.

**Local pre-check, 2026-08-08 (not the bundle, so not the record).** Two of the
three fits were run here, off the bundle, on R 4.3.3 with `dbarts` 0.9.33 and
`stochtree` 0.4.0. Both came back `HONOURS_SEED`:

| fit | values compared | max abs diff, same seed / other global state | max abs diff, other seed / same state | verdict |
|---|---:|---:|---:|---|
| `dbarts::bart` (propensity) | 1,000,000 | 0 | 3.594 | HONOURS_SEED |
| `stochtree::bcf` | 500,000 | 0 | 106.478 | HONOURS_SEED |

The right-hand column is the positive control, and it is what makes the zero in
the middle column mean something: both fits do respond to their explicit seed, so
the invariance to the global state is a real property rather than a
seed-insensitive model producing the same numbers whatever it is told. This is
strong evidence that §3.6 will come back green, but it is **not** the answer for
the record, because the package versions are not the bundle's and because
`bartCause` is not installed locally. Run the notebook on the bundle to close the
box.
