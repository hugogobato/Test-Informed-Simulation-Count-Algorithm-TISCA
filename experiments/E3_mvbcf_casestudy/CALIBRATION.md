# E3 — MVBCF case study: calibration and execution log

**Status:** artifacts staged; **execution in progress on the user's Colab sessions**.

## Prepared artifacts (this task, 2026-08-04)

| Artifact | Notes |
|---|---|
| `run_cell.R` | the P3-T5(a) driver, derived from `GitHub_DGP1.R` (see its header for the deviation list). Validated locally: parses; the DGP formula blocks, the L'Ecuyer-CMRG stream construction, and the 170-column schema pass `ALL_VALIDATION_PASSED`. |
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
