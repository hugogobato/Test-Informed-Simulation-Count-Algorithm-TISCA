# E3 — MVBCF case study: calibration log

**Status:** placeholder for the execution records that Phase 0/1 fill in.

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
- [ ] P0-T4 library bundle: exact `sessionInfo()`, the published tarball URL,
      and any version-pin deviations.

Relevant plan sections: §P0-T3, §P0-T4, §P3-T5(e,f). Acceptances: the
`stochtree` gate passes before any confirmatory replication runs.
