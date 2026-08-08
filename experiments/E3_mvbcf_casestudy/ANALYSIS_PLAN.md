# E3 MVBCF case-study: pre-registered analysis plan

**Status: PRE-REGISTERED.** This document is committed to the repository **before**
the confirmatory replications are analysed (plan P3-T5[b]). Analysis decisions are
frozen as written here; a change after results are seen invalidates the
pre-registration and must be logged as a deviation in `CALIBRATION.md` and the
response letter. It answers the reviewer's question (IJDA minor #10) about whether
the comparisons were pre-specified.

**Scope of the plan:** the `run_cell.R` re-run of DGP1/2/3 at n = 500 and DGP1 at
n = 100, 1000 replications per cell, plus the 50-replication pilots. The recorded
columns and their definitions are fixed in `run_cell.R` (see the header and the
`COL_NAMES` schema, 171 columns, including per-replication `replication_seconds`).

---

## AMENDMENT 1 — 2026-08-06: pilot block and pilot size

**This amendment is made after the confirmatory replications were collected and
after their column means were inspected. It is logged here, in `CALIBRATION.md`, and
must be disclosed in the response letter.** It changes only how `J*` is *computed*;
it changes no hypothesis, no metric, no family, no correction and no margin.

**Implemented in code, and no cell is re-run.** The amendment is a partition of rows
that already exist, not a new data-generation instruction, so `run_cell.R` is
unchanged: it remains the generator that produced seeds `0…999` per cell from cell
master 1. The split is applied at analysis time by
`collect_shards.py --pilot-size 100`, which writes
`results/E3/DGP{d}_n{n}_amended_pilot_replications.csv` and
`…_amended_confirmatory_replications.csv` and labels every row with an
`analysis_phase` column alongside the generator's own `seq_phase`; the reanalysis
itself is `analyse_e3.py`. The four Round-0 rows in `shard_table.csv` now carry
`status = superseded`, so the collector no longer requires the reserved master-2
pilot files that item 1 below rules out.

**What changed, and why.**

1. **The dedicated pilot seed block was not run for three of the four cells.** §5 as
   originally written reserved seeds `1_000_001…1_000_050` (`seed_cell_master = 2`)
   for the pilot. A pilot at those seeds exists only for DGP1 n = 500, and it was
   produced by a **superseded** `run_cell.R`: it reports BCF PEHE 10.92 and 95%
   coverage 0.759, against 9.70 and 0.970 from the driver that produced the
   confirmatory data. It is therefore not a valid pilot for this analysis and is not
   used.

2. **The pilot is instead the first block of the confirmatory seeds.** Seeds
   `0…J0-1` (`seed_cell_master = 1`) form the pilot; seeds `J0…999` form the
   confirmatory set. The two blocks are disjoint replications drawn from disjoint
   L'Ecuyer-CMRG substreams, so they are independent samples, which is the property
   Algorithm 1 actually requires — the spec asks for a pilot independent of the
   confirmatory data, not for a pilot in a particular seed block. **The pilot rows
   are discarded from all final inference**, so the reported contrasts are computed
   on `1000 - J0` replications. Retaining them would convert the design into the
   internal-pilot variant D2, which the operating-characteristic study measures at
   +0.0033 on the unconditional level (`docs/phase3_findings.md`, F2).

   Because the pilot is used only to size `J` and is then discarded, this deviation
   cannot bias the reported estimates or inflate the level of the reported tests.

3. **`J0` is raised from 50 to 100.** Justified on three grounds declared here, none
   of which refers to the resulting `J*`:

   * **Stability of the selected `J`.** Module C of the operating-characteristic
     study (5,000 outer repetitions per cell, design D4) measures
     `sd(J)` = 20.4 / 12.9 / **8.6** and `q95(J)` = 98 / 79 / **68** at
     `J0` = 25 / 50 / 100, with achieved power converging onto the 0.80 target from
     above (0.878 / 0.869 / 0.858) rather than overshooting it. Larger pilots do not
     buy power; they buy a budget that is predictable, which is what a
     replication-planning method should optimise.
   * **The family is six-dimensional.** §3 plans Romano-Wolf and joint power from
     `Σ̂_D`, a 6×6 covariance over C1–C6. Estimating a 6×6 covariance from 50
     observations is thin; 100 is defensible.
   * **Relative cost.** The pilot is discarded, so `J0` is a direct cost. At the
     `J*` scale this case study operates at (order 200–300), 50 extra discarded
     replications are roughly 14% of total spend, and all 1000 are already in hand.

   **Declared before computing `J*`:** the choice of `J0 = 100` is fixed by the
   three criteria above and will not be revisited after seeing `J*`. For
   transparency, `J*` at `J0 = 50` and at `J0 = 25` is reported alongside as a
   sensitivity row, and the direction of any difference between them is reported
   whether or not it favours the paper.

**What did not change:** the six primary contrasts C1–C6, mode M1 two-sided, MVBCF
as control, the Bonferroni-adjusted planning `α = 0.05/6`, Romano-Wolf as the family
procedure with Bonferroni and Holm as sensitivity, `δ = 0.25 · sd(τ(X_1,test))`, the
variance-uncertainty inflation at `γ = 0.20`, the coverage framing of §2, and the
MCS layer of §3.

---

## 1. Primary family and hypotheses

The proposed method is **MVBCF**. The re-run fits four models on identical data in
every replication: `mvbcf` (proposed), `bcf` (single calibrated benchmark),
`bart` (univariate BART), `mvbart` (multivariate skew-BART).

The **primary family** is the set of paired, two-sided comparisons of MVBCF against
each benchmark on the headline loss metrics, defined on paired differences
`D_j = L_MVBCF,j − L_benchmark,j` across common replications `j` (paired structure,
plan §2.1/§2.3, reviewer IJDA #2):

| # | Contrast | Metric | Family role |
|---|---|---|---|
| C1 | MVBCF vs BCF | PEHE Y1 | primary (headline) |
| C2 | MVBCF vs BCF | PEHE Y2 | primary |
| C3 | MVBCF vs BART | PEHE Y1 | primary |
| C4 | MVBCF vs BART | PEHE Y2 | primary |
| C5 | MVBCF vs MVBART | PEHE Y1 | primary |
| C6 | MVBCF vs MVBART | PEHE Y2 | primary |
| C7 | MVBCF vs BCF | 95% CATE coverage (deviation from 0.95), Y1/Y2 | secondary, calibration |

**Hypothesis mode:** M1, two-sided equality `H0: θ = 0` vs `H1: θ ≠ 0`, where `θ = E[D_j]`
is the population mean paired difference. One-sample paired t on `D_j` (and, for robustness,
the studentized paired bootstrap) as the test statistic; two-sided critical value at the
multiplicity-adjusted `α` (below). Sidedness of the test equals sidedness used in planning.

**Control method:** MVBCF. There is no prior expectation of direction, so two-sided is used
throughout; the plan does not pre-select a "better" direction.

## 2. Coverage hypothesis (calibration framing)

Coverage is **not** analysed as "higher is better" (plan C6, reviewer IJDA #10). The primary
coverage estimand is the deviation from the nominal level, `|cov_j − 0.95|`, compared between
MVBCF and each benchmark with a two-sided / minimum-effect framing and an equivalence band
`Δ_prac` set below. Reported alongside are the 95% coverage itself (descriptive) and the
averaged credible-interval width (sharpness). The **interval score** (Winkler) is the scalar
uncertainty loss used in the Model Confidence Set layer, not raw coverage.

## 3. Multiplicity procedure

- **Family definition:** the six primary comparisons C1–C6 form one family. Secondary
  coverage comparisons (C7) form a second family and are not familyed with the primary tests.
- **Correction:** the primary family is evaluated with **Romano–Wolf stepdown** (FWER control,
  exploiting the joint dependence of the paired differences), with Bonferroni and Holm reported
  as a sensitivity table. The planning `α` for the precision/power criterion uses the
  Bonferroni-adjusted level `α/K` for the six primary contrasts per plan 3.5.
- **Family-level success criterion:** conjunctive (all six primary comparisons reported with
  family-adjusted p-values; the headline is MVBCF vs BCF on PEHE). Marginal and family-level
  power are both reported, per plan §2.2/3.5.
- **Model Confidence Set (MCS):** reported over the four models using the **interval score**
  and (separately) PEHE as scalar losses, with the §2.3 caveats stated verbatim in the
  manuscript. MCS is an **inference layer, not a stopping rule**; no `J` is chosen from MCS
  cardinality.

## 4. Design targets

J is chosen from the two-layer protocol (plan §2.1).

- **Precision layer (default, always active):** smallest `J` such that the Monte Carlo
  standard error of the paired contrast, `σ̂_D/√J`, is at or below a target `m`, and/or the
  CI half-width `t_{1−α/2,J−1} σ̂_D/√J ≤ h`. For the headline PEHE comparisons, an absolute MCSE
  target is set relative to the scale of the marginal SD of `τ(X_i)` (below), not as a bare unit.
- **Decision layer (optional, active because a confirmatory comparative claim is made):**
  smallest `J` with paired two-sided power `≥ 1 − β = 0.80` at the planning alternative `δ`,
  computed with the **noncentral t with J−1 dof** at the **Bonferroni-adjusted** `α = 0.05/6`.

**Planning alternative δ (and its domain justification, IJDA #8a / minor #11):** a PEHE
contrast `δ` is *not* justified by the bare fact that the true ATE ≈ 20. The meaningful
scale is the **dispersion of the unit-level treatment effect `τ(X_i)`**, because PEHE is an
error on the CATE surface, not on the marginal ATE. Accordingly `δ` is expressed relative to
`sd(τ(X_i))` (computed from the DGP over the test set). The planning alternative for the
headline comparison is **`δ = 0.25 · sd(τ(X_1,test))`** (a quarter of the treatment-effect
cross-sectional SD), declared here. The equal-marginal-variances Cohen's d on the *unpaired*
marginal distributions is rejected as the planning metric because the design is paired
(IJDA minor #12).

**Variance-uncertainty inflation (IJDA #8b):** `σ̂_D` is estimated from the independent pilot
of size J0; the planned J uses the upper-confidence-bound inflation
`σ̂_D × sqrt((J0−1)/χ²_{γ,J0−1})` with `γ = 0.20` so that planning is not systematically
under-sized.

## 5. Pilot and confirmatory protocol (two-stage default)

- **Pilot:** ~~the first 50 replications of each cell use the **independent pilot seed block**
  (`seed_cell_master = 2`, emitted seeds `1_000_001…1_000_050`)~~ **superseded by
  AMENDMENT 1 (2026-08-06):** the pilot is seeds `0…99` at `seed_cell_master = 1`
  (`J0 = 100`). Pilot rows are used **only** to estimate `σ̂_D`, `Σ̂_D` and the DGP's
  `sd(τ(X))`; they are **discarded from the final inference** (plan P1-T2 §5).
- **Confirmatory:** ~~`J` replications at `seed_cell_master = 1`, seeds `0…J−1`.~~
  **superseded:** `J` replications at seeds `100…99+J`, i.e. the pilot block is
  excluded, so the confirmatory set is drawn from the remaining 900. `J_final = max`
  over the six primary contrasts of `max(J_precision, J_power)`, up to the declared budget cap
  `J_max = 1000` (the run on every cell is 1000 by author decision; J* is reported, not imposed).
- **Report order:** every contrast reports the **estimate and its Monte Carlo CI first, the
  p-value second** (IJDA #1e).

## 6. Why PEHE-vs-MVBART was omitted in v1 (or included now)

In the v1 paper the headline uncertainty and precision narrative compared MVBCF to MVBART on
**coverage**, and the PEHE comparisons were drawn against `bcf_*` / `wsbcf_*` (the benchmark
that reproduced the original PEHE). PEHE-vs-MVBART was *not* a headline v1 claim; it was a
supporting comparison. In this revision the PEHE-vs-MVBART contrast is **included and promoted
to a primary contrast (C5, C6)** for three reasons: (1) consistency (with the single
calibrated-BCF benchmark fix (§1.3) the four-model family is compared on every metric, and
omitting one model on one metric revives the very "benchmark selection" pathology reviewer #7
objects to; (2) MVBART, like MVBCF, is a multivariate model, so PEHE-vs-MVBART is the most
direct test of the multivariate-sharing claim; and (3) the model count is small (4), so the
multiplicity cost of including it is minor and is already absorbed by the protected family.

## 7. Fixed analysis code and outputs

The analysis reads the committed shard CSVs, asserts seed completeness (no gaps, no
duplicates, expected row counts per cell), and produces `results/E3/<cell>_replications.csv`
per cell plus:
- paired-contrast estimates, Monte Carlo CIs and p-values (`results/E3/paired_contrasts.csv`);
- planning table: J_precision, J_power, J_final per contrast, with the pilot `σ̂_D` (CI-inflated);
- family inference (Romano–Wolf, Bonferroni, Holm) and the MCS;
- calibration (deviation from 0.95, interval score) analysis;
- the `stochtree::bcf` calibration-gate comparison to McJames et al. Table 2 (P3-T5[e]).

## 8. Deviations register

Any post-hoc change to this plan (metric, family, α, δ, correction, or analysis method) is
recorded here with date, the change, and the reason. A change made after seeing results is
reported as a deviation, not silently absorbed.

| Date | Section | Change | Reason |
|---|---|---|---|
| _2026-08-04_ | — | plan written and committed before confirmatory analysis | — |

## 9. P3-T5(e) calibration gate (runs inside the pilot)

The `stochtree::bcf` benchmark is validated against McJames et al.'s published DGP1/n=500
Table 2 on the **50-pilot replications** before any confirmatory shard is launched:

| Metric | Target (published) | Pass band |
|---|---:|---:|
| BCF PEHE, Y1 | 9.63 ± 0.16 | 9.3–10.0 |
| BCF PEHE, Y2 | 9.96 ± 0.16 | 9.6–10.3 |
| BCF τ 95% coverage, Y1 | 0.97 | 0.95–0.98 |
| BCF τ 95% coverage, Y2 | 0.96 | 0.94–0.98 |

Coverage in band with PEHE still ~10% high triggers the diagnosis path (prior scaling,
`include_pi`/propensity handling, muscale/tauscale analogues); coverage near 0.85 means the
chain is not converging and `num_gfr`/`num_mcmc` must be increased. **No confirmatory
replications start until this passes** (gate acceptance).

---

_Authors of the analysis:_ H. G. Souto, F. Louzada Neto. This plan commits the analytical
decisions; it is not a statement about the results, which have not been run yet at the time
of writing.
