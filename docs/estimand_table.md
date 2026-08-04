# P1-T1 — Estimand table for TISCA v2

**Status:** final for review (L-SPEC, task P1-T1)
**Revision:** 2026-08-04 (rev. 2). Corrections in rev. 2, all backed by `notebooks/P1_math_verification.ipynb`: row 3b's delta-method MCSE was missing the inner square root; row 8's calibration estimand was defined per replication, which confounds miscalibration with Monte Carlo dispersion (§3.1); rows 4/5/6/8/9 are now parameterised by nominal level so the 50% columns the case study records are not discarded (§3.2).
**Files it governs:** the estimates/inference layer of the `tisca/` package (`estimands.py`, and the parallel R port), the manuscript equations (1)–(7), and every planned paired contrast in P4-T1 and the E3 case study.
**Acceptance (restated from the plan):** one row per supported performance measure; every equation in Section 2 of the manuscript maps to exactly one row; a second reader can implement any row from the table alone.

This document is deliberately notation-heavy and self-contained. Conventions established here are reused verbatim by `tisca_v2_spec.md` (P1-T2) and by the manuscript Section 3.

---

## 1. Notation and the four explicit conventions the table relies on

1. **Unit index `i`, replication index `j`.** `i` runs over the units in a single simulated/test dataset (typically `i = 1, …, n`); `j` runs over replications or Monte Carlo datasets (`j = 1, …, J`). The per-unit CATE estimates are written `\hat{τ}_i` and the true per-unit CATE `τ_i`. The per-replication ATE estimates are written `\hat{τ}_j` (or `\overline{τ}^{(j)}` where the superscript is needed to avoid collision) and the true per-replication ATE `τ_j` (equivalently `\overline{τ}^{(j)}`).

   This disambiguates the current manuscript, which reuses `τ_i` and `τ_j` across equations (1), (3), (5), (7) in an inconsistent way (SoutoNeto §1, IJDA minor #1). Adopt: **`τ_i` = unit-level CATE, `τ_j` = replication-level ATE.**

2. **Rooted versus unrooted root-mean-square error (IJDA minor #2).** Let `Q_j = mean_i (\hat{τ}_{j,i} − τ_{j,i})²` be the replication-level mean squared CATE error. Then
   - the **rooted** quantity `PEHE_j = sqrt(Q_j)` is a per-replication loss in the natural units of `τ`; its estimand is `E[PEHE_j]`;
   - the **unrooted** mean squared error has as its natural estimand `E[Q_j]`.
   These two are **different estimands**. In general **`E[sqrt(Q)] ≠ sqrt(E[Q])`** (Jensen, strict for a non-degenerate `Q`). MVBART / BCF / MVBART coverage and PEHE tables in the paper and in the original authors' code report the rooted value `PEHE_j`; the revision records **both** (see §3 rows 1–2) and reports the estimator for each separately. Whenever a headline is written as a single number per method, it is the rooted `PEHE` estimand `E[PEHE_j]` unless stated otherwise.

3. **Paired design is the default (IJDA #2, R2 ¶1, SoutoNeto §3).** All contrasts are formed on common replications: replication `j` must produce a value for every metric of every model being compared to `j`'s peers, so that any two methods `A`, `B` yield a paired, same-`j` pair of column values. Inference is on `D_j = L_{A,j} − L_{B,j}`. NAs are dropped **listwise across the pair** (a replication `j` is dropped from the contrast `(A,B)` if either `L_{A,j}` or `L_{B,j}` is missing), and the number of dropped replications is reported. Per-column NA dropping is forbidden (paper §1.6 defect).

4. **Estimator notation.** `\widehat{\theta}` = sample estimator of the estimand `θ`. `\mbox{MCSE}(\widehat{\theta})` = Monte Carlo standard error of that estimator across the `J` replications actually run.

---

## 2. Generic across-replication estimand and its MCSE

Let a replication-level quantity be `L_j`, `j = 1, …, J`, i.i.d. with unknown mean `θ = E[L_j]` and unknown variance `σ_L² = Var(L_j)`. The natural estimator is the sample mean
`\widehat{\theta} = \overline{L} = (1/J) ∑_{j=1}^{J} L_j`,
with Monte Carlo standard error
```
MCSE( L̄ ) = s_L / sqrt(J),   s_L² = (1/(J−1)) ∑_j (L_j − L̄)².
```
A `(1−α)` normal-approximation CI for `θ` is `L̄ ± t_{1−α/2, J−1} · s_L / sqrt(J)`. This generic row anchors every "expectation of a replicate-level loss" row below. Use the **studentized paired bootstrap** rather than the CLT CI whenever the third absolute moment of `L_j` is large relative to `s_L³` (see P3-T4 / the Berry–Esseen-driven operational rule); the table still reports the analytic MCSE because it is cheap and exact for the variance regardless of distribution.

For every row, a paired contrast `D_j = L_{A,j} − L_{B,j}` inherits all of the above with `θ_D = E[D_j]`, `σ_D² = Var(D_j)`, and `MCSE = s_D/sqrt(J)`.

---

## 3. The table

Legend for the `Admissible test` column: `paired-t` = one-sample t on `D_j` (with `J−1` df); `paired-boot` = studentized paired bootstrap on `D_j`; `McNemar` = exact/binomial McNemar on the paired binary indicators (only for 0/1 `L_j`); `describe` = reported as a point estimate + CI, no hypothesis test (the estimand is not a scalar to be tested). `δ_j` = per-replication deviation used in the calibration row.

| # | Measure | Unit-level quantity | Replicate-level `L_j` | Across-replication estimand `θ` | Estimator `\widehat{θ}` | MCSE formula | Paired contrast `D_j` | Admissible test | Notes / dependencies |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **PEHE (rooted)** | `(\hat{τ}_{j,i} − τ_{j,i})²`, summed by `i` | `PEHE_j = sqrt( mean_i (\hat{τ}_{j,i} − τ_{j,i})² )` | `θ_1 = E[PEHE_j]` | `mean(PEHE_j)` | `s_{PEHE}/sqrt(J)` | `D_j = PEHE_{A,j} − PEHE_{B,j}` | paired-t, paired-boot | grounded in the manuscript Eq. (2) with `τ_i` (unit CATE). Declared rooted convention (row 1 of §1). |
| 2 | **CATE MSE (unrooted)** | `(\hat{τ}_{j,i} − τ_{j,i})²` | `Q_j = mean_i (\hat{τ}_{j,i} − τ_{j,i})²` | `θ_2 = E[Q_j]` | `mean(Q_j)` | `s_Q/sqrt(J)` | `D_j = Q_{A,j} − Q_{B,j}` | paired-t, paired-boot | `E[sqrt(Q)] ≠ sqrt(E[Q])` (see §1.2). Recorded to make that distinction explicit. |
| 3a | **ATE squared error** | `( \hat{τ}_j − τ_j )²` for the ATE `τ_j` of replication `j` | `Q_{ATE,j} = ( \hat{τ}_j − τ_j )²` | `θ_{3a} = E[Q_{ATE,j}]` (MSE) | `mean(Q_{ATE,j})` | `s_{Q_{ATE}}/sqrt(J)` | `D_j = Q_{ATE,A,j} − Q_{ATE,B,j}` | paired-t, paired-boot | IJDA #11 core point: a single-replication "per-replication RMSE_ATE" is **not** a well-defined scalar estimand; report the MSE `E[Q_{ATE}]` and separately the RMSE `sqrt(E[Q_{ATE}])` (see row 3b). |
| 3b | **RMSE_ATE** | (derived from row 3a) | (none; derived statistic) | `θ_{3b} = sqrt( E[ Q_{ATE,j} ] )` | `sqrt( mean(Q_{ATE,j}) )` | **delta method:** `MCSE ≈ s_{Q_{ATE}} / ( 2·sqrt(J)·sqrt(mean(Q_{ATE})) )`; or **paired bootstrap** on `Q_{ATE,j}` (recommended) | on the squared scale as row 3a | paired-boot (preferred), delta-method CI only as a cross-check | Contrasting RMSE_ATE between methods is a contrast **on the squared scale**, then re-rooted; do not pair/re-root per-replication raw values. `g(x)=sqrt(x)`, `g'(x)=1/(2 sqrt(x))`, so the denominator carries `sqrt(mean(Q))`, **not** `mean(Q)`. Verified numerically (`P1_math_verification.ipynb` §5): dropping the inner root inflates the MCSE by 2× at `E[Q]=0.25` and understates it by 5× at `E[Q]=25`, i.e. exactly the PEHE-scale regime of this case study. |
| 4 | **ATE coverage (95%)** | indicator `1(τ_j ∈ CI_j)` for ATE `τ_j` | `Cov_{ATE,j} = 1(τ_j ∈ CI_j) ∈ {0,1}` | `θ_4 = P(cover)`, nominal target `0.95` | `mean(Cov_{ATE,j})` | `sqrt( p̂(1−p̂)/J )`, where `p̂ = mean(Cov_{ATE,j})` | binary: `D_j = Cov_{ATE,A,j} − Cov_{ATE,B,j}` | **McNemar** (paired binary) | nominal target makes this a calibration object, not "lower is better" (IJDA #10). Compare deviation from nominal (row 8) or the interval score (row 9). |
| 5 | **CATE coverage (95%)** | indicator `1(τ_{j,i} ∈ CI_{j,i})` per unit `i` | `Cov_{CATE,j} = mean_i 1(τ_{j,i} ∈ CI_{j,i})` | `θ_5 = E[ proportion covered ]`, nominal target `0.95` | `mean(Cov_{CATE,j})` | `s_{Cov_{CATE}}/sqrt(J)` | `D_j = Cov_{CATE,A,j} − Cov_{CATE,B,j}` | paired-t; paired-boot strongly preferred (column is bounded and strongly skewed, see §1.1 of the plan) | bounded `[0,1]`; skew shown empirically (pairing r up to 0.80). Report the deviation-from-nominal row separately for calibration (row 8). |
| 6 | **Mean credible-interval width (CIL), CATE and ATE** | per-unit width `( \hat{τ}^{up} − \hat{τ}^{lo} )`; at ATE level `( \hat{τ}_j^{up} − \hat{τ}_j^{lo} )` | `CIL_j = mean_i (\hat{τ}_{j,i}^{up} − \hat{τ}_{j,i}^{lo})` (CATE) or `CIL_j^{ATE} = ( \hat{τ}_j^{up} − \hat{τ}_j^{lo} )` | `θ_6 = E[CIL_j]` | `mean(CIL_j)` | `s_{CIL}/sqrt(J)` | `D_j = CIL_{A,j} − CIL_{B,j}` | describe (point + CI) by default; paired-t if a test is (rarely) wanted | corresponds to manuscript Eq. (6), (7). Eq. (7) needs an index: write `CIL_j^{ATE}`, average over `j` (see §4 text). Sharpness object for the C6 analysis. |
| 7 | **Bias (CATE)** | per-unit `(\hat{τ}_{j,i} − τ_{j,i})` | `Bias_j = mean_i (\hat{τ}_{j,i} − τ_{j,i})` | `θ_7 = E[Bias_j]` | `mean(Bias_j)` | `s_{Bias}/sqrt(J)` | `D_j = Bias_{A,j} − Bias_{B,j}` | paired-t, paired-boot | reported for diagnostics and the pre-registered analysis; not a headline decision metric. |
| 8 | **Calibration deviation (coverage)** | – | *(no per-replication `L_j`; see the note)* `L_j = Cov_j` is carried, and the deviation is formed **after** averaging | `θ_8 = \|E[Cov_j] − 0.95\|`, target `≈ 0` with a declared **equivalence band** `[0, Δ_cal]` | `\|mean(Cov_j) − 0.95\|` | `s_{Cov}/sqrt(J)` (the MCSE of `mean(Cov_j)`; the absolute value is a monotone transform of a scalar, so the band test is run on `mean(Cov_j)` itself) | `D_j = Cov_{A,j} − Cov_{B,j}`, tested for equivalence against `±Δ_cal` | **equivalence (TOST)** on `E[Cov_j]` against `[0.95−Δ_cal, 0.95+Δ_cal]`; for a two-method comparison, TOST on `E[D_j]` | This is how an over-/under-covering method is judged; it is **not** "lower coverage is better" (IJDA #10). The band `Δ_cal` is fixed before analysis (e.g. `±0.01`); see `tisca_v2_spec.md` §1 mode M5. |
| 9 | **Interval score (Gneiting-Raftery; Winkler at the two levels)** | per-unit interval score | `IS_j^{(0.95)} = mean_i [ (u_i−l_i) + (2/α)(l_i − x_i) 1(x_i<l_i) + (2/α)(x_i − u_i) 1(x_i>u_i) ]` for the 95% central interval `(l_i,u_i)` and observed `x_i`, `α=0.05` | `θ_9 = E[IS_j]` (lower is better) | `mean(IS_j)` | `s_{IS}/sqrt(J)` | `D_j = IS_{A,j} − IS_{B,j}` | paired-t, paired-boot; **primary MCS loss for uncertainty quantification** (lower better, scalar, one value per replication) | proper scoring rule jointly penalising width and non-coverage (C6); used as the MCS loss. The original authors already record CRPS; `CRPS_j = mean_i` CRPS over units is a valid MCS loss too and is preferred where it applies (§1.7b). |
| 10 | **CRPS** | per-unit CRPS | `CRPS_j = mean_i CRPS( \hat{F}_{j,i}, x_{j,i} )` | `θ_{10} = E[CRPS_j]`, lower better | `mean(CRPS_j)` | `s_{CRPS}/sqrt(J)` | `D_j = CRPS_{A,j} − CRPS_{B,j}` | paired-t, paired-boot; **MCS loss** | already recorded by the original authors (`scoringRules`); the natural proper-scoring-rule loss for the uncertainty comparison (§2.3). |
| 11 | **Runtime** | seconds to fit one model on one replication | `T_j` (seconds) | `θ_{11} = E[T_j]`; **report distribution** | `mean(T_j)`, plus `q05/q50/q95` and `sd` | `s_T/sqrt(J)` | `D_j = T_{A,j} − T_{B,j}` (if compared) | describe (distribution); do not hypothesis-test speed as a decision | satisfies IJDA #14a (timing distributions), and feeds the hardware-specific cost illustration in P5-T6. |
| 12 | **Convergence / failure flag** | – | `F_j = 1` if replication `j` of this method failed to converge or errored, `0` otherwise | `θ_{12} = P(fail)` | `mean(F_j)` | `sqrt( p̂_f(1−p̂_f)/J )` | – | report count and proportion; separate diagnostic, never pooled silently | `converged_flag`/`error_message` per model; governs whether a replication is listwise-eligible for any other row. |

### 3.1 Why row 8 is `|E[Cov] − 0.95|` and **not** `E[|Cov_j − 0.95|]`

An earlier draft of this table defined the calibration estimand as the mean per-replication absolute deviation `E[|Cov_j − 0.95|]`. That is the wrong object and it must not be used, because it **confounds miscalibration with replication-level dispersion**. `Cov_j` is itself a Monte Carlo quantity: even a perfectly calibrated method has `Cov_j` scattered around 0.95, and `E|Cov_j − 0.95| > 0` purely because of that scatter. Verified numerically (`P1_math_verification.ipynb` §6), with `J = 500`:

| method | true coverage | `E[\|Cov_j − 0.95\|]` | verdict under the wrong estimand |
|---|---|---|---|
| A: calibrated, `sd(Cov_j) = 0.010` | 0.950 | 0.008 | "best calibrated" |
| B: calibrated, `sd(Cov_j) = 0.040` | 0.950 | 0.032 | "worst calibrated" — but it is exactly as calibrated as A |
| C: **mis**calibrated by 0.02, `sd(Cov_j) = 0.005` | 0.930 | 0.020 | scores *better* than the perfectly calibrated B |

So the per-replication absolute deviation would rank a genuinely miscalibrated method above an unbiased one and would reward whichever method happens to have the tighter across-replication spread — precisely the kind of metric-direction error IJDA #10 raises. The deviation must be taken **after** the across-replication expectation.

Consequence for the software: `Cov_j` is the column that is carried and paired; the absolute value is applied to the estimate, never per replication. Sharpness (row 6) and the interval score (row 9) remain the separate, and jointly complete, half of the C6 story. The dispersion `sd(Cov_j)` is still worth reporting, but as a **stability diagnostic**, explicitly labelled as such, not as calibration.

### 3.2 Coverage and interval score are reported at **both** nominal levels

The original authors record 50% *and* 95% coverage and interval widths (revision plan §1.7b), and the E3 superset keeps `cov_50`, `cov_95`, `width_50`, `width_95`, `ate_cov_50`, `ate_cov_95`, `ate_width_50`, `ate_width_95`. Rows 4, 5, 6, 8 and 9 are therefore **parameterised by the nominal level `1−c ∈ {0.50, 0.95}`**, giving one instantiated row per level, with the nominal target in rows 4, 5 and 8 set to `1−c` rather than hard-coded to 0.95. Two levels is what turns "coverage was near nominal" into an actual calibration curve, and it is the stronger answer to IJDA #10; reporting only the 95% level discards half the evidence that was already collected.

Note the level parameter `c` in the interval-score formula (row 9) is the **interval's** nominal miscoverage, and is a different object from the **test** level `α` used everywhere in `tisca_v2_spec.md`. The table uses `c` for the former to keep them apart; row 9's `(2/α)` factors are to be read as `(2/c)`, with `c = 0.05` at the 95% level and `c = 0.50` at the 50% level.

> **Equation mapping (manuscript Section 2):** Eq. (1) → defines `τ_i`, unit CATE, and the DGP (see §4); Eq. (2) → row 1 (PEHE, rooted) and row 2 (unrooted); Eq. (3) → rows 3a/3b with the index fix `∑_{i=j}^{J} → ∑_{j=1}^{J}`; Eq. (4) → row 5 (CATE coverage) with `τ_i`; Eq. (5) → row 4 (ATE coverage) with `τ_j`, with the index fix; Eq. (6) → row 6 CIL_CATE; Eq. (7) → row 6 CIL_ATE, with the addition of an explicit index and an average over `j`.

---

## 4. The text to paste into `main.tex` (edits to equations (1), (3), (5), (7))

These are the concrete substitutions for P1-T1's "also fix here" requirement (they also get applied in P5-T2). Inline math is given as LaTeX; the maintainer drops these into the respective `\begin{equation}` blocks.

**Eq. (1) — distinguish potentials from observed, and state the shared-disturbance assumption.**
Current: `Y_i(Z_i) = f(\mathbf{X}_i) + \tau(\mathbf{X}_i) \cdot Z_i + \epsilon_i`.
Replace with the two-line version that separates `Y_i(0), Y_i(1)` from the observed `Y_i`:
```
Y_i(1) = f(\mathbf{X}_i) + \tau(\mathbf{X}_i) + \epsilon_i, \quad
Y_i(0) = f(\mathbf{X}_i) + \epsilon_i,
\qquad Z_i \perp Y_i(0), Y_i(1) \mid \mathbf{X}_i,
\qquad Y_i = Y_i(Z_i).
```
Text to accompany: because a single shared disturbance `ε_i` affects both potential outcomes, the unit-level treatment effect is exactly `τ(X_i) = Y_i(1) − Y_i(0)`; there is no additional individual-level effect to estimate. This is what makes `τ(X_i)` both the estimand and the "true CATE" the PEHE is taken against.

**Eq. (3) — RMSE_ATE index fix (IJDA minor #3).**
Current: `\text{RMSE}_{ATE} = \sqrt{ \frac{1}{J} \sum_{i=j}^{J} (\hat{\tau}_j − \tau_j)^2 }`.
Replace lower limit and the running index:
```
\text{RMSE}_{ATE} = \sqrt{ \frac{1}{J} \sum_{j=1}^{J} (\hat{\tau}_j − \tau_j)^2 }.
```
Keep the convention that `τ_j` is the replication-`j` ATE (not to be confused with the unit CATE `τ_i`). In the new representation this is `\sqrt{E[Q_{ATE}]}` where `Q_{ATE,j} = (\hat{τ}_j − τ_j)^2`; the per-replication squared error is row 3a, the rooted ACROSS-replication quantity is row 3b. Note the paper replaces `\sum_{i=j}^{J}` (which equals the single addend at `i=j`) with the correct average.

**Eq. (5) — Coverage_ATE index fix (IJDA minor #3).**
Current: `\text{Coverage}_{ATE} = \frac{1}{J} \sum_{i=j}^{J} \mathbf{1}( \tau_j \in [\hat{\tau}_j^{lo}, \hat{\tau}_j^{up}] )`.
Replace:
```
\text{Coverage}_{ATE} = \frac{1}{J} \sum_{j=1}^{J} \mathbf{1}\Big( \tau_j \in \big[ \hat{\tau}_j^{lo},\; \hat{\tau}_j^{up} \big] \Big).
```

**Eq. (7) — CIL_ATE needs an index and an average over `j` (IJDA minor #4).**
Current: `\text{CIL}_{ATE} = ( \hat{\tau}_j^{up} − \hat{\tau}_j^{lo} )` (no `j` averaging).
Replace with a per-replication width, then the reporting estimand:
```
\text{CIL}^{(j)}_{ATE} = \hat{\tau}_j^{up} − \hat{\tau}_j^{lo}, \qquad
\text{CIL}_{ATE} = \frac{1}{J} \sum_{j=1}^{J} \text{CIL}^{(j)}_{ATE}.
```

**`τ_i` versus `τ_j` (SoutoNeto §1).** Throughout Section 2, reserve `τ_i` for the per-unit CATE (rows 1, 2, 5, 6, 7) and `τ_j` (equivalently `\overline{τ}^{(j)}`) for the replication-level ATE (rows 3a, 3b, 4, 6). Never mix them in one formula.

---

## 5. Implementation checklist (what each row must become in code)

For every performance metric in columns of the E3 CSVs (which carry the original `GitHub_DGP1.R` columns plus the planned additions `cate_mse`, `ate_sq_err`, `fit_seconds`, `converged_flag`, per the P3-T5 superset list):

1. Identify the replicate-level column(s) that produce `L_j`. If the column is already replication-level (e.g. `mvbcf_pehe1`), `L_j` is that column directly. If the CSV only stores a per-unit downstream quantity, first aggregate to `L_j` **before** any across-replication statistic.
2. Confirm `L_j` is usable for a paired contrast: the same `j` must have values for both methods. If `L_j` is missing for either method, drop replication `j` from **that** contrast (listwise) and record the drop count; never drop per-column.
3. Compute `L̄`, `s_L`, `MCSE`, and the CI per §2.
4. For binary rows (4, 12) verify `L_j ∈ {0,1}`; for bounded rows (5) verify `0 ≤ L_j ≤ 1`; otherwise raise a validation error in `validate.py` (P2-T1) before reporting.
5. Refuse (IJDA #11e): a request to contrast a column that is an *across-replication aggregate* rather than a replicate-level value; a per-replication "RMSE_ATE" scalar where only the squared error is a valid per-replication object; or any metric with a nominal target (rows 4, 5, 8) that the caller tries to treat as "lower is better" (route those to the interval score row instead).
6. For row 8, apply the absolute value to the **estimate**, never to `Cov_j` (§3.1). `validate.py` must reject a caller-supplied column of per-replication absolute deviations.
7. Instantiate rows 4, 5, 6, 8 and 9 once per nominal level in `{0.50, 0.95}` (§3.2), and carry the level in the output row's metric name so the two levels are never silently pooled.

**End of P1-T1.**
