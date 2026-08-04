# P1-T2 — TISCA v2 formal specification

**Status:** draft for review (L-SPEC, task P1-T2)
**Revision:** 2026-08-04 (rev. 2). Corrections in rev. 2, all backed by `notebooks/P1_math_verification.ipynb`: M5 feasibility rule was wrong (§1); the §3 variance-inflation rationale was stated as a bias argument when it is an assurance argument; §4.4's Romano–Wolf planning covariance was off by a factor of `J`; §4.3's BH planning level was anti-conservative; §4.5 (new) works out what the v2 defaults do to the case study's `J`.
**Dependencies:** `docs/estimand_table.md` (P1-T1) — this document inherits its notation verbatim.
**Acceptance (restated from the plan):** a competent statistician can implement v2 from this document alone with no reference to v1 code. Every one of C1–C8 is addressed by a named subsection.
**Scope:** the design/planning algorithm only. The paired-inference and Model-Confidence-Set machinery invoked here is specified at the level needed for planning; the full resampling spec lives with P2-T1/P2-T2 (`inference.py`, `mcs.py`) and Section 3.6 of the manuscript.

The strategic intent is unchanged from the plan's §2.1: TISCA v2 is a **two-layer simulation-design protocol**. A precision layer (default) is always meaningful; an optional decision (power) layer is added only when a confirmatory comparative claim is being made. `J_final = max` over comparisons × active layers, bounded by a declared `J_max`.

---

## 0. What is preserved from v1 and what is removed

Preserved: the aim of selecting a replication count `J` on a principled basis; the per-pair `(A, B)` structure; the plugin approach (use a small pilot to estimate unknown variance, then solve for `J`).

Removed and replaced:

| v1 element | v2 replacement | Rationale / task |
|---|---|---|
| Welch's unpaired t as the test of choice | **paired** one-sample inference on `D_j = L_{A,j} − L_{B,j}` (estimand table §1, convention 3) | IJDA #2, R2 ¶1, SoutoNeto §3 |
| Iterative accumulate-and-recheck loop as the *default* | **two-stage default**: independent pilot then one closed-form sized confirmatory run (Algorithm 1 below) | IJDA #5, R2 ¶4 |
| Analytical power with `alpha/2` critical value, no multiplicity | power at the **adjusted** α and same sidedness as the final test, with explicit multiplicity handling (§4.3, §4) | IJDA #4, #6 |
| MDE as the sole user input | "planning alternative" δ plus an optional separate "minimally important difference"; both with explicit mode + sidedness | IJDA #4 |
| `correction_method = "none"` default | **multiplicity-aware planning is on by default** whenever more than one primary contrast is declared; `"none"` a conscious, declared choice | R2 ¶2 |

---

## 1. Hypothesis modes M1–M5

For a scalar loss where **lower is better** (PEHE, interval score, CRPS, ATE squared error), define the paired contrast `D_j = L_{A,j} − L_{B,j}`, `j = 1, …, J`, i.i.d. with `θ = E[D_j]` and `σ_D² = Var(D_j)`. A negative `θ` means method A beats B. For a bounded target metric (coverage) the comparison is re-expressed via the calibration deviation `δ_j == |Cov_j − nominal|` (estimand table row 8) or the interval score (row 9) so that **lower is better** holds everywhere a test is applied. This one rule ("express every tested difference so that lower is better, then apply the modes below") removes the coverage-is-not-monotone problem (IJDA #10) and lets every mode below be stated once.

Test statistic for the mean of `D_j`:
```
D̄ = (1/J) Σ_j D_j,   s_D² = Σ_j (D_j − D̄)²/(J−1),   T = sqrt(J)·(D̄ − θ0)/s_D ~ t_{J−1} under H0.
```
Planning uses the **noncentral t with `J−1` df** (the R `pt(q, df, ncp)`; in Python `scipy.stats.nct`). **Do not** use `scipy.stats.t.cdf(crit, df, ncp)` with `ncp` in the third (loc) slot — that is the v1 bug (paper §1.6); it happens to be small here but must not be propagated.

`Δ` denotes a **margin** (a fixed, externally justified boundary, not estimated); `δ` denotes the **planning alternative** (the value of `θ` at which we want a given power); `α` is the (possibly adjusted) level; `J` is the running count. Power functions below give `Pr(reject | θ = δ)`.

### M1 — Two-sided equality (`H0: θ = 0` vs `H1: θ ≠ 0`)
- Acceptance region: `|T| ≤ t_{1−α/2, J−1}`.
- Power:
  ```
  β1(J) = 1 − [  F_{t,J−1,ncp}( t_{1−α/2,J−1} ) − F_{t,J−1,ncp}( −t_{1−α/2,J−1} ) ],
  ncp = sqrt(J)·δ/s_D,   (F_{t,·,·} = noncentral t CDF)
  ```
- This reproduces the v1 two-sided target and the plan's §1.2 J values (verified in `notebooks/P1_math_verification.ipynb`: at δ/σ_D = 0.5, 80% power needs J = 34; the plan's six §1.2 cells reproduce as 114/97/205/185/49/66 — coverage-vs-wsBCF, PEHE-vs-BCF and coverage-vs-MVBART for Y1/Y2 respectively — when `σ_D` is recomputed from `DGP1_500_results.csv` rather than from the 4-decimal values quoted in the plan).

### M2 — Directional superiority, lower is better (`H0: θ ≥ 0` vs `H1: θ < 0`)
- Reject iff `T < t_{α, J−1}`.
- Power: `β2(J) = F_{t,J−1,ncp}( t_{α, J−1} )`, `ncp = sqrt(J)·δ/s_D` with `δ < 0`.

### M3 — Minimum-effect superiority (`H0: θ ≥ −Δ` vs `H1: θ < −Δ`)
- Reject iff `T < t_{α, J−1}` where `T` is centered on the inner boundary `−Δ`.
- Power at `δ < −Δ`: `β3(J) = F_{t,J−1,ncp}( t_{α, J−1} )`, `ncp = sqrt(J)·(δ + Δ)/s_D`.
- Used when a "practically meaningful" beat (at least Δ) is the claim; if `Δ` is not known, prefer M1/M2 rather than guessing a margin.

### M4 — Non-inferiority (here: `H0: θ ≥ Δ` vs `H1: θ < Δ`)
- Reject iff `T' < t_{α, J−1}` where `T' = sqrt(J)·(D̄ − Δ)/s_D` (test the difference against the upper margin).
- Power at `δ < Δ`: `β4(J) = F_{t,J−1,ncp'}( t_{α, J−1} )`, `ncp' = sqrt(J)·(δ − Δ)/s_D` (`F_{t,·,ncp}` = noncentral t CDF; the rejection region is the lower tail).
- Interpretation in the loss framing: "A is not worse than B by more than Δ" is established by rejecting that `θ ≥ Δ`. (Verified by Monte Carlo against the analytic form.)

### M5 — Equivalence, two one-sided tests (TOST) (`H0: |θ| ≥ Δ` vs `H1: |θ| < Δ`)
- Reject iff `T_lo > t_{1−α, J−1}` AND `T_up < −t_{1−α, J−1}`, where `T_lo = sqrt(J)·(D̄ + Δ)/s_D`, `T_up = sqrt(J)·(D̄ − Δ)/s_D`.
- Power: both TOST statistics must fall in their combined rejection triangle. With `T_lo = T_up + 2·sqrt(J)·Δ/s_D`, the power at `δ` with `|δ| < Δ` is
  ```
  β5(J) = F_{t,J−1}( b_hi ) − F_{t,J−1}( b_lo ),
  b_lo = ( −Δ + t_{1−α,J−1}·s_D/sqrt(J) − δ ) · sqrt(J)/s_D,
  b_hi = (  Δ − t_{1−α,J−1}·s_D/sqrt(J) − δ ) · sqrt(J)/s_D ,
  F_{t,J−1} = central t CDF (the signed boundary of D̄).
  ```
  This is the **known-σ approximation**: it fixes the acceptance boundary at `m = Δ − t_{1−α,J−1}·σ_D/sqrt(J)`, whereas the test actually run compares against `Δ − t_{1−α,J−1}·s_D/sqrt(J)` with the *sample* `s_D`. Verified against the exact TOST decision rule (`P1_math_verification.ipynb` §3, both by 4·10⁵-draw Monte Carlo of the real studentized rule and by exact quadrature over `s_D ~ σ_D·sqrt(χ²_{J−1}/(J−1))`): the approximation is accurate to **≤ 0.009 in power** over `J ∈ [30, 200]`, `|δ| ≤ 0.4Δ`, and returns the **same planned `J`** at the 0.80 target. It is therefore fit for planning; it is **not** the object to use if an exact TOST power is ever needed for reporting, where the quadrature form
  ```
  β5_exact(J) = E_{s_D} [ Φ( (m(s_D) − δ)·sqrt(J)/σ_D ) − Φ( (−m(s_D) − δ)·sqrt(J)/σ_D ) ] ,
  m(s) = max( Δ − t_{1−α,J−1}·s/sqrt(J), 0 ),  s_D ~ σ_D·sqrt(χ²_{J−1}/(J−1))
  ```
  should be used instead. `planning.py` implements both and defaults to the approximation with the exact form available as `exact=True`.
- **Feasibility (corrected).** If `t_{1−α, J−1}·s_D / sqrt(J) ≥ Δ` then `m ≤ 0`, the two acceptance boundaries cross, and TOST power is 0 **at that `J`**. This is a small-`J` artifact, **not** a statement about `Δ`: `m` increases monotonically toward `Δ` as `J` grows, so the design becomes feasible at larger `J`. Verified: at `Δ = 0.5`, `σ_D = 1`, `α = 0.05`, `J = 10` gives `m = −0.08` (power 0) while `J = 40` gives power 0.86. **The planner must therefore keep increasing `J` past a zero-power point, and must never report global infeasibility from it.**

  The genuine infeasibility condition is `|δ| ≥ Δ`: the planning alternative lies on or outside the equivalence margin, and no `J` reaches the target power. At `|δ| = Δ` the power converges to `α` from below (verified: 0.050 at `J = 50, 200, 1000, 5000`); at `|δ| > Δ` it converges to 0 (0.010 → 0.000). **Rejection rule:** return "infeasible: |δ| ≥ Δ, the planning alternative is not inside the equivalence margin" when `|δ| ≥ Δ`, and otherwise solve normally up to `J_max`.

**Mode choice rule:** the mode must be chosen once, before any confirmatory data, and must be the same one used for the final test and for planning (C4). The default is M1 two-sided for confirmatory comparative claims; M2 is the smallest practical choice when a strict one-sided claim is intended; M3/M5 require a defensible margin.

---

## 2. Design targets

Two targets; `J_final` is the maximum over comparisons and over active targets, bounded by the budget `J_max`.

### 2.1 Precision target (default, always available) — C2
Select the smallest `J` such that the Monte Carlo standard error meets a target:
- **MCSE target:** `s_D / sqrt(J) ≤ m` for a declared absolute target `m` on `MCSE(θ̂)`.
- **CI half-width target:** `t_{1−α/2, J−1} · s_D / sqrt(J) ≤ h` for a declared half-width `h`. Solve by scanning `J` until the inequality holds (it is monotone in `J` for fixed `s_D`).
Precision is meaningful whether or not any null is rejected; this is the layer Morris–White–Crowther / Burton et al. / Koehler et al. govern, and it is the recommended default (IJDA #1). For each row of the estimand table, `m` should be set as a fraction of the estimand's plausible scale or, for bounded metrics, of the target.

### 2.2 Power target (optional decision layer) — C4
Select the smallest `J` with `Pr(reject | θ = δ) ≥ 1 − β`, where the probability is computed with the noncentral t of **J−1** df and the **adjusted** level `α_adj` (not the unadjusted `α`). Each mode M1–M5 from §1 supplies its `β(J)`. `α_adj` is set by the multiplicity rule (§4.3) when a family of contrasts is planned. The planning alternative `δ` is the *value of θ we want to be able to detect*, distinct from any minimally important difference; it must be given a domain justification (see P3-T5(b) and §7 of `power_target_guidance.md`).

### 2.3 Combination and budget cap
`J_final = min( max_{active targets, comparisons} J_detail , J_max )`. If the closed-form `J_final` would exceed `J_max`, stop at `J_max` and report the achieved (reduced) precision/power rather than silently continuing (removes the v1 infinite-loop-on-sd=0 behaviour; see P3-T7).

---

## 3. Variance-uncertainty propagation (IJDA #8b) — C4/C8

A pilot of size `J0` yields `s_D`, an estimate of `σ_D`. Because `s_D` is noisy, roughly half of all studies planned from it end up **under-powered**, even though the *average* `J` is about right. Inflate the pilot estimate to a one-sided upper confidence bound before solving:
```
σ̂_D_UB = s_D · sqrt( (J0 − 1) / χ²_{γ, J0−1} ),   γ = 0.20 default.
```
`χ²_{γ, J0−1}` is the γ-quantile of a chi-square with `J0−1` df; `γ = 0.20` means "with 80% confidence the true `σ_D` is at most `σ̂_D_UB`". Verified inflation factors: J0 = 25 → 1.153, J0 = 50 → 1.099, J0 = 100 → 1.067.

**State the rationale correctly — it is an *assurance* argument, not a bias argument.** Planning from the raw `s_D` is *not* systematically under-sized in `E[J]`; it is under-sized in the **probability of hitting the target**. Measured (`P1_math_verification.ipynb` §4; true `σ_D = 1`, `δ = 0.5`, target power 0.80, oracle `J = 34`):

| pilot `J0` | plan from | `E[J]` | `E[achieved power]` | `Pr(achieved power ≥ 0.80)` | 5th pct of achieved power |
|---:|---|---:|---:|---:|---:|
| 25 | raw `s_D` | 33.8 | 0.781 | **0.48** | 0.587 |
| 25 | `σ̂_D_UB` | 44.1 | 0.873 | **0.81** | 0.688 |
| 50 | raw `s_D` | 33.8 | 0.793 | **0.49** | 0.650 |
| 50 | `σ̂_D_UB` | 40.4 | 0.859 | **0.81** | 0.739 |
| 100 | raw `s_D` | 33.8 | 0.799 | **0.51** | 0.706 |
| 100 | `σ̂_D_UB` | 38.2 | 0.846 | **0.82** | 0.754 |

`E[J]` is within 0.6% of the oracle in every raw-`s_D` row, so any text claiming that the point estimate "under-sizes on average" is wrong and must not be written into the manuscript. What the inflation buys is `Pr(achieved power ≥ target)` rising from ≈ 0.5 to ≈ 1 − γ = 0.8, at a cost of 12–30% more replications. That is the tradeoff to state, and `γ` is the dial: it is the assurance level, and it should be reported alongside `J`.

Use `σ̂_D_UB` everywhere `σ_D`/`s_D` enters the planning solve (§2.1, §2.2, §4.4), **not** the bare pilot `s_D`. The pilot rows are then **discarded** from final inference, so they cannot introduce optimism into the final test.

The same inflation is applied to the pilot covariance `Σ̂_D` (each diagonal element), which is what the multiplicity-aware planning (§4.4) and any bootstrap-family planning use.

---

## 4. Multiplicity (C5, IJDA #6, R2 ¶2)

### 4.1 Family definition
The user declares a **family** = the set of primary contrasts that in aggregate answer the study's confirmatory question. Everything OUTSIDE the declared family is exploratory (no multiplicity control; reported with a flag). Family size `K` = number of primary contrasts.

### 4.2 Success criterion
The user must declare which of the following counts as "success", and planning/power reporting must match it:
- **Conjunctive** ("all K claims must hold"): plan so that each marginal claim holds with enough margin that the joint requirement is met; report conjunctive power (explicitly below-1 each cell).
- **Disjunctive** ("at least one holds"): the multiplicity burden is the lightest; report disjunctive power.
- **Named subset**: exactly the reader's list. Default when unspecified is conjunctive on all K primary claims.

### 4.3 Correction mapping
| Correction | How to adjust the planning level |
|---|---|
| Bonferroni | plan every contrast at `α/K`; simplest, always valid |
| Holm | asymptotically no worse than Bonferroni; plan at `α/K` is conservative and valid; report the achieved (actual) level |
| Benjamini–Hochberg (FDR) | **do not plan at the marginal `α`.** BH's per-claim threshold ranges from `α/K` (the smallest p-value) to `α` (the largest), and which one binds depends on how many of the `K` nulls are false — an unknown at design time. Plan at `α·r/K` where `r` is the *pre-declared, conservative* number of claims expected to be rejected (`r = 1` reduces to Bonferroni and is the safe default for a conjunctive criterion), or simulate the BH decision from `Σ̂_D` as for Romano–Wolf. Report family-level FDR, not per-claim T1E |
| Romano–Wolf stepdown | **plan by simulating power from the pilot covariance `Σ̂_D`** (§4.4): the joint distribution of `(T_1,…,T_K)` under the planning alternative is approximated as multivariate normal with covariance `Σ̂_D`, and power = fraction of draws in which the stepdown set includes the true positives. This is the strictly better answer to IJDA #6 because it exploits the true correlation instead of assuming the worst case. |

`"none"` is **not** the default. It is allowed only when K = 1 or when the user explicitly declares "I am not adjusting for multiplicity and accept the implications", which the software must record in the output.

### 4.4 Planning with `Σ̂_D`
For Bonferroni/Holm this reduces to the per-contrast solves at `α_adj = α/K`.

For Romano–Wolf the planning is a nested Monte Carlo over candidate `J`. For each candidate `J`:

1. Draw `R_plan` datasets of `J` rows each from `MVN( (δ_1,…,δ_K), Σ̂_D_UB )` — that is, simulate the **replication-level contrast vectors** `D_j`, not their means.
2. From each dataset form the `K` studentized statistics `T_k = sqrt(J)·D̄_k / s_{D,k}` exactly as the final analysis will, run the Romano–Wolf stepdown at level `α` with `B_plan` bootstrap resamples, and record which claims were rejected.
3. Power = the fraction of the `R_plan` datasets meeting the success criterion of §4.2.

Take `J` as the smallest candidate reaching the target.

**The covariance argument is `Σ̂_D` (the covariance of one replication's contrast vector), not `Σ̂_D/J`, precisely because step 1 draws rows rather than means.** Drawing `(D̄_1,…,D̄_K)` directly with covariance `Σ̂_D` — as an earlier draft of this section specified — inflates the sampling variability of the means by a factor of `J` and would make the planner return wildly oversized `J`. If a fast approximation is wanted, draw the means from `MVN(δ, Σ̂_D_UB / J)` and treat the variances as known, but then the stepdown must be run on `z`-statistics, and the approximation must be checked against the row-level version at the chosen `J`.

This is a nested-MC step; its own quantile noise must be reported (report an MCSE on the planned `J`, or plan on a grid and report the grid resolution).

**Report both** marginal (per-contrast) power and the family-level (conjunctive/disjunctive) power, as the plan requires.

### 4.5 Worked example — what the v2 defaults do to the case study's `J` (read this before writing the abstract)

The revision plan's §1.2 computes `J = 205` for the MVBCF family and calls the resulting ~5× saving "the single most important strategic fact in the revision". **That number is computed at an unadjusted `α = 0.05` and from the bare pilot `s_D`, i.e. with two of this specification's own defaults switched off.** Applying the spec as written to the same six contrasts and the same `DGP1_500_results.csv` (verified numerically; the family is `K = 6`, Bonferroni, `J0 = 50`, `γ = 0.20`):

| contrast | δ | `J` at α = 0.05, raw `s_D` | `J` at α/6 | `J` at α/6 **and** `σ̂_D_UB` |
|---|---|---:|---:|---:|
| coverage Y1, MVBCF vs wsBCF | 0.015 | 114 | 177 | 213 |
| coverage Y2, MVBCF vs wsBCF | 0.015 | 97 | 149 | 180 |
| PEHE Y1, MVBCF vs BCF | 0.5 | 205 | 316 | **382** |
| PEHE Y2, MVBCF vs BCF | 0.5 | 185 | 285 | 344 |
| coverage Y1, MVBCF vs MVBART | 0.015 | 49 | 75 | 90 |
| coverage Y2, MVBCF vs MVBART | 0.015 | 66 | 101 | 122 |
| **`J_final` = max over the family** | | **205** | **316** | **382** |

So the spec-compliant answer for this case study is `J ≈ 382`, a **2.6× saving against the original 1000**, not 5×. Three consequences that must be carried into the manuscript and the response letter:

1. **Do not publish the 5× figure.** A referee who applies §4 to the paper's own case study lands on 382, and the discrepancy would read as the paper not following its own procedure. The honest headline is "≈ 2.6× fewer replications, with multiplicity and pilot-variance uncertainty both accounted for" — still a substantial saving, and far more defensible than a number obtained by turning off the two corrections the revision exists to add.
2. **Pairing is still worth ~2.3×, and that is the claim to make.** Under the *same* multiplicity and inflation settings, the unpaired v1 machinery needs `J_final = 882` against the paired 382. The like-for-like gain from fixing the reviewer's #2 objection is 2.31×; it should be reported as a comparison at matched settings, never as 1000/205.
3. **`J_final` here is driven entirely by the PEHE-vs-BCF contrasts, and those are the contrasts §1.3 of the plan is about to change.** The BCF benchmark is being recalibrated (P3-T5(e)); when the cold-start `bcf` column is replaced by the calibrated warm-start one, `σ_D` for the PEHE contrasts changes and `J_final` moves with it. Every number in this table is therefore **provisional until the P3-T5(e) gate passes**, and no version of it belongs in the abstract before then.

---

## 5. The two-stage default procedure — Algorithm 1 (C7, C8)

This is the new default. It directly answers R2 ¶4 ("the iterative algorithm seems unnecessarily complicated"). Seeds follow `docs/seed_rng_protocol.md` (P0-T2): pilot uses a disjoint seed block `[1_000_001, 1_000_000 + J0]`, confirmatory uses `[1, J]`.

```
Algorithm 1 — two-stage design (default)
Inputs: sim_func(seed) -> per-replication metric vectors;
        primary_contrasts (family, each with mode M, δ, target type);
        J0 (pilot size), γ=0.20, J_max, α, multiplicity rule, success criterion.
1. Run pilot: for j in 1..J0 (pilot seed block): L_j <- sim_func(seed_pilot[j]).
2. For each contrast (A,B): build D_j = L_{A,j} − L_{B,j} (listwise over the pair; record how many j dropped).
   Estimate D̄, s_D; form σ̂_D_UB per §3; (for RW) Σ̂_D_UB.
3. Solve for J per active target (§2.1, §2.2) at α_adj (§4.3) and σ̂_D_UB (§3).
   J_final <- clamp(max over contrasts and targets, J_max).
4. Discard pilot rows. Run confirmatory: for j in 1..J_final (confirmatory seed block): L_j <- sim_func(seed[j]).
5. Final inference: for each contrast, run the pre-specified test ONCE at α_adj (§6);
   report estimate + MCSE/CI first, p-value second (IJDA #1e).
6. (If more than two models) run the bootstrap family inference (§6.2).
7. Report operating characteristics of the design (achieved precision/power, E[J], sd(J), P(J=J_max)), per C8, and attach MCSEs. Report any J_final = J_max cap in a flag.
```

Properties to be validated empirically in P3-T2: unconditional Type I error, power, CI coverage of θ, bias of θ̂, and distribution of `J_final`.

Because the pilot is discarded, `J_final` is a function of the pilot data only, and the confirmatory sample is independent of it. Conditional on `J_final` the one-sample t-test is exact, so the **unconditional** Type I error is exactly `α_adj` by iterated expectation — the two-stage design has none of the internal-pilot inflation of an adaptive design. P3-T2 measures this rather than assuming it, and the expected result is agreement with `α_adj` to within Monte Carlo error (±0.0031 at `R = 5000` near 0.05). Note the corresponding statement for **power** is different and weaker: unconditional power is *not* the nominal `1 − β`, because `J_final` is random. It averages below target when planning from the raw `s_D` and lands near the `1 − γ` assurance level when planning from `σ̂_D_UB` (§3).

---

## 6. Adaptive mode (optional) — Algorithm 2 (C7, IJDA #5, R2 ¶4)

Retained as an **optional** algorithm, **off by default**. The default recommendation is **do not use it** unless the user has a documented operational need and accepts that its unconditional error rates must be validated (they are, by P3-T2, before it can be recommended).

```
Algorithm 2 — adaptive / internal-pilot re-estimation (optional, validated only)
Inputs: as Algorithm 1, plus: nmax_looks (max re-estimation points), J_max (hard cap).
1. J <- J0; run J pilot/accumulating replications (seed block per protocol).
2. At each of up to nmax_looks planned look times (or when J crosses them):
     re-estimate s_D and σ̂_D_UB (§3); recompute the required J_uncond (closed-form);
     if a futility/stop rule is used, it must be pre-registered here.
3. J <- min(cap, max(J, required)); if more looks remain and J < cap, goto 2.
4. Final test ONCE at α_adj on the FULL accumulated set INCLUDING the pilot rows. Reusing the pilot rows is
   exactly what distinguishes Algorithm 2 from Algorithm 1, and it is also the source of its error-rate risk.
   Flag in the output that the pilot rows were reused and were NOT discarded.
5. Report the design type and that its unconditional error rates come from P3-T2 (or P3-T3's verdict), not assumed.
```

Caveats that must be printed near Algorithm 2 (verbatim intent of plan §2.3): an internal-pilot re-estimation that **reuses** the pilot rows in the final test has well-known Type I error behaviour that is *not* equal to nominal (Kieser & Friede 2000; Proschan 2005), and stopping on a data-dependent criterion such as "when the MCS becomes a singleton" is optimistically biased. Do **not** stop on MCS cardinality (plan §2.3). The two-stage Algorithm 1 avoids all of these by discarding the pilot; that is why it is the default.

---

## 7. Inference layer (C1, C6, §2.3 of the plan)

### 7.1 Paired inference on a single contrast
For a scalar lower-is-better contrast `D_j` (rows 1, 2, 3a, 3b, 9, 10 of the estimand table): paired one-sample t on `D_j` (df `J_final − 1`) as the fast default; the **studentized paired bootstrap** (B = 4999, report its MCSE) as the recommended check when the skewness of `D_j` is non-trivial, with the Berry–Esseen-based operational rule (P3-T4) deciding when they may differ materially. For bounded/binary rows (4, 12) use McNemar's exact/binomial test on the paired indicators.

### 7.2 Family inference on a set of contrasts
Three bootstrap routines (specified fully in P2-T1 `mcs.py` / R port; unit-test against CRAN `MCS`):
- **Romano–Wolf stepdown** (bootstrap, FWER-controlling): replaces the `p.adjust` menu for the declared family; uses the actual correlation structure.
- **Hansen SPA / White's Reality Check**: tests "the proposed model beats all benchmarks" for a nominated champion.
- **Model Confidence Set (Hansen, Lunde & Nason 2011)**: the set of models indistinguishable from the best at level `α`, using **one scalar loss per replication** (lower better) such as the interval score (row 9) or CRPS (row 10), or PEHE (row 1). Report both the elimination path and the final set, and report the bootstrap-p-value MCSE (B = 4999).

Every p-value/CI/estimate carries a Monte Carlo standard error as demanded by C8.

---

## 8. Input validation and rejection rules (IJDA #11e) — `validate.py`

The implementation must refuse metric structures it cannot validly test. Concretely, at the point a caller requests a contrast, reject with an explicit message:

1. **Across-replication aggregate passed as a per-replication value.** If a caller hands a column that is already an across-replication mean of something (so there is no per-replication `L_j`), it cannot form `D_j`; require the per-replication column.
2. **Per-replication "RMSE_ATE".** A single squared ATE error is a valid per-replication object (row 3a); "the RMSE of one replication" is not. Refuse to accept a column named like a re-rooted per-replication RMSE as a mean to be tested; route through row 3b.
3. **Coverage-as-loss.** A metric with a nominal target (rows 4, 5, 8) must not be tested with the "lower is better" machinery directly. Either use the calibration deviation (row 8) or the interval score (row 9); refuse a bare-coverage contrast that the caller tags `lower-is-better`.
4. **Mismatched mode/level.** If the caller's planned mode differs from the test mode for the same contrast (e.g. planning two-sided, final test one-sided), refuse — C4 requires identity.
5. **s_D = 0 degenerate.** If `s_D < epsilon` (machine scale), M1/M2/M3 targets are vacuous and the v1 `sd > 0` guard forced power to 0 → infinite loop. v2 must terminate: treat power targets as satisfied (no more replications needed) and report the degenerate cell (P3-T7 required a graceful path and, when DGP is the same model, `θ = 0, D_j ≡ 0`).
6. **TOST infeasibility.** Refuse an M5 design with `|δ| ≥ Δ`: the planning alternative is on or outside the equivalence margin and no `J` reaches the target (§1 M5). **Do not** infer infeasibility from a zero-power evaluation at a small `J` — that is a `J`-artifact and the solver must keep increasing `J`. Beyond `J_max`, return "target not reached within budget", which is a different message from "infeasible".

---

## 9. Acceptance self-check against C1–C8

| Change | Addressed by |
|---|---|
| C1 paired contrasts, no per-column NA | §0, §1, §7.1; estimand table §1.3 / rows |
| C2 precision as default | §2.1 |
| C3 estimand table | P1-T1 (this task's sibling doc) |
| C4 aligned hypotheses; planning alternative vs MID | §1 (modes), §2.2, §4; P1-T3 |
| C5 multiplicity-aware planning | §4 |
| C6 coverage as calibration + sharpness | estimand table rows 8–10; §1 (lower-is-better rule); §7.2 |
| C7 two-stage default, adaptive optional | §5, §6 |
| C8 validated operating characteristics | §5 step 7; §7.1; P3-T2/P3-T3; P1-T3 |

Implementation targets (P2-T1/P2-T2): modules `planning.py` (modes M1–M5 × targets), `multiplicity.py` (Bonferroni, Holm, BH, Romano–Wolf), `mcs.py` (SPA, Reality Check, MCS with T_R and T_max statistics), `procedure.py` (Algorithm 1 + 2), `validate.py` (§8), `estimands.py` (estimand table rows).

**End of P1-T2.**
