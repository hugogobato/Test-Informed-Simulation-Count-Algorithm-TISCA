# P1-T3 — When is a power target appropriate? Decision guidance

**Status:** draft for review (L-SPEC, task P1-T3)
**Revision:** 2026-08-04
**Depends on:** `docs/tisca_v2_spec.md` (P1-T2), whose hypothesis modes (§1), precision target (§2.1) and power target (§2.2) this document qualifies.
**Audience:** the scientist running TISCA. The goal is a two-page, decision-oriented guide: when precision planning alone suffices, when to add a decision (power) layer, why large J is a *feature* of testing not a bug, and why a conditional-on-one-DGP result does not license "generally superior" claims.
**Acceptance (restated):** contains an explicit "do not use TISCA for X" list.

---

## 1. The two layers in one sentence

TISCA v2 plans the replication count from a **precision target** (how precisely do I want each estimand known?) by default, and adds a **power target** (how often must my pre-specified confirmatory test reject when an effect of size δ holds?) only when the study is making a confirmatory comparative claim. Precision is always meaningful; power is only meaningful when there is a well-defined hypothesis to be answered at a controlled error rate.

---

## 2. Use precision planning alone when ...

- **The deliverable is an estimate or an uncertainty band**, not a pass/fail. Example: "report the PEHE of each method with a Monte Carlo CI of half-width h". This is the commonest case and the default (C2, IJDA #1).
- **The comparison is exploratory or descriptive.** "We want to show the four methods' loss distributions" is a precision question about each estimand, and only about the estimands you pre-specify. Pairing still applies; power does not.
- **There is no null to reject.** Coverage/calibration (how close to nominal, with a band) and interval score are estimands of this type. Testing "is coverage exactly 0.95?" is usually the wrong question; asking "is coverage within ±0.01 of nominal, tightly estimated?" is a precision/calibration question (IJDA #10, C6).
- **The user cannot justify a single δ or margin.** If no defensible planning alternative or margin exists, a power target is a gesture, not a design decision; precision planning is honest where a made-up δ is not (see §4).
- **You only need to detect "some" difference with no operational consequence tied to its size.** Then the whole decision layer adds nothing.

The guidance rule: **if the analysis would be identical (modulo resolution) under the best- and worst-case value of the test statistic, you only needed precision.** Power buys the ability to *stand behind* a binary comparison, and only then.

---

## 3. Add a decision (power) layer when ...

- **A confirmatory comparative claim is the headline.** E.g. "MVBCF has lower PEHE than BCF on this DGP". Now a Type I / Type II error budget is meaningful, and the power layer sizes J so the study can actually detect the planned difference if present (R2 ¶2, IJDA #4).
- **You must quote an error rate.** The reviewer or the protocol asks "what was the power to detect the effect you claim?" Only a power-planned design can answer this; precision-only planning cannot supply a power number.
- **A decision rule with cost lives downstream.** If choosing method A over B triggers a purchase, a label, or a regulatory/clinical step at a specific effect size, you want the error rates controlled at the size you can stand behind.
- **A minimum practically important effect is known.** You have a defensible `Δ` (or a δ that is not bare units — see §4) that is the *smallest effect that would change a decision*. That is a signal to add the decision layer; the planning alternative must then be set and justified (P3-T5(b)).

The guidance rule: **add power only when the report will include a designed-for binary decision.** Otherwise the decision layer is a self-inflicted multiple-testing burden (every added target raises J_final via the max).

---

## 4. "Arbitrarily large J makes tiny differences significant" is a feature, not a bug (IJDA #1)

A common objection is that a large J will make an irrelevant difference "significant". Two separations resolve it:

1. **Testing and estimation are different gates.** A test at a fixed α only rejects if the estimate is far from the null *on the scale of its own standard error*. With large J a small true effect will be detected. That is exactly what a well-powered test should do: detect the true effect, however small. The problem the objector wants to avoid is not "significant", it is "significant *and I didn't plan to act on that size*". The remedy is the **planning alternative / minimally-important-difference**, not a smaller J.
2. **Precision planning keeps the estimate honest.** Even when a difference is significant at large J, the paired contrast still reports the estimate and its MCSE/CI first (IJDA #1e). A significant-but-tiny θ̂ is a correctly-reported scientific finding; the reader's reaction to its *size* is not a property of J.

So: never choose J to keep a nuisance difference from being significant. Choose J to reach a **target precision** or to reach a **power at a justified δ**, then report the effect with its uncertainty. It is fine (and correct) that TISCA's output is "the estimate, its MCSE, and, if designed for it, the pass/fail at the planned error rate".

---

## 5. One DGP ≠ "generally superior" (IJDA #1c)

Power and precision are **per-DGP, per-cell** quantities. A p-value or CI computed conditional on one data generating process, one training size n, and one family of models says nothing about any other DGP. Consequences:

- Never write "Method A is generally better than B" from a single-DGP power result; write "A beat B on this DGP at this n with planned error rates".
- Extrapolation claims require the cross-DGP evidence (Tiers B/C of the case study, IJDA #12), and the pre-registered analysis plan must say so before results are seen (P3-T5(b)).
- When several DGPs are studied, plan each cell's J from its **own** pilot variance (the correct J may differ a lot across cells — that is a finding, not a bug; see plan P3-T5(c) point 2). Do not force one J on all cells to "be fair": that conflates equal sample size with equal precision.

---

## 6. Two helpings of honesty that must accompany any power claim

- **The planning alternative is an input, not a guarantee.** Power = "Pr(detect | θ = δ **and** σ_D = σ̂_D_UB)". The realised (unconditional) power across the whole procedure — including the pilot uncertainty in §3 of the spec — is a joint property that must be **measured**, not assumed (C8, P3-T2). Report it with its MCSE, and plan conservatively with the σ̂_D_UB inflation.
- **The adaptive loop changes the error rate.** If you do not use the two-stage Algorithm 1, you must have the P3-T2/P3-T3 evidence that the design you did use holds its unconditional error rates. Do not assume the internal-pilot design is at nominal (Kieser & Friede 2000; Proschan 2005).

---

## 7. The planning alternative δ and the minimally-important-difference

- The **planning alternative δ** is the value of θ at which you want power. It must not be chosen to make a desired J come out small — that is circular MDE-selection (IJDA #8).
- Give δ a **domain justification** (IJDA #8a): e.g. a PEHE margin expressed relative to the scale/SD of the true `τ(X)`, not as a bare number. The fact that "the true ATE ≈ 20" does **not** by itself justify "δ = 0.5" for a PEHE (IJDA minor #11); the margin belongs on a quantity the contrast really is (a per-unit CATE error), not on the aggregate ATE.
- If a *practical significance* claim is made (“δ = 0.5 PEHE units matters”), there must be a separate **minimally important difference** documented independently of the planning exercise. The two names must not be conflated.

---

## 8. The "do not use TISCA for X" list (mandatory content)

TISCA v2 is a *design* tool for the replication count of a Monte Carlo simulation study that compares estimands via paired, replicate-level losses. **Do not use TISCA for:**

1. **Choosing data size n or covariate design for a single real data set.** TISCA sizes the number of *simulation replications J*. It says nothing about experimental design, sample size n for a field study, or how many units to recruit. (A related but separate sample-size problem needs its own machinery.)
2. **Making "generally superior" causal claims from one DGP.** See §5. TISCA reports per-cell results; extrapolation needs multiple DGPs, not more replications of one.
3. **Testing a single run / per-replication "significance".** There is no paired contrast without a set of `j`. A single fitted model yields one `L_j`, not a distribution.
4. **As a stopping rule to get a desired headline.** Do not iterate J until "some" test turns significant. J is chosen by target; the test runs once. "Keep sampling till p < α" is a different, invalid procedure (and it is not TISCA).
5. **Stopping on MCS-cardinality** ("stop when the MCS is a singleton"). Optimistically biased; default is never. Stop on precision, then report the MCS.
6. **Selecting the best benchmark post-hoc.** MCS/SPA/Romano–Wolf are for inference *after* a declared comparison family, not a licence to search benchmarks and report the winner; benchmark identity is pre-specified once (P3-T5(b)).
7. **Comparing coverage as "lower is better".** Coverage has a target. Route through calibration deviation or interval score (C6).
8. **Any claim whose error behaviour you refuse to validate.** If you insist on the adaptive loop or "none" multiplicity, you must present the P3-T2/P3-T3 operating-characteristics evidence, not assume nominal rates.
9. **When the metric is an across-replication aggregate, not a replicate-level `L_j`.** TISCA cannot pair what was never measured per replication (spec §8 validation rejection).

When in doubt, run the precision layer alone; if you cannot say what you would do differently given a “reject” versus a “not reject”, you did not need the decision layer.

**End of P1-T3.**
