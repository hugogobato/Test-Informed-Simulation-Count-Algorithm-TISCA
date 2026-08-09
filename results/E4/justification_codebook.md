# P5.5-T2 replication-number justification codebook

## Scope and denominator

This is an availability-stratified census of the pre-specified 100-paper corpus
used in P3-T1. The primary bibliometric analysis retains its original
denominator of 99 records with a numeric `J`. The optional justification study
is restricted to records satisfying both conditions: a numeric `J` is present,
and the paper has accessible full text as either an arXiv preprint or a journal
article classified as gold or diamond open access in the existing DOAJ-backed
venue table. This produces 17 eligible records. The other 82 numeric-`J`
records are retained in the row-level output with
`eligibility_status=not_coded_inaccessible`; they are not coded as
unjustified. The one screened record without a numeric `J` has
`eligibility_status=not_applicable_J_not_reported`.

The accessible records were read at the simulation-methods and appendix level
where available. A design or data-generating-process rationale was not counted
as a replication-count rationale unless it addressed the choice of `J` itself.
Likewise, varying sample sizes, data-generating processes, priors, matching
criteria, or MCMC iterations was not counted as sensitivity to `J`.

## Fields and permitted values

`corpus_row` is the one-based row in
`results/E4/bibliometric_coded.csv`. It is the merge key because the mechanical
corpus contains duplicate cite keys. `access_basis` records why a paper was or
was not in the accessible stratum. `eligibility_status` records whether the
paper entered the justification analysis. `justification_eligible` is `Y` only
for the 17 readable numeric-`J` records.

`replication_unit` records the unit repeated by the paper, such as outer
simulated data sets, simulation seeds, matched data sets, or collections of
data sets. `stated_count` is the corpus `J` value. `reported_count_scope`
records whether the value applies per scenario, per DGP, per experiment, or
another scope. `other_counts_reported` records other repeated quantities that
could otherwise be confused with `J`, including test samples, posterior draws,
tuning runs, MCMC iterations, or a second scenario-specific count.

`justification_class` uses the following mutually exclusive categories:

1. `explicit`: the paper states a count-specific criterion such as power,
   precision, pilot variance, precedent, or computational budget.
2. `implicit_or_convention`: the paper gives a recognizable convention or
   cited precedent for the count without a formal count-specific calculation.
3. `unjustified`: the count is reported, but no count-specific reason is
   stated or cited in the accessible source.
4. `unclear_report`: the source or corpus mapping does not establish which
   count should be coded, or the report is internally/version-wise inconsistent.
5. `not_applicable`: the paper was not eligible for this reading analysis,
   either because full text was inaccessible or because numeric `J` was not
   reported. The associated `eligibility_status` identifies which case holds.

`justification_status` keeps a missing rationale separate from missing or
unavailable reporting. `justification_given` is a coarser field with values
`yes`, `no`, `unclear`, and `not_applicable`; the more informative
`justification_class` is used for the reported percentages. `justification_criterion`
uses `power`, `precision`, `pilot_variance`, `precedent`,
`computational_budget`, `other`, `none`, `unclear`, or `not_applicable`.
`source_cited` records the cited source supporting the count rationale, or
`none` when no such source was found.

`sensitivity_to_J` uses `yes`, `no`, `unclear`, and `not_applicable`. It is
`yes` only when the paper varies `J` itself or directly studies Monte Carlo
error as a function of `J`. Sensitivity to a DGP, prior, ICC, sample size,
matching rule, or estimator is not sufficient. `sensitivity_evidence` records
the basis for the code.

`source_url`, `source_location`, and `source_version_checked` make the reading
auditable. `justification_quote` is a short count-report quotation, not a
substitute for the classification rule. `coding_notes` records ambiguities,
secondary repetition, and reasons for not treating a design rationale as a
count rationale. `double_coded` is `N` for every row because no independent
second coder was available for this pass; this is reported as a limitation,
not hidden as an agreement check.

## Analysis and uncertainty intervals

Run:

```text
python experiments/E4_bibliometrics/code_justifications.py
```

The script reads the mechanical base, joins the manual coding by
`corpus_row`, and writes `bibliometric_justifications_coded.csv`. It then
writes `justification_summary.csv`. Every percentage in the summary is
computed from the row-level output. The intervals are two-sided 95% Wilson
intervals for a binomial proportion, reported descriptively. They are not
adjustments for article-selection or access bias and should not be interpreted
as estimates for the 82 inaccessible papers.

The accessible-stratum results are 0/17 explicit, 0/17 implicit or
convention-based, 16/17 unjustified, and 1/17 unclear. The unclear record is
the arXiv paper whose corpus row records `J=50` but whose accessible arXiv v3
reports 100 collections per scenario. No accessible paper reports a direct
sensitivity analysis varying `J` in this coding pass. These results are
descriptive evidence about reporting practice in the readable subset, not a
claim that a universal replication number exists.
