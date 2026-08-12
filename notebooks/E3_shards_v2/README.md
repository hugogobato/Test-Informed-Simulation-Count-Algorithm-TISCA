# E3 schema-`e3-v2.0` rerun

> **STATUS: COMPLETE (2026-08-10).** All 33 shards were run and downloaded, and
> every one returned the verdict `BIT_IDENTICAL`. Across 4000 seeds and 580,000
> recorded cells the v1 and v2 campaigns agree exactly; the worst relative
> deviation is 0 and the width audit is 0. The collected replications are in
> `results/E3/v2/`, the audit log is `results/E3/v1_v2_parity.json`, and the
> exact interval-score analysis is in `results/E3/mcs_interval_score.csv`,
> `interval_score_decomposition.csv` and `interval_score_contrasts.csv`. The run
> order and descoping guidance below are retained as the record of how the
> campaign was planned.


33 pre-filled Colab notebooks that rerun the confirmatory E3 campaign with
`experiments/E3_mvbcf_casestudy/run_cell_v2.R`, plus the manifest
`experiments/E3_mvbcf_casestudy/shard_table_v2.csv`. They are the deliverable of
plan task P5.5-T1. The v1 directory `notebooks/E3_shards/` is not touched and
its CSVs remain the record of the first campaign.

## Why the rerun exists

The v1 driver stored, per replication and per model, the mean credible-interval
width and the mean coverage proportion. It did not store the miss distances, so
the Winkler interval score

```
IS_a(l, u; y) = (u - l) + (2/a)(l - y)1{y < l} + (2/a)(y - u)1{y > u}
```

cannot be reconstructed from the committed shards at any level. That score is
what the revision needs in two places. Reviewer 2's third point is that coverage
has a target rather than a direction, and that interval quality should be judged
by deviation from the nominal level together with interval length; the interval
score is exactly that comparison collapsed into one number. And
`REVISION_PLAN.md` §2.3 requires a scalar, lower-is-better, per-replication loss
for the model confidence set, which raw coverage cannot supply and which the
interval score can, because it is a proper scoring rule (Gneiting & Raftery
2007, §6.2). Without it the uncertainty comparison falls back to CRPS, which is
proper but scores the whole predictive distribution and never refers to the
declared 95% level at all.

The plan's own superset rule (§4.5, "record a superset; extra columns cost
nothing, a missing column costs a full round") should have caught this before
the campaign launched. It did not, and this directory is the cost.

## What the rerun also buys

The interval score is a deterministic function of posterior draws that v1
already computed. `run_cell_v2.R` therefore draws no random numbers that v1 did
not, in the same order, from the same L'Ecuyer substreams, under the same fixed
seeds `0..999`. Every v1 column must come back **identical**.

That turns a metric patch into a full replication of the v1 campaign, run months
later, on different Colab machines, in different sessions. The final cell of
every notebook checks it: it downloads the matching v1 CSV from this repository
and compares, seed by seed, the RNG state hashes, the model seeds and every v1
metric. The original study seeded from `as.numeric(Sys.time())` and cannot make
a statement of this kind about its own results, so the audit is worth reporting
in its own right rather than being kept as an internal check.

Timing and provenance columns (`hostname`, `git_sha`, `session_hash`,
`fit_seconds_*`, `replication_seconds`) are expected to differ and are reported
but never asserted.

## What changed in the schema

Every v1 column keeps its name, its value and its position; the 53 new columns
are appended on the right, so the v1 header is a strict prefix of the v2 header.
Per model in {`mvbcf`, `bcf`, `bart`, `mvbart`}, level L in {50, 95} and outcome
k in {1, 2}:

```
<model>_is<L><k>       mean unit-level CATE interval score
<model>_pen<L><k>      its non-coverage penalty; score = wid + pen, and
                       pen/(2/a) is the mean miss distance, the one quantity
                       v1 failed to record
<model>_ate_is<L><k>   interval score of the ATE credible interval
```

plus `tau_true_mean1`, `tau_true_mean2` (the true test-set ATEs, which make the
ATE scores auditable without the `ate - bias` identity),
`width_audit_max_abs_dev`, `schema_version` and `interval_levels`.

`width_audit_max_abs_dev` is the internal consistency check. The score's width
component and the v1 `wid<L><k>` column are two independent code paths to the
same empirical quantile difference. The driver records their largest
disagreement over all model, outcome and level combinations in the replication,
and fails the replication if it exceeds `1e-8` times the width scale. It should
print at the 1e-15 level.

The 50% level is carried alongside the 95% level throughout. Calibration
evidence at two nominal levels is considerably harder to argue with than at one,
and it costs nothing here because v1 already recorded both.

## Scope: confirmatory only

Only the 33 confirmatory shards are regenerated. Under `ANALYSIS_PLAN.md`
AMENDMENT 1 the pilot is the first `J0 = 100` seeds of the confirmatory block,
applied as a split at analysis time, and the four Round 0 pilot notebooks are
marked `superseded` in the v1 manifest. Rerunning them would spend 200
replications that no analysis reads.

Projected cost is about 206 session-hours, the same as v1, distributed as 3
sessions across 11 accounts by the `account_slot` and `session_slot` columns of
`shard_table_v2.csv`.

**If that budget is not available**, rerun one cell rather than none. The model
confidence set is computed per cell, so a completed DGP1 n=500 (10 shards, ~62
session-hours) yields a genuine exact interval-score MCS for the primary cell,
with the CRPS substitution and a stated scope note for the other three. Prefer a
complete cell to a partial spread across cells.

## Regenerating the notebooks

```bash
python notebooks/_generators/build_e3_v2_notebooks.py \
  --bundle-folder-url 'https://drive.google.com/drive/folders/1w3quuskj25CBOFCGG0mTRGUHcufPpdb3?usp=sharing' \
  --bundle-sha256 '12d223bc0fcef624c1ff4cc35c5d7ecc1b1f9b05aa84ecd9d9e4a5a3382bae3c'
```

The generator imports the environment cells (R install, bundle restore, MVBCF
compile) from `build_e3_notebooks.py` unchanged and replaces only the
driver-download cell, so the v2 sessions differ from the v1 sessions in the
driver and nothing else. It asserts that it found exactly one such cell, so a
future edit to the v1 generator cannot silently desynchronise the two.

## Before uploading anything

`run_cell_v2.R` must be committed and pushed to `main` first. The notebooks
download the driver from the raw GitHub URL and assert that it contains both the
v1 invariants and the v2 interval-score code; a stale `main` stops the notebook
in that cell rather than spending six hours producing another v1 shard under a
v2 filename.

## Run order

There is no calibration gate this time. It was passed in v1, the model
configuration is unchanged, and the v1 parity check is a strictly stronger
statement than the calibration bands were.

1. Push `run_cell_v2.R`.
2. Run one shard end to end, `E3v2_DGP1_n500_confirmatory_shard01_seeds000-099`,
   and read its parity verdict before distributing anything else. A verdict of
   `BIT_IDENTICAL` or `IDENTICAL_TO_FP_NOISE` clears the remaining 32. A
   `MISMATCH` means the rerun is not a replication and must be diagnosed before
   more compute is spent: check the reported column first, then whether the R
   bundle SHA still matches the one v1 used.
3. Distribute the remaining 32 by `account_slot` / `session_slot`.

Each notebook is idempotent on re-upload: it keeps successful seed rows, backs
up and drops failed or malformed rows, and reruns only missing contiguous seed
ranges. A dead session loses that session's `/content` checkpoint, so the final
download is the durable copy.

## Collect and audit

Copy the downloaded `E3v2_*.csv` files into one directory, then:

```bash
python experiments/E3_mvbcf_casestudy/check_v1_v2_parity.py \
  --v2-dir /path/to/copied/E3_v2_csvs \
  --report results/E3/v1_v2_parity.json
```

This is the audit log the revision cites. It checks, in order: the v1 header is
a prefix of the v2 header; the RNG state hashes and model seeds match exactly;
every v1 metric agrees to `1e-12` relative; and the new columns are finite,
have non-negative penalties, and satisfy `is == wid + pen`. Exit status is
non-zero if any check fails, so it can gate the analysis. `--allow-missing-shards`
audits a partial rerun.

The pure-mathematics checks on the score itself are separate and need no
campaign data:

```bash
Rscript experiments/E3_mvbcf_casestudy/test_interval_score_v2.R
```

It validates the driver's helpers against a naive per-observation loop, checks
that the two matrix orientations agree, that the width component reproduces the
v1 `cred_width()` path exactly, that the penalty recovers the mean miss
distance, that a wrong orientation argument raises rather than recycling, that a
calibrated predictive distribution scores better than an over- or
under-dispersed one, and one hand-computed known answer.

## Reading the result

Lower interval score is better, and it is deliberately **not** the same ranking
as coverage. A model can cover at 0.98 by being wide and still score badly,
which is the whole reason for recording it. Report the score as the primary
uncertainty-quantification loss and the MCS loss; keep raw coverage and mean
width alongside as descriptive diagnostics, and keep the calibration estimand as
the deviation of *mean* coverage from nominal, computed across replications, not
as the mean of per-replication absolute deviations. Those two are not the same
quantity.
