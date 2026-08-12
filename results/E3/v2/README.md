# E3 schema-`e3-v2.0` collected replications

The four confirmatory cells of the case study, 1000 rows each, collected from
`notebooks/E3_shards_v2/` by `experiments/E3_mvbcf_casestudy/collect_shards.py`
with `--table shard_table_v2.csv`.

These are the analysis inputs. `analyse_e3.py` reads this directory by default;
pass `--replications-dir ../` (that is, `results/E3`) to reproduce the pre-P5.5
analysis from the v1 collection, which is left in place unchanged.

## Relationship to the v1 files one directory up

The v1 header is a strict prefix of the v2 header: every v1 column keeps its
name, its position and its value, and 53 interval-score columns are appended.
`run_cell_v2.R` draws no random numbers `run_cell.R` did not, in the same order,
from the same L'Ecuyer substreams, under the same fixed seeds `0..999`, so the
shared columns had to return identical values.

They did. `check_v1_v2_parity.py` compared all 33 shards, 4000 seeds and 580,000
recorded cells and found every one bit-identical, with a worst relative
deviation of exactly zero. The log is `results/E3/v1_v2_parity.json`. Every
analysis table except the new interval-score tables is byte-identical whether it
is computed from the v1 or the v2 files, which is checked by recomputing both.

The rerun happened months after the original, on different machines, in
different sessions, so this is a genuine independent replication of the campaign
rather than a re-read of the same output.

## What the appended columns are

Per model in {`mvbcf`, `bcf`, `bart`, `mvbart`}, level L in {50, 95} and outcome
k in {1, 2}:

```
<model>_is<L><k>       mean unit-level CATE Winkler interval score
<model>_pen<L><k>      its non-coverage penalty; is == wid + pen exactly, and
                       pen/(2/a) is the mean miss distance
<model>_ate_is<L><k>   interval score of the ATE credible interval
```

plus `tau_true_mean1`, `tau_true_mean2`, `width_audit_max_abs_dev`,
`schema_version` and `interval_levels`. `width_audit_max_abs_dev` is the
driver's internal check that the score's width component and the separately
recorded `wid<L><k>` agree; it prints at the 1e-15 level.

Only the four full 1000-row files are kept here. The pilot/confirmatory split at
`J0 = 100` is applied at analysis time by `analyse_e3.py`, so the split files are
derived rather than stored.
