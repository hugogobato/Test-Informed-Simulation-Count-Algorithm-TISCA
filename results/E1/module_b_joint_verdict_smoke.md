# Joint Module B verdict (smoke)

Status: complete (32 of 32 cells).

Module B now reports marginal and family-level operating characteristics. FWER, conjunctive power, disjunctive power, and FDR are distinct estimands. BH FDR is stored separately even under the global null, where it equals FWER.

## D4 global-null summary

```text
 K  correction  cells  marginal  mean_fwer  min_fwer  max_fwer
 1          bh      1  0.015000   0.015000     0.015     0.015
 1  bonferroni      1  0.050000   0.050000     0.050     0.050
 1        holm      1  0.035000   0.035000     0.035     0.035
 1        none      1  0.055000   0.055000     0.055     0.055
 1 romano_wolf      1  0.050000   0.050000     0.050     0.050
 3          bh      1  0.023333   0.060000     0.060     0.060
 3  bonferroni      1  0.016667   0.050000     0.050     0.050
 3        holm      1  0.011667   0.035000     0.035     0.035
 3        none      1  0.065000   0.175000     0.175     0.175
 3 romano_wolf      1  0.025000   0.050000     0.050     0.050
 6          bh      1  0.009167   0.045000     0.045     0.045
 6  bonferroni      1  0.009167   0.055000     0.055     0.055
 6        holm      1  0.006667   0.035000     0.035     0.035
 6        none      1  0.036667   0.190000     0.190     0.190
 6 romano_wolf      3  0.007222   0.041667     0.040     0.045
```

For D4, unadjusted FWER rises from about 0.055 at K=1 to 0.175 at K=3 and 0.190 at K=6. In these global-null cells, Bonferroni, Holm, BH, and Romano-Wolf remain near the family level across K.

## D4 alternative summary

```text
 K  correction  marginal  conjunctive  disjunctive
 1          bh  0.865000        0.865        0.865
 1  bonferroni  0.850000        0.850        0.850
 1        holm  0.805000        0.805        0.805
 1        none  0.890000        0.890        0.890
 1 romano_wolf  0.810000        0.810        0.810
 3          bh  0.983333        0.960        0.995
 3  bonferroni  0.920000        0.780        1.000
 3        holm  0.961667        0.920        0.995
 3        none  0.916667        0.780        0.995
 3 romano_wolf  0.865000        0.725        0.975
 6          bh  0.986667        0.925        1.000
 6  bonferroni  0.948333        0.735        1.000
 6        holm  0.987500        0.950        1.000
 6        none  0.929167        0.690        1.000
 6 romano_wolf  0.836667        0.605        0.995
```

## Construction and limitations

Every repetition contains one shared method-A loss and K separately generated benchmark losses. The paired vector is computed by the declared linear map C = [1, -I_K]. Romano-Wolf resamples the full K-column contrast block with one common row-index draw, preserving cross-contrast dependence.

The two-column empirical source cannot identify a K-benchmark joint law. Accordingly, 2 saved covariance rows are labeled synthetic Gaussian-copula extensions. K=1 row-bootstrap cells retain the observed pair exactly.

All 32 saved covariance rows pass the exact D = C L mapping; the minimum contrast-covariance eigenvalue is 0.64, and no diagnostic found repeated benchmark columns.

The canonical theta grid contains either K nulls or K alternatives, not a mixed null/alternative family. BH FDR is therefore informative at the global null and equals zero in the all-alternative cells; partial-null FDR is not identified by this factorial.

The canonical Romano-Wolf cells retain the existing planning schedule at family alpha. Their reported marginal and family power is the measured power of the final joint stepdown decisions, not a claim that a nested Romano-Wolf planner guarantees the nominal conjunctive target.

E3 Romano-Wolf, Holm, and Bonferroni results remain the case-study inference layer. No E3 model fit or E3 seed-verification run is part of this task. G3 bibliometric justification coding is out of scope.
