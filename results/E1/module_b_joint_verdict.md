# Joint Module B verdict (full grid)

Status: complete (660 of 660 cells).

Module B now reports marginal and family-level operating characteristics. FWER, conjunctive power, disjunctive power, and FDR are distinct estimands. BH FDR is stored separately even under the global null, where it equals FWER.

## D4 global-null summary

```text
 K  correction  cells  marginal  mean_fwer  min_fwer  max_fwer
 1          bh     11  0.053227   0.053227    0.0400    0.0670
 1  bonferroni     11  0.052727   0.052727    0.0445    0.0665
 1        holm     11  0.052182   0.052182    0.0415    0.0685
 1        none     11  0.052909   0.052909    0.0390    0.0735
 1 romano_wolf     11  0.054455   0.054455    0.0440    0.0680
 3          bh     11  0.020030   0.051045    0.0370    0.0780
 3  bonferroni     11  0.017909   0.051091    0.0400    0.0725
 3        holm     11  0.017621   0.048182    0.0400    0.0730
 3        none     11  0.053636   0.144273    0.1165    0.1750
 3 romano_wolf     11  0.017091   0.045500    0.0365    0.0525
 6          bh     11  0.011356   0.051500    0.0395    0.0860
 6  bonferroni     11  0.009348   0.050818    0.0350    0.0770
 6        holm     11  0.010038   0.051864    0.0350    0.0875
 6        none     11  0.051386   0.243864    0.1835    0.2965
 6 romano_wolf     11  0.009402   0.047864    0.0345    0.0575
```

For D4, unadjusted FWER rises from about 0.053 at K=1 to 0.144 at K=3 and 0.244 at K=6. In these global-null cells, Bonferroni, Holm, BH, and Romano-Wolf remain near the family level across K.

## D4 alternative summary

```text
 K  correction  marginal  conjunctive  disjunctive
 1          bh  0.846091     0.846091     0.846091
 1  bonferroni  0.839318     0.839318     0.839318
 1        holm  0.845682     0.845682     0.845682
 1        none  0.840591     0.840591     0.840591
 1 romano_wolf  0.844091     0.844091     0.844091
 3          bh  0.956894     0.897182     0.993591
 3  bonferroni  0.906621     0.775000     0.992091
 3        holm  0.950985     0.891682     0.992000
 3        none  0.899364     0.756864     0.992136
 3 romano_wolf  0.855121     0.722955     0.962227
 6          bh  0.984121     0.927182     0.999273
 6  bonferroni  0.937871     0.741682     0.998727
 6        holm  0.978697     0.923545     0.998864
 6        none  0.920182     0.681045     0.998955
 6 romano_wolf  0.853455     0.628545     0.986409
```

## Construction and limitations

Every repetition contains one shared method-A loss and K separately generated benchmark losses. The paired vector is computed by the declared linear map C = [1, -I_K]. Romano-Wolf resamples the full K-column contrast block with one common row-index draw, preserving cross-contrast dependence.

The two-column empirical source cannot identify a K-benchmark joint law. Accordingly, 240 saved covariance rows are labeled synthetic Gaussian-copula extensions. K=1 row-bootstrap cells retain the observed pair exactly.

All 660 saved covariance rows pass the exact D = C L mapping; the minimum contrast-covariance eigenvalue is 0.19, and no diagnostic found repeated benchmark columns.

The canonical theta grid contains either K nulls or K alternatives, not a mixed null/alternative family. BH FDR is therefore informative at the global null and equals zero in the all-alternative cells; partial-null FDR is not identified by this factorial.

The canonical Romano-Wolf cells retain the existing planning schedule at family alpha. Their reported marginal and family power is the measured power of the final joint stepdown decisions, not a claim that a nested Romano-Wolf planner guarantees the nominal conjunctive target.

E3 Romano-Wolf, Holm, and Bonferroni results remain the case-study inference layer. No E3 model fit or E3 seed-verification run is part of this task. G3 bibliometric justification coding is out of scope.
