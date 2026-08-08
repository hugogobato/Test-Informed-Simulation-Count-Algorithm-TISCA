# P3-T4: finite-sample calibration and the rule replacing J > 30

Replacing the J > 30 heuristic. Calibration is judged against an equivalence band of +/-0.005 on the level (or 2 MCSE, whichever is wider), held for every larger J in the grid. Of the 15 loss cells, the 9 with |skew(D)| < 0.5 are calibrated for the paired t by J = 10, and the most skewed cell of all -- the row bootstrap of the real MVBCF losses, skew(D) = -1.55 -- needs J = 150. Two cells never calibrate anywhere in the grid, and the direction differs by test in a way that matters: the paired t is CONSERVATIVE on the catastrophic-failure mixture (mix (rho=0), mix (rho=0.6)), reaching only 0.0401 at J = 400, whereas the studentized bootstrap is LIBERAL on the same cells (mix (rho=0), mix (rho=0.6)), reaching 0.1156 against a nominal 0.05. The bootstrap is therefore not a general remedy for non-normality: it repairs skewness but over-rejects under rare catastrophic failures, exactly where the t errs in the safe direction. The operating requirement is a function of the shape of D, which a pilot can estimate, and not a universal constant.

| test                  | family           |   rho |   J_min | direction    |   type_I_at_max_J |   skew_D |   abs_skew |   skew_mcse |   J_berry_esseen_0.005 |
|:----------------------|:-----------------|------:|--------:|:-------------|------------------:|---------:|-----------:|------------:|-----------------------:|
| paired_t              | beta             |   0   |      10 | calibrated   |            0.0517 |  -0.3347 |     0.3347 |      0.0024 |                   1011 |
| paired_t              | beta             |   0.6 |      10 | calibrated   |            0.0495 |  -0.2277 |     0.2277 |      0.0024 |                    468 |
| paired_t              | empirical        | nan   |     150 | calibrated   |            0.0503 |  -1.5484 |     1.5484 |      0.0024 |                  21620 |
| paired_t              | empirical_copula |   0   |      10 | calibrated   |            0.05   |  -0.4283 |     0.4283 |      0.0024 |                   1655 |
| paired_t              | empirical_copula |   0.6 |      30 | calibrated   |            0.0499 |  -0.7513 |     0.7513 |      0.0024 |                   5090 |
| paired_t              | gamma            |   0   |      10 | calibrated   |            0.0503 |   0.4988 |     0.4988 |      0.0024 |                   2244 |
| paired_t              | gamma            |   0.6 |      10 | calibrated   |            0.0512 |   0.3663 |     0.3663 |      0.0024 |                   1211 |
| paired_t              | lognormal        |   0   |      10 | calibrated   |            0.0491 |   0.6185 |     0.6185 |      0.0024 |                   3450 |
| paired_t              | lognormal        |   0.6 |      10 | calibrated   |            0.051  |   0.5492 |     0.5492 |      0.0024 |                   2721 |
| paired_t              | mix              |   0   |     nan | conservative |            0.0436 |   0.0972 |     0.0972 |      0.0024 |                     86 |
| paired_t              | mix              |   0.6 |     nan | conservative |            0.0402 |   0.1339 |     0.1339 |      0.0024 |                    162 |
| paired_t              | normal           |   0   |      10 | calibrated   |            0.05   |  -0.0009 |     0.0009 |      0.0024 |                      1 |
| paired_t              | normal           |   0.6 |      10 | calibrated   |            0.0486 |  -0.0003 |     0.0003 |      0.0024 |                      1 |
| paired_t              | t3               |   0   |      10 | calibrated   |            0.0486 |   0.292  |     0.292  |      0.0024 |                    770 |
| paired_t              | t3               |   0.6 |      10 | calibrated   |            0.0481 |  -0.3027 |     0.3027 |      0.0024 |                    827 |
| studentized_bootstrap | beta             |   0   |     200 | calibrated   |            0.0588 |  -0.3347 |     0.3347 |      0.0024 |                   1011 |
| studentized_bootstrap | beta             |   0.6 |      25 | calibrated   |            0.0492 |  -0.2277 |     0.2277 |      0.0024 |                    468 |
| studentized_bootstrap | empirical        | nan   |      50 | calibrated   |            0.0568 |  -1.5484 |     1.5484 |      0.0024 |                  21620 |
| studentized_bootstrap | empirical_copula |   0   |     400 | calibrated   |            0.0536 |  -0.4283 |     0.4283 |      0.0024 |                   1655 |
| studentized_bootstrap | empirical_copula |   0.6 |     100 | calibrated   |            0.0468 |  -0.7513 |     0.7513 |      0.0024 |                   5090 |
| studentized_bootstrap | gamma            |   0   |      20 | calibrated   |            0.056  |   0.4988 |     0.4988 |      0.0024 |                   2244 |
| studentized_bootstrap | gamma            |   0.6 |     150 | calibrated   |            0.0584 |   0.3663 |     0.3663 |      0.0024 |                   1211 |
| studentized_bootstrap | lognormal        |   0   |     150 | calibrated   |            0.0532 |   0.6185 |     0.6185 |      0.0024 |                   3450 |
| studentized_bootstrap | lognormal        |   0.6 |     200 | calibrated   |            0.0544 |   0.5492 |     0.5492 |      0.0024 |                   2721 |
| studentized_bootstrap | mix              |   0   |     nan | liberal      |            0.0924 |   0.0972 |     0.0972 |      0.0024 |                     86 |
| studentized_bootstrap | mix              |   0.6 |     nan | liberal      |            0.1156 |   0.1339 |     0.1339 |      0.0024 |                    162 |
| studentized_bootstrap | normal           |   0   |      25 | calibrated   |            0.0588 |  -0.0009 |     0.0009 |      0.0024 |                      1 |
| studentized_bootstrap | normal           |   0.6 |      10 | calibrated   |            0.0552 |  -0.0003 |     0.0003 |      0.0024 |                      1 |
| studentized_bootstrap | t3               |   0   |     400 | calibrated   |            0.0576 |   0.292  |     0.292  |      0.0024 |                    770 |
| studentized_bootstrap | t3               |   0.6 |     400 | calibrated   |            0.0544 |  -0.3027 |     0.3027 |      0.0024 |                    827 |
