# P4-T1 precision gain from 500 to 1000 replications

This ancillary diagnostic reads the complete 1000-row E3 confirmatory blocks already collected. It performs no model fits and does not replace the amended final inference, which discards seeds 0..99 and uses seeds 100..999.

For each of the six primary paired PEHE contrasts, `J = 500` uses the first 500 seed-sorted rows and `J = 1000` uses all 1000 rows. The expected MCSE ratio (`MCSE_500 / MCSE_1000`) is `sqrt(1000/500) = 1.4142`, equivalent to a 29.2893% reduction in MCSE.

| cell | mean MCSE ratio | MCSE reduction range | mean reduction | max absolute estimate change |
|---|---:|---:|---:|---:|
| DGP1 n=500 | 1.4029 | 24.9657%--31.6811% | 28.6354% | 0.2311 |
| DGP2 n=500 | 1.4293 | 27.1914%--32.1019% | 30.0014% | 0.1254 |
| DGP3 n=500 | 1.4347 | 26.9138%--33.0414% | 30.2422% | 0.1540 |
| DGP1 n=100 | 1.4250 | 27.5241%--32.0497% | 29.7834% | 0.1179 |

Across all 24 cell-by-contrast rows, the observed MCSE ratio ranges from `1.3327` to `1.4935`, with mean `1.4230`. The observed MCSE reduction ranges from `24.9657%` to `33.0414%`, with mean `29.6656%`. This is consistent with the square-root precision law, up to finite-sample variation.

The row-level details, including estimates, MCSEs, 95% precision half-widths, and the expected benchmark, are in `precision_gain.csv`.
