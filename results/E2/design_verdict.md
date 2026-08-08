# P3-T3 verdict: is the iterative loop worth it?

Use the two-stage design (D4) unless a larger budget is acceptable in itself: the loop's extra power (+0.069) is bought with 2.06x the replications, so it is a bigger design rather than a better one. At the planning alternative D4 achieves power 0.840 (target 0.8) at E[J] = 73, against the oracle's 0.816 at E[J] = 64; its unconditional Type I error is 0.0494 (nominal 0.05, MCSE 0.0031). The adaptive loop's effect on the level is +0.0033 (within 2 MCSE), and it spends 2.06x the replications for +0.069 power.

TISCA v1 (D3) is CONSERVATIVE, not liberal: its level falls to 0.0004 at rho = 0.9 against a nominal 0.05, and it spends 2.39x the replications of v2. Fixing the paired-design error therefore buys efficiency, not validity.

The precision target (D5) is far more demanding than the power target at these settings: 8.5x the replications of D4, hitting the budget cap in 17.1% of repetitions. Which target binds is a reporting obligation, not a detail.

## Design table (Module A)

| design   | label                                           |   n_cells |   type_I_mean |   type_I_max |   type_I_mcse |   n_cells_over_level |   ci_cover_mean |   power_mean |   power_min |   n_cells_under_power |   E_J_null |   E_J_alt |   sd_J_alt |   q95_J_alt |   P_J_at_cap |   bias_alt |   rmse_alt | scope                   |
|:---------|:------------------------------------------------|----------:|--------------:|-------------:|--------------:|---------------------:|----------------:|-------------:|------------:|----------------------:|-----------:|----------:|-----------:|------------:|-------------:|-----------:|-----------:|:------------------------|
| D1       | D1 fixed-J, independent pilot (raw pilot sigma) |       144 |        0.0497 |       0.0736 |        0.0031 |                    5 |          0.9503 |       0.782  |      0.6954 |                    20 |    61.7388 |   62.1478 |    42.0299 |    152.596  |       0.0034 |     0.0004 |     0.1968 | Module A (all families) |
| D2       | D2 internal-pilot re-estimation (adaptive)      |       144 |        0.0523 |       0.0592 |        0.0031 |                    5 |          0.9477 |       0.9258 |      0.8624 |                     0 |    89.2697 |   89.4356 |    46.2281 |    197.211  |       0.0056 |    -0      |     0.1362 | Module A (all families) |
| D3       | D3 TISCA v1: unpaired Welch, iterative          |       144 |        0.0314 |       0.0882 |        0.0025 |                    7 |          0.8966 |       0.8684 |      0.7612 |                     4 |    80.7025 |   80.9368 |    42.8725 |    175.597  |       0.0036 |     0.0005 |     0.1453 | Module A (all families) |
| D4       | D4 TISCA v2: paired, two-stage (Algorithm 1)    |       144 |        0.0494 |       0.0762 |        0.0031 |                    2 |          0.9506 |       0.8404 |      0.7286 |                     4 |    73.3357 |   73.2875 |    47.0905 |    180.683  |       0.0048 |     0.001  |     0.1814 | Module A (all families) |
| D5       | D5 paired fixed-precision (MCSE target)         |       144 |        0.0485 |       0.0554 |        0.003  |                    0 |          0.9515 |       0.9958 |      0.8936 |                     0 |   563.973  |  564.424  |   126.462  |    740.706  |       0.1714 |     0.0001 |     0.0595 | Module A (all families) |
| D6       | D6 oracle fixed-J (true sigma known)            |       144 |        0.0497 |       0.071  |        0.0031 |                    3 |          0.9503 |       0.816  |      0.7498 |                     3 |    64.2222 |   64.2222 |     0      |     64.2222 |       0      |     0.0007 |     0.1715 | Module A (all families) |

## Matched comparisons

| comparison   | question                                         |   d_type_I_mean |   d_type_I_max |   d_power_mean |   d_power_worst |   d_E_J_mean |   E_J_ratio |   n_matched |
|:-------------|:-------------------------------------------------|----------------:|---------------:|---------------:|----------------:|-------------:|------------:|------------:|
| D2 - D4      | adaptive loop vs two-stage (Reviewer 2 par. 4)   |          0.0033 |         0.0174 |         0.0687 |         -0.0076 |      19.6231 |      2.0608 |         135 |
| D3 - D4      | TISCA v1 (unpaired, iterative) vs TISCA v2       |         -0.0165 |         0.0762 |         0.022  |         -0.1178 |      15.2177 |      2.3864 |         135 |
| D1 - D4      | raw-pilot fixed-J vs variance-inflated two-stage |          0.0006 |         0.01   |        -0.0605 |         -0.0808 |     -13.2487 |      0.841  |          54 |
| D5 - D4      | precision target vs power target                 |         -0.0007 |         0.0218 |         0.1488 |          0.0828 |     525.598  |      8.4714 |          54 |
| D4 - D6      | two-stage vs the oracle it approximates          |         -0.0007 |         0.019  |         0.0282 |         -0.0914 |      11.8614 |      1.1795 |          54 |
