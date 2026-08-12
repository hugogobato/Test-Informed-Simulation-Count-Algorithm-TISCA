# E6: calibration of the declared paired tests on the real case-study losses

Protocol identical to E1b: equivalence band max(+/-0.005, 2 MCSE) on the level, held for every larger J on the grid [10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 400]; R = 40000 for the paired t and 10000 for the studentized bootstrap with B = 499. The null is imposed by the row bootstrap of the real (M, 2) loss matrix, which preserves the real joint dependence, marginal shapes and variance ratio.

|skew(D)| across the 24 real contrasts ranges from 0.03 to 1.47; the within-replication Pearson correlation ranges from 0.74 to 0.98.

Paired t: J_min ranges from 10 to 100 (median 28); 11 of 24 contrasts require more than 30 replications and 0 never calibrate on the grid.
Studentized bootstrap: J_min ranges from 10 to 400 (median 200); 13 never calibrate.

The bootstrap's residual excess level averages +0.0048 over J >= 200. Raising the resample count from B = 499 to 1999 removes 0.0021 of it, so part of the excess is the finite-resample noise in the bootstrap quantiles rather than the bootstrap principle; the remainder is not removed by more resamples. See bootstrap_B_sensitivity.csv.

|   dgp |   n | contrast   | test                  |   J_min | direction   |   type_I_at_max_J |   skew_D |   abs_skew |   rho_pearson |   rho_spearman |   sd_ratio_A_over_B |   observed_mean_D |   observed_sd_D |
|------:|----:|:-----------|:----------------------|--------:|:------------|------------------:|---------:|-----------:|--------------:|---------------:|--------------------:|------------------:|----------------:|
|     3 | 500 | C2         | paired_t              |      10 | calibrated  |            0.05   |  -0.026  |     0.026  |        0.8896 |         0.8612 |              0.945  |           -0.1289 |          1.1378 |
|     3 | 500 | C1         | paired_t              |      10 | calibrated  |            0.0496 |   0.0327 |     0.0327 |        0.8738 |         0.8354 |              0.9264 |           -0.3432 |          1.3366 |
|     2 | 500 | C1         | paired_t              |      10 | calibrated  |            0.0518 |   0.0784 |     0.0784 |        0.9798 |         0.9774 |              0.9876 |           -1.348  |          1.0861 |
|     1 | 500 | C1         | paired_t              |      10 | calibrated  |            0.0501 |  -0.1999 |     0.1999 |        0.8728 |         0.8373 |              0.9681 |           -0.5042 |          1.216  |
|     1 | 500 | C2         | paired_t              |      10 | calibrated  |            0.0509 |  -0.203  |     0.203  |        0.887  |         0.8456 |              0.9604 |           -0.3691 |          1.1513 |
|     2 | 500 | C2         | paired_t              |      10 | calibrated  |            0.0489 |  -0.242  |     0.242  |        0.918  |         0.893  |              0.917  |           -0.3259 |          0.7874 |
|     2 | 500 | C5         | paired_t              |      10 | calibrated  |            0.0496 |  -0.2645 |     0.2645 |        0.9799 |         0.9773 |              0.9908 |           -1.436  |          1.08   |
|     2 | 500 | C3         | paired_t              |      10 | calibrated  |            0.0496 |  -0.4376 |     0.4376 |        0.9735 |         0.9713 |              0.9861 |           -1.4381 |          1.2461 |
|     1 | 100 | C6         | paired_t              |      10 | calibrated  |            0.0502 |  -0.4699 |     0.4699 |        0.9115 |         0.8514 |              0.8508 |           -2.4374 |          2.6042 |
|     1 | 100 | C4         | paired_t              |      10 | calibrated  |            0.051  |  -0.4741 |     0.4741 |        0.9386 |         0.8884 |              0.8885 |           -1.4471 |          2.0907 |
|     1 | 100 | C1         | paired_t              |      20 | calibrated  |            0.051  |  -0.5626 |     0.5626 |        0.8844 |         0.8122 |              0.7767 |           -1.8031 |          3.3298 |
|     1 | 100 | C2         | paired_t              |      25 | calibrated  |            0.0504 |  -0.5632 |     0.5632 |        0.8794 |         0.8026 |              0.743  |           -2.0938 |          3.5524 |
|     1 | 100 | C3         | paired_t              |      30 | calibrated  |            0.0499 |  -0.7572 |     0.7572 |        0.9218 |         0.8748 |              0.8842 |           -1.8808 |          2.3794 |
|     1 | 100 | C5         | paired_t              |      40 | calibrated  |            0.0506 |  -0.9193 |     0.9193 |        0.8849 |         0.8271 |              0.8147 |           -3.1826 |          3.122  |
|     2 | 500 | C4         | paired_t              |      75 | calibrated  |            0.0491 |  -1.0709 |     1.0709 |        0.7623 |         0.7435 |              0.7684 |           -1.1178 |          1.5336 |
|     3 | 500 | C5         | paired_t              |      50 | calibrated  |            0.0511 |  -1.0942 |     1.0942 |        0.7586 |         0.7202 |              0.7836 |           -1.5102 |          2.1061 |
|     3 | 500 | C6         | paired_t              |      50 | calibrated  |            0.0519 |  -1.1937 |     1.1937 |        0.8024 |         0.7595 |              0.8757 |           -1.0202 |          1.605  |
|     2 | 500 | C6         | paired_t              |      50 | calibrated  |            0.051  |  -1.2103 |     1.2103 |        0.7981 |         0.7826 |              0.8046 |           -1.0536 |          1.3634 |
|     1 | 500 | C5         | paired_t              |      40 | calibrated  |            0.0503 |  -1.2302 |     1.2302 |        0.7792 |         0.7556 |              0.8521 |           -1.1912 |          1.7529 |
|     1 | 500 | C4         | paired_t              |      75 | calibrated  |            0.05   |  -1.2528 |     1.2528 |        0.7416 |         0.699  |              0.8322 |           -1.5444 |          1.9234 |
|     1 | 500 | C6         | paired_t              |     100 | calibrated  |            0.0517 |  -1.2635 |     1.2635 |        0.7994 |         0.7668 |              0.9155 |           -1.3242 |          1.5804 |
|     3 | 500 | C4         | paired_t              |      50 | calibrated  |            0.0498 |  -1.2748 |     1.2748 |        0.7438 |         0.708  |              0.8141 |           -1.2123 |          1.9296 |
|     1 | 500 | C3         | paired_t              |      75 | calibrated  |            0.051  |  -1.4475 |     1.4475 |        0.748  |         0.7407 |              0.8119 |           -1.578  |          1.9439 |
|     3 | 500 | C3         | paired_t              |      75 | calibrated  |            0.052  |  -1.4684 |     1.4684 |        0.7424 |         0.7345 |              0.7791 |           -1.5391 |          2.1796 |
|     3 | 500 | C2         | studentized_bootstrap |     nan | liberal     |            0.0553 |  -0.026  |     0.026  |        0.8896 |         0.8612 |              0.945  |           -0.1289 |          1.1378 |
|     3 | 500 | C1         | studentized_bootstrap |     100 | calibrated  |            0.0539 |   0.0327 |     0.0327 |        0.8738 |         0.8354 |              0.9264 |           -0.3432 |          1.3366 |
|     2 | 500 | C1         | studentized_bootstrap |     nan | liberal     |            0.0552 |   0.0784 |     0.0784 |        0.9798 |         0.9774 |              0.9876 |           -1.348  |          1.0861 |
|     1 | 500 | C1         | studentized_bootstrap |     400 | calibrated  |            0.0525 |  -0.1999 |     0.1999 |        0.8728 |         0.8373 |              0.9681 |           -0.5042 |          1.216  |
|     1 | 500 | C2         | studentized_bootstrap |     nan | liberal     |            0.0564 |  -0.203  |     0.203  |        0.887  |         0.8456 |              0.9604 |           -0.3691 |          1.1513 |
|     2 | 500 | C2         | studentized_bootstrap |     200 | calibrated  |            0.0539 |  -0.242  |     0.242  |        0.918  |         0.893  |              0.917  |           -0.3259 |          0.7874 |
|     2 | 500 | C5         | studentized_bootstrap |     nan | liberal     |            0.0565 |  -0.2645 |     0.2645 |        0.9799 |         0.9773 |              0.9908 |           -1.436  |          1.08   |
|     2 | 500 | C3         | studentized_bootstrap |     400 | calibrated  |            0.0549 |  -0.4376 |     0.4376 |        0.9735 |         0.9713 |              0.9861 |           -1.4381 |          1.2461 |
|     1 | 100 | C6         | studentized_bootstrap |     nan | liberal     |            0.0558 |  -0.4699 |     0.4699 |        0.9115 |         0.8514 |              0.8508 |           -2.4374 |          2.6042 |
|     1 | 100 | C4         | studentized_bootstrap |     nan | liberal     |            0.0562 |  -0.4741 |     0.4741 |        0.9386 |         0.8884 |              0.8885 |           -1.4471 |          2.0907 |
|     1 | 100 | C1         | studentized_bootstrap |     nan | liberal     |            0.0561 |  -0.5626 |     0.5626 |        0.8844 |         0.8122 |              0.7767 |           -1.8031 |          3.3298 |
|     1 | 100 | C2         | studentized_bootstrap |      10 | calibrated  |            0.0548 |  -0.5632 |     0.5632 |        0.8794 |         0.8026 |              0.743  |           -2.0938 |          3.5524 |
|     1 | 100 | C3         | studentized_bootstrap |     150 | calibrated  |            0.0526 |  -0.7572 |     0.7572 |        0.9218 |         0.8748 |              0.8842 |           -1.8808 |          2.3794 |
|     1 | 100 | C5         | studentized_bootstrap |     nan | liberal     |            0.0556 |  -0.9193 |     0.9193 |        0.8849 |         0.8271 |              0.8147 |           -3.1826 |          3.122  |
|     2 | 500 | C4         | studentized_bootstrap |     200 | calibrated  |            0.0532 |  -1.0709 |     1.0709 |        0.7623 |         0.7435 |              0.7684 |           -1.1178 |          1.5336 |
|     3 | 500 | C5         | studentized_bootstrap |      75 | calibrated  |            0.0527 |  -1.0942 |     1.0942 |        0.7586 |         0.7202 |              0.7836 |           -1.5102 |          2.1061 |
|     3 | 500 | C6         | studentized_bootstrap |     nan | liberal     |            0.0551 |  -1.1937 |     1.1937 |        0.8024 |         0.7595 |              0.8757 |           -1.0202 |          1.605  |
|     2 | 500 | C6         | studentized_bootstrap |     400 | calibrated  |            0.0534 |  -1.2103 |     1.2103 |        0.7981 |         0.7826 |              0.8046 |           -1.0536 |          1.3634 |
|     1 | 500 | C5         | studentized_bootstrap |     nan | liberal     |            0.0581 |  -1.2302 |     1.2302 |        0.7792 |         0.7556 |              0.8521 |           -1.1912 |          1.7529 |
|     1 | 500 | C4         | studentized_bootstrap |     nan | liberal     |            0.0578 |  -1.2528 |     1.2528 |        0.7416 |         0.699  |              0.8322 |           -1.5444 |          1.9234 |
|     1 | 500 | C6         | studentized_bootstrap |     nan | liberal     |            0.0581 |  -1.2635 |     1.2635 |        0.7994 |         0.7668 |              0.9155 |           -1.3242 |          1.5804 |
|     3 | 500 | C4         | studentized_bootstrap |     200 | calibrated  |            0.0533 |  -1.2748 |     1.2748 |        0.7438 |         0.708  |              0.8141 |           -1.2123 |          1.9296 |
|     1 | 500 | C3         | studentized_bootstrap |     400 | calibrated  |            0.0539 |  -1.4475 |     1.4475 |        0.748  |         0.7407 |              0.8119 |           -1.578  |          1.9439 |
|     3 | 500 | C3         | studentized_bootstrap |     nan | liberal     |            0.0568 |  -1.4684 |     1.4684 |        0.7424 |         0.7345 |              0.7791 |           -1.5391 |          2.1796 |
