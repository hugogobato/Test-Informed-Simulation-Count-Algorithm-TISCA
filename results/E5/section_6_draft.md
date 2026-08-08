
### Draft for Section 6: a demonstration outside causal inference

Nothing in the design layer of Section 3 refers to treatment effects. To make that
concrete, and to make it checkable, we apply the identical workflow to one-step-ahead
forecasting of a stationary AR(1) process, where the competing forecasters can be
ranked in closed form. With y_t = 0.7 y_(t-1) + e_t and e_t ~ N(0, 1.0), the
population one-step-ahead mean squared errors are 1.000 for an
AR(1) forecast at the true parameter, 1.176 for the random walk
and 1.961 for the sample-mean forecast; an AR(2) fitted to AR(1)
data is asymptotically tied with the AR(1) but pays one further estimation penalty of
approximately sigma^2 / n at any finite training length, here
0.0025. Before seeing any data we declared a practical
margin of 0.05 on the loss scale, five per cent of the innovation variance, and used
it both as the planning alternative and as the margin of the final claim.

The pilot of 50 independent replications, after variance-uncertainty inflation and
Bonferroni planning at alpha/3, required J = 1162 replications, set by the
precision target rather than the power target and driven by the widest contrast
(AR(1) against the sample mean). On the confirmatory block the ranking of all four
forecasters was recovered exactly, every estimate agreed with its finite-sample
reference to within 0.0027, and every Monte Carlo
standard error met its target.

The AR(1)-versus-AR(2) contrast is the instructive one. Its estimate is
-0.0027 with a Monte Carlo interval of
[-0.0033, -0.0021], and the two-sided test of
equality rejects overwhelmingly (p = 1.5e-24). It would be
wrong to report that as evidence that AR(1) is the better forecaster in any useful
sense: the difference is 5% of the margin we
declared to matter, and the minimum-effect test against that margin does not reject
(p = 1.000). The two genuinely inferior forecasters are
rejected at the margin (p < 0e+00),
and the Model Confidence Set eliminates both. This is the distinction of Section 3.4
made visible on a problem where the truth is known: a sufficiently large J makes any
non-zero difference significant, so the design must be sized for precision, and the
claim must be tested against a margin declared in advance.

The point is not that forecasting is a novel application. It is that the estimand
table, the paired contrasts, the two-stage sizing and the multiple-comparison layer
transfer without modification, and that on a problem whose answer is known they return
it -- including the part of the answer that says a statistically certain difference is
not worth acting on.
