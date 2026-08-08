
### Draft for Section 3.10 (Limitations): the no-difference case

Three situations must be distinguished, and they do not behave alike. When two
distinct methods happen to have equal expected loss, the design is in its ordinary
null state: across the 216 null cells of the operating-characteristic study
the two-stage procedure rejected at 0.0494 against a
nominal 0.05 (worst cell over all designs 0.0882), and the confidence
interval for the contrast covered zero at 0.9506. The
number of replications is unaffected, because it is chosen from the planning
alternative and the pilot variance, neither of which depends on the true effect.

When the two methods are the same method, the paired contrast is identically zero,
its variance is zero, and no amount of replication can change that. The procedure
must recognise the situation rather than plan for it. TISCA v2 returns an estimate
of exactly zero, a p-value of 1, no rejection, and stops at J = 2;
the studentized bootstrap reports the sample as degenerate instead of studentizing
by zero. This is the case in which v1 misbehaved: its power was estimated from the
two marginal standard deviations and set to zero whenever either was zero, so the
target was never met and the loop consumed the entire budget
(J = 5050 here) before stopping at the cap.

The third case is the limit of the second. As the two methods become perfectly
correlated, sigma_D falls to zero continuously and the required J falls with it,
which is the pairing gain in its extreme form; at correlation exactly one with
matched marginals both design targets become vacuous, since any J attains any
precision. Here the current implementation is serviceable but not yet ideal: the
planning functions return the smallest admissible J rather than declining to plan,
and the two-stage result reports a marginal power of 1 for a contrast that is
identically zero. Neither harms the procedure, which terminates and does not
reject, but a zero pilot standard deviation is far more often a duplicated metric
column than a genuine tie, so it should be surfaced as a refusal rather than
absorbed silently. That is a recommendation on the software, and it is recorded
here as one.
