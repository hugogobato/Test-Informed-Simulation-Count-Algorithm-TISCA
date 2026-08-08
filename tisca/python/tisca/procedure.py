"""Two-stage default and optional adaptive procedure — `procedure.py`.

Implements **Algorithm 1** (two-stage default, §5) and **Algorithm 2** (adaptive /
internal-pilot re-estimation, optional, §6) of `docs/tisca_v2_spec.md`.

Algorithm 1 is the default: an *independent-seed* pilot of size ``J0`` is used only
to estimate the contrast variances; the planned ``J`` is solved in closed form;
the pilot rows are **discarded** from final inference; the confirmatory run of
``J`` balanced replications is tested **once** at the adjusted ``alpha_adj``.
Because ``J`` is a function of the pilot alone and the confirmatory sample is
independent, the unconditional Type I error is exactly ``alpha_adj`` by iterated
expectation.

Algorithm 2 reuses the pilot rows in the final test (its distinguishing and risky
feature). The spec and plan flag that its error rates are *not* nominal and must
be measured by the P3-T2 outer-MC study; this module merely orchestrates the
single run and reports the design caveats.

Both algorithms are pure orchestrators: they call an injected ``sim_func(seed) ->
array`` (shape ``(replications, n_metrics)`` or a list) for every replication,
so they are trivially grid-searchable by the outer-MC harness (P2-T3).
"""

from __future__ import annotations

import numpy as np

from . import multiplicity as _mply
from . import planning as _plan
from . import validate

__all__ = [
    "TwoStageDesign",
    "AdaptiveDesign",
    "plan_and_run",
    "two_stage",
    "adaptive",
    "run_confirmatory",
]


class TwoStageDesign:
    """Algorithm 1: independent pilot, closed-form planning, discarding pilot.

    Parameters
    ----------
    sim_func : callable
        ``sim_func(seed) -> (n_rows, n_metrics)`` array of per-replication metric
        vectors, or a list of row vectors. Row ``j`` must contain every metric of
        every method needed by the declared contrasts.
    primary_contrasts : list of dict
        Each dict declares one primary contrast with keys ``A``, ``B`` (column
        indices or names into the metric rows), ``mode`` (M1..M5), ``delta``,
        and any target kwargs (``target_mcse``, ``halfwidth``, ``target_power``).
    J0 : int
        Pilot (stage-1) size.
    gamma : float
        Variance-uncertainty inflation tail probability (§3), default 0.20.
    alpha : float
        Nominal test level; divided by ``K`` for the adjusted level unless an
        explicit ``alpha_adj`` is given.
    Correction : str or None
        Multiplicity correction used to set ``alpha_adj``.
    J_max : int
        Budget cap on the confirmatory size.
    seed_pilot_base, seed_conf_base : int
        Seed block bases (pilot and confirmatory blocks are disjoint per the
        seed/RNG protocol).
    n_metrics : int
        Number of columns that ``sim_func`` returns per replication.
    """

    def __init__(
        self,
        sim_func,
        primary_contrasts,
        *,
        J0: int,
        gamma: float = 0.20,
        alpha: float = 0.05,
        correction="bonferroni",
        J_max: int = _plan._DEFAULT_J_MAX,
        seed_pilot_base: int = 1_000_001,
        seed_conf_base: int = 1,
        n_metrics: int | None = None,
        success_criterion: str = "conjunctive",
        r_bh: int | None = None,
    ) -> None:
        self.sim_func = sim_func
        self.contrasts = primary_contrasts
        self.J0 = int(J0)
        self.gamma = float(gamma)
        self.alpha = float(alpha)
        self.correction = correction
        self.J_max = int(J_max)
        self.seed_pilot_base = seed_pilot_base
        self.seed_conf_base = seed_conf_base
        self.n_metrics = n_metrics
        self.success_criterion = success_criterion
        self.r_bh = r_bh
        self.K = len(primary_contrasts)
        if self.K < 1:
            raise validate.ValidationError("At least one primary contrast is required.")
        self.alpha_adj, self._note = _mply.planning_alpha(
            correction, self.K, alpha=self.alpha, r=r_bh
        )

    def _run_rows(self, base, n):
        """Run ``sim_func`` over ``n`` consecutive seeds and stack the rows."""
        rows = []
        for j in range(n):
            row = self.sim_func(base + j)
            arr = np.asarray(row, dtype=float)
            if arr.ndim == 0:
                arr = arr.reshape(1, 1)
            elif arr.ndim == 1:
                arr = arr.reshape(1, -1)
            rows.append(arr)
        return np.concatenate(rows, axis=0) if rows else np.empty((0, self.n_metrics or 0))

    def _metric(self, row, key):
        """Get a per-replication column from a metric row dict or index."""
        idx = row[key]
        return idx

    def run(self, verbose: bool = True) -> dict:
        """Execute the two-stage design end to end."""
        if self.n_metrics is None:
            raise validate.ValidationError("n_metrics must be provided (columns in sim rows).")

        # Stage 1: independent pilot (discarded from final inference).
        pilot = self._run_rows(self.seed_pilot_base, self.J0)

        Js = []
        plan_rows = []
        for c in self.contrasts:
            a = pilot[:, self._metric(c, "A")]
            b = pilot[:, self._metric(c, "B")]
            # Use the LISTWISE-DELETED pair (C1). Calling the validator and then
            # discarding its return value left any NaN in place, so s_D came
            # back NaN and the planner silently mis-sized J.
            a, b, dropped = validate.validate_contrast_pair(a, b, name=c.get("name", "contrast"))
            D = a - b
            sd = float(np.std(D, ddof=1))
            mode = c.get("mode", "M1")
            delta = c.get("delta", 0.0)
            margin = c.get("margin")
            J, sigma_ub = _plan.required_J(
                sd,
                self.J0,
                gamma=self.gamma,
                mode=mode,
                delta=delta,
                target_mcse=c.get("target_mcse"),
                halfwidth=c.get("halfwidth"),
                alpha=self.alpha_adj,
                target_power=c.get("target_power"),
                margin=margin,
                J_max=self.J_max,
                exact=c.get("exact", False),
            )
            Js.append(J)
            plan_rows.append(
                dict(name=c.get("name", "contrast"), sd=sd, sigma_ub=sigma_ub, J=J, n_dropped=dropped)
            )

        J_final = _plan.combine_J(Js, self.J_max)
        capped = J_final >= self.J_max

        # Stage 2: confirmatory run (the pilot is NOT reused).
        conf = self._run_rows(self.seed_conf_base, J_final)

        # Final inference: run each pre-specified test once at alpha_adj, in the
        # SAME mode it was planned in (C4 / spec §8.4).
        results = []
        rejected = []
        for c in self.contrasts:
            a = conf[:, self._metric(c, "A")]
            b = conf[:, self._metric(c, "B")]
            a, b, dropped = validate.validate_contrast_pair(a, b, name=c.get("name", "contrast"))
            D = a - b
            mode = c.get("mode", "M1")
            alpha_c = self.alpha_adj
            t_res, reject = self._final_test(D, mode, c.get("margin"), alpha_c)
            results.append({**t_res, "mode": mode, "n_dropped": dropped})
            rejected.append(bool(reject))

        # Family-level reporting. Marginal power is the PLANNED power at delta
        # and the planned J, evaluated in the contrast's own mode -- not the
        # realised rejections (an earlier version returned a vector of 1.0s,
        # which read as "100% power" in the output).
        #
        # The `if sigma_ub > 0 else 1.0` short circuit that used to guard this call
        # reintroduced the same defect for degenerate contrasts: a pilot with zero
        # observed contrast variance was reported at power 1.0 whatever delta was,
        # including delta = 0, where nothing is detectable at all. The mode-specific
        # power functions now handle sigma = 0 themselves, so the call is
        # unconditional and each degenerate contrast is FLAGGED instead
        # (spec §8.5 asks for a graceful path *and* for the cell to be reported).
        marginal_power = np.array(
            [
                _plan.power_function(
                    c.get("mode", "M1"),
                    J_final,
                    c.get("delta", 0.0),
                    pr["sigma_ub"],
                    margin=c.get("margin"),
                    alpha=self.alpha_adj,
                )
                for c, pr in zip(self.contrasts, plan_rows)
            ],
            dtype=float,
        )
        for pr in plan_rows:
            pr["degenerate"] = bool(pr["sigma_ub"] <= 0.0)
        degenerate = [pr["name"] for pr in plan_rows if pr["degenerate"]]
        conjunctive = bool(all(rejected)) if self.success_criterion == "conjunctive" else None
        disjunctive = bool(any(rejected)) if self.success_criterion == "disjunctive" else None

        return {
            "J_final": J_final,
            "Js_by_contrast": plan_rows,
            "J_max": self.J_max,
            "capped": capped,
            "alpha_adj": self.alpha_adj,
            "correction": self.correction,
            "correction_note": self._note,
            "contrast_results": results,
            "rejected": rejected,
            "family_rejected_conjunctive": conjunctive,
            "family_rejected_disjunctive": disjunctive,
            "marginal_power": marginal_power,
            "degenerate_contrasts": degenerate,
            "design": "D4",
            "gamma": self.gamma,
            "success_criterion": self.success_criterion,
        }

    @staticmethod
    def _final_test(D, mode, margin, alpha):
        """Run the pre-specified final test for one contrast, in its own mode.

        Each mode has its own null VALUE, not just its own sidedness. Testing
        every mode against ``mu0 = 0`` (as an earlier version did) silently
        turns M3 into M2, M4 into a superiority test, and M5 into a two-sided
        equality test -- which is precisely the planning/testing misalignment
        C4 exists to forbid.

        * M1 two-sided,        ``H0: theta = 0``
        * M2 lower one-sided,  ``H0: theta >= 0``
        * M3 lower one-sided,  ``H0: theta >= -Delta``  (centre on ``-Delta``)
        * M4 lower one-sided,  ``H0: theta >=  Delta``  (centre on ``+Delta``)
        * M5 TOST,             ``H0: |theta| >= Delta`` (both arms at alpha)
        """
        from .inference import paired_t

        mode = (mode or "M1").upper()
        if mode == "M1":
            res = paired_t(D, 0.0, alternative="two-sided")
            return res, res["p_value"] <= alpha
        if mode == "M2":
            res = paired_t(D, 0.0, alternative="less")
            return res, res["p_value"] <= alpha
        if mode in ("M3", "M4"):
            if margin is None:
                raise validate.ValidationError(f"Mode {mode} requires a margin.")
            mu0 = -float(margin) if mode == "M3" else float(margin)
            res = paired_t(D, mu0, alternative="less")
            res["mu0"] = mu0
            return res, res["p_value"] <= alpha
        if mode == "M5":
            if margin is None:
                raise validate.ValidationError("Mode M5 requires a margin.")
            m = float(margin)
            lo = paired_t(D, -m, alternative="greater")   # H0: theta <= -Delta
            up = paired_t(D, m, alternative="less")       # H0: theta >=  Delta
            reject = (lo["p_value"] <= alpha) and (up["p_value"] <= alpha)
            res = dict(up)
            res.update(
                p_value=max(lo["p_value"], up["p_value"]),  # TOST p = max of the arms
                p_lower=lo["p_value"], p_upper=up["p_value"],
                margin=m, alternative="TOST",
            )
            return res, reject
        raise validate.ValidationError(f"Unknown mode {mode!r} in the final test.")


class AdaptiveDesign:
    """Algorithm 2: internal-pilot re-estimation (optional, validation required).

    Distinctive and risky: the pilot rows are **reused** in the final test. This
    function reports that reuse explicitly and does not cache the adaptive error
    rates as nominal. It is provided only for the P3-T2 operating-characteristics
    study and for callers who accept that validation burden.
    """

    def __init__(
        self,
        sim_func,
        primary_contrasts,
        *,
        J0: int,
        gamma: float = 0.20,
        alpha: float = 0.05,
        correction="bonferroni",
        J_max: int = _plan._DEFAULT_J_MAX,
        nmax_looks: int = 5,
        seed_base: int = 1,
        n_metrics: int | None = None,
        r_bh: int | None = None,
    ) -> None:
        self.sim_func = sim_func
        self.contrasts = primary_contrasts
        self.J0 = int(J0)
        self.gamma = float(gamma)
        self.alpha = float(alpha)
        self.correction = correction
        self.J_max = int(J_max)
        self.nmax_looks = int(nmax_looks)
        self.seed_base = seed_base
        self.n_metrics = n_metrics
        self.r_bh = r_bh
        self.K = len(primary_contrasts)
        self.alpha_adj, _ = _mply.planning_alpha(correction, self.K, alpha=alpha, r=r_bh)

    def run(self) -> dict:
        """Run Algorithm 2: accumulate, re-estimate up to ``nmax_looks`` times."""
        if self.n_metrics is None:
            raise validate.ValidationError("n_metrics must be provided.")
        current = self.J0
        look = 0
        while look < self.nmax_looks:
            accumulated = self._run_rows(self.seed_base, current)
            # Re-estimate the required closed-form J from the current sample.
            req = []
            for c in self.contrasts:
                a = accumulated[:, self._metric(c, "A")]
                b = accumulated[:, self._metric(c, "B")]
                a, b, _ = validate.validate_contrast_pair(a, b, name=c.get("name", "contrast"))
                D = a - b
                sd = float(np.std(D, ddof=1))
                J, _ = _plan.required_J(
                    sd,
                    current,
                    gamma=self.gamma,
                    mode=c.get("mode", "M1"),
                    delta=c.get("delta", 0.0),
                    margin=c.get("margin"),
                    target_mcse=c.get("target_mcse"),
                    halfwidth=c.get("halfwidth"),
                    alpha=self.alpha_adj,
                    target_power=c.get("target_power"),
                    J_max=self.J_max,
                )
                req.append(J)
            required = _plan.combine_J(req, self.J_max)
            if required <= current or current >= self.J_max:
                break
            current = min(required, self.J_max)  # grow toward the closed-form need
            look += 1

        # Final test ONCE on the FULL accumulated set, INCLUDING the pilot rows.
        accumulated = self._run_rows(self.seed_base, current)
        results = []
        rejected = []
        for c in self.contrasts:
            a = accumulated[:, self._metric(c, "A")]
            b = accumulated[:, self._metric(c, "B")]
            a, b, dropped = validate.validate_contrast_pair(a, b, name=c.get("name", "contrast"))
            D = a - b
            mode = c.get("mode", "M1")
            t_res, reject = TwoStageDesign._final_test(D, mode, c.get("margin"), self.alpha_adj)
            results.append({**t_res, "mode": mode, "n_dropped": dropped})
            rejected.append(bool(reject))

        return {
            "J_final": current,
            "n_looks": look + 1,
            "alpha_adj": self.alpha_adj,
            "correction": self.correction,
            "contrast_results": results,
            "rejected": rejected,
            "pilot_reused_in_final_test": True,
            "warning": (
                "Algorithm 2 reuses the pilot rows in the final test; its unconditional "
                "Type I error is NOT nominal and must be measured by the P3-T2 outer-MC "
                "study before any decision is based on it."
            ),
            "design": "D2",
        }

    def _run_rows(self, base, n):
        rows = []
        for j in range(n):
            arr = np.asarray(self.sim_func(base + j), dtype=float)
            if arr.ndim == 0:
                arr = arr.reshape(1, 1)
            elif arr.ndim == 1:
                arr = arr.reshape(1, -1)
            rows.append(arr)
        return np.concatenate(rows, axis=0) if rows else np.empty((0, self.n_metrics or 0))

    def _metric(self, row, key):
        return row[key]


# Functional conveniences --------------------------------------------------- #

def two_stage(sim_func, contrasts, *, J0, n_metrics, **kwargs):
    """Conv predicate wrapper around :class:`TwoStageDesign`."""
    return TwoStageDesign(sim_func, contrasts, J0=J0, n_metrics=n_metrics, **kwargs).run()


def adaptive(sim_func, contrasts, *, J0, n_metrics, **kwargs):
    """Conv predicate wrapper around :class:`AdaptiveDesign`."""
    return AdaptiveDesign(sim_func, contrasts, J0=J0, n_metrics=n_metrics, **kwargs).run()


def plan_and_run(design, sim_func, contrasts, *, J0, n_metrics, **kwargs):
    """Dispatch to the two-stage or the adaptive procedure by ``design`` string."""
    d = (design or "").lower()
    if d in ("d4", "two_stage", "two-stage", "twostage"):
        return TwoStageDesign(sim_func, contrasts, J0=J0, n_metrics=n_metrics, **kwargs).run()
    if d in ("d2", "adaptive", "internal_pilot"):
        return AdaptiveDesign(sim_func, contrasts, J0=J0, n_metrics=n_metrics, **kwargs).run()
    raise validate.ValidationError(f"plan_and_run: unknown design {design!r}.")


def run_confirmatory(sim_func, base_seed, J, contrasts, n_metrics, alpha, *, correction="bonferroni"):
    """Run a fixed-``J`` confirmatory stage only (used by outermc D1/D5/D6)."""
    ds = TwoStageDesign(sim_func, contrasts, J0=1, n_metrics=n_metrics, alpha=alpha,
                        correction=correction)
    ds.J0 = 1
    rows = ds._run_rows(base_seed, J)
    return rows
