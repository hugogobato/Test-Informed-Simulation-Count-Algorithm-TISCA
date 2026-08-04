"""Tests for `procedure.py` — Algorithm 1 (two-stage default) and Algorithm 2
(optional adaptive), including the degenerate D~N(0,0) case (P3-T7)."""

from __future__ import annotations

import numpy as np
import pytest

from tisca import procedure as PROC
from tisca import validate


def _make_sim(delta_true, sigma=1.0, seed=0):
    def f(seed):
        rng = np.random.default_rng(seed)
        D = rng.normal(delta_true, sigma)
        e1, e2 = rng.normal(0, 1, 2)
        return [e1 + D / 2, e2 - D / 2]
    return f


def test_two_stage_design_discards_pilot():
    """Algorithm 1: independent pilot, closed-form J, pilot NOT reused."""
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.3, target_power=0.80, name="A-vs-B")]
    res = PROC.two_stage(_make_sim(-0.3), contrasts, J0=50, n_metrics=2,
                         alpha=0.05, J_max=5000)
    assert res["design"] == "D4"
    assert res["J_final"] > 0
    assert res["capped"] is False
    # pilot sd used; confirmatory estimate close to truth and significant
    est = res["contrast_results"][0]["estimate"]
    assert est < 0 and res["contrast_results"][0]["p_value"] < 0.05
    assert res["rejected"][0]
    # alpha_adj correct for K=1 bonferroni (which is just alpha)
    assert res["alpha_adj"] == pytest.approx(0.05)
    assert "Js_by_contrast" in res and res["Js_by_contrast"][0]["J"] == res["J_final"]


def test_two_stage_precision_only_target():
    contrasts = [dict(A=0, B=1, mode="M1", delta=0.0, target_mcse=0.05, name="A-vs-B")]
    res = PROC.two_stage(_make_sim(0.0), contrasts, J0=50, n_metrics=2,
                         alpha=0.05, J_max=5000)
    assert res["J_final"] >= 2
    # precision target dominates (delta=0 -> power target vacuous if not set)
    assert res["Js_by_contrast"][0]["J"] >= 2


def test_two_stage_family_multiple_contrasts():
    """K=6 family drives alpha_adj down and J up (multiplicity-aware planning)."""
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.5, target_power=0.80, name=f"c{i}")
                 for i in range(6)]
    res = PROC.two_stage(_make_sim(-0.5), contrasts, J0=50, n_metrics=2,
                         alpha=0.05, correction="bonferroni", J_max=5000)
    assert res["alpha_adj"] == pytest.approx(0.05 / 6)
    assert len(res["contrast_results"]) == 6
    assert res["family_rejected_conjunctive"] or res["family_rejected_disjunctive"] is not None


def test_two_stage_capped_at_J_max():
    # A small planning alternative over a large sigma needs far more than 40 reps.
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.05, target_power=0.90, name="A-vs-B")]
    res = PROC.two_stage(_make_sim(-0.05), contrasts, J0=25, n_metrics=2,
                         alpha=0.05, J_max=40)
    assert res["J_final"] <= 40
    assert res["capped"] is True


def test_two_stage_requires_n_metrics():
    with pytest.raises(validate.ValidationError):
        PROC.two_stage(_make_sim(0.0), [dict(A=0, B=1)], J0=10, n_metrics=None)


def test_degenerate_same_model_terminates_gracefully():
    """P3-T7: both models identical => D ~ N(0, 0); must terminate, no infinite loop."""
    def identical(seed):
        rng = np.random.default_rng(seed)
        v = rng.normal()
        return [v, v]

    contrasts = [dict(A=0, B=1, mode="M1", delta=0.5, target_power=0.80, name="same")]
    res = PROC.two_stage(identical, contrasts, J0=50, n_metrics=2,
                         alpha=0.05, J_max=5000)
    assert res["J_final"] >= 2
    assert res["contrast_results"][0]["estimate"] == 0.0
    assert res["contrast_results"][0]["p_value"] == 1.0
    assert res["rejected"] == [False]


def test_adaptive_flags_pilot_reuse():
    """Algorithm 2 reuses the pilot and warns that error rates must be validated."""
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.3, target_power=0.80, name="A-vs-B")]
    res = PROC.adaptive(_make_sim(-0.3), contrasts, J0=30, n_metrics=2,
                        nmax_looks=5, J_max=500)
    assert res["design"] == "D2"
    assert res["pilot_reused_in_final_test"] is True
    assert "Type I error" in res["warning"]


def test_plan_and_run_dispatch():
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.3, target_power=0.80, name="c")]
    a = PROC.plan_and_run("D4", _make_sim(-0.3), contrasts, J0=30, n_metrics=2)
    b = PROC.plan_and_run("D2", _make_sim(-0.3), contrasts, J0=30, n_metrics=2)
    assert a["design"] == "D4" and b["design"] == "D2"
    with pytest.raises(validate.ValidationError):
        PROC.plan_and_run("D9", _make_sim(-0.3), contrasts, J0=30, n_metrics=2)


def test_two_stage_marginal_family_power_present():
    """The marginal power vector and family flags are reported (C8)."""
    contrasts = [dict(A=0, B=1, mode="M1", delta=-0.5, target_power=0.80, name="c")]
    res = PROC.two_stage(_make_sim(-0.5), contrasts, J0=50, n_metrics=2,
                         alpha=0.05, success_criterion="conjunctive")
    assert res["family_rejected_conjunctive"] is True
    assert res["family_rejected_disjunctive"] is None
