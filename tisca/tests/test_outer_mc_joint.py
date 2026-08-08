"""Validation tests for the genuine joint Module-B path (G1)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from tisca import multiplicity
from tisca.outermc import e1_grid, engine
from tisca.outermc.joint_engine import run_joint_cell, run_joint_results
from tisca.outermc.joint_families import (
    contrast_map,
    contrasts_from_losses,
    covariance_record,
    covariance_spec,
    draw_joint_losses,
)


def _empirical_matrix():
    rng = np.random.default_rng(5)
    z = rng.multivariate_normal([0, 0], [[1, 0.62], [0.62, 1]], size=500)
    return np.column_stack([9 + np.exp(0.45 * z[:, 0]),
                            11 + 1.4 * np.exp(0.35 * z[:, 1])])


def _fixed(**updates):
    config = {
        "family": "normal", "rho": 0.9, "K": 3, "theta": 0.0,
        "R": 600, "J0": 20, "Jmax": 80, "fixed_J": 80,
        "design": "D4", "correction": "none", "alpha": 0.05,
        "alpha_adj": 0.05, "delta": 0.5, "power_target": 0.80,
        "gamma": 0.20, "seed": 991, "B": 99,
    }
    config.update(updates)
    return config


def test_joint_k1_reduces_to_scalar_module_b_within_monte_carlo_error():
    """K=1 keeps the scalar estimand, planner and test."""
    cfg = dict(engine.DEFAULT_CONFIG)
    cfg.update(
        family="normal", rho=0.3, K=1, theta=0.5, R=2500, J0=50,
        Jmax=400, design="D4", correction="none", alpha=0.05,
        alpha_adj=0.05, delta=0.5, power_target=0.80, gamma=0.20, seed=123,
    )
    scalar = engine.run_e1(cfg)[0]
    joint = run_joint_cell(cfg)[0]
    # The addressable joint RNG is intentionally different from the old monolithic
    # block RNG, so equality is statistical rather than bitwise.
    combined_mcse = np.sqrt(
        scalar["reject_rate"] * (1 - scalar["reject_rate"]) / cfg["R"]
        + joint["marginal_level_or_power"] * (1 - joint["marginal_level_or_power"]) / cfg["R"]
    )
    assert abs(scalar["reject_rate"] - joint["marginal_level_or_power"]) <= 4 * combined_mcse
    assert abs(scalar["E_J"] - joint["E_J"]) < 2.0
    assert abs(joint["bias"]) < 0.03


def test_k_plus_one_block_maps_exactly_to_declared_contrasts():
    losses = draw_joint_losses("normal", 4, 30, 6, rho=0.3, theta=0.5, seed=7)
    assert losses.shape == (4, 30, 7)
    direct = contrasts_from_losses(losses)
    mapped = np.einsum("kj,rtj->rtk", contrast_map(6), losses)
    np.testing.assert_array_equal(direct, mapped)
    for k in range(6):
        np.testing.assert_array_equal(direct[..., k], losses[..., 0] - losses[..., k + 1])


def test_implied_sigma_d_is_psd_for_every_joint_grid_configuration():
    matrix = _empirical_matrix()
    grid = e1_grid.module_b(matrix=matrix, planning_alpha=multiplicity.planning_alpha)
    seen = set()
    for cell in grid:
        cfg = cell["config"]
        key = (cfg["family"], cfg["K"], cfg["rho"])
        if key in seen:
            continue
        seen.add(key)
        spec = covariance_spec(
            cfg["family"], cfg["K"], rho=cfg["rho"], matrix=matrix,
            numeric_n=12_000,
        )
        assert spec.min_eigenvalue >= -1e-8, key
        assert spec.loss_min_eigenvalue >= -1e-8, key
        np.testing.assert_allclose(spec.contrast_covariance,
                                   spec.contrast_covariance.T, atol=1e-12)
    assert len(seen) == 33


def test_unadjusted_fwer_increases_with_family_size_under_global_null():
    estimates = [run_joint_cell(_fixed(K=K, R=1800))[0]["fwer"] for K in (1, 3, 6)]
    assert estimates[0] < estimates[1] < estimates[2], estimates
    assert estimates[2] > 3 * estimates[0], estimates


@pytest.mark.parametrize("correction", ["bonferroni", "holm", "romano_wolf"])
def test_joint_fwer_control_near_alpha(correction):
    alpha_adj = 0.05 if correction == "romano_wolf" else 0.05 / 6
    R = 600
    result = run_joint_cell(_fixed(
        K=6, R=R, correction=correction, alpha_adj=alpha_adj,
        B=149, seed=992,
    ))[0]
    mcse = np.sqrt(0.05 * 0.95 / R)
    assert result["fwer"] <= 0.05 + 4 * mcse, (correction, result["fwer"])
    assert result["fwer"] >= max(0.0, 0.05 - 4 * mcse), (correction, result["fwer"])


def test_bh_reports_fdr_as_a_distinct_field():
    null = run_joint_cell(_fixed(K=6, R=900, correction="bh", alpha_adj=0.05 / 6))[0]
    assert np.isfinite(null["fdr"])
    assert null["fdr"] == pytest.approx(null["fwer"])
    assert null["fdr_mcse"] == pytest.approx(null["fwer_mcse"])
    non_bh = run_joint_cell(_fixed(K=6, R=100, correction="holm", alpha_adj=0.05 / 6))[0]
    assert np.isnan(non_bh["fdr"])


def test_cross_contrast_dependence_changes_joint_not_marginal_estimand():
    base = _fixed(K=3, R=2200, rho=0.3, theta=0.25, fixed_J=40,
                  Jmax=80, seed=993)
    weak = run_joint_cell({**base, "benchmark_residual_corr": 0.0})[0]
    strong = run_joint_cell({**base, "benchmark_residual_corr": 0.8})[0]
    mcse = np.sqrt(0.25 * 0.75 / base["R"])
    assert abs(weak["marginal_level_or_power"] - strong["marginal_level_or_power"]) < 4 * mcse
    assert strong["conjunctive_power"] > weak["conjunctive_power"] + 0.05
    assert strong["disjunctive_power"] < weak["disjunctive_power"] - 0.05


def test_romano_wolf_joint_resampling_preserves_dependence():
    rng = np.random.default_rng(1)
    cov = np.full((3, 3), 0.75)
    np.fill_diagonal(cov, 1.0)
    D = rng.multivariate_normal(np.zeros(3), cov, size=300)
    rw = multiplicity.romano_wolf_stepdown(D, B=799, seed=2)
    observed = rw["observed_contrast_correlation"]
    bootstrap = rw["bootstrap_stat_correlation"]
    assert np.mean(observed[np.triu_indices(3, 1)]) > 0.65
    assert np.mean(bootstrap[np.triu_indices(3, 1)]) > 0.60
    np.testing.assert_allclose(bootstrap, observed, atol=0.12)


def test_fixed_seed_and_outer_shard_offset_are_bit_reproducible():
    cfg = _fixed(K=3, R=64, seed=994)
    first = run_joint_results(cfg)
    second = run_joint_results(cfg)
    for name in first.__dataclass_fields__:
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))

    offset = run_joint_results(cfg, outer_start=17, outer_stop=41)
    keep = (first.outer_index >= 17) & (first.outer_index < 41)
    for name in first.__dataclass_fields__:
        np.testing.assert_array_equal(getattr(first, name)[keep], getattr(offset, name))


def test_k_greater_than_one_never_repeats_one_benchmark_column():
    for family, rho, matrix in [
        ("normal", 0.6, None),
        ("empirical_copula", 0.6, _empirical_matrix()),
        ("empirical", None, _empirical_matrix()),
    ]:
        losses = draw_joint_losses(
            family, 2, 200, 6, rho=rho, matrix=matrix, seed=44,
        )
        assert all(
            not np.array_equal(losses[..., left], losses[..., right])
            for left in range(1, 7)
            for right in range(left + 1, 7)
        )
        assert np.max(np.abs(np.corrcoef(
            contrasts_from_losses(losses).reshape(-1, 6), rowvar=False
        )[np.triu_indices(6, 1)])) < 0.999


def test_covariance_audit_record_contains_required_matrices():
    cfg = _fixed(K=3, R=20)
    record = covariance_record(cfg, cell_id="B_0000", diagnostic_n=800, numeric_n=2000)
    assert record["exact_contrast_mapping"]
    assert record["shared_control_column"]
    assert not record["repeated_benchmark_columns"]
    for key in ("declared_loss_covariance", "implied_sigma_D", "implied_corr_D",
                "achieved_empirical_corr_D"):
        value = np.asarray(json.loads(record[key]))
        expected = 4 if key == "declared_loss_covariance" else 3
        assert value.shape == (expected, expected)
