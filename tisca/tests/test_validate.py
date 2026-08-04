"""Tests for `validate.py` — the input rejection rules (spec §8, estimand table §5)."""

from __future__ import annotations

import numpy as np
import pytest

from tisca import validate


def test_reject_aggregate_column_by_name():
    with pytest.raises(validate.ValidationError):
        validate.reject_aggregate_column(np.arange(10.0), "coverage_mean")
    with pytest.raises(validate.ValidationError):
        validate.reject_aggregate_column(np.array([1.234]), "pehe_j_mean")
    # A single value is an aggregate regardless of name.
    with pytest.raises(validate.ValidationError):
        validate.reject_aggregate_column(np.array([5.0]), "pehe")


def test_reject_per_replication_rmse():
    with pytest.raises(validate.ValidationError):
        validate.reject_per_replication_rmse("rmse_ate")
    with pytest.raises(validate.ValidationError):
        validate.reject_per_replication_rmse("RMSE_ATE")
    # Squared error names are fine.
    validate.reject_per_replication_rmse("ate_sq_err")


def test_reject_bare_coverage_as_loss():
    with pytest.raises(validate.ValidationError):
        validate.reject_bare_coverage_as_loss("cov_95", nominal=0.95)
    with pytest.raises(validate.ValidationError):
        validate.reject_bare_coverage_as_loss("cate_coverage", nominal=None)
    # Interval score / CRPS are valid lower-is-better losses.
    validate.reject_bare_coverage_as_loss("interval_score", nominal=None)
    validate.reject_bare_coverage_as_loss("crps", nominal=None)


def test_validate_mode_targets():
    validate.validate_mode_targets("M1", "power")
    with pytest.raises(validate.ValidationError):
        validate.validate_mode_targets("M9", "power")
    with pytest.raises(validate.ValidationError):
        validate.validate_mode_targets("M1", "bogus")


def test_require_matching_mode():
    validate.require_matching_mode("M1", "M1", "c1")
    with pytest.raises(validate.ValidationError):
        validate.require_matching_mode("M1", "M2", "c1")


def test_validate_loss_array():
    good = validate.validate_loss_array(np.arange(30.0).reshape(10, 3))
    assert good.shape == (10, 3)
    with pytest.raises(validate.ValidationError):
        validate.validate_loss_array(np.array([1.0, 2.0, 3.0]))
    with pytest.raises(validate.ValidationError):
        validate.validate_loss_array(np.array([[1.0, np.nan], [2.0, 3.0]]))
    with pytest.raises(validate.ValidationError):
        validate.validate_loss_array(np.ones((5, 1)))
    with pytest.raises(validate.ValidationError):
        validate.validate_loss_array(np.array([[-1.0, 2.0], [3.0, 4.0]]), positive=True)  # negative
    validate.validate_loss_array(np.array([[1.0, 2.0], [3.0, 4.0]]), positive=True)


def test_validate_contrast_pair_listwise():
    a, b, dropped = validate.validate_contrast_pair(
        np.array([1.0, np.nan, 3.0]), np.array([0.0, 0.0, 0.0])
    )
    assert dropped == 1 and a.shape == (2,) and b.shape == (2,)
    with pytest.raises(validate.ValidationError):
        validate.validate_contrast_pair(np.array([1.0, 2.0]), np.array([1.0]))


def test_validate_binary_bounded():
    validate.validate_binary(np.array([0, 1, 0, 1]), name="cov")
    with pytest.raises(validate.ValidationError):
        validate.validate_binary(np.array([0, 1, 2]), name="cov")
    validate.validate_bounded(np.array([0.5, 0.9]), lo=0, hi=1, name="p")
    with pytest.raises(validate.ValidationError):
        validate.validate_bounded(np.array([1.2]), lo=0, hi=1, name="p")


def test_validate_pilot_samples():
    assert validate.validate_pilot_samples(1.0, 50, gamma=0.20) == pytest.approx(1.099, abs=1e-3)
    assert validate.validate_pilot_samples(0.0, 50, gamma=0.20) == 0.0
    with pytest.raises(validate.ValidationError):
        validate.validate_pilot_samples(1.0, 1, gamma=0.2)
    with pytest.raises(validate.ValidationError):
        validate.validate_pilot_samples(1.0, 50, gamma=1.5)


def test_validate_margin_feasibility():
    validate.validate_margin_feasibility(0.2, 0.5, mode="M5")  # inside
    with pytest.raises(validate.ValidationError):
        validate.validate_margin_feasibility(0.5, 0.5, mode="M5")  # |delta|>=Delta
    with pytest.raises(validate.ValidationError):
        validate.validate_margin_feasibility(-0.7, 0.5, mode="M5")
    with pytest.raises(validate.ValidationError):
        validate.validate_margin_feasibility(0.2, -0.5, mode="M5")  # nonpositive margin


def test_require_per_replication():
    arr = validate.require_per_replication(np.arange(10.0), "x", inference_method="paired-t")
    assert arr.shape == (10,)
    with pytest.raises(validate.ValidationError):
        validate.require_per_replication(np.array([[1.0, 2.0], [3.0, 4.0]]), "x", inference_method="paired-t")
