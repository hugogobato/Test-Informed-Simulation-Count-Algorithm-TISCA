"""Input validation and rejection rules (IJDA #11e) — `validate.py`.

Implements the eight rejection rules enumerated in `docs/tisca_v2_spec.md` §8
plus the estimand-level structural checks in `docs/estimand_table.md` §5.
The package must *refuse* metric structures it cannot validly test instead of
silently producing a number.

Every rejector raises `ValidationError` (a `ValueError` subclass) with an
explicit message.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "ValidationError",
    "require_per_replication",
    "reject_aggregate_column",
    "reject_per_replication_rmse",
    "reject_bare_coverage_as_loss",
    "validate_mode_targets",
    "require_matching_mode",
    "validate_loss_array",
    "validate_contrast_pair",
    "validate_binary",
    "validate_bounded",
    "validate_pilot_samples",
    "validate_margin_feasibility",
]

# Machine-scale floor used to decide whether a contrast variance is (numerically) zero.
_EPS = np.finfo(float).eps


class ValidationError(ValueError):
    """Raised when a metric structure or configuration is invalid for TISCA v2."""


_META_SUFFIX_PATTERNS = (
    # Names a caller might use for an across-replication aggregate.
    "mean",
    ".mean",
    "_mean",
    "avg",
    ".avg",
    "_avg",
    "aggregate",
)


def require_per_replication(values, name: str, *, inference_method: str) -> np.ndarray:
    """Return `values` as a 1-D array after checking it is per-replication.

    Rejects rule 8.1 — an across-replication aggregate passed where a per
    replication ``L_j`` is required. A column whose length is 1 cannot carry
    ``J`` replications and is rejected outright.
    """
    arr = np.asarray(values, dtype=float)
    if arr.ndim > 1:
        raise ValidationError(
            f"Column '{name}' has shape {arr.shape}; a per-replication column must "
            f"be 1-D with one value per replication."
        )
    if arr.size < 2:
        raise ValidationError(
            f"Column '{name}' has {arr.size} value(s). A per-replication column must "
            f"carry at least 2 replications to estimate a variance; "
            f"the across-replication estimand requires replication-level L_j, "
            f"not a single aggregate."
        )
    inference_method = (inference_method or "").lower()
    if inference_method in {"paired-t", "paired-boot", "mcnemar"}:
        # For these methods we need the individual L_j; nothing more to check here.
        pass
    return arr


def reject_aggregate_column(values, name: str) -> None:
    """Refuse a column that is an across-replication aggregate.

    Rule 8.1. If the caller passes a 0-D / length-1 array or a name that plainly
    encodes an aggregation (``.mean``, ``_avg``, ``aggregate``), there is no
    per-replication ``L_j`` to contrast: reject.
    """
    low_name = (name or "").lower()
    if any(p in low_name for p in _META_SUFFIX_PATTERNS) and "per" not in low_name:
        raise ValidationError(
            f"Column '{name}' is named like an across-replication aggregate, so it "
            f"cannot be treated as a per-replication value L_j. Provide the "
            f"replication-level column (one value per replication) instead."
        )
    arr = np.asarray(values)
    if arr.ndim == 0 or arr.size <= 1:
        raise ValidationError(
            f"Column '{name}' carries {arr.size} value(s): a single across-replication "
            f"aggregate cannot be contrasted. Require the per-replication column."
        )


def reject_per_replication_rmse(name: str) -> None:
    """Refuse a per-replication RMSE name (rule 8.2).

    A single squared error ``(tau_hat_j - tau_j)^2`` is a valid per-replication
    object (estimand row 3a); "the RMSE of one replication" is not. RMSE is
    defined only AFTER averaging the squared errors across replications
    (row 3b). If a caller points a contrast at a column called like a
    re-rooted per-replication RMSE, reject and route to the squared-error row.
    """
    low = (name or "").lower()
    if "rmse" in low or "root_mean_sq" in low or low == "rmse_ate":
        raise ValidationError(
            f"Column '{name}' looks like a per-replication RMSE, which is not a "
            f"well-defined scalar estimand (estimand row 3b). Contrast the per-replication "
            f"squared error Q_ATE,j = (tau_hat_j - tau_j)^2 (row 3a) and re-root the "
            f"averaged value, or use the paired bootstrap on the squared scale."
        )


def reject_bare_coverage_as_loss(name: str, *, nominal: float | None) -> None:
    """Refuse raw coverage tested with the lower-is-better machinery (rule 8.3).

    Coverage has a nominal target; ``lower is better`` does not apply (IJDA #10,
    estimand rows 4, 5, 8). The caller must route coverage through the interval
    score (row 9) or, for calibration, through the deviation-from-nominal
    equivalence test (row 8).
    """
    low = (name or "").lower()
    if nominal is not None or any(k in low for k in ("cov", "coverage", "cover")):
        raise ValidationError(
            f"Column '{name}' is a coverage/nominal-target metric and must NOT be tested "
            f"with the lower-is-better machinery. Use the calibration deviation "
            f"('deviation from nominal' equivalence row 8) or the interval score "
            f"(row 9 / CRPS row 10) instead of raw coverage."
        )


def validate_mode_targets(mode: str, target: str) -> None:
    """Check that the mode string and the target string are known (rule 8.4)."""
    modes = {"M1", "M2", "M3", "M4", "M5", None}
    targets = {"power", "precision_mcse", "precision_halfwidth", None}
    if mode not in modes:
        raise ValidationError(
            f"Unknown hypothesis mode '{mode}'. Choose one of M1 (two-sided), "
            f"M2 (directional superiority), M3 (minimum-effect), M4 (non-inferiority), "
            f"M5 (equivalence TOST)."
        )
    if target not in targets:
        raise ValidationError(
            f"Unknown design target '{target}'. Choose one of 'power', "
            f"'precision_mcse', 'precision_halfwidth'."
        )


def require_matching_mode(plan_mode: str, test_mode: str, contrast_name: str) -> None:
    """Enforce C4: the mode used to plan must equal the mode used to test (rule 8.4)."""
    if plan_mode != test_mode:
        raise ValidationError(
            f"Contrast '{contrast_name}': planning mode '{plan_mode}' differs from "
            f"test mode '{test_mode}'. TISCA v2 requires the same sidedness/level and "
            f"hypothesis mode for planning and the final test (change C4)."
        )


def validate_loss_array(loss, *, name: str | None = None, positive: bool = False) -> np.ndarray:
    """Validate a scalar-loss matrix/array used for MCS/SPA (2-D: T x m, lower better)."""
    name = name or "loss"
    arr = np.asarray(loss, dtype=float)
    if arr.ndim != 2:
        raise ValidationError(
            f"Loss matrix '{name}' must be 2-D (rows = replications, columns = models); "
            f"got shape {arr.shape}."
        )
    if arr.shape[0] < 2 or arr.shape[1] < 2:
        raise ValidationError(
            f"Loss matrix '{name}' has shape {arr.shape}; need at least 2 replications "
            f"(rows) and 2 models (columns) for a model comparison."
        )
    if not np.all(np.isfinite(arr)):
        raise ValidationError(
            f"Loss matrix '{name}' contains NaN or infinite entries; these are not "
            f"allowed in a model comparison."
        )
    if positive and np.any(arr < 0):
        raise ValidationError(
            f"Loss matrix '{name}' contains negative entries but must be non-negative."
        )
    return arr


def validate_contrast_pair(a, b, *, name: str | None = None) -> tuple[np.ndarray, np.ndarray, int]:
    """Listwise-delete a contrast pair over common replications.

    Implements convention 3 of the estimand table: replication ``j`` is dropped
    from the contrast ``(A, B)`` if either ``L_A,j`` or ``L_B,j`` is missing.
    Returns the aligned pair and the number of dropped replications. Per-column
    NA dropping is forbidden; this is the only valid NA policy for a paired
    contrast.
    """
    name = name or "contrast"
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.ndim != 1 or b.ndim != 1:
        raise ValidationError(f"Contrast '{name}' requires 1-D per-replication arrays.")
    if a.size != b.size:
        raise ValidationError(
            f"Contrast '{name}': lengths differ ({a.size} vs {b.size}). A paired "
            f"contrast requires the same replication index set for both methods "
            f"before any NA handling."
        )
    mask = np.isfinite(a) & np.isfinite(b)
    n_dropped = int(a.size - int(mask.sum()))
    return a[mask], b[mask], n_dropped


def validate_binary(values, *, name: str) -> None:
    """Require a 0/1 column (estimand rows 4, 12) before McNemar-type analysis."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValidationError(f"Column '{name}' must be 1-D.")
    uniq = np.unique(arr)
    if not np.all(np.isin(uniq, [0.0, 1.0])):
        raise ValidationError(
            f"Column '{name}' is not binary: found values {uniq}. McNemar inference "
            f"requires paired indicators in {{0, 1}}."
        )


def validate_bounded(values, *, lo: float, hi: float, name: str) -> None:
    """Require a column within ``[lo, hi]`` (bounded rows such as coverage)."""
    arr = np.asarray(values, dtype=float)
    if np.any(arr < lo - _EPS) or np.any(arr > hi + _EPS):
        raise ValidationError(
            f"Column '{name}' violates the bounded range [{lo}, {hi}]: min/max "
            f"= {arr.min():.6g}/{arr.max():.6g}."
        )


def validate_pilot_samples(sd_pilot, J0: int, *, gamma: float) -> float:
    """Return an upper-confidence-bound-inflated standard deviation (§3).

    ``γ`` is the tail probability such that with confidence ``1 - γ`` the true
    ``σ_D`` is at most the returned value. The inflation factor is
    ``sqrt((J0-1) / χ²_{γ, J0-1})``. Requires ``J0 >= 2`` and ``0 < γ < 1``.
    """
    if J0 < 2:
        raise ValidationError(f"Pilot size J0 = {J0} < 2: a variance cannot be estimated.")
    if not 0.0 < gamma < 1.0:
        raise ValidationError(f"gamma = {gamma} must lie in (0, 1).")
    if sd_pilot < 0.0 or not math.isfinite(float(sd_pilot)):
        raise ValidationError(f"sd_pilot = {sd_pilot} is not a finite non-negative value.")
    if sd_pilot == 0.0:
        return 0.0
    from scipy.stats import chi2  # deferred import keeps module import light

    q = chi2.ppf(gamma, df=J0 - 1)
    if not (q > 0):
        raise ValidationError(
            f"chi2({gamma}, {J0 - 1}) quantile is not positive; cannot inflate."
        )
    factor = math.sqrt((J0 - 1) / q)
    return float(sd_pilot) * factor


def validate_margin_feasibility(delta: float, margin: float, *, mode: str) -> None:
    """Enforce the genuine infeasibility condition for margin-based modes.

    For M5 (equivalence) a planning alternative ``|δ| >= Δ`` is genuinely
    infeasible: no ``J`` reaches the target power. This is distinct from the
    small-``J`` zero-power artifact, which the planner must simply step past.
    For M3/M4 a margin below 0 is nonsensical.
    """
    if margin is None:
        return
    if margin <= 0:
        raise ValidationError(
            f"Mode {mode} requires a strictly positive margin delta; got {margin}."
        )
    if mode == "M5" and abs(delta) >= margin - np.finfo(float).eps * 1e6:
        raise ValidationError(
            f"M5 (equivalence) planning is infeasible: |delta| = {abs(delta):.6g} >= "
            f"margin {margin:.6g}. The planning alternative must lie strictly inside "
            f"the equivalence margin for any J to reach the target power."
        )
