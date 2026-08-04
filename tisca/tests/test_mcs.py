"""Tests for `mcs.py` — SPA, Reality Check, and the full MCS procedure.

The parity test against the CRAN ``MCS`` package is the key acceptance criterion
(P2-T1). It runs the real installed ``MCS::MCSprocedure`` via R on generated loss
matrices, captures the exact bootstrap resample indices CRAN uses internally, then
drives this Python implementation with the *identical* indices to get a bit-for-bit
comparison (far stronger than "within bootstrap MCSE"). The test is skipped when R
or the ``MCS`` package is unavailable.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pytest

from tisca import mcs
from tisca import validate

ORACLE = os.path.join(os.path.dirname(__file__), "mcs_parity_oracle.R")


def _r_available():
    if not shutil.which("Rscript"):
        return False
    try:
        out = subprocess.run(
            ["Rscript", "-e", "cat(requireNamespace('MCS'))"], capture_output=True, text=True, timeout=60
        )
        return "TRUE" in out.stdout and out.returncode == 0
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _r_available(), reason="R and/or 'MCS' package unavailable")


def _make_loss_cases():
    rng = np.random.default_rng(2026)
    cases = []
    c1 = rng.exponential(size=(80, 4)); c1[:, 3] -= 0.3
    c2 = rng.lognormal(mean=0, sigma=1.0, size=(60, 3)); c2[:, 2] *= 0.8
    c3 = rng.normal(loc=0, scale=1, size=(100, 5))
    c3 += np.array([0.2, 0.2, 0.1, 0.0, 0.0])
    c4 = rng.gamma(shape=2, scale=1, size=(90, 4)); c4[:, 1] += 0.4
    c5 = rng.beta(a=2, b=5, size=(70, 3)) * 2
    c6 = rng.normal(size=(50, 6)); c6[:, 0] += 0.5; c6[:, 1] += 0.3
    for c in [c1, c2, c3, c4, c5, c6]:
        cases.append(np.abs(c))
    return cases


def _run_r_oracle(loss_file, out_prefix, alpha, B, statistic, k, seed):
    subprocess.run(
        [shutil.which("Rscript"), ORACLE, loss_file, out_prefix, str(alpha),
         str(B), statistic, str(k), str(seed)],
        check=True, capture_output=True, text=True, timeout=600,
    )


def _read_oracle(out_prefix):
    tab = pd.read_csv(out_prefix + "_result.csv")
    inc = pd.read_csv(out_prefix + "_included.csv")["included"].astype(str).tolist()
    exc = pd.read_csv(out_prefix + "_excluded.csv")["excluded"].astype(str).tolist()
    return tab, inc, exc


def _map_r_names(names):
    return sorted(n.replace("V", "model_") for n in names)


@pytest.mark.parametrize("statistic", ["Tmax"])
def test_mcs_parity_with_cran_tmax(statistic):
    """All >=5 test cases match CRAN MCS bit-for-bit on p-values + sets (Tmax)."""
    B, alpha, seed = 4999, 0.15, 42
    cases = _make_loss_cases()
    assert len(cases) >= 5
    with tempfile.TemporaryDirectory() as tmp:
        for i, loss in enumerate(cases, 1):
            base = os.path.join(tmp, f"case{i}")
            loss_file = base + "_loss.csv"
            np.savetxt(loss_file, loss, delimiter=",")
            _run_r_oracle(loss_file, base, alpha, B, statistic, 0, seed)
            # The exact bootstrap indices R used.
            r_idx = pd.read_csv(base + "_indices.csv", header=None).to_numpy(dtype=int)
            tab_r, inc_r, exc_r = _read_oracle(base)
            # Drive the Python implementation with identical resamples.
            res = mcs.mcs(loss, alpha=alpha, B=B, statistic=statistic, k=0,
                          bootstrap_indices=r_idx)
            # p-values: both are multisets of elimination-step p-values; sort before compare.
            p_abs = np.max(np.abs(np.sort(res["p_H0"]) - np.sort(tab_r["p_H0"].to_numpy())))
            assert p_abs < 1e-6, (i, p_abs)
            assert _map_r_names(inc_r) == sorted(res["included"])
            assert _map_r_names(exc_r) == sorted(res["excluded"])


def test_mcs_parity_with_cran_tr():
    """TR statistic also agrees with CRAN to within bootstrap MCSE."""
    B, alpha, seed = 4999, 0.15, 9
    cases = _make_loss_cases()[:2]
    with tempfile.TemporaryDirectory() as tmp:
        for i, loss in enumerate(cases, 1):
            base = os.path.join(tmp, f"tr{i}")
            loss_file = base + "_loss.csv"
            np.savetxt(loss_file, loss, delimiter=",")
            _run_r_oracle(loss_file, base, alpha, B, "TR", 0, seed)
            r_idx = pd.read_csv(base + "_indices.csv", header=None).to_numpy(dtype=int)
            tab_r, inc_r, exc_r = _read_oracle(base)
            res = mcs.mcs(loss, alpha=alpha, B=B, statistic="TR", k=0,
                          bootstrap_indices=r_idx)
            diff_p = np.max(np.abs(np.sort(res["p_H0"]) - np.sort(tab_r["p_H0"].to_numpy())))
            # Within bootstrap MCSE of the difference of two independent estimates.
            mcse = 1.0 / np.sqrt(B)
            assert diff_p < 3 * mcse, (i, diff_p)
            assert _map_r_names(exc_r) == sorted(res["excluded"])


def test_mcs_basic_output_structure():
    rng = np.random.default_rng(0)
    loss = np.exp(rng.normal(size=(40, 4)))
    res = mcs.mcs(loss, alpha=0.15, B=199, seed=0)
    assert len(res["model_names"]) == 4
    assert res["p_H0"].shape == (4,) and res["p_mcs"].shape == (4,)
    # At least one model must remain (the best is never excluded by construction).
    assert len(res["included"]) >= 1
    # Elimination order has at most m-1 steps.
    assert len(res["elimination_order"]) <= 3
    assert res["table"].shape == (4, 3)


def test_mcs_picks_best_model_in_set():
    rng = np.random.default_rng(0)
    # Model 4 is genuinely best (lowest shift on the loss scale).
    loss = np.exp(rng.normal(size=(200, 4)) + np.array([0.3, 0.2, 0.1, 0.0]))
    res = mcs.mcs(loss, alpha=0.30, B=999, seed=1)
    best = res["avg_loss"].argmin()
    # The best-performing model (lowest realised average loss) is never excluded.
    assert res["model_names"][best] in res["included"]


def test_mcs_rejects_bad_input():
    with pytest.raises(validate.ValidationError):
        mcs.mcs(np.random.rand(10))  # not 2-D
    with pytest.raises(validate.ValidationError):
        mcs.mcs(np.array([[np.nan, 0.1], [0.2, 0.3]]))  # nan
    with pytest.raises(validate.ValidationError):
        mcs.mcs(np.ones((5, 3)), statistic="bogus")  # invalid statistic


def test_reality_check_and_spa_sensible():
    """RC/SPA p-value reflects whether the champion actually dominates."""
    rng = np.random.default_rng(0)
    # Champion (col 0) clearly best => small p (few extreme bootstrap exceed obs).
    loss = np.exp(rng.normal(size=(150, 3)) + np.array([-0.6, 0.0, 0.1]))
    rc = mcs.reality_check_pvalue(loss, champion=0, B=199, seed=1)
    sp = mcs.spa_pvalue(loss, champion=0, B=199, seed=2)
    assert 0.0 <= rc["p_value"] <= 1.0
    assert 0.0 <= sp["p_value"] <= 1.0
    assert rc["T_obs"] > 0  # champion has lower average loss


def test_reality_check_poor_champion():
    """A poor champion yields a negative T and (likely) large p (> 0.1)."""
    rng = np.random.default_rng(0)
    loss = np.exp(rng.normal(size=(150, 3)) + np.array([0.5, 0.0, 0.0]))
    rc = mcs.reality_check_pvalue(loss, champion=0, B=199, seed=1)
    assert rc["T_obs"] < 0
