"""Canonical E1 cell grid (REVISION_PLAN.md P3-T2).

Single source of truth for the four modules, imported by both the Colab runners
(``notebooks/_generators/build_e1_modules_nb.py``) and the local parallel runner
(``experiments/E1_operating_characteristics/run_e1_grid.py``). The grid used to be
duplicated as a string literal inside the notebook generator; keeping it here is
what stops the two from drifting.

Module sizes::

    A  operating characteristics        864 cells  R = 5000
    B  multiplicity                      660 cells  R = 2000, bootstrap B = 999
    C  tuning sensitivity                243 cells  R = 5000
    D  variance ratio                    216 cells  R = 5000
                                       -----------
                                        1983 cells

A and B differ from the plan's 840 and 600 because the empirical family is now run
in two variants (``families.py``): the row bootstrap, which keeps the real joint
dependence and therefore admits no design ``rho``, and the copula variant, in which
``rho`` is a free factor. The plan's single "row-bootstrap of the real matrix" entry
could not do both, and in v1 it did neither -- it ignored every factor it was given.
"""

from __future__ import annotations

import itertools

ALPHA = 0.05
DELTA = 0.5
JMAX = 1000

PARAMETRIC = ["normal", "lognormal", "gamma", "mix", "beta", "t3"]
RHOS = [-0.3, 0.0, 0.3, 0.6, 0.9]
THETA_MULTS = [0.0, 0.5, 1.0, 2.0]
DESIGNS = ["D1", "D2", "D3", "D4", "D5", "D6"]

_SEED_BASE = {"A": 100_000, "B": 400_000, "C": 200_000, "D": 300_000}


def _cell(module, index, factors, **config):
    cfg = {
        "design": "D4", "family": "normal", "rho": 0.0, "sigma_a": 1.0,
        "sigma_b": 1.0, "theta": 0.0, "sigma_D": None, "R": 5000, "J0": 50,
        "Jmax": JMAX, "alpha": ALPHA, "alpha_adj": ALPHA, "mode": 1,
        "delta": DELTA, "power_target": 0.80, "gamma": 0.20,
        "correction": "none", "K": 1, "matrix": None,
        "seed": _SEED_BASE[module] + index, "fixed_J": None, "mcse": 0.05,
        "batch": 50, "B": 50,
    }
    cfg.update(config)
    return {"cell_id": f"{module}_{index:04d}", "module": module,
            "factors": {"cell_index": index, **factors}, "config": cfg}


def module_a(matrix=None):
    """F x rho x theta x Design, K = 1, J0 = 50: the core operating-characteristic grid."""
    cells = []
    for family, rho, tm, design in itertools.product(PARAMETRIC, RHOS, THETA_MULTS, DESIGNS):
        cells.append(_cell("A", len(cells),
                           dict(module_cell="A", family=family, rho=rho, theta_mult=tm,
                                design=design, J0=50, B=50, sigma_ratio=1.0),
                           design=design, family=family, rho=rho, theta=DELTA * tm))
    # family (g'): real marginals, design-controlled dependence
    for rho, tm, design in itertools.product(RHOS, THETA_MULTS, DESIGNS):
        cells.append(_cell("A", len(cells),
                           dict(module_cell="A", family="empirical_copula", rho=rho,
                                theta_mult=tm, design=design, J0=50, B=50, sigma_ratio=1.0),
                           design=design, family="empirical_copula", rho=rho,
                           theta=DELTA * tm, matrix=matrix))
    # family (g): row bootstrap. rho is the data's own, so it is not a factor.
    for tm, design in itertools.product(THETA_MULTS, DESIGNS):
        cells.append(_cell("A", len(cells),
                           dict(module_cell="A", family="empirical", rho=None,
                                theta_mult=tm, design=design, J0=50, B=50, sigma_ratio=1.0),
                           design=design, family="empirical", rho=None,
                           theta=DELTA * tm, matrix=matrix))
    assert len(cells) == 864, len(cells)
    return cells


def module_b(matrix=None, planning_alpha=None):
    """K x correction x rho x F x theta x Design at R = 2000, bootstrap B = 999."""
    cells = []

    def _add(K, correction, rho, family, theta, design):
        alpha_plan, note = (planning_alpha(correction, K, alpha=ALPHA, r=1)
                            if planning_alpha else (ALPHA / max(1, K), "bonferroni-like"))
        cells.append(_cell("B", len(cells),
                           dict(module_cell="B", family=family, rho=rho,
                                theta_mult=round(theta / DELTA, 3), design=design,
                                J0=50, B=999, K=K, correction=correction,
                                alpha_plan=alpha_plan, alpha_note=note),
                           design=design, family=family, rho=rho, theta=theta,
                           R=2000, alpha_adj=alpha_plan, correction=correction, K=K,
                           matrix=matrix if family.startswith("empirical") else None,
                           batch=50, B=999))

    for K, corr, rho, family, theta, design in itertools.product(
            [1, 3, 6], ["none", "bonferroni", "holm", "bh", "romano_wolf"], RHOS,
            ["normal", "empirical_copula"], [0.0, DELTA], ["D3", "D4"]):
        _add(K, corr, rho, family, theta, design)
    for K, corr, theta, design in itertools.product(
            [1, 3, 6], ["none", "bonferroni", "holm", "bh", "romano_wolf"],
            [0.0, DELTA], ["D3", "D4"]):
        _add(K, corr, None, "empirical", theta, design)
    assert len(cells) == 660, len(cells)
    return cells


def module_c(matrix=None):
    """J0 x B x F x rho x Design at the planning alternative: tuning sensitivity."""
    cells = []
    for J0, B, family, rho, design in itertools.product(
            [25, 50, 100], [25, 50, 100], ["normal", "lognormal", "gamma"],
            [-0.3, 0.3, 0.9], ["D2", "D3", "D4"]):
        cells.append(_cell("C", len(cells),
                           dict(module_cell="C", family=family, rho=rho, theta_mult=1.0,
                                design=design, J0=J0, B=B, sigma_ratio=1.0),
                           design=design, family=family, rho=rho, theta=DELTA,
                           J0=J0, batch=B, B=B))
    assert len(cells) == 243, len(cells)
    return cells


def module_d(matrix=None):
    """sigma-ratio x F x rho x theta x Design: the variance-ratio module."""
    cells = []
    for sr, family, rho, tm, design in itertools.product(
            [1.0, 2.0], ["normal", "lognormal", "gamma"], [-0.3, 0.3, 0.9],
            [0.0, 1.0], DESIGNS):
        cells.append(_cell("D", len(cells),
                           dict(module_cell="D", family=family, rho=rho, theta_mult=tm,
                                design=design, J0=50, B=50, sigma_ratio=sr),
                           design=design, family=family, rho=rho, theta=DELTA * tm,
                           sigma_a=sr))
    assert len(cells) == 216, len(cells)
    return cells


def make_grid(modules="ABCD", matrix=None, planning_alpha=None):
    """Concatenate the requested modules in A, B, C, D order."""
    out = []
    for m in "ABCD":
        if m not in modules:
            continue
        if m == "B":
            out += module_b(matrix=matrix, planning_alpha=planning_alpha)
        else:
            out += {"A": module_a, "C": module_c, "D": module_d}[m](matrix=matrix)
    return out


EXPECTED = {"A": 864, "B": 660, "C": 243, "D": 216}
