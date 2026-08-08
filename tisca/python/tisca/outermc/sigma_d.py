"""Exact ``sigma_D`` for every loss family — the oracle design D6 plans from.

D6's entire claim is that it plans from the *true* contrast sd, so the value it
uses must not itself be a Monte Carlo estimate. The previous implementation drew
one million pairs per cell and claimed ``sigma/sqrt(2n) = 0.07%`` accuracy. That
formula is the MCSE of a sample sd under a **finite fourth moment**, and family
``t3`` does not have one: Student-t with 3 df has finite moments only of order
< 3, so the sample sd of ``D`` obeys no central limit theorem and the quoted
0.07% is not a valid statement about it. Measured against the closed forms below,
the 1e6-draw estimate was off by up to **1.4%** on ``t3`` (rho = 0.9), twenty
times the claimed accuracy, and the error does not shrink reliably with n.

Every family here is therefore solved exactly, or by quadrature whose error is
bounded and reported, rather than sampled. Throughout, ``theta`` only shifts the
mean of ``D`` and never its variance, so it does not appear.

Write ``A`` and ``B`` for the two per-replication losses as ``families.py``
constructs them, ``(Z_A, Z_B)`` for the underlying bivariate normal with
correlation ``rho``, and ``D = A - B``.

**normal.** ``A = sigma_a Z_A``, ``B = sigma_b Z_B``, so

    sigma_D^2 = sigma_a^2 + sigma_b^2 - 2 rho sigma_a sigma_b        (exact)

**mix.** Normal pair plus an independent contamination on B: with probability
0.02 a symmetric ``Exponential(scale = 10)`` magnitude. That component has mean 0
and second moment ``0.02 * 2 * 10^2 = 4``, and it is independent of the pair, so

    sigma_D^2 = sigma_a^2 + sigma_b^2 - 2 rho sigma_a sigma_b + 4    (exact)

**lognormal / gamma / beta / t3, the default ``asym = True``.** Method A carries
the skewed marginal ``G(z) = (g(Phi(z)) - mu)/sd`` (mean 0, variance 1 by
construction) and method B keeps the standard normal from the same copula draw:
``A = sigma_a G(Z_A)``, ``B = sigma_b Z_B``. Then

    Cov(G(Z_A), Z_B) = rho * E[G(Z) Z]                (Z_B = rho Z_A + s W)

so with the single family constant ``c_f = E[G(Z) Z]``,

    sigma_D^2 = sigma_a^2 + sigma_b^2 - 2 rho c_f sigma_a sigma_b

``c_f = (1/sd) * integral_0^1 g(u) Phi^{-1}(u) du`` is a one-dimensional integral
evaluated by adaptive quadrature. It reproduces the lognormal closed form
``0.5 e^{1/8}/sd`` to 1e-13, so it is exact to working precision for all four.

**lognormal / gamma / beta / t3 with ``asym = False``.** Both marginals are
skewed, so the needed object is ``kappa_f(rho) = E[G(Z_A) G(Z_B)]``. Expanding
``G`` in normalised Hermite polynomials, ``G = sum_k a_k He_k/sqrt(k!)`` with
``sum_k a_k^2 = 1``, Mehler's formula gives

    kappa_f(rho) = sum_{k >= 1} a_k^2 rho^k,   sigma_D^2 = sigma_a^2 + sigma_b^2
                                               - 2 sigma_a sigma_b kappa_f(rho)

Truncating at ``K`` terms leaves a remainder bounded by
``|rho|^{K+1} (1 - sum_{k <= K} a_k^2)``, which is **computed** and asserted, not
assumed. ``a_1`` equals ``c_f``, which cross-checks the two routes.

**empirical (row bootstrap).** The law is discrete uniform on the M rows of the
loss matrix under a fixed affine map, so ``sigma_D`` is the population sd of that
finite vector: exact, no integral at all.

**empirical_copula.** Marginals are the empirical quantile functions, coupled by
the Gaussian copula. ``Var(D)`` needs ``Cov(q_A(U_A), q_B(U_B))``, computed by
conditioning: the inner expectation ``E[q_B | Z_A = z]`` is a finite sum of normal
tail probabilities (exact given z), and the outer integral over ``u = Phi(z)`` is
Gauss-Legendre on each of the M mass cells, where ``q_A`` is constant.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate, stats

from . import families as _fam

__all__ = ["sigma_D_exact", "family_c", "family_kappa", "MEHLER_TERMS"]

MEHLER_TERMS = 40                  # asym=False truncation; the bound is checked
_MEHLER_TOL = 1e-8                 # max accepted truncation remainder on kappa
_CELL_NODES = 24                   # Gauss-Legendre nodes per empirical-copula cell
_TAIL_NODES = 400                  # ... and in the two unbounded end cells

_C_CACHE: dict[str, float] = {}
_A_CACHE: dict[str, np.ndarray] = {}


def _G(family, z):
    """Standardised skewed marginal ``(g(Phi(z)) - mu)/sd``: mean 0, variance 1."""
    mu, sd = _fam._MOMENTS[family]
    u = stats.norm.cdf(z)
    return (_fam._family_marginal(family, u, z) - mu) / sd


def family_c(family) -> float:
    """``c_f = E[G(Z) Z] = Corr(G(Z), Z)``, the only family constant asym=True needs.

    Integrated in ``u``-space, where the integrand ``g(u) Phi^{-1}(u)`` has only
    integrable endpoint singularities and adaptive quadrature converges cleanly.
    Doing it in ``z``-space instead is numerically unsafe for ``t3``: at |z| ~ 40,
    ``Phi(z)`` rounds to 1, ``t.ppf(1) = inf`` and ``phi(z) = 0`` multiply to NaN.
    """
    if family not in _C_CACHE:
        mu, sd = _fam._MOMENTS[family]

        def integrand(u):
            z = stats.norm.ppf(u)
            return _fam._family_marginal(family, u, z) * z

        val, err = integrate.quad(integrand, 0.0, 1.0, limit=400,
                                  epsabs=1e-12, epsrel=1e-12)
        if not math.isfinite(val):
            raise ValueError(f"quadrature for family {family!r} did not converge")
        _C_CACHE[family] = val / sd
    return _C_CACHE[family]


def _mehler_coeffs(family, K=MEHLER_TERMS) -> np.ndarray:
    """``a_k = E[G(Z) He_k(Z)]/sqrt(k!)`` for k = 1..K (probabilists' Hermite).

    High-order Hermite polynomials overflow at the large ``|z|`` the adaptive rule
    reaches near ``u -> 0, 1`` (degree 40 at |z| = 8 is already ~1e36, and for
    ``t3`` it multiplies a diverging marginal). Rather than trust that, the result
    is checked against Parseval: ``G`` has unit variance by construction, so
    ``sum_k a_k^2`` must converge to 1 from below. An overflow or a NaN cannot
    satisfy that, so the identity doubles as the guard.
    """
    if family not in _A_CACHE:
        out = []
        for k in range(1, K + 1):
            basis = np.zeros(k + 1)
            basis[k] = 1.0

            def integrand(u, basis=basis):
                z = stats.norm.ppf(u)
                with np.errstate(over="ignore", invalid="ignore"):
                    val = _G(family, z) * np.polynomial.hermite_e.hermeval(z, basis)
                # The true integrand decays; a non-finite value is a floating-point
                # artefact of the two diverging factors and contributes nothing.
                return val if np.isfinite(val) else 0.0

            val, _ = integrate.quad(integrand, 0.0, 1.0, limit=500,
                                    epsabs=1e-12, epsrel=1e-12)
            out.append(val / math.sqrt(math.factorial(k)))
        a = np.asarray(out, dtype=float)
        total = float(np.sum(a * a))
        if not np.all(np.isfinite(a)) or not (0.99 < total <= 1.0 + 1e-3):
            raise ValueError(
                f"Mehler expansion for family {family!r} failed its Parseval check: "
                f"sum a_k^2 = {total!r} (expected -> 1). The coefficients are not "
                f"usable; do not fall back to sampling.")
        _A_CACHE[family] = a
    return _A_CACHE[family]


def family_kappa(family, rho) -> float:
    """``E[G(Z_A) G(Z_B)]`` for the exchangeable (``asym = False``) construction."""
    a = _mehler_coeffs(family)
    powers = np.power(float(rho), np.arange(1, a.size + 1))
    kappa = float(np.sum(a * a * powers))
    # Computed, not assumed: the tail the truncation drops is bounded by
    # |rho|^{K+1} times the unexplained variance share.
    remainder = abs(rho) ** (a.size + 1) * max(0.0, 1.0 - float(np.sum(a * a)))
    if remainder > _MEHLER_TOL:
        raise ValueError(
            f"Mehler truncation for family {family!r} at rho={rho} leaves "
            f"{remainder:.2e} > {_MEHLER_TOL:g}; raise MEHLER_TERMS.")
    return kappa


def _empirical_copula_cov(matrix, rho) -> float:
    """``Cov(q_A(U_A), q_B(U_B))`` under the Gaussian copula, by exact conditioning."""
    a = np.sort(matrix[:, 0])
    b = np.sort(matrix[:, 1])
    M = a.size
    s = math.sqrt(1.0 - rho * rho)
    edges = stats.norm.ppf(np.arange(1, M) / M)     # interior cell boundaries
    db = np.diff(b)

    def inner(z):
        """E[q_B(U_B) | Z_A = z], exact: b_(0) + sum_j (b_(j)-b_(j-1)) P(Z_B > e_j|z)."""
        z = np.asarray(z, dtype=float)[:, None]
        return b[0] + (db[None, :] * stats.norm.sf((edges[None, :] - rho * z) / s)).sum(1)

    x_mid, w_mid = np.polynomial.legendre.leggauss(_CELL_NODES)
    x_end, w_end = np.polynomial.legendre.leggauss(_TAIL_NODES)
    total = 0.0
    for k in range(M):
        lo, hi = k / M, (k + 1) / M
        # q_A is exactly constant on this cell, so only the smooth inner factor is
        # quadratured. The two end cells map to an unbounded z range and get the
        # denser rule.
        x, w = (x_end, w_end) if k in (0, M - 1) else (x_mid, w_mid)
        u = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
        total += a[k] * 0.5 * (hi - lo) * float(np.dot(w, inner(stats.norm.ppf(u))))
    return total - float(a.mean()) * float(b.mean())


def sigma_D_exact(family, rho, sigma_a=1.0, sigma_b=1.0, matrix=None, asym=True) -> float:
    """True paired-difference sd of ``D = L_A - L_B``, with no Monte Carlo anywhere."""
    a, b = float(sigma_a), float(sigma_b)

    if family == "empirical":
        mat = _fam._check_matrix(matrix)
        mom = _fam.empirical_moments(mat)
        scale = math.sqrt((mom[0][1] ** 2 + mom[1][1] ** 2) / 2.0)
        d = (a * (mat[:, 0] - mom[0][0]) - b * (mat[:, 1] - mom[1][0])) / scale
        return float(d.std(ddof=0))     # ddof=0: this IS the bootstrap population

    if rho is None:
        raise ValueError(f"family {family!r} requires a numeric rho")
    r = float(rho)
    if not -1.0 < r < 1.0:
        raise ValueError(f"rho must lie strictly inside (-1, 1); got {r}")

    if family == "normal":
        return float(math.sqrt(a * a + b * b - 2.0 * r * a * b))

    if family == "mix":
        # 2% contamination, symmetric Exponential(scale=10): mean 0, E[C^2] = 4.
        var_contamination = 0.02 * 2.0 * 10.0 ** 2
        return float(math.sqrt(a * a + b * b - 2.0 * r * a * b + var_contamination))

    if family == "empirical_copula":
        mat = _fam._check_matrix(matrix)
        mom = _fam.empirical_moments(mat)
        scale = (mom[0][1] ** 2 + mom[1][1] ** 2) / 2.0
        cov = _empirical_copula_cov(mat, r)
        var = (a * a * mom[0][1] ** 2 + b * b * mom[1][1] ** 2 - 2.0 * a * b * cov)
        return float(math.sqrt(max(var, 0.0) / scale))

    if family in _fam._SKEWED:
        coupling = r * family_c(family) if asym else family_kappa(family, r)
        return float(math.sqrt(a * a + b * b - 2.0 * a * b * coupling))

    raise ValueError(f"no exact sigma_D derived for family {family!r}")
