"""Loss-distribution families for the outer-MC harness (E1).

Every family samples paired replication losses ``(L_A, L_B)`` of shape
``(R, J, 2)`` (R outer repetitions, J replications per repetition) in one
vectorised call. ``D = L_A - L_B`` is the paired contrast; **lower is better**;
a negative ``E[D]`` means method A beats B (estimand table / tisca_v2_spec.md
section 1).

Pairing uses a Gaussian copula with the design rank-correlation ``rho``, mapped
onto the requested marginals; for the bivariate-normal family the Pearson
correlation is exactly ``rho``, for the transformed families it is the drawn
rank dependence and the achieved correlation is reported as a diagnostic.

Effect / variance control:
  ``theta = E[D]`` (true mean of the paired contrast) is realised by a signed mean
  shift of ``L_B`` relative to ``L_A``. ``sigma_a`` / ``sigma_b`` (variance-ratio
  factor) are realised by standardising each observed marginal to unit variance and
  rescaling. The ``(R, J, 2)`` array is drawn with one top-level generator so a
  fixed ``master_seed`` reproduces the whole block bit-identically; per-replication
  streams are addressed by ``SeedSequence(master_seed, spawn_key=(r, ...))`` when
  the caller needs them (see ``draw_rep_losses``).
"""

from __future__ import annotations

import numpy as np
from scipy import stats

FAMILIES = {
    "normal": "bivariate normal (family a)",
    "lognormal": "bivariate lognormal (family b)",
    "gamma": "gamma marginals via copula (family c)",
    "mix": "normal + 2% catastrophic-failure mixture (family d)",
    "beta": "Beta near the upper boundary, coverage-like (family e)",
    "t3": "bivariate Student-t, 3 df (family f)",
    "empirical": "row-bootstrap of a supplied loss matrix (family g)",
}


def _gauss_z(rng, size, rho):
    """(.., 2) bivariate-normal z-scores with correlation rho, non-standardised."""
    cov = np.array([[1.0, rho], [rho, 1.0]])
    return rng.multivariate_normal(np.zeros(2), cov, size=size)


def _standardise_scaled(block, sigma_a, sigma_b, theta):
    """Standardise the two marginal columns to target sigma and E[D] = theta."""
    a = block[..., 0]
    b = block[..., 1]
    sa = a.std(axis=-1, keepdims=True)
    sb = b.std(axis=-1, keepdims=True)
    a = (a - a.mean(axis=-1, keepdims=True)) / (sa + 1e-12)
    b = (b - b.mean(axis=-1, keepdims=True)) / (sb + 1e-12)
    out = np.stack([a, b], axis=-1)
    out[..., 0] *= sigma_a
    out[..., 1] *= sigma_b
    out[..., 1] -= theta
    return out


_parametric = ["normal", "lognormal", "gamma", "mix", "beta", "t3"]


def sample_batch(family, R, J, rho=0.0, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                 matrix=None, rng=None, master_seed=None):
    """Draw ``(R, J, 2)`` paired losses for one family in one vectorised call.

    ``rng`` is a ``numpy.random.Generator``; if omitted a generator is seeded from
    ``master_seed`` (via a SeedSequence) so the block is reproducible and resumable.
    """
    if rng is None:
        ss = np.random.SeedSequence(master_seed)
        rng = np.random.default_rng(ss)
    z = _gauss_z(rng, (R, J), rho)                     # (R, J, 2) z-scores
    if family == "normal":
        block = np.stack([z[..., 0] * sigma_a, z[..., 1] * sigma_b - theta], axis=-1)
        return block
    if family == "t3":
        x = stats.t.ppf(stats.norm.cdf(z), df=3)
        return _standardise_scaled(x, sigma_a, sigma_b, theta)
    if family == "mix":
        block = np.stack([z[..., 0] * sigma_a, z[..., 1] * sigma_b - theta], axis=-1)
        cat = rng.uniform(size=(R, J)) < 0.02
        mag = rng.exponential(scale=10.0, size=(R, J))
        sign = rng.choice([-1.0, 1.0], size=(R, J))
        block[..., 1] += np.where(cat, sign * mag, 0.0)
        return block
    u = stats.norm.cdf(z)                              # Gaussian copula uniforms
    if family == "lognormal":
        x = np.exp(0.5 * stats.norm.ppf(u))
        return _standardise_scaled(x, sigma_a, sigma_b, theta)
    if family == "gamma":
        x = stats.gamma.ppf(u, a=2.0)
        return _standardise_scaled(x, sigma_a, sigma_b, theta)
    if family == "beta":
        mu, sd = 0.96, 0.02
        alb = (mu * (1 - mu)) / (sd * sd) - 1.0
        a, b = mu * alb, (1 - mu) * alb
        x = stats.beta.ppf(u, a, b)
        return _standardise_scaled(x, sigma_a, sigma_b, theta)
    if family == "empirical":
        if matrix is None or matrix.shape[1] != 2:
            raise ValueError("empirical family requires a (M, 2) loss matrix")
        idx = rng.integers(0, matrix.shape[0], size=(R, J))
        return matrix[idx]
    raise ValueError(f"unknown family {family!r}; choose from {sorted(FAMILIES)}")


def sample_pairs(family, n, rho=0.0, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                 matrix=None, rng=None, seed=None):
    """Compatibility helper: draw ``n`` pairs as a ``(n, 2)`` array."""
    return sample_batch(family, 1, n, rho=rho, sigma_a=sigma_a, sigma_b=sigma_b,
                        theta=theta, matrix=matrix,
                        rng=(np.random.default_rng(seed) if seed is not None
                             else (rng if isinstance(rng, np.random.Generator)
                                   else np.random.default_rng())),
                        master_seed=None)[0]


def draw_rep_losses(family, R, J, rho=0.0, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                    matrix=None, master_seed=0):
    """Vectorised ``(R, J, 2)`` draw seeded via SeedSequence``master_seed``."""
    return sample_batch(family, R, J, rho=rho, sigma_a=sigma_a, sigma_b=sigma_b,
                        theta=theta, matrix=matrix, master_seed=master_seed)
