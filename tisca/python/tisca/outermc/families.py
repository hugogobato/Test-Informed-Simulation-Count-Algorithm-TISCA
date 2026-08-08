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
    "empirical": ("row bootstrap of a supplied loss matrix, real joint dependence "
                  "(family g; rho is fixed by the data, so omit it or pass None)"),
    "empirical_copula": ("empirical marginals of a supplied loss matrix coupled by "
                         "the design's Gaussian copula (family g', rho is free)"),
}


class _Unset:
    """Sentinel for "the caller named no ``rho``", which is NOT the same as ``None``.

    ``rho`` used to default to ``0.0``, so ``sample_batch("empirical", matrix=...)``
    raised on its own default: the row bootstrap refuses an imposed rho (see the
    family note below), and a defaulted ``0.0`` is indistinguishable from a caller
    who explicitly asked for independence. Every call site then had to spell out
    ``rho=None`` or crash. With this sentinel the three states are separate:

      unnamed (``UNSET``)   -> ``0.0`` for the copula families, accepted by the
                               row bootstrap, which has no rho to set;
      explicit ``None``     -> "this family's rho is fixed by its data"; still an
                               error for the copula families, which need a number;
      explicit number       -> honoured, and still an error for ``empirical``.
    """

    __slots__ = ()

    def __repr__(self):                                   # pragma: no cover - debug aid
        return "<rho unset>"


UNSET = _Unset()


def _gauss_z(rng, size, rho):
    """(.., 2) bivariate-normal z-scores with correlation rho, non-standardised."""
    cov = np.array([[1.0, rho], [rho, 1.0]])
    return rng.multivariate_normal(np.zeros(2), cov, size=size)


def _standardise_scaled(block, sigma_a, sigma_b, theta, moments=None):
    """Standardise the two marginal columns to target sigma and ``E[D] = theta``.

    **The centring and scaling constants must be POPULATION constants, not the
    per-repetition sample mean and sd.** Standardising each repetition by its own
    ``mean(axis=-1)`` / ``std(axis=-1)`` forces ``mean(L_A) == mean(L_B) == 0``
    *within every repetition*, hence ``D_bar == theta`` exactly in every
    repetition: the sampling variability of the estimator is removed. The paired
    t statistic then collapses to ``sqrt(J) theta / s_D``, so at ``theta = 0``
    the measured Type I error is exactly 0 and CI coverage exactly 1 -- for the
    lognormal, gamma, beta and t3 families, i.e. four of the seven. Those are
    precisely the families the P3-T4 skewness sub-study reads.

    ``moments`` is the ``(mean, sd)`` pair of the *marginal law*, computed once in
    ``_MOMENTS`` by high-accuracy quadrature/closed form.
    """
    if moments is None:
        raise ValueError("population (mean, sd) required; per-sample standardisation is invalid")
    mu, sd = moments
    a = (block[..., 0] - mu) / sd
    b = (block[..., 1] - mu) / sd
    out = np.stack([a, b], axis=-1)
    out[..., 0] *= sigma_a
    out[..., 1] *= sigma_b
    out[..., 1] -= theta
    return out


def _marginal_moments(family):
    """Population ``(mean, sd)`` of each family's marginal, on the raw transform scale."""
    if family == "lognormal":            # exp(0.5 Z), Z ~ N(0,1) -> lognormal(0, 0.5)
        s2 = 0.25
        m = np.exp(s2 / 2.0)
        v = (np.exp(s2) - 1.0) * np.exp(s2)
        return float(m), float(np.sqrt(v))
    if family == "gamma":                # Gamma(shape=2, scale=1)
        return 2.0, float(np.sqrt(2.0))
    if family == "t3":                   # Student t, 3 df: mean 0, var = 3/(3-2)
        return 0.0, float(np.sqrt(3.0))
    if family == "beta":                 # Beta with mean .96, sd .02 (coverage-like)
        return 0.96, 0.02
    raise ValueError(f"no population moments registered for family {family!r}")


_MOMENTS = {f: _marginal_moments(f) for f in ("lognormal", "gamma", "t3", "beta")}


# --- empirical family (g) -----------------------------------------------------
#
# The v1 implementation was ``matrix[idx]``: a row bootstrap that silently ignored
# ``theta``, ``rho``, ``sigma_a`` and ``sigma_b``. That made every design factor
# inert. Its ``theta = 0`` cells were not nulls at all -- the real MVBCF-minus-BCF
# PEHE matrix carries E[D] = -2.20 -- so the measured "Type I error" was 1.000 for
# D1-D5 and 0.948 for the oracle D6, and the five rho levels were exact duplicates.
#
# The repair splits family (g) in two, because no single construction can both
# preserve the real joint dependence and let the design set rho:
#
#   ``empirical``          row bootstrap. Both losses come from the SAME observed
#                          replication, so the real joint structure survives intact
#                          and so does skew(D) = -1.55, which is the realism this
#                          family contributes to the P3-T4 skewness sub-study. rho
#                          is therefore not a free factor: requesting one raises.
#   ``empirical_copula``   the real marginals (empirical quantile functions), coupled
#                          by the same Gaussian copula the other families use, so
#                          rho IS a free factor. Shape is real, dependence is by
#                          design, and skew(D) is close to 0 because differencing two
#                          similarly right-skewed marginals under an exchangeable
#                          copula is nearly symmetric.
#
# Both centre and scale with the empirical column moments, which ARE the population
# constants of the bootstrap law, so the warning in ``_standardise_scaled`` about
# per-sample standardisation is respected. ``empirical_natural_rho`` reports the
# matrix's own dependence (the level ``empirical_copula`` has to match to be
# comparable) and ``empirical_effect`` the real standardised effect.


def _check_matrix(matrix):
    """Validate and return the (M, 2) loss matrix both empirical families need."""
    if matrix is None:
        raise ValueError("empirical families require a (M, 2) loss matrix")
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2 or mat.shape[1] != 2 or mat.shape[0] < 2:
        raise ValueError(f"loss matrix must have shape (M >= 2, 2), got {mat.shape}")
    if not np.all(np.isfinite(mat)):
        raise ValueError("loss matrix contains non-finite values")
    return mat


def _empirical_standardise(pair, matrix, sigma_a, sigma_b, theta):
    """Centre and scale drawn pairs by the matrix's own moments; shift to ``E[D] = theta``.

    Each column is centred on its own mean (a location shift changes neither shape
    nor dependence, and it is what makes ``E[D] = theta`` exact), but BOTH are
    divided by the SAME constant, the root-mean-square of the two column sds. A
    common divisor is an affine map of the pair, so the row bootstrap keeps the real
    joint shape exactly: skew(D) stays at the data's -1.55 and the real marginal
    variance ratio (2.17 vs 3.21 on the PEHE scale) is retained rather than being
    forced to 1. Dividing each column by its own sd would instead reweight the
    contrast and pull skew(D) to -0.98, discarding part of the realism this family
    is here to supply. ``sigma_a`` / ``sigma_b`` then multiply on top of the real
    ratio, so the design's variance-ratio factor still works.
    """
    mom = empirical_moments(matrix)
    scale = float(np.sqrt((mom[0][1] ** 2 + mom[1][1] ** 2) / 2.0))
    out = np.stack([(pair[..., 0] - mom[0][0]) / scale * sigma_a,
                    (pair[..., 1] - mom[1][0]) / scale * sigma_b], axis=-1)
    out[..., 1] -= theta
    return out


def _empirical_quantile(column, u):
    """Empirical inverse CDF of ``column`` evaluated at uniforms ``u``.

    Equal mass on each of the M observed values, so the mean and sd of the result
    are exactly the column's own mean and (population) sd -- which is what makes
    the centring in ``sample_batch`` exact rather than approximate.
    """
    order = np.sort(np.asarray(column, dtype=float))
    m = order.shape[0]
    idx = np.minimum((np.asarray(u) * m).astype(np.intp), m - 1)
    return order[idx]


def empirical_moments(matrix):
    """``((mean_A, sd_A), (mean_B, sd_B))`` of the two columns, population (ddof=0)."""
    matrix = np.asarray(matrix, dtype=float)
    return tuple((float(matrix[:, k].mean()), float(matrix[:, k].std(ddof=0)))
                 for k in (0, 1))


def empirical_natural_rho(matrix, kind="spearman"):
    """Dependence actually present in the supplied loss matrix.

    ``kind='spearman'`` is the one comparable with the design's ``rho``, because
    ``rho`` enters as a Gaussian-copula RANK dependence; ``kind='pearson'`` is
    reported alongside it for the manuscript.
    """
    matrix = np.asarray(matrix, dtype=float)
    if kind == "pearson":
        return float(np.corrcoef(matrix[:, 0], matrix[:, 1])[0, 1])
    return float(stats.spearmanr(matrix[:, 0], matrix[:, 1]).statistic)


def empirical_effect(matrix):
    """Real standardised paired effect ``E[D] / sd(D)`` of the supplied matrix.

    Reported so the manuscript can say where the actual case-study contrast sits
    on the ``theta`` axis of the E1 grid instead of asserting a round number.
    """
    matrix = np.asarray(matrix, dtype=float)
    d = matrix[:, 0] - matrix[:, 1]
    return float(d.mean() / d.std(ddof=1))

# Families whose *marginal* is skewed/heavy-tailed. Applying the same marginal to
# both methods and coupling them with an exchangeable Gaussian copula makes the
# pair (L_A, L_B) exchangeable when sigma_a == sigma_b, and an exchangeable pair
# has a SYMMETRIC difference: skew(D) is then identically 0, whatever the
# marginal. The P3-T4 sub-study measures Type I error AGAINST the standardised
# third moment of D, so it needs cells where that moment is non-zero.
#
# `asym = True` (the default) therefore applies the skewed marginal to method A
# only and gives method B a normal marginal of matching variance -- the same
# device the `mix` family already used, and a realistic scenario ("one method
# occasionally does much worse than the other"). Set `asym = False` to recover
# the exchangeable construction, which is a useful null case: it should show
# essentially nominal Type I error at every J.
_parametric = ["normal", "lognormal", "gamma", "mix", "beta", "t3"]
_SKEWED = ("lognormal", "gamma", "beta", "t3")


def _family_marginal(family, u, z):
    """Map copula uniforms ``u`` (and z-scores ``z``) onto the family's marginal."""
    if family == "t3":
        return stats.t.ppf(stats.norm.cdf(z), df=3)
    if family == "lognormal":
        return np.exp(0.5 * stats.norm.ppf(u))
    if family == "gamma":
        return stats.gamma.ppf(u, a=2.0)
    if family == "beta":
        mu, sd = 0.96, 0.02
        alb = (mu * (1 - mu)) / (sd * sd) - 1.0
        return stats.beta.ppf(u, mu * alb, (1 - mu) * alb)
    raise ValueError(f"no marginal registered for family {family!r}")


def sample_batch(family, R, J, rho=UNSET, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                 matrix=None, rng=None, master_seed=None, asym=True):
    """Draw ``(R, J, 2)`` paired losses for one family in one vectorised call.

    ``rng`` is a ``numpy.random.Generator``; if omitted a generator is seeded from
    ``master_seed`` (via a SeedSequence) so the block is reproducible and resumable.
    ``rho`` defaults to the ``UNSET`` sentinel, which means 0.0 for the copula
    families and "not applicable" for the ``empirical`` row bootstrap; see
    ``_Unset``. ``asym`` controls whether a skewed family is applied to method A
    only (the default, which gives a skewed contrast ``D``) or to both marginals
    (which makes the pair exchangeable and forces ``skew(D) = 0``); see ``_SKEWED``.
    """
    if rng is None:
        ss = np.random.SeedSequence(master_seed)
        rng = np.random.default_rng(ss)
    if family == "empirical":
        # Row bootstrap: the pair is taken from ONE observed replication, so the
        # real joint dependence -- and with it the real skew(D) = -1.55, which is
        # what this family exists to carry into the P3-T4 sub-study -- is retained
        # exactly. rho cannot be imposed on top of that, so requesting one is an
        # error rather than a silently ignored argument (see the note above).
        mat = _check_matrix(matrix)
        if not (rho is None or isinstance(rho, _Unset)):
            raise ValueError(
                "family 'empirical' is a row bootstrap: it reproduces the loss "
                "matrix's own dependence and cannot be set to rho=%r. Omit rho or "
                "pass rho=None, or use family 'empirical_copula' to control rho."
                % (rho,))
        idx = rng.integers(0, mat.shape[0], size=(R, J))
        return _empirical_standardise(mat[idx], mat, sigma_a, sigma_b, theta)
    if isinstance(rho, _Unset):
        rho = 0.0
    if rho is None:
        raise ValueError(f"family {family!r} requires a numeric rho")
    z = _gauss_z(rng, (R, J), rho)                     # (R, J, 2) z-scores
    if family == "normal":
        block = np.stack([z[..., 0] * sigma_a, z[..., 1] * sigma_b - theta], axis=-1)
        return block
    if family == "mix":
        block = np.stack([z[..., 0] * sigma_a, z[..., 1] * sigma_b - theta], axis=-1)
        cat = rng.uniform(size=(R, J)) < 0.02
        mag = rng.exponential(scale=10.0, size=(R, J))
        sign = rng.choice([-1.0, 1.0], size=(R, J))
        block[..., 1] += np.where(cat, sign * mag, 0.0)
        return block
    if family == "empirical_copula":
        mat = _check_matrix(matrix)
        u = stats.norm.cdf(z)                          # Gaussian copula uniforms
        pair = np.stack([_empirical_quantile(mat[:, 0], u[..., 0]),
                         _empirical_quantile(mat[:, 1], u[..., 1])], axis=-1)
        return _empirical_standardise(pair, mat, sigma_a, sigma_b, theta)
    if family in _SKEWED:
        u = stats.norm.cdf(z)                          # Gaussian copula uniforms
        x = _family_marginal(family, u, z)
        block = _standardise_scaled(x, sigma_a, sigma_b, theta, moments=_MOMENTS[family])
        if asym:
            # Method B keeps the standard-normal marginal implied by the same
            # copula draw, so it has the identical variance and the identical
            # rank dependence with A -- only the SHAPE differs, and D inherits
            # A's skewness. Without this the difference of two exchangeable
            # marginals is symmetric and the skewness factor is inert.
            block[..., 1] = z[..., 1] * sigma_b - theta
        return block
    raise ValueError(f"unknown family {family!r}; choose from {sorted(FAMILIES)}")


def contrast_skewness(family, rho=UNSET, sigma_a=1.0, sigma_b=1.0, asym=True,
                      n=400_000, seed=0, matrix=None):
    """Standardised third moment of ``D = L_A - L_B`` for a family/rho cell.

    This is the x-axis of the P3-T4 figure ("Type I error vs skewness of D_j"),
    so it has to be reported per cell rather than assumed. Estimated once by a
    large independent draw; its own MCSE is ``~sqrt(6/n)``.
    """
    if family == "empirical" and matrix is None:
        raise ValueError("empirical family needs its loss matrix to report skewness")
    b = sample_batch(family, 1, n, rho=rho, sigma_a=sigma_a, sigma_b=sigma_b,
                     theta=0.0, master_seed=seed, asym=asym, matrix=matrix)
    d = b[0, :, 0] - b[0, :, 1]
    d = d - d.mean()
    s = d.std()
    return float(np.mean(d ** 3) / s ** 3) if s > 0 else 0.0


def sample_pairs(family, n, rho=UNSET, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                 matrix=None, rng=None, seed=None):
    """Compatibility helper: draw ``n`` pairs as a ``(n, 2)`` array."""
    return sample_batch(family, 1, n, rho=rho, sigma_a=sigma_a, sigma_b=sigma_b,
                        theta=theta, matrix=matrix,
                        rng=(np.random.default_rng(seed) if seed is not None
                             else (rng if isinstance(rng, np.random.Generator)
                                   else np.random.default_rng())),
                        master_seed=None)[0]


def draw_rep_losses(family, R, J, rho=UNSET, sigma_a=1.0, sigma_b=1.0, theta=0.0,
                    matrix=None, master_seed=0, asym=True):
    """Vectorised ``(R, J, 2)`` draw seeded via SeedSequence``master_seed``."""
    return sample_batch(family, R, J, rho=rho, sigma_a=sigma_a, sigma_b=sigma_b,
                        theta=theta, matrix=matrix, master_seed=master_seed,
                        asym=asym)
