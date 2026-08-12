# TISCA: Test-Informed Simulation Count Algorithm

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Repository for:** Souto, H. G. and Louzada Neto, F. "Beyond Arbitrary
Replications: A Principled Approach to Simulation Design in Causal Inference."
International Journal of Data Science and Analytics (major revision, 2026).

This repository hosts the **TISCA v2** code, experiments, results, and figures
that accompany the revised manuscript, together with an auditable environment
specification. The **original v1** submission code is preserved, untouched, under
[`legacy/`](legacy/README.md) for audit.

> **Status: revised-analysis package.** The repository contains the TISCA
> implementation, the revised experiment specifications, the fixed-seed
> results, and the notebooks used for the case study. The original submission
> code remains under [`legacy/`](legacy/README.md) for audit and is not used by
> the revised analysis. Before a public release, the repository snapshot and
> the manuscript source must be synchronized.

## What TISCA v2 is

TISCA v2 is a **two-layer simulation-design protocol** for choosing the number
of Monte Carlo replications `J` in a simulation study:

- **Design layer (default):** choose `J` from a Monte Carlo precision target
  (MCSE or CI half-width) on pre-specified paired estimands. This follows
  Morris, White & Crowther (2019), Burton et al. (2006) and Koehler, Brown &
  Haneuse (2009).
- **Decision layer (optional):** if the study makes a confirmatory comparative
  claim, add a power target for the pre-specified hypothesis, computed at the
  same sidedness and the same multiplicity-adjusted `α` as the final test.

`J_final = max` over comparisons and over whichever layers are active. The
default procedure is **two-stage**: an independent-seed pilot sizes the
confirmatory run, whose replications are then analysed with **paired contrasts**
(common replications, common random numbers). When more than two models are
compared, a bootstrap **Model Confidence Set** (Hansen, Lunde & Nason, 2011)
reports the set of models indistinguishable from the best. v1's unpaired Welch
power-search loop is superseded in v2; see the revision rationale in
`docs/` (and `REVISION_PLAN.md` at the paper directory).

## Repository structure

```
tisca/                  # installable package (R and Python reference)
  R/                    # tisca v2 R functions
  python/tisca/         # tisca v2 Python functions (reference implementation)
  tests/                # parity tests R <-> Python
experiments/
  E1_operating_characteristics/   # outer-MC study of the whole procedure
  E2_design_comparison/           # two-stage vs adaptive vs fixed-J verdict
  E3_mvbcf_casestudy/             # MVBCF case-study re-run (DGP1-DGP3)
  E4_bibliometrics/               # bibliometric re-coding and recount
  E5_generality_demo/             # non-causal generality demonstration
results/                # committed CSVs, one per experiment, versioned
figures/                # figures emitted by experiments
docs/                   # specs: estimand table, TISCA v2 spec, seed/RNG protocol
env/                    # R dependency installer (renv.lock arrives with P0-T4)
notebooks/              # Colab notebooks (library bundle, pilots, shards)
legacy/                 # v1 submission code, preserved for audit (read-only mindset)
LICENSE                 # MIT
```

## Environment

- **R:** `env/install_R_dependencies.R` is the current source of truth (R 4.3.3
  baseline). A genuine `renv.lock` with pinned versions is emitted by the P0-T4
  bundle build (`renv::snapshot()` over the built library) and committed once
  the bundle exists — a lockfile without resolved versions would fail
  `renv::restore()`, so none is committed before then. On Google Colab the heavy packages
  (`stochtree`, `dbarts`, `bartCause`, `skewBART`, `mvbcf`, `mvtnorm`,
  `bcf`, `scoringRules`, `matrixStats`, `progress`, `MCS`) are precompiled into a
  single library bundle and restored in seconds: see the
  `notebooks/P0T4_build_rlib_bundle.ipynb` notebook and the seed/RNG protocol.
- **Python:** see `environment.yml` and `requirements.txt`.

Every experiment records its seeds per replication and keeps its model-fitting
stream separate from its data-generation stream, per
[`docs/seed_rng_protocol.md`](docs/seed_rng_protocol.md). Shards assert seed
completeness (no gaps, no duplicates) when concatenated.

## Reproducibility

`./run_all.sh` checks the repository layout, imports the Python reference
implementation, and runs the E1 acceptance gate. The E4 coding scripts, E3
collection and analysis scripts, and E5 verification outputs are run through
their experiment-specific entry points documented alongside the corresponding
artifacts. The fixed seeds and released results permit the manuscript numbers to
be checked without rerunning the full model-fitting campaign.

The `stochtree::bcf` benchmark diagnostic is recorded in
[`experiments/E3_mvbcf_casestudy/CALIBRATION.md`](experiments/E3_mvbcf_casestudy/CALIBRATION.md),
and the real-driver seed checks are recorded in `results/E3/`. The case-study
analysis is retrospective: the final pilot designation and benchmark diagnostic
were completed alongside or after generation of the full fixed-seed block. The
original authors' code is linked from the paper and the methods, not copied
here.

## Citation

If you use TISCA in your research, please cite the paper (and the v1 arXiv
preprint if you rely on the v1 algorithm):

```
@misc{https://doi.org/10.48550/arxiv.2409.05161,
  doi = {10.48550/ARXIV.2409.05161},
  url = {https://arxiv.org/abs/2409.05161},
  author = {Souto,  Hugo Gobato and Neto,  Francisco Louzada},
  title = {Beyond Arbitrary Replications: A Principled Approach to Simulation Design in Causal Inference},
  publisher = {arXiv},
  year = {2024}
}
```

## License

MIT. See the [LICENSE](LICENSE) file. The MVBCF case study reproduces the
design of McJames et al.; their code is attributed and linked from the paper
rather than distributed here.
