# TISCA: Test-Informed Simulation Count Algorithm

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Repository for:** Souto, H. G. and Louzada Neto, F. "Beyond Arbitrary
Replications: A Principled Approach to Simulation Design in Causal Inference."
The current manuscript is under review.

This repository hosts the current **TISCA v2** code, experiments, results, and
figures, together with an auditable environment specification. The **original
v1** code is preserved under [`legacy/`](legacy/README.md) for historical
reference and audit.

> **Status: current analysis package, manuscript under review.** The repository
> contains the current TISCA v2 implementation, experiment specifications,
> fixed-seed results, and case-study notebooks. The 2024 arXiv preprint is an
> older version of the work and describes TISCA v1. Use the current v2 code and
> manuscript associated with this repository for new analyses, reproduction,
> and scientific references. The arXiv version is retained only for historical
> provenance.

## Installation

The current TISCA v2 libraries are distributed from the [TISCA GitHub
repository](https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA).
The installation commands below install the current library, not the historical
code under `legacy/`.

### Python

Python 3.9 or later is required. Install the Python reference implementation
directly from GitHub with:

```bash
python -m pip install "git+https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
```

For development of the core library, clone the repository and install the
development dependencies in editable mode:

```bash
git clone https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git
cd Test-Informed-Simulation-Count-Algorithm-TISCA
python -m pip install -e ".[dev]"
```

For the full repository experiments, also install the pinned broader runtime
requirements:

```bash
python -m pip install -r requirements.txt
```

Then import the v2 modules, for example:

```python
from tisca import inference, planning
```

### R

R 4.0 or later is required. Install the R package from its `tisca/` subdirectory
using `remotes`:

```r
install.packages("remotes")
remotes::install_github(
  "hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA",
  subdir = "tisca"
)
library(tisca)
```

To install from a local clone instead, run `R CMD INSTALL tisca` from the
repository root. The heavy packages used by the MVBCF case study are separate
from the core TISCA library. See [`env/install_R_dependencies.R`](env/install_R_dependencies.R)
and the Colab bundle instructions above when reproducing that case study.

## Basic usage

Researchers should first define the replication-level estimand and the paired
contrast between methods. For a loss where lower values are better, let
`D_j = L_{A,j} - L_{B,j}`. A negative planning alternative therefore represents
an expected advantage for method A. The default workflow is to run an
independent pilot, estimate the standard deviation of `D_j`, use TISCA to plan
the confirmatory replication count, discard the pilot rows, and then run the
confirmatory block with common random numbers for both methods.

The following Python example plans for both an MCSE target and 80% power for a
one-sided superiority claim. Replace `pilot_results` and `confirm_results` with
the losses produced by the researcher's own simulation code:

```python
import numpy as np
from tisca import inference, planning

# Replace `pilot_results` with the researcher's pilot output. It must contain
# one paired loss for each of J0 independent pilot replications.
pilot_A = np.asarray(pilot_results["method_A"], dtype=float)
pilot_B = np.asarray(pilot_results["method_B"], dtype=float)
pilot_D = pilot_A - pilot_B

J0 = pilot_D.size
J_final, sigma_ub = planning.required_J(
    np.std(pilot_D, ddof=1),
    J0,
    gamma=0.20,
    mode="M2",             # lower-is-better directional superiority
    delta=-0.20,            # planned advantage for method A
    target_mcse=0.05,
    target_power=0.80,
    alpha=0.05,             # use 0.05 / K for K Bonferroni-planned contrasts
    J_max=10000,
)
print({"J_final": J_final, "sigma_upper_bound": sigma_ub})

# After planning, run J_final new, independent confirmatory replications using
# common seeds and store the two loss vectors in `confirm_results`.
confirm_A = np.asarray(confirm_results["method_A"], dtype=float)
confirm_B = np.asarray(confirm_results["method_B"], dtype=float)
confirm_D = confirm_A - confirm_B

final = inference.paired_t(confirm_D, alternative="less")
print({"estimate": final["estimate"], "p_value": final["p_value"]})
```

The equivalent R workflow is:

```r
library(tisca)

# Replace `pilot_results` with the researcher's pilot output. It must contain
# one paired loss for each of J0 independent pilot replications.
pilot_A <- pilot_results$method_A
pilot_B <- pilot_results$method_B
pilot_D <- pilot_A - pilot_B

J0 <- length(pilot_D)
sigma_upper <- sigma_ub(sd(pilot_D), J0 = J0, gamma = 0.20)
power_plan <- solve_power_J(
  mode = "M2", delta = -0.20, sigma_D = sigma_upper,
  alpha_adj = 0.05, target_power = 0.80, J_max = 10000
)
mcse_plan <- solve_mcse_J(sigma_upper, m = 0.05, J_max = 10000)
J_final <- combine_J(c(power_plan$J, mcse_plan$J), J_max = 10000)$J_final
print(list(J_final = J_final, sigma_upper_bound = sigma_upper))

# After planning, run J_final new, independent confirmatory replications using
# common seeds and store the two loss vectors in `confirm_results`.
confirm_A <- confirm_results$method_A
confirm_B <- confirm_results$method_B
contrast <- contrast_from_columns(confirm_A, confirm_B)
final <- paired_t(contrast$D, alternative = "less")
print(list(estimate = final$estimate, p_value = final$p_value))
```

For two-sided equality use `mode = "M1"` and `alternative = "two-sided"`.
Modes `M3`, `M4`, and `M5` support minimum-effect, non-inferiority, and
equivalence claims and require a `margin`/`Delta`. For multiple pre-specified
contrasts, apply the same multiplicity-adjusted level in planning and final
inference, and take the maximum required `J` across contrasts. See
[`docs/tisca_v2_spec.md`](docs/tisca_v2_spec.md) for the full protocol and
[`docs/power_target_guidance.md`](docs/power_target_guidance.md) for guidance on
when to use the precision layer, the power layer, or both.

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
power-search loop is superseded in v2. The formal v2 specification and planning
guidance are in [`docs/tisca_v2_spec.md`](docs/tisca_v2_spec.md) and
[`docs/power_target_guidance.md`](docs/power_target_guidance.md).

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
env/                    # R dependency installer and library-bundle digest
notebooks/              # Colab notebooks (library bundle, pilots, shards)
legacy/                 # v1 submission code, preserved for audit (read-only mindset)
LICENSE                 # MIT
```

## Environment

- **R:** `env/install_R_dependencies.R` is the local installation entry point
  (R 4.3.3 baseline). On Google Colab, the heavy packages (`stochtree`,
  `dbarts`, `bartCause`, `skewBART`, `mvbcf`, `mvtnorm`, `bcf`, `scoringRules`,
  `matrixStats`, `progress`, `MCS`) are compiled into a library bundle and
  restored by the [`P0T4_build_rlib_bundle.ipynb`](notebooks/P0T4_build_rlib_bundle.ipynb)
  workflow. The expected bundle digest is recorded in
  [`env/tisca_rlib.sha256`](env/tisca_rlib.sha256).
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

## Version and citation

The current manuscript is under review, so publication and venue information is
intentionally not included here. If you use TISCA, cite the current manuscript
version associated with this repository and use the v2 implementation. The
arXiv record below is provided only as a link to the older v1 version. It should
not be used as the basis for new analyses or for documenting the current TISCA
v2 method.

Historical arXiv reference (older v1 version): [arXiv:2409.05161](https://arxiv.org/abs/2409.05161)

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
