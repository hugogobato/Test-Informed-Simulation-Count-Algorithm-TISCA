# `notebooks/` — Google Colab notebooks

Notebooks that run on Google Colab (Python runtime, R via `Rscript` calls) to
build and execute the heavy parts of the project across the 33-session execution
model described in `REVISION_PLAN.md` §5.

| Notebook | Task | Purpose |
|---|---|---|
| `P0T4_build_rlib_bundle.ipynb` | **P0-T4** | compile the heavy R packages **once** into `tisca_rlib.tar.gz`, record `sessionInfo()`, publish, and restore-and-verify in a fresh session |
| `P1_math_verification.ipynb` | **P1-T2/T1 (light)** | numerically verify the v2 planning formulas (power modes M1-M5, the plan's §1.2 J table, variance-inflation, rooted-vs-unrooted estimand, the `scipy.stats.t` bug magnitude). Pure NumPy/SciPy, seconds to minutes; no R bundle needed. |
| `E1_harness_selfcheck.ipynb` | **P2-T3** | self-check the outer-MC harness: run the acceptance gate (bivariate normal, rho=0, oracle sigma -> Type I .05 and power .80 within 2 MCSE), smoke all designs D1..D6 and all loss families, emit a small results CSV. ~1 min on a free CPU runtime. |
| `E1_modA_C_D.ipynb` | **P3-T2 Modules A/C/D** | combined 1,323-cell runner: Module A (864), Module C (243), and Module D (216), with per-cell local checkpointing under `/content/TISCA_E1` and final download. The grid is **imported** from `tisca/python/tisca/outermc/e1_grid.py`, not restated in the notebook. |
| `E1_modB_joint_shard1.ipynb` | **P3-T2 Module B, joint family** | first 330 cells of the canonical 660-cell multiplicity grid. Every repetition generates one shared-control `(J, K+1)` loss block and tests all `K` contrasts jointly. The notebook runs the deterministic validation gate first, checkpoints each completed cell to Drive, records a manifest, and downloads every generated file. |
| `E1_modB_joint_shard2.ipynb` | **P3-T2 Module B, joint family** | second, disjoint set of 330 cells, with the same joint simulation, validation, checkpoint, manifest, resume, and download protocol. |
| `E1_modB_shard1.ipynb` | **superseded marginal Module B runner** | provenance artifact for the earlier two-column simulation. It changes the adjusted alpha with `K` but does not generate a joint family, so it must not be used for G1 family-level results. |
| `E1_modB_shard2.ipynb` | **superseded marginal Module B runner** | second provenance shard of the earlier two-column simulation. Use `E1_modB_joint_shard2.ipynb` for the revision results. |
| `E1_moduleA.ipynb` | **superseded P3-T2 runner** | earlier standalone Module-A runner. The notebook itself is no longer checked in; only its generator survives, marked superseded, because it carries a private 840-cell grid that predates the split of the empirical family. Use `E1_modA_C_D.ipynb`. |
| `E3_mvbcf_shard.ipynb` | **P3-T5(d)** | historical E3 worker entry point, regenerated from the first confirmatory shard; the ready-made per-shard notebooks are in `E3_shards/`. It restores the P0-T4 bundle, compiles the upstream `MVBCF_Code.cpp`, checkpoints each replication under `/content`, and verifies seed completeness. |
| `E3_round0_pilot_calibration.ipynb` | **P3-T5(e)** | Round 0 pilot (50 independent-seed reps of DGP1/n=500) plus the **P3-T5(e) gate**: compare the calibrated `stochtree::bcf` means to McJames et al. DGP1 Table 2. No confirmatory shard runs until this passes. |
| `E3_round0_dual_bcf_pilot.ipynb` | **BCF diagnostic Round 0** | Paired DGP1/n=500 pilot using the corrected `stochtree::bcf` translation and the original `bcf` package settings. Reports separate PEHE and coverage bands in one CSV. |
| `E4_bibliometric_analysis.ipynb` | **P3-T1(b)** | clones the repo when uploaded alone, runs `code_bibliometrics.py` if `bibliometric_coded.csv` is absent, then recounts the section **from code** (no percentage literals), corrects denominators (99/88), and emits `Fig_bib_*.png`. |
| `E2_design_comparison.ipynb` | **P3-T3** | analysis-only over the E1 output: the six designs D1-D6 on achieved level, achieved power and E[J], the matched design-versus-design table, and the generated verdict sentence answering Reviewer 2 paragraph 4. Writes `results/E2/`. |
| `E1b_skewness_calibration.ipynb` | **P3-T4** | the evidence replacing the deleted "J > 30" claim: Type I error of the paired t and of the studentized paired bootstrap against J, panelled by the standardised third moment of D, with Berry-Esseen as the explanatory frame. Runs its own fixed-J sweep (the E1 designs each choose their own J). Writes `results/E1b/`. |
| `E1c_no_difference.ipynb` | **P3-T7** | the "what if there is no difference" case, separated into equal-expected-loss, same-model-twice, and perfectly-correlated. Shows where v1 runs to the budget cap, and drafts the Section 3.10 paragraph. Writes `results/E1c/`. |
| `E5_generality_demo.ipynb` | **P3-T6** | the complete v2 workflow on AR(1) forecasting, where the ranking of the four competitors is known in closed form, so every output is checked rather than reported. Writes `results/E5/`. About a minute. |
| `E3_seed_verification.ipynb` | **P3-T5 / P0-T2 acceptance** | re-runs 3 randomly chosen seeds and compares them to the stored rows exactly, plus the four-way `mc.cores` x shard-offset identity test. **Run 2026-08-08, all green**, 0 of 159 columns differing everywhere. Reports per-model **run-level reproducibility**; it does **not** settle whether each model honours its explicit seed (`CALIBRATION.md` §3.6), because `run_cell.R` positions the global RNG stream before every fit, so a model that ignored its seed argument would reproduce identically. ~50 min, one session. Writes `results/E3/seed_verification.csv` and the block to paste into `CALIBRATION.md`. |
| `E3_seed_honouring.ipynb` | **P0-T2 §3.6** | the discriminating comparison the verification notebook could not make: same model, same data, same explicit seed, **two different global RNG states**, with a different-seed positive control so a seed-insensitive fit cannot masquerade as a pass. Covers `dbarts` (propensity, and underneath `bartCause::bartc`) and `stochtree::bcf`; `fast_bart()` and `MultiskewBART` are out of scope because neither exposes a seed argument. **Re-runs no replication and reads no stored row**: one DGP1 n=100 dataset, nine fits, ~20 min. Writes `results/E3/seed_honouring.csv`. |

`P0T4_build_rlib_bundle.ipynb` is the main notebook to run once during Phase 0.
The E3 shard set is generated with the public Drive folder containing the
tarball and `tisca_rlib.sha256`; each E3 notebook downloads that folder and
verifies the embedded digest before restoring `.libPaths(<dest>/tisca_rlib/rlib)`.

`P1_math_verification.ipynb` is the Phase-1 notebook you run yourself (single
run, on the laptop or a free Colab CPU runtime, no R bundle). Every cell prints
`[PASS]`/`[FAIL]`; a clean run is the audit record for the Phase-1 specs in
`docs/`.

`E1_harness_selfcheck.ipynb` is the Phase-2 notebook to run yourself (P2-T3):
it pulls the current harness (local checkout or the public GitHub repo), runs the
acceptance gate, and smoke-tests all six designs and all loss families.

The complete Phase-3 E1 package is `E1_modA_C_D.ipynb` plus
`E1_modB_joint_shard1.ipynb` and `E1_modB_joint_shard2.ipynb`. Module B now
reports marginal rejection probability, FWER, conjunctive family power,
disjunctive family power, and BH FDR as distinct estimands. Its `K+1` loss block
shares method A across all contrasts. The empirical `K=1` case uses the committed
`legacy/Paper_Experiments/DGP1_500_results.csv`; empirical configurations with
`K>1` use the documented synthetic Gaussian-copula extension because the source
file contains only one observed benchmark.

Run `E1_harness_selfcheck.ipynb` first. Then the A/C/D notebook and the two joint
Module B shards may run concurrently in separate Colab sessions. Each joint shard
also runs `selfcheck_module_b_joint.py` before the expensive grid. Its checkpoints
live under the configured Drive directory and resume by canonical `cell_id`, so a
restart does not append duplicates. The old Module B notebooks and
`E1_moduleA.ipynb` remain available as provenance artifacts but are not part of
the complete three-notebook launch.

`E3_mvbcf_shard.ipynb` is the E3 case-study worker: it restores the P0-T4 bundle,
downloads `run_cell.R` (this repo) and the upstream `MVBCF_Code.cpp`, compiles
the cpp with `sourceCpp`, and runs one shard, checkpointing each replication under
`/content`. `E3_round0_pilot_calibration.ipynb` runs the independent pilot and
the P3-T5(e) `stochtree::bcf` calibration gate that must pass before any
confirmatory shard launches. `E3_round0_dual_bcf_pilot.ipynb` is a separate
paired diagnostic: it runs the corrected stochtree configuration and the paper
`bcf` implementation on the same 50 seeds. Its diagnostic target failures do
not stop the notebook, so both method results remain available for comparison.

The notebooks are generated from `notebooks/_generators/` (`build_p0t4_bundle_nb.py`,
`build_p1_verification_nb.py`, `build_e1_selfcheck_nb.py`,
`build_e1_moduleA_nb.py`, `build_e1_modules_nb.py`,
`build_e1_modB_joint_notebooks.py`, `build_e3_notebooks.py`,
`build_e3_dual_round0_notebook.py`, `build_e4_analysis_nb.py`). Regenerate the
dual BCF notebook with
`python notebooks/_generators/build_e3_dual_round0_notebook.py` and regenerate the
three legacy E1 runners with
`python notebooks/_generators/build_e1_modules_nb.py`. Regenerate the two joint
Module B shards with
`python notebooks/_generators/build_e1_modB_joint_notebooks.py`. Keep each
generator and its `.ipynb` files in sync when you edit either.


## Running E1 without Colab

`E1_modA_C_D.ipynb` and the two Module-B shards remain the Colab path, but the whole
E1 grid is a few minutes of pure NumPy on a workstation. Prefer:

```bash
python experiments/E1_operating_characteristics/run_e1_grid.py        # all 1,983 cells
python experiments/E1_operating_characteristics/run_e1_grid.py -j 6   # fewer workers
```

It writes `results/E1/operating_characteristics.csv` directly, is resumable by
`cell_id`, and asserts grid completeness before finishing. **Worker count is chosen
from available RAM, not from core count**: each worker holds an `(R, J_max, 2)` block
and peaks around 0.5 GB, so the default takes at most half of currently available
memory. Pass `-j` explicitly only if you know the machine is otherwise idle.

The canonical grid lives in `tisca/python/tisca/outermc/e1_grid.py` and is imported by
both the local runner and the notebook generator, so the two cannot drift.

The genuine joint Module B has a separate resumable runner. Its bare invocation
runs the canonical grid, while the other two required entry points inspect the
grid and resume validated per-cell checkpoints:

```bash
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --dry-run
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --resume
```

For Colab, generate the notebooks locally, open each generated notebook in a
separate runtime, and run all cells. After both Drive shard directories have been
downloaded or mounted, combine them with:

```bash
python notebooks/_generators/build_e1_modB_joint_notebooks.py
python experiments/E1_operating_characteristics/combine_e1_module_b_joint.py \
  results/E1/module_b_joint_shard1 \
  results/E1/module_b_joint_shard2 \
  --outdir results/E1
```

The combiner rejects missing cells, unexpected seeds, and duplicate cell IDs. A
successful combination contains exactly the 660 canonical Module B rows and
writes `module_b_joint_operating_characteristics.csv`,
`module_b_joint_covariance.csv`, and `module_b_joint_manifest.csv`.

## The four analysis notebooks

`E2_design_comparison.ipynb` and `E1c_no_difference.ipynb` read
`results/E1/operating_characteristics.csv`, so run the E1 grid first and re-run
these two whenever it changes. `E1b_skewness_calibration.ipynb` does **not** read it
(it runs its own fixed-J sweep, because each E1 design chooses its own J), and
`E5_generality_demo.ipynb` is self-contained; both still depend on the `tisca`
package, so re-run them when `planning.py` or `outermc/` changes.

None of the four needs Colab or R: they are pure NumPy, pandas and matplotlib, they
detect an enclosing checkout before falling back to cloning, and they write to
absolute `results/` and `figures/` paths. Executing their code cells in order from
`notebooks/` reproduces every output, which is how they were last refreshed. All
four keep the `google.colab.files.download` fallback. Regenerate them from
`notebooks/_generators/build_e2_design_comparison_nb.py`,
`build_e1b_skewness_nb.py`, `build_e1c_no_difference_nb.py` and
`build_e5_generality_nb.py`; shared boilerplate is in `_generators/_nbcommon.py`.
