# `notebooks/` — Google Colab notebooks

Notebooks that run on Google Colab (Python runtime, R via `Rscript` calls) to
build and execute the heavy parts of the project across the 33-session execution
model described in `REVISION_PLAN.md` §5.

| Notebook | Task | Purpose |
|---|---|---|
| `P0T4_build_rlib_bundle.ipynb` | **P0-T4** | compile the heavy R packages **once** into `tisca_rlib.tar.gz`, record `sessionInfo()`, publish, and restore-and-verify in a fresh session |
| `P1_math_verification.ipynb` | **P1-T2/T1 (light)** | numerically verify the v2 planning formulas (power modes M1-M5, the plan's §1.2 J table, variance-inflation, rooted-vs-unrooted estimand, the `scipy.stats.t` bug magnitude). Pure NumPy/SciPy, seconds to minutes; no R bundle needed. |
| `E1_harness_selfcheck.ipynb` | **P2-T3** | self-check the outer-MC harness: run the acceptance gate (bivariate normal, rho=0, oracle sigma -> Type I .05 and power .80 within 2 MCSE), smoke all designs D1..D6 and all loss families, emit a small results CSV. ~1 min on a free CPU runtime. |
| `E1_modA_C_D.ipynb` | **P3-T2 Modules A/C/D** | combined 1,299-cell runner: Module A (840), Module C (243), and Module D (216), with per-cell local checkpointing under `/content/TISCA_E1` and final download. |
| `E1_modB_shard1.ipynb` | **P3-T2 Module B** | first 300 cells of the 600-cell multiplicity grid, with `R=2000`, `B=999`, local checkpointing, and final download. |
| `E1_modB_shard2.ipynb` | **P3-T2 Module B** | second 300 cells of the 600-cell multiplicity grid, with `R=2000`, `B=999`, local checkpointing, and final download. |
| `E1_moduleA.ipynb` | **legacy P3-T2 runner** | earlier standalone Module-A runner retained for provenance; use `E1_modA_C_D.ipynb` for the complete planned A/C/D sweep. |
| `E3_mvbcf_shard.ipynb` | **P3-T5(d)** | historical E3 worker entry point, regenerated from the first confirmatory shard; the ready-made per-shard notebooks are in `E3_shards/`. It restores the P0-T4 bundle, compiles the upstream `MVBCF_Code.cpp`, checkpoints each replication under `/content`, and verifies seed completeness. |
| `E3_round0_pilot_calibration.ipynb` | **P3-T5(e)** | Round 0 pilot (50 independent-seed reps of DGP1/n=500) plus the **P3-T5(e) gate**: compare the calibrated `stochtree::bcf` means to McJames et al. DGP1 Table 2. No confirmatory shard runs until this passes. |
| `E4_bibliometric_analysis.ipynb` | **P3-T1(b)** | clones the repo when uploaded alone, runs `code_bibliometrics.py` if `bibliometric_coded.csv` is absent, then recounts the section **from code** (no percentage literals), corrects denominators (99/88), and emits `Fig_bib_*.png`. |

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
`E1_modB_shard1.ipynb` and `E1_modB_shard2.ipynb`. They use the same P2-T3
NumPy/SciPy harness. The empirical family is derived automatically from the
committed `legacy/Paper_Experiments/DGP1_500_results.csv`, using the
pre-specified `mvbcf_pehe1` versus `bcf_pehe1` pair.

Run `E1_harness_selfcheck.ipynb` first. Then the A/C/D notebook and the two Module
B shards may run concurrently in separate Colab sessions. Their CSVs are written
under `/content` and downloaded at the end. The old `E1_moduleA.ipynb` remains available
as a provenance artifact but is not part of the complete three-notebook launch.

`E3_mvbcf_shard.ipynb` is the E3 case-study worker: it restores the P0-T4 bundle,
downloads `run_cell.R` (this repo) and the upstream `MVBCF_Code.cpp`, compiles
the cpp with `sourceCpp`, and runs one shard, checkpointing each replication under
`/content`. `E3_round0_pilot_calibration.ipynb` runs the independent pilot and
the P3-T5(e) `stochtree::bcf` calibration gate that must pass before any
confirmatory shard launches.

The notebooks are generated from `notebooks/_generators/` (`build_p0t4_bundle_nb.py`,
`build_p1_verification_nb.py`, `build_e1_selfcheck_nb.py`,
`build_e1_moduleA_nb.py`, `build_e1_modules_nb.py`, `build_e3_notebooks.py`,
`build_e4_analysis_nb.py`). Regenerate the three complete E1 runners with
`python notebooks/_generators/build_e1_modules_nb.py`. Keep the generator and the
`.ipynb` files in sync when you edit either.
