# `notebooks/` — Google Colab notebooks

Notebooks that run on Google Colab (Python runtime, R via `Rscript` calls) to
build and execute the heavy parts of the project across the 33-session execution
model described in `REVISION_PLAN.md` §5.

| Notebook | Task | Purpose |
|---|---|---|
| `P0T4_build_rlib_bundle.ipynb` | **P0-T4** | compile the heavy R packages **once** into `tisca_rlib.tar.gz`, record `sessionInfo()`, publish, and restore-and-verify in a fresh session |
| `P1_math_verification.ipynb` | **P1-T2/T1 (light)** | numerically verify the v2 planning formulas (power modes M1-M5, the plan's §1.2 J table, variance-inflation, rooted-vs-unrooted estimand, the `scipy.stats.t` bug magnitude). Pure NumPy/SciPy, seconds to minutes; no R bundle needed. |
| `E1_harness_selfcheck.ipynb` | **P2-T3** | self-check the outer-MC harness: run the acceptance gate (bivariate normal, rho=0, oracle sigma -> Type I .05 and power .80 within 2 MCSE), smoke all designs D1..D6 and all loss families, emit a small results CSV. ~1 min on a free CPU runtime. |
| `E1_modA_C_D.ipynb` | **P3-T2 Modules A/C/D** | combined 1,299-cell runner: Module A (840), Module C (243), and Module D (216), with per-cell Drive checkpointing to `E1_modA_C_D_results.csv`. |
| `E1_modB_shard1.ipynb` | **P3-T2 Module B** | first 300 cells of the 600-cell multiplicity grid, with `R=2000`, `B=999`, and per-cell Drive checkpointing. |
| `E1_modB_shard2.ipynb` | **P3-T2 Module B** | second 300 cells of the 600-cell multiplicity grid, with `R=2000`, `B=999`, and per-cell Drive checkpointing. |
| `E1_moduleA.ipynb` | **legacy P3-T2 runner** | earlier standalone Module-A runner retained for provenance; use `E1_modA_C_D.ipynb` for the complete planned A/C/D sweep. |
| `E3_mvbcf_shard.ipynb` | **P3-T5(d)** | the parameterised E3 worker: edit Cell 2 (`SHARD_ID, DGP, N, SEED_START, SEED_END, MC_CORES, MODE`), restore the P0-T4 bundle, compile the upstream `MVBCF_Code.cpp`, run `Rscript run_cell.R`, checkpoint each replication to a Drive CSV, verify seed completeness, and mirror a copy via `files.download`. One notebook = one 2-vCPU session. |
| `E3_round0_pilot_calibration.ipynb` | **P3-T5(e)** | Round 0 pilot (50 independent-seed reps of DGP1/n=500) plus the **P3-T5(e) gate**: compare the calibrated `stochtree::bcf` means to McJames et al. DGP1 Table 2. No confirmatory shard runs until this passes. |
| `E4_bibliometric_analysis.ipynb` | **P3-T1(b)** | recount the bibliometric section **from code** (no percentage literals), correct denominators (99/88), regenerate Fig 1/2/3a/3b with larger fonts, and add the publisher/OA stratification. Runs on laptop or a free Colab runtime; emits `Fig_bib_*.png`. |

`P0T4_build_rlib_bundle.ipynb` is the main notebook to run once during Phase 0.
Once it has produced and published the tarball, paste the direct-download URL
into Cell 9; every later worker notebook will restore the bundle with a short
`wget <URL> && tar xzf && .libPaths(<dest>/tisca_rlib/rlib)`.

`P1_math_verification.ipynb` is the Phase-1 notebook you run yourself (single
run, on the laptop or a free Colab CPU runtime, no R bundle). Every cell prints
`[PASS]`/`[FAIL]`; a clean run is the audit record for the Phase-1 specs in
`docs/`.

`E1_harness_selfcheck.ipynb` is the Phase-2 notebook to run yourself (P2-T3):
it pulls the current harness (local checkout or the public GitHub repo), runs the
acceptance gate, and smoke-tests all six designs and all loss families.

The complete Phase-3 E1 package is `E1_modA_C_D.ipynb` plus
`E1_modB_shard1.ipynb` and `E1_modB_shard2.ipynb`. They use the same P2-T3
NumPy/SciPy harness. Each notebook finds `E1_empirical_loss_matrix.npy` (or
`.csv`, with two columns) in the `TISCA_E1` Drive folder, or offers the standard
Colab upload dialog, because Modules A and B include the empirical family.

Run `E1_harness_selfcheck.ipynb` first. Then the A/C/D notebook and the two Module
B shards may run concurrently in separate Colab sessions. Their CSVs are written
to Drive and are resumable by cell. The old `E1_moduleA.ipynb` remains available
as a provenance artifact but is not part of the complete three-notebook launch.

`E3_mvbcf_shard.ipynb` is the E3 case-study worker: it restores the P0-T4 bundle,
downloads `run_cell.R` (this repo) and the upstream `MVBCF_Code.cpp`, compiles
the cpp with `sourceCpp`, and runs one shard, checkpointing each replication to
a Drive CSV. `E3_round0_pilot_calibration.ipynb` runs the independent pilot and
the P3-T5(e) `stochtree::bcf` calibration gate that must pass before any
confirmatory shard launches.

The notebooks are generated from `notebooks/_generators/` (`build_p0t4_bundle_nb.py`,
`build_p1_verification_nb.py`, `build_e1_selfcheck_nb.py`,
`build_e1_moduleA_nb.py`, `build_e1_modules_nb.py`, `build_e3_notebooks.py`,
`build_e4_analysis_nb.py`). Regenerate the three complete E1 runners with
`python notebooks/_generators/build_e1_modules_nb.py`. Keep the generator and the
`.ipynb` files in sync when you edit either.
