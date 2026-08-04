# P0-T2 — Seed and RNG protocol

**Status:** draft for review  
**Revision:** 2026-08-04  
**Applies to:** every experiment in the revised `Test-Informed-Simulation-Count-Algorithm-TISCA` repository, but most critically to `E3_mvbcf_casestudy` (the MVBCF re-run) and `E1_operating_characteristics` (the outer Monte Carlo study).

## 1. Why this exists

Reviewer IJDA #14 and the final report note that the v1 replication used
`seed_val <- as.numeric(Sys.time())` and never stored the seed, so the published
replication is not reproducible at all (§1.4 of the revision plan). This document
is the fix. It specifies, once, how seeds are allocated, how the data-generation
and model-fitting streams are kept apart, and how independence of the pilot from
the confirmatory run is guaranteed.

## 2. Terminology

- **Replication** `j` = one complete run of the pipeline for one DGP/training-size/method cell. In the case study a replication produces one row of the output CSV.
- **Data-generation stream** = the RNG draws used to create the covariates `X`, the treatment `Z` and the potential outcomes `Y(0), Y(1)` for one replication.
- **Model-fitting stream** = the RNG draws consumed by each stochastic model (BCF/BART/MVBART), including MCMC sampling, warm starts, and Gibbs steps.
- **Pilot** = the first (independent-seed) `J0` replications used only to size the confirmatory run (`σ̂_D`, `Σ̂_D`). Pilot rows are **discarded from final inference** (plan P1-T2 §5).
- **Confirmatory run** = the `J` replications whose rows are used for the final reported inference.

## 3. The rule (R)

### 3.1 Stream separation
Every replication `j` is initialised with **two** independent seeds derived from a single base seed `S_j`:

```
data_seed_j = S_j                 (used without any model draws before it)
fit_seed_j  = S_j + 2_000_000_000 (a fixed, documented offset far outside the pilot block)
```

Why an offset and not a second `set.seed`: it is trivially auditable and keeps the two streams disjoint even if a model library resets or advances the RNG mid-run. The offset is a constant; it never collides with the pilot block (below) because that block lives in `[1_000_001, 1_000_000 + J0]`.

### 3.2 Stream in the field names
Each row of the output CSV records exactly the stream seeds actually used:

| column | meaning |
|---|---|
| `seed_data` | the data-generation seed for this replication |
| `seed_fit` | the model-fitting seed for this replication |
| `rng_kind_data` / `rng_kind_fit` | `"L'Ecuyer-CMRG"` (R) or `"numpy.SeedSequence"` (Python) |

Storing the seed *in every row* is the minimum reproducibility fix and matches the original authors' own pattern (`GitHub_DGP1.R` stores `seed_val`).

### 3.3 L'Ecuyer-CMRG with parallel streams (the standard recommended path)
For a fixed `J` cell on a single node we use `parallel::nextRNGStream`:

```r
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(S_cell)                       # master seed for the cell
data_streams <- vector("list", J)
fit_streams  <- vector("list", J)
data_streams[[1]] <- .Random.seed
fit_streams[[1]]  <- .Random.seed + 2e9
for (j in seq_len(J)[-1]) {
  data_streams[[j]] <- nextRNGStream(data_streams[[j - 1]])
  fit_streams[[j]]  <- nextRNGStream(fit_streams[[j - 1]])
}
# each worker sets: RNGkind("L'Ecuyer-CMRG"); assign(".Random.seed", data_streams[[j]], envir=.GlobalEnv)
# then draws data; then assign(".Random.seed", fit_streams[[j]], ...) and fits the model.
```

`nextRNGStream` guarantees the streams for different `j` are (stochastically) disjoint. Workers consume their pre-assigned stream; they never call `set.seed(j)` themselves, because `set.seed(j)` with a fixed `j` uses the *current* RNG kind and does not give the disjoint-stream guarantee across parallel workers.

### 3.4 The idiomatic `set.seed(base + j)` rule, and when it is acceptable
The plan text (`P0-T2`) phrases the rule as "`set.seed(base_seed + j)` for data generation". That is correct and sufficient **when the data-generation stream and the model-fitting stream are explicitly separated**, i.e. when the code literally does

```r
set.seed(base + j)          # data draws
<draw data; store data_seed = base + j>
set.seed(base + j + 2e9)    # separate, documented fit stream; store fit_seed
<fit model>
```

If instead a simulation reuses `.Random.seed` from the last data draw to start the MCMC, the two streams are entangled and `P0-T2`'s acceptance test (bit-identical reruns across machines/worker counts) will not hold. The `set.seed(base + j)` notation in the plan is therefore to be read as shorthand for the two-seed pattern in §3.1–3.3, **not** as licence to share one stream.

### 3.5 Reproducibility acceptance (P0-T2)
A given cell, run as `run_cell(dgp, n, seed_start, seed_end)`, must be bit-identical whether it runs on a 2-core or 20-core machine, and whether `mc.cores` is 1 or 2. This holds because (a) every replication draws only from its own designated stream, and (b) each worker's model is fitted single-threaded (`n.threads = 1`, per the colab calibration note). If a model library quantises floating-point results, the columns affected (e.g. captured time) are exempted and documented; the count/indicator columns must still match exactly.

## 4. The rule (Python)

Use `numpy.random.SeedSequence().spawn()` for disjoint streams:

```python
import numpy as np
base_seed = 12345
ss = np.random.SeedSequence(base_seed)
data_child = ss.spawn(J)   # one SeedSequence per replication, for data
fit_child  = ss.spawn(J)   # distinct spawn for model fits
for j in range(J):
    rng_data = np.random.default_rng(data_child[j])
    ...                       # draw data with rng_data
    rng_fit  = np.random.default_rng(fit_child[j])
    ...                       # fit stochastic model with rng_fit
```

(The two `spawn(J)` calls are backed by different branches of the same root `SeedSequence`, so they are independent.) Log `rng_data[j].bit_generator.state['seed']`-style identifiers per row where cheap.

## 5. Pilot/confirmatory seed allocation

The plan reserves disjoint blocks so the pilot is genuinely independent:

| block | range | purpose |
|---|---|---|
| Confirmatory | `1 … J` | final inference |
| Pilot | `1_000_001 … 1_000_000 + J0` | sizing only; **discarded** |

Concretely, `base` for a confirmatory cell is the small integer `0` (or the CLI `SEED_START`) and for the pilot cell it is `SEED_START = 1_000_001`. The two blocks never overlap, so the pilot's variance estimate `σ̂_D` is not biased by having been "spent" on tuning. This is what makes the pilot genuinely independent (IJDA #5d, #8b).

**Decision:** use fixed integer seeds `0 … 999` for confirmatory cells (author decision in plan P3-T5(a)), and `1_000_001 … 1_000_000 + J0` for pilots. `as.numeric(Sys.time())` is deleted from all driver code.

## 6. Sharding and seed completeness

In the Colab execution model (§5 of the plan), a cell is sharded into contiguous seed ranges. Each shard takes its `[SEED_START, SEED_END]` as CLI arguments and appends rows to a per-shard CSV **on Drive**, one row per committed replication. After concatenation the driver asserts **no gaps, no duplicates, expected row count per cell** using the seed column; the assertion output is published (plan P3-T5(d)§7). This is why every row must carry its seeds.

## 7. What the code must do (checklist)

- [ ] Delete every `as.numeric(Sys.time())` seed.
- [ ] Add `seed_data`, `seed_fit` (and RNG-kind) columns to every output row.
- [ ] Separate data and fit streams per §3 (R) / §4 (Python).
- [ ] Use disjoint pilot and confirmatory blocks per §5.
- [ ] Workers take seeds from pre-assigned streams, never `set.seed(j)` mid-cell.
- [ ] Model fits run single-threaded (`nthread = 1` / `n.threads = 1`).
- [ ] Shard concatenation asserts seed completeness (§6).

## 8. Fidelity to the original authors' design

We do **not** try to reproduce McJames et al.'s RNG for the DGP: their runs were likewise seeded by `as.numeric(Sys.time())` and the exact draws are unrecoverable. What we reproduce is the *stochastic model* (same model settings, tree counts, priors) so that distributional statistics (means, SEs, coverage) match their published Table 2 within Monte Carlo tolerance. Exact per-replication identity with their runs is neither possible nor required.
