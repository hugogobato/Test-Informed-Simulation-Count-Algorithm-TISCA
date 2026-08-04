# P0-T2 — Seed and RNG protocol

**Status:** specified; **not yet accepted** — acceptance requires the §3.5 four-way identity test to pass on the real E3 driver  
**Revision:** 2026-08-04 (rev. 2: §3.1/§3.3 stream construction corrected and verified in R 4.3.3; §3.6 model-owned RNG added)  
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
Every replication `j` consumes **two** disjoint RNG streams:

| stream | used for |
|---|---|
| `data_stream_j` | the covariates `X`, the treatment `Z`, and the potential outcomes `Y(0), Y(1)` |
| `fit_stream_j` | every stochastic model fit on that dataset (MCMC, warm start, Gibbs) |

The two are generated as consecutive L'Ecuyer-CMRG substreams from one cell master seed (§3.3), so replication `j`'s data draws can never advance or contaminate its model draws, and no two replications share state. The pair `(2j, 2j+1)` indexing means a shard can reconstruct exactly the streams of its own seed range without running any earlier replication.

**Do not derive one stream from the other by arithmetic on `.Random.seed`.** An earlier draft of this protocol used `fit_streams[[1]] <- .Random.seed + 2e9`. That is invalid and was verified to fail: under `L'Ecuyer-CMRG`, `.Random.seed` is a length-7 **integer** vector whose first element is the RNG-kind code (`10407`), so adding a constant corrupts the kind code and pushes the six state words outside the CMRG's valid moduli. In R 4.3.3 this produces `Error in nextRNGStream: invalid value of 'seed'`, and, worse, `assign(".Random.seed", <that vector>, ...)` only emits `'.Random.seed' is not an integer vector but of type 'double', so ignored` — a **warning, not an error**, so the fit stream is silently *not* set and the run is quietly irreproducible. Use §3.3 verbatim.

### 3.2 Stream identifiers in the field names
Each row of the output CSV records exactly the streams actually used:

| column | meaning |
|---|---|
| `seed` | the replication index / cell seed (`0 … 999` confirmatory, `1_000_001 …` pilot) |
| `seed_cell_master` | the cell master seed the streams were generated from |
| `seed_data_hash` / `seed_fit_hash` | a short digest of the 7-word `.Random.seed` vector actually installed for each stream (`digest::digest()` or `paste(collapse="_")`) |
| `model_seed_<m>` | the integer seed passed explicitly to model `m` (§3.6), where the model takes one |
| `rng_kind` | `"L'Ecuyer-CMRG"` (R) or `"numpy.SeedSequence"` (Python) |

Storing the seed *in every row* is the minimum reproducibility fix and matches the original authors' own pattern (`GitHub_DGP1.R` stores `seed_val`). Storing the stream digests as well is what makes a re-run verifiable rather than merely re-runnable.

### 3.3 L'Ecuyer-CMRG with parallel streams (the required construction)
Streams are generated in interleaved pairs from one cell master seed, so replication `j` owns substreams `2j` (data) and `2j+1` (fit). Verified in R 4.3.3: identical results at `mc.cores = 1` and `mc.cores = 2`, identical across re-runs, and a shard built for `[j_from, j_to]` produces byte-identical streams to the corresponding slice of a full-cell build.

```r
library(parallel)
RNGkind("L'Ecuyer-CMRG")

# Streams for replications j_from..j_to of one cell, reconstructible by any shard.
make_streams <- function(cell_master, j_from, j_to) {
  set.seed(cell_master)
  s <- .Random.seed
  data_s <- fit_s <- vector("list", j_to - j_from + 1L)
  idx <- 0L
  for (j in 0:j_to) {
    d <- s; s <- nextRNGStream(s)      # substream 2j   -> data
    f <- s; s <- nextRNGStream(s)      # substream 2j+1 -> model fits
    if (j >= j_from) { idx <- idx + 1L; data_s[[idx]] <- d; fit_s[[idx]] <- f }
  }
  list(data = data_s, fit = fit_s)
}

run_rep <- function(d_stream, f_stream, ...) {
  assign(".Random.seed", d_stream, envir = .GlobalEnv)
  # <draw X, Z, Y(0), Y(1)>
  assign(".Random.seed", f_stream, envir = .GlobalEnv)
  # <fit every model; see 3.6 for models that carry their own RNG>
}
```

`nextRNGStream` guarantees the substreams are (stochastically) disjoint. Workers consume their pre-assigned stream; they never call `set.seed(j)` themselves, because `set.seed(j)` with a fixed `j` does not give the disjoint-stream guarantee across parallel workers. Advancing the master `2 * j_to` times costs microseconds even at `j_to = 999`, so every shard rebuilds its own streams independently and no state has to be shipped between sessions.

### 3.4 The idiomatic `set.seed(base + j)` rule, and when it is acceptable
The plan text (`P0-T2`) phrases the rule as "`set.seed(base_seed + j)` for data generation". That is a *shorthand*. It is acceptable only in the weaker form

```r
RNGkind("L'Ecuyer-CMRG")
set.seed(base + j)                       # per-replication master
d <- .Random.seed                        # data substream
f <- nextRNGStream(d)                    # fit substream (NEVER d + constant)
```

which keeps the two streams separated and valid, at the cost of the cross-replication disjointness guarantee that §3.3 provides (two different `set.seed` values are not *provably* non-overlapping, they are merely overwhelmingly unlikely to overlap). **§3.3 is the required construction for E3 and E1; §3.4 is documented only so that the plan's shorthand is not misread as licence to share one stream, or to derive a stream by adding a constant to `.Random.seed` (§3.1).**

If a simulation instead reuses `.Random.seed` from the last data draw to start the MCMC, the two streams are entangled and `P0-T2`'s acceptance test will not hold.

### 3.5 Reproducibility acceptance (P0-T2)
A given cell, run as `run_cell(dgp, n, seed_start, seed_end)`, must be bit-identical whether it runs on a 2-core or 20-core machine, and whether `mc.cores` is 1 or 2. This holds because (a) every replication draws only from its own designated stream, (b) each worker's model is fitted single-threaded (`n.threads = 1`, per the colab calibration note), and (c) models with their own RNG receive an explicit seed (§3.6). If a model library quantises floating-point results, the columns affected (e.g. captured time) are exempted and documented; the count/indicator columns must still match exactly.

**The acceptance test is a script, not an intention.** `tisca/tests/test_seed_protocol.R` must run replication 37 of DGP1/n=500 four ways — `mc.cores = 1`, `mc.cores = 2`, in a shard whose range starts at 37, and in a shard whose range starts at 0 — and assert that all four rows are identical on every non-timing column. This test is what signs off P0-T2; until it exists and passes on the real driver, P0-T2 is specified but not accepted.

### 3.6 Models that carry their own RNG (the failure mode §3.3 alone does not cover)
Setting `.Random.seed` only governs draws that go through R's RNG. Several models in this study do **not**:

* `stochtree::bcf()` and `stochtree::bart()` take an explicit `random_seed` argument and use their own C++ generator when it is supplied. McJames et al.'s own script already does this (`bcf(..., random_seed = seed_val)`).
* `dbarts::bart()` / `bartCause::bartc()` accept `seed` and, with `n.threads > 1`, are not reproducible at all — which is a second, independent reason for `n.threads = 1`.
* `fast_bart()` from `MVBCF_Code.cpp` must be inspected: if it draws through `Rcpp`/R's RNG (`GetRNGstate`/`PutRNGstate`), §3.3 covers it; if it seeds its own generator, it needs an explicit seed argument like the others.

**Rule:** every model that exposes a seed argument is passed a deterministic integer derived from that replication's fit stream, and that integer is written to the output row as `model_seed_<m>`:

```r
assign(".Random.seed", f_stream, envir = .GlobalEnv)
model_seed <- sample.int(.Machine$integer.max, 1L)   # drawn FROM the fit stream, so reproducible
fit <- stochtree::bcf(..., random_seed = model_seed)
```

Verifying which of the four models actually honour their seed is a Round 0 pilot task, not an assumption: fit each model twice on the same replication and assert identical output. Record the outcome in `experiments/E3_mvbcf_casestudy/CALIBRATION.md`.

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

(The two `spawn(J)` calls are backed by different branches of the same root `SeedSequence`, so they are independent — `spawn` is stateful and advances the spawn key.) Log `rng_data[j].bit_generator.state['seed']`-style identifiers per row where cheap.

For **shard reconstruction** (E1 shards a cell block without running earlier repetitions), prefer the addressable form, which needs no prior spawning and is the exact analogue of §3.3:

```python
ss_data = np.random.SeedSequence(base_seed, spawn_key=(2 * j,))
ss_fit  = np.random.SeedSequence(base_seed, spawn_key=(2 * j + 1,))
```

## 5. Pilot/confirmatory seed allocation

The plan reserves disjoint blocks so the pilot is genuinely independent:

| block | range | purpose |
|---|---|---|
| Confirmatory | `1 … J` | final inference |
| Pilot | `1_000_001 … 1_000_000 + J0` | sizing only; **discarded** |

Concretely, `base` for a confirmatory cell is the small integer `0` (or the CLI `SEED_START`) and for the pilot cell it is `SEED_START = 1_000_001`. The two blocks never overlap, so the pilot's variance estimate `σ̂_D` is not biased by having been "spent" on tuning. This is what makes the pilot genuinely independent (IJDA #5d, #8b).

**Independence requires a different cell master, not just a different index range.** Under §3.3 the substream index is `2j`/`2j+1` from one master, so a pilot at `j = 1_000_001` and a confirmatory run at `j = 0…999` already occupy non-overlapping substreams of the same master and are independent in the L'Ecuyer sense. Use *the same* `cell_master` per cell and let the disjoint index blocks do the work; do **not** re-seed the pilot from a different master, or the shard-reconstruction property of §3.3 is lost. Advancing to substream `2_000_002` is still an O(1)-per-step loop of ~2 million `nextRNGStream` calls, which takes a few seconds — so pilot shards should build their streams once and cache them, or use a dedicated pilot master `cell_master + 1` with the offset documented in the output row. **Decision: use `cell_master_pilot = cell_master + 1` with pilot indices `0 … J0−1`, recorded in `seed_cell_master`.** The blocks are then disjoint by construction and cheap to rebuild, and the `seed` column still carries the `1_000_001 …` labels for human auditing.

**Decision:** use fixed integer seeds `0 … 999` for confirmatory cells (author decision in plan P3-T5(a)), and `1_000_001 … 1_000_000 + J0` for pilots. `as.numeric(Sys.time())` is deleted from all driver code.

## 6. Sharding and seed completeness

In the Colab execution model (§5 of the plan), a cell is sharded into contiguous seed ranges. Each shard takes its `[SEED_START, SEED_END]` as CLI arguments and appends rows to a per-shard CSV **on Drive**, one row per committed replication. After concatenation the driver asserts **no gaps, no duplicates, expected row count per cell** using the seed column; the assertion output is published (plan P3-T5(d)§7). This is why every row must carry its seeds.

## 7. What the code must do (checklist)

- [ ] Delete every `as.numeric(Sys.time())` seed.
- [ ] Add `seed`, `seed_cell_master`, `seed_data_hash`, `seed_fit_hash` and `rng_kind` columns to every output row.
- [ ] Separate data and fit streams per §3.3 (R) / §4 (Python) — **never** by adding a constant to `.Random.seed` (§3.1).
- [ ] Pass an explicit `model_seed_<m>`, drawn from the fit stream, to every model that exposes a seed argument; record it (§3.6).
- [ ] Use disjoint pilot and confirmatory blocks per §5.
- [ ] Workers take seeds from pre-assigned streams, never `set.seed(j)` mid-cell.
- [ ] Model fits run single-threaded (`nthread = 1` / `n.threads = 1`).
- [ ] Shard concatenation asserts seed completeness (§6).
- [ ] `tisca/tests/test_seed_protocol.R` exists and passes the four-way identity test (§3.5). **P0-T2 is not accepted until this runs green on the real driver.**
- [ ] Round 0 records which models actually honour their seed (§3.6) in `CALIBRATION.md`.

## 8. Fidelity to the original authors' design

We do **not** try to reproduce McJames et al.'s RNG for the DGP: their runs were likewise seeded by `as.numeric(Sys.time())` and the exact draws are unrecoverable. What we reproduce is the *stochastic model* (same model settings, tree counts, priors) so that distributional statistics (means, SEs, coverage) match their published Table 2 within Monte Carlo tolerance. Exact per-replication identity with their runs is neither possible nor required.
