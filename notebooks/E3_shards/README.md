# E3 shard launch procedure

The directory contains 37 pre-filled E3 notebooks and the committed manifest
`experiments/E3_mvbcf_casestudy/shard_table.csv`. There are 4 Round 0 pilot
notebooks, one for each cell (DGP1/2/3 at `n=500` and DGP1 at `n=100`), with 50
pilot replications each. There are 33 Round 1 confirmatory notebooks: 10 each
for DGP1/2/3 at `n=500`, and 3 for DGP1 at `n=100`, covering 1000 replications
per cell. The 33 Round 1 notebooks fit the 11 accounts times 3 concurrent
sessions arrangement. Round 0 uses four of those sessions and is run before
Round 1.

Before uploading, build the P0-T4 R bundle and regenerate the notebooks with the
shared Drive folder and the SHA256 recorded in `env/tisca_rlib.sha256`:

```bash
python notebooks/_generators/build_e3_notebooks.py \
  --bundle-folder-url 'https://drive.google.com/drive/folders/1w3quuskj25CBOFCGG0mTRGUHcufPpdb3?usp=sharing' \
  --bundle-sha256 '12d223bc0fcef624c1ff4cc35c5d7ecc1b1f9b05aa84ecd9d9e4a5a3382bae3c'
```

The generated notebooks download the folder with `gdown`, locate both
`tisca_rlib.tar.gz` and `tisca_rlib.sha256`, verify the downloaded checksum
against the embedded digest and the tarball bytes, and then restore the R
library. No cell needs editing after that generation step. Each notebook
downloads the upstream `MVBCF_Code.cpp` at run time, sets `TISCA_MVBCF_CPP` to
that downloaded copy, restores the R bundle, uses `mc.cores = 2`, and runs
model fits with one thread. The driver appends every completed replication,
including its `replication_seconds`, to a CSV under `/content/TISCA_E3`; the
final cell includes the `google.colab.files.download` fallback.

## Run order

1. Upload and run `E3_DGP1_n500_pilot_shard01_seeds1000001-1000050.ipynb`
   alone first. This is the recommended Round 0 smoke test and carries the
   P3-T5(e) `stochtree::bcf` calibration gate.
2. Inspect the gate verdict and record it in
   `experiments/E3_mvbcf_casestudy/CALIBRATION.md`. The pass bands and decision
   rule are also frozen in `ANALYSIS_PLAN.md` §9. Do not start a Round 1
   notebook unless this gate passes.
3. Run the other three Round 0 notebooks, completing the pilots for DGP2 n=500,
   DGP3 n=500, and DGP1 n=100. Round 0 must be complete before Round 1.
4. After the gate and Round 0 completion, distribute the 33 Round 1 notebooks,
   one per account/session, according to the `account_slot` and `session_slot`
   columns in `shard_table.csv`. The deterministic assignment is the manifest
   order: the first three confirmatory rows go to account01 sessions 1, 2, 3;
   the next three go to account02 sessions 1, 2, 3; continue this pattern
   through account11. The ranges are contiguous and cover `0..999` for each
   confirmatory cell. The four pilot rows have `account_slot=round0` and can
   use any available session after the first pilot and gate check.

These notebooks intentionally assume that each Colab session runs to
completion. If a session dies, its `/content` CSV is lost and that shard must
be rerun; the final download is the durable copy. Do not edit the seed range or
rename the checkpoint CSV during a successful run. A duplicate or out-of-range
checkpoint is still rejected loudly before collection.

## Collect the results

Copy the downloaded per-shard CSVs into one local directory, preserving the
filenames from `shard_table.csv`, then run:

```bash
python experiments/E3_mvbcf_casestudy/collect_shards.py \
  --drive-dir /path/to/copied/E3_csvs
```

The collector asserts every manifest file exists, every shard has the expected
row count, seeds have no gaps or duplicates, metadata and schemas agree, and
each cell is complete. It reports the `converged_flag == 0` count and writes
mode-specific files such as
`results/E3/DGP1_n500_confirmatory_replications.csv`, together with combined
pilot plus confirmatory files such as `results/E3/DGP1_n500_replications.csv`.

The local machine does not currently have `bartCause`, `skewBART`, or
`scoringRules` installed, so `run_cell.R` cannot be executed end to end here.
Its structure, seed protocol, schema, checkpointing, and dimension guards are
verified locally; the first real execution is the DGP1 n=500 Round 0 pilot.
