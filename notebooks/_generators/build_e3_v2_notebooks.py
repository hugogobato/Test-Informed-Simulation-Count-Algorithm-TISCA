#!/usr/bin/env python3
"""Generate the E3 schema-v2 Colab shard set (plan P5.5-T1).

The v1 campaign recorded mean credible-interval width and mean coverage but not
the miss distances, so the Winkler interval score cannot be reconstructed from
the committed shards.  ``run_cell_v2.R`` adds it.  Because the addition is pure
post-processing of posterior draws and draws no random numbers, a v2 shard must
reproduce its v1 counterpart exactly on every shared column, and each generated
notebook ends by checking that against the v1 CSV published in this repository.
The rerun is therefore two things at once: the missing metric, and an
independent replication of the whole v1 campaign on different machines in a
different session.

Only the 33 confirmatory shards are generated.  Under ``ANALYSIS_PLAN.md``
AMENDMENT 1 the pilot is the first 100 seeds of the confirmatory block, applied
as a split at analysis time, and the four Round 0 pilot notebooks are marked
``superseded`` in the v1 manifest.  Regenerating them would rerun 200
replications that no analysis reads.

The environment cells (R install, bundle restore, MVBCF compile) are taken from
``build_e3_notebooks.py`` unchanged and only the driver-download cell is
replaced, so the v2 sessions differ from v1 in the driver alone.

    python notebooks/_generators/build_e3_v2_notebooks.py \\
      --bundle-folder-url '<drive folder>' --bundle-sha256 '<64 hex>'
"""

from __future__ import annotations

import argparse
import csv
import json
import textwrap
from pathlib import Path

from build_e3_notebooks import (
    N100_MINUTES,
    N500_MINUTES,
    SPEEDUP,
    _bundle_args,
    _ranges,
    _sanity,
    _shared_cells,
    _validated_bundle,
    code,
    md,
    notebook,
)

SCHEMA_VERSION = "e3-v2.0"
V1_RAW_BASE = (
    "https://raw.githubusercontent.com/hugogobato/"
    "Test-Informed-Simulation-Count-Algorithm-TISCA/main/notebooks/E3_shards/"
)


def build_shards_v2():
    """The v1 confirmatory grid, re-labelled.  Seed ranges are unchanged.

    The seed ranges must stay identical to v1 or the parity check has nothing
    to compare against; this function deliberately mirrors ``build_shards()``
    rather than re-deriving the partition.
    """
    rows = []

    def add_cell(dgp, n, total, shard_size, shard_offset):
        minutes = N500_MINUTES if n == 500 else N100_MINUTES
        if n == 100 and total == 1000:
            ranges = [(0, 332), (333, 665), (666, 999)]
        else:
            ranges = _ranges(total, shard_size)
        for local_id, (cli_start, cli_end) in enumerate(ranges, start=shard_offset):
            count = cli_end - cli_start + 1
            seed_label = f"{cli_start:03d}-{cli_end:03d}"
            v1_stem = (f"E3_DGP{dgp}_n{n}_confirmatory_shard{local_id:02d}_"
                       f"seeds{seed_label}")
            stem = (f"E3v2_DGP{dgp}_n{n}_confirmatory_shard{local_id:02d}_"
                    f"seeds{seed_label}")
            rows.append({
                "round": "Round1",
                "schema_version": SCHEMA_VERSION,
                "dgp": dgp,
                "n": n,
                "mode": "confirmatory",
                "seed_start": cli_start,
                "seed_end": cli_end,
                "cli_seed_start": cli_start,
                "cli_seed_end": cli_end,
                "replication_count": count,
                "projected_hours": f"{count * minutes / SPEEDUP / 60.0:.3f}",
                "notebook_filename": stem + ".ipynb",
                "output_filename": stem + ".csv",
                "v1_output_filename": v1_stem + ".csv",
            })

    add_cell(1, 500, 1000, 100, 1)
    add_cell(2, 500, 1000, 100, 1)
    add_cell(3, 500, 1000, 100, 1)
    add_cell(1, 100, 1000, 333, 1)

    for key in {(r["dgp"], r["n"]) for r in rows}:
        cell = sorted((r for r in rows if (r["dgp"], r["n"]) == key),
                      key=lambda r: r["seed_start"])
        assert cell[0]["seed_start"] == 0, key
        for left, right in zip(cell, cell[1:]):
            assert left["seed_end"] + 1 == right["seed_start"], key
        assert sum(r["replication_count"] for r in cell) == 1000, key
        assert cell[-1]["seed_end"] == 999, key
    assert len(rows) == 33, len(rows)
    assert len({r["notebook_filename"] for r in rows}) == len(rows)
    assert len({r["output_filename"] for r in rows}) == len(rows)
    for index, row in enumerate(rows):
        row["account_slot"] = f"account{index // 3 + 1:02d}"
        row["session_slot"] = f"session{index % 3 + 1}"
    return rows


def _v2_driver_cell():
    """Replaces the v1 driver-download cell.

    The guard list is the v2 contract: if GitHub main is still serving a driver
    without the interval-score code, the notebook must stop here rather than
    spend six hours producing another v1 shard under a v2 filename.
    """
    return textwrap.dedent(f"""\
        import os, urllib.request

        RUNCELL_URL = (
            "https://raw.githubusercontent.com/hugogobato/"
            "Test-Informed-Simulation-Count-Algorithm-TISCA/main/"
            "experiments/E3_mvbcf_casestudy/run_cell_v2.R"
        )
        MVBCF_CPP_URL = (
            "https://raw.githubusercontent.com/Nathan-McJames/MVBCF_Paper/"
            "main/MVBCF_Code.cpp"
        )
        os.makedirs("/content/e3", exist_ok=True)
        urllib.request.urlretrieve(RUNCELL_URL, "/content/e3/run_cell_v2.R")
        # The upstream C++ is downloaded at runtime and is never committed here.
        urllib.request.urlretrieve(MVBCF_CPP_URL, "/content/e3/MVBCF_Code.cpp")
        assert os.path.getsize("/content/e3/run_cell_v2.R") > 1000
        assert os.path.getsize("/content/e3/MVBCF_Code.cpp") > 10000
        with open("/content/e3/run_cell_v2.R") as f:
            run_cell_source = f.read()
        required_fixes = [
            # v1 invariants: these must still be exactly as they were, because
            # the parity check below assumes the fitting code is unchanged.
            "nthread = nthread_global",
            "num_threads = nthread_global",
            "acquired <- dir.create(lk",
            "num_gfr = 0",
            "sigma2_leaf_init = 1^2 / n_tree_mu",
            "sigma2_leaf_init = 0.375^2 / n_tree_tau",
            'propensity_covariate = "prognostic"',
            "sample_sigma2_leaf = FALSE",
            # v2 additions (P5.5-T1).
            "E3_SCHEMA_VERSION",
            '"{SCHEMA_VERSION}"',
            "interval_score_mat <- function",
            "interval_score_vec <- function",
            "fill_interval_scores <- function",
            "width_audit_max_abs_dev",
        ]
        missing_fixes = [item for item in required_fixes if item not in run_cell_source]
        assert not missing_fixes, (
            "GitHub main is serving a driver without the schema-{SCHEMA_VERSION} "
            "interval-score code. Commit and push run_cell_v2.R before running "
            f"this notebook; missing: {{missing_fixes}}"
        )
        print("downloaded run_cell_v2.R (schema {SCHEMA_VERSION}) and upstream MVBCF_Code.cpp")
        """)


def _shared_cells_v2(bundle_source, bundle_sha, row, repo_root):
    cells = _shared_cells(bundle_source, bundle_sha, row, repo_root)
    targets = [i for i, c in enumerate(cells)
               if "RUNCELL_URL" in "".join(c["source"])]
    assert len(targets) == 1, (
        "expected exactly one driver-download cell in the shared v1 preamble, "
        f"found {len(targets)}; build_e3_notebooks.py has changed shape"
    )
    cells[targets[0]] = {"cell_type": "code", "execution_count": None,
                         "metadata": {}, "outputs": [],
                         "source": _v2_driver_cell().splitlines(keepends=True)}
    cells[0]["source"] = textwrap.dedent(f"""\
        # E3 schema-{SCHEMA_VERSION} shard (rerun for the interval score)

        Pre-filled for **DGP{row['dgp']}**, training **n={row['n']}**, mode
        **confirmatory**, seeds **{row['seed_start']}..{row['seed_end']}**. Run
        cells from top to bottom; there is nothing to edit.

        This is the same cell, the same seeds and the same models as the v1 run.
        The driver is `run_cell_v2.R`, which adds the Winkler interval score at
        the 50% and 95% levels (unit-level CATE and ATE, per model and outcome)
        plus its non-coverage penalty component. Everything else is byte-
        identical v1 code, and the added quantities consume no random numbers.

        Two products, therefore. First, the interval score, which is the proper
        scoring rule the model confidence set needs for the uncertainty
        comparison and which answers the reviewer's point that coverage has a
        target rather than a direction. Second, and free: the final cell checks
        this shard against the published v1 CSV
        `{row['v1_output_filename']}`, seed by seed and column by column. The
        RNG state hashes and every v1 metric must match. A different machine,
        a different session and a different day reproducing the numbers exactly
        is the replication evidence the original study, which seeded from
        `as.numeric(Sys.time())`, could not provide.
        """).splitlines(keepends=True)
    return cells


def _run_cells_v2(row):
    """The v1 run cells with v2 filenames and the v2 driver name."""
    from build_e3_notebooks import _run_cells

    cells = _run_cells(row)
    for cell in cells:
        source = "".join(cell["source"])
        if "run_cell.R" in source:
            # Global, so the invocation, the error message and the RNG-stream
            # comment all name the driver this notebook actually runs. v1 cells
            # never contain "run_cell_v2.R", so this cannot double-substitute.
            cell["source"] = source.replace(
                "run_cell.R", "run_cell_v2.R").splitlines(keepends=True)
    return cells


def _parity_cells(row):
    cells = []
    md(cells, textwrap.dedent("""\
    ## v1 parity check (this is the replication evidence)

    Downloads the v1 CSV for these exact seeds and compares it with what this
    session just produced. Three groups of columns are treated differently.

    **RNG state and model seeds** must be identical strings. These are the
    L'Ecuyer substream states and the integer seeds drawn from them, so a match
    proves both runs sampled from the same streams in the same order.

    **v1 metrics** must be identical to floating-point noise. They are computed
    by unchanged code from those streams.

    **Provenance and timing** (`hostname`, `git_sha`, `session_hash`,
    `fit_seconds_*`, `replication_seconds`) are *expected* to differ: this is a
    different machine on a different day. They are reported, never asserted.

    The verdict is printed and written to a JSON file, and it is not fatal. A
    long fit should not be discarded by a Colab cell; whether a deviation
    matters is an analysis decision, made once against all 33 shards.
    """))
    code(cells, textwrap.dedent(f"""\
        import csv, json, math, os, urllib.request

        V1_URL = {V1_RAW_BASE + row['v1_output_filename']!r}
        V1_LOCAL = "/content/TISCA_E3/v1_" + {row['v1_output_filename']!r}
        PARITY_JSON = OUTPUT_CSV.replace(".csv", "_parity.json")

        STREAM_EXACT = ["rng_kind", "seed_data_hash", "seed_fit_hash",
                        "seed_cell_master", "seq_phase", "n", "dgp",
                        "model_seed_propensity", "model_seed_mvbcf",
                        "model_seed_bcf1", "model_seed_bcf2",
                        "model_seed_bart1", "model_seed_bart2",
                        "model_seed_mvbart"]
        PROVENANCE = ["hostname", "git_sha", "session_hash", "error_message",
                      "replication_seconds", "fit_seconds_mvbcf",
                      "fit_seconds_bcf1", "fit_seconds_bcf2",
                      "fit_seconds_bart1", "fit_seconds_bart2",
                      "fit_seconds_mvbart"]

        urllib.request.urlretrieve(V1_URL, V1_LOCAL)
        with open(V1_LOCAL, newline="") as f:
            v1_rows = {{int(r["seed"]): r for r in csv.DictReader(f)}}
        with open(OUTPUT_CSV, newline="") as f:
            v2_reader = csv.DictReader(f)
            v2_fields = v2_reader.fieldnames
            v2_rows = {{int(r["seed"]): r for r in v2_reader}}

        v1_fields = list(next(iter(v1_rows.values())).keys())
        new_columns = [c for c in v2_fields if c not in v1_fields]
        dropped_columns = [c for c in v1_fields if c not in v2_fields]
        shared_seeds = sorted(set(v1_rows) & set(v2_rows))

        # v2 must be a strict superset of the v1 schema, in v1's order.
        prefix_ok = v2_fields[:len(v1_fields)] == v1_fields

        stream_mismatch = {{}}
        for column in STREAM_EXACT:
            if column not in v1_fields or column not in v2_fields:
                continue
            bad = [s for s in shared_seeds
                   if v1_rows[s][column] != v2_rows[s][column]]
            if bad:
                stream_mismatch[column] = len(bad)

        compare = [c for c in v1_fields
                   if c in v2_fields and c not in STREAM_EXACT
                   and c not in PROVENANCE and c != "seed"]
        worst = {{"column": None, "seed": None, "abs": 0.0, "rel": 0.0}}
        exact = 0
        nonnumeric = []
        for column in compare:
            for seed in shared_seeds:
                a, b = v1_rows[seed][column], v2_rows[seed][column]
                if a == b:
                    exact += 1
                    continue
                try:
                    x, y = float(a), float(b)
                except (TypeError, ValueError):
                    nonnumeric.append((column, seed))
                    continue
                d = abs(x - y)
                r = d / max(abs(x), 1e-12)
                if r > worst["rel"]:
                    worst = {{"column": column, "seed": seed, "abs": d, "rel": r}}

        total = len(compare) * len(shared_seeds)
        if stream_mismatch or nonnumeric or not prefix_ok:
            verdict = "MISMATCH"
        elif exact == total:
            verdict = "BIT_IDENTICAL"
        elif worst["rel"] < 1e-12:
            verdict = "IDENTICAL_TO_FP_NOISE"
        elif worst["rel"] < 1e-6:
            verdict = "AGREES_LOOSELY"
        else:
            verdict = "MISMATCH"

        report = {{
            "shard": os.path.basename(OUTPUT_CSV),
            "v1_shard": {row['v1_output_filename']!r},
            "schema_version": {SCHEMA_VERSION!r},
            "verdict": verdict,
            "seeds_compared": len(shared_seeds),
            "v1_only_seeds": sorted(set(v1_rows) - set(v2_rows)),
            "v2_only_seeds": sorted(set(v2_rows) - set(v1_rows)),
            "v1_schema_is_prefix_of_v2": prefix_ok,
            "new_columns": new_columns,
            "dropped_columns": dropped_columns,
            "compared_cells": total,
            "bit_identical_cells": exact,
            "stream_mismatch": stream_mismatch,
            "nonnumeric_mismatch": nonnumeric[:20],
            "worst_relative_deviation": worst,
        }}
        with open(PARITY_JSON, "w") as f:
            json.dump(report, f, indent=2)

        print("v1 parity verdict:", verdict)
        print("seeds compared:", len(shared_seeds),
              "| cells compared:", total,
              "| bit-identical:", exact)
        print("new v2 columns:", len(new_columns))
        print("v1 schema is a prefix of v2:", prefix_ok)
        if dropped_columns:
            print("WARNING dropped v1 columns:", dropped_columns)
        if stream_mismatch:
            print("WARNING RNG/seed columns differ:", stream_mismatch)
        if worst["column"]:
            print("worst deviation: %s seed %s abs=%.3g rel=%.3g"
                  % (worst["column"], worst["seed"], worst["abs"], worst["rel"]))
        try:
            from google.colab import files
            files.download(PARITY_JSON)
            print("Downloaded:", PARITY_JSON)
        except Exception as e:
            print("(Not on Colab / download skipped):", e)
        """))
    md(cells, textwrap.dedent("""\
    ## Interval-score sanity summary

    A cheap look at what the rerun bought, before any of it reaches an
    analysis. The width audit compares the interval score's own width component
    against the v1 `wid` columns; the driver already fails a replication whose
    disagreement exceeds `1e-8` times the width scale, so this should print
    zero. The score identity `is = wid + pen` is checked here in the collected
    file rather than trusted from the driver.
    """))
    code(cells, textwrap.dedent("""\
        import csv, statistics

        with open(OUTPUT_CSV, newline="") as f:
            rows = list(csv.DictReader(f))

        audit = [float(r["width_audit_max_abs_dev"]) for r in rows]
        print("width audit, max over replications: %.3g" % max(audit))

        worst_identity = 0.0
        for model in ("mvbcf", "bcf", "bart", "mvbart"):
            for level in ("50", "95"):
                for k in ("1", "2"):
                    for r in rows:
                        score = float(r[f"{model}_is{level}{k}"])
                        parts = (float(r[f"{model}_wid{level}{k}"])
                                 + float(r[f"{model}_pen{level}{k}"]))
                        worst_identity = max(
                            worst_identity,
                            abs(score - parts) / max(abs(score), 1e-12))
        print("max relative violation of is == wid + pen: %.3g" % worst_identity)

        print()
        print("%-8s %-6s %-10s %-10s %-10s %-10s"
              % ("model", "level", "IS(mean)", "width", "penalty", "coverage"))
        for model in ("mvbcf", "bcf", "bart", "mvbart"):
            for level in ("50", "95"):
                for k in ("1", "2"):
                    def mean_of(prefix):
                        return statistics.mean(
                            float(r[f"{model}_{prefix}{level}{k}"]) for r in rows)
                    print("%-8s %-6s %-10.4f %-10.4f %-10.4f %-10.4f"
                          % (f"{model}.y{k}", level + "%", mean_of("is"),
                             mean_of("wid"), mean_of("pen"), mean_of("cov")))
        print()
        print("Lower interval score is better. It is NOT the same ranking as "
              "coverage: a model can cover well by being wide and still score "
              "badly, which is the point of recording it.")
        """))
    return cells


def write_outputs_v2(repo_root, rows, bundle_source, bundle_sha):
    notebook_dir = repo_root / "notebooks" / "E3_shards_v2"
    notebook_dir.mkdir(parents=True, exist_ok=True)

    for row in rows:
        cells = _shared_cells_v2(bundle_source, bundle_sha, row, repo_root)
        cells.extend(_run_cells_v2(row))
        cells.extend(_parity_cells(row))
        nb = notebook(cells)
        _sanity(nb, row["notebook_filename"])
        source = json.dumps(nb)
        assert "run_cell.R" not in source.replace("run_cell_v2.R", ""), \
            f"{row['notebook_filename']} still references the v1 driver"
        assert row["output_filename"] in source
        assert row["v1_output_filename"] in source
        (notebook_dir / row["notebook_filename"]).write_text(
            json.dumps(nb, indent=1) + "\n")

    table_path = (repo_root / "experiments" / "E3_mvbcf_casestudy"
                  / "shard_table_v2.csv")
    with table_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    total_hours = sum(float(r["projected_hours"]) for r in rows)
    print(f"wrote {len(rows)} schema-{SCHEMA_VERSION} shard notebooks to {notebook_dir}")
    print(f"wrote {table_path}")
    print(f"projected compute: {total_hours:.1f} session-hours across "
          f"{len({r['account_slot'] for r in rows})} accounts")


def main():
    parser = argparse.ArgumentParser()
    _bundle_args(parser)
    args = parser.parse_args()
    bundle_url, bundle_sha = _validated_bundle(args)
    repo_root = Path(__file__).resolve().parents[2]
    write_outputs_v2(repo_root, build_shards_v2(), bundle_url, bundle_sha)


if __name__ == "__main__":
    main()
