#!/usr/bin/env python3
"""Audit the schema-e3-v2.0 rerun against the v1 campaign (plan P5.5-T1).

Each v2 shard was produced by ``run_cell_v2.R``, which differs from
``run_cell.R`` only by adding the Winkler interval score.  The addition draws no
random numbers, so a v2 shard must reproduce its v1 counterpart on every shared
column, for every seed.  This script checks that claim across all 33
confirmatory shards and writes the audit log that the revision cites as
reproducibility evidence.

It checks four things, in order of what a failure would mean.

1.  **Schema.**  The v1 header must be a prefix of the v2 header.  A dropped or
    reordered v1 column means the rerun is not a superset and the two campaigns
    cannot be pooled or compared column-wise.

2.  **RNG state.**  ``seed_data_hash``, ``seed_fit_hash`` and every
    ``model_seed_*`` must match exactly.  These are L'Ecuyer substream states
    and the integers drawn from them.  If they match, both runs sampled the same
    streams in the same order, which is the strongest available statement that
    the seed protocol is machine-independent.

3.  **Metrics.**  Every v1 metric must agree.  Timing and provenance columns are
    excluded because they are expected to differ across machines.

4.  **The new columns.**  Finite, and internally consistent: the interval score
    must equal width plus penalty, and the score's width component must match
    the separately recorded v1 width (``width_audit_max_abs_dev``).

Usage:

    python experiments/E3_mvbcf_casestudy/check_v1_v2_parity.py \\
        --v1-dir notebooks/E3_shards \\
        --v2-dir /path/to/downloaded/E3_v2_csvs \\
        --report results/E3/v1_v2_parity.json

Exit status is non-zero if any hard check fails, so it can gate a rerun.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

MODELS = ("mvbcf", "bcf", "bart", "mvbart")
LEVELS = ("50", "95")
OUTCOMES = ("1", "2")

# Must be identical strings: RNG state and everything derived from it.
STREAM_EXACT = [
    "rng_kind", "seed_data_hash", "seed_fit_hash", "seed_cell_master",
    "seq_phase", "n", "dgp", "model_seed_propensity", "model_seed_mvbcf",
    "model_seed_bcf1", "model_seed_bcf2", "model_seed_bart1",
    "model_seed_bart2", "model_seed_mvbart",
]

# Expected to differ: this is a different machine on a different day.
PROVENANCE = [
    "hostname", "git_sha", "session_hash", "error_message",
    "replication_seconds", "fit_seconds_mvbcf", "fit_seconds_bcf1",
    "fit_seconds_bcf2", "fit_seconds_bart1", "fit_seconds_bart2",
    "fit_seconds_mvbart",
]

# Relative agreement at or below this is floating-point noise, not a difference
# in the computation. Identical code on identical inputs should give 0.0; the
# margin exists only for a different BLAS or CPU in the Colab pool.
METRIC_TOLERANCE = 1e-12


def read_shard(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return reader.fieldnames, {int(r["seed"]): r for r in rows}


def compare_shard(v1_path, v2_path):
    v1_fields, v1_rows = read_shard(v1_path)
    v2_fields, v2_rows = read_shard(v2_path)

    result = {
        "v1_shard": v1_path.name,
        "v2_shard": v2_path.name,
        "failures": [],
        "warnings": [],
    }

    # (1) schema
    prefix_ok = v2_fields[:len(v1_fields)] == v1_fields
    result["v1_schema_is_prefix_of_v2"] = prefix_ok
    result["new_columns"] = [c for c in v2_fields if c not in v1_fields]
    if not prefix_ok:
        result["failures"].append(
            "v1 header is not a prefix of the v2 header (columns dropped or "
            "reordered)")

    missing = sorted(set(v1_rows) - set(v2_rows))
    extra = sorted(set(v2_rows) - set(v1_rows))
    result["seeds_v1"] = len(v1_rows)
    result["seeds_v2"] = len(v2_rows)
    result["missing_seeds"] = missing[:20]
    result["unexpected_seeds"] = extra[:20]
    if missing:
        result["failures"].append(f"{len(missing)} v1 seeds absent from v2")
    if extra:
        result["failures"].append(f"{len(extra)} v2 seeds not present in v1")

    unconverged = [s for s, r in v2_rows.items() if r.get("converged_flag") != "1"]
    result["unconverged"] = len(unconverged)
    if unconverged:
        result["failures"].append(
            f"{len(unconverged)} v2 replications have converged_flag != 1")

    seeds = sorted(set(v1_rows) & set(v2_rows))

    # (2) RNG state
    stream_mismatch = {}
    for column in STREAM_EXACT:
        if column not in v1_fields or column not in v2_fields:
            continue
        bad = [s for s in seeds if v1_rows[s][column] != v2_rows[s][column]]
        if bad:
            stream_mismatch[column] = {"count": len(bad), "first_seeds": bad[:5]}
    result["stream_mismatch"] = stream_mismatch
    if stream_mismatch:
        result["failures"].append(
            "RNG state or model seeds differ; the two runs did not sample the "
            "same streams")

    # (3) v1 metrics
    compared = [c for c in v1_fields
                if c in v2_fields and c not in STREAM_EXACT
                and c not in PROVENANCE and c != "seed"]
    exact = 0
    total = 0
    worst = {"column": None, "seed": None, "abs": 0.0, "rel": 0.0}
    unparseable = []
    for column in compared:
        for seed in seeds:
            total += 1
            a, b = v1_rows[seed][column], v2_rows[seed][column]
            if a == b:
                exact += 1
                continue
            try:
                x, y = float(a), float(b)
            except (TypeError, ValueError):
                unparseable.append({"column": column, "seed": seed,
                                    "v1": a, "v2": b})
                continue
            deviation = abs(x - y)
            relative = deviation / max(abs(x), 1e-12)
            if relative > worst["rel"]:
                worst = {"column": column, "seed": seed,
                         "abs": deviation, "rel": relative}
    result["metric_columns_compared"] = len(compared)
    result["cells_compared"] = total
    result["cells_bit_identical"] = exact
    result["worst_relative_deviation"] = worst
    result["unparseable_mismatch"] = unparseable[:20]
    if unparseable:
        result["failures"].append(
            f"{len(unparseable)} shared cells differ and are not numeric")
    if worst["rel"] > METRIC_TOLERANCE:
        result["failures"].append(
            f"metric {worst['column']} differs by {worst['rel']:.3g} relative "
            f"at seed {worst['seed']}, above the {METRIC_TOLERANCE:.0g} tolerance")

    if exact == total and not result["failures"]:
        result["verdict"] = "BIT_IDENTICAL"
    elif not result["failures"]:
        result["verdict"] = "IDENTICAL_TO_FP_NOISE"
    else:
        result["verdict"] = "MISMATCH"

    # (4) the new columns
    audit_max = 0.0
    identity_max = 0.0
    nonfinite = []
    for seed in seeds:
        row = v2_rows[seed]
        if "width_audit_max_abs_dev" in row:
            audit_max = max(audit_max, float(row["width_audit_max_abs_dev"]))
        for model in MODELS:
            for level in LEVELS:
                for k in OUTCOMES:
                    keys = [f"{model}_is{level}{k}", f"{model}_pen{level}{k}",
                            f"{model}_wid{level}{k}", f"{model}_ate_is{level}{k}"]
                    try:
                        score, penalty, width, ate = (float(row[key]) for key in keys)
                    except (KeyError, TypeError, ValueError):
                        nonfinite.append({"seed": seed, "keys": keys})
                        continue
                    if not all(math.isfinite(v) for v in (score, penalty, width, ate)):
                        nonfinite.append({"seed": seed, "keys": keys})
                        continue
                    if penalty < 0:
                        result["failures"].append(
                            f"negative interval-score penalty at seed {seed}, "
                            f"{model} level {level} outcome {k}")
                    identity_max = max(
                        identity_max,
                        abs(score - (width + penalty)) / max(abs(score), 1e-12))
    result["width_audit_max_abs_dev"] = audit_max
    result["max_relative_violation_of_score_identity"] = identity_max
    result["nonfinite_new_columns"] = len(nonfinite)
    if nonfinite:
        result["failures"].append(
            f"{len(nonfinite)} interval-score entries are missing or non-finite")
    if identity_max > 1e-12:
        result["failures"].append(
            f"interval score != width + penalty (max relative {identity_max:.3g})")

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v1-dir", type=Path, default=Path("notebooks/E3_shards"))
    parser.add_argument("--v2-dir", type=Path, required=True,
                        help="directory holding the downloaded E3v2_*.csv files")
    parser.add_argument("--manifest", type=Path,
                        default=Path("experiments/E3_mvbcf_casestudy/shard_table_v2.csv"))
    parser.add_argument("--report", type=Path,
                        default=Path("results/E3/v1_v2_parity.json"))
    parser.add_argument("--allow-missing-shards", action="store_true",
                        help="audit whatever is present instead of requiring all 33")
    args = parser.parse_args()

    with args.manifest.open(newline="") as f:
        manifest = list(csv.DictReader(f))

    shards = []
    absent = []
    for row in manifest:
        v2_path = args.v2_dir / row["output_filename"]
        v1_path = args.v1_dir / row["v1_output_filename"]
        if not v1_path.exists():
            raise SystemExit(f"missing v1 shard: {v1_path}")
        if not v2_path.exists():
            absent.append(row["output_filename"])
            continue
        shards.append(compare_shard(v1_path, v2_path))

    if absent and not args.allow_missing_shards:
        raise SystemExit(
            f"{len(absent)} v2 shards have not been collected: {absent[:5]}"
            " (pass --allow-missing-shards to audit a partial rerun)")

    verdicts = {}
    for shard in shards:
        verdicts[shard["verdict"]] = verdicts.get(shard["verdict"], 0) + 1
    failed = [s for s in shards if s["failures"]]
    worst = max((s["worst_relative_deviation"] for s in shards),
                key=lambda w: w["rel"], default={"rel": 0.0})

    report = {
        "schema_version": "e3-v2.0",
        "shards_audited": len(shards),
        "shards_absent": absent,
        "verdicts": verdicts,
        "shards_with_failures": len(failed),
        "worst_relative_metric_deviation": worst,
        "max_width_audit": max((s["width_audit_max_abs_dev"] for s in shards),
                               default=0.0),
        "total_cells_compared": sum(s["cells_compared"] for s in shards),
        "total_cells_bit_identical": sum(s["cells_bit_identical"] for s in shards),
        "shards": shards,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")

    print(f"audited {len(shards)} shards -> {args.report}")
    for verdict, count in sorted(verdicts.items()):
        print(f"  {verdict}: {count}")
    print(f"  cells compared: {report['total_cells_compared']}, "
          f"bit-identical: {report['total_cells_bit_identical']}")
    print(f"  worst relative metric deviation: {worst['rel']:.3g}"
          + (f" ({worst['column']}, seed {worst['seed']})" if worst.get("column") else ""))
    print(f"  max width audit: {report['max_width_audit']:.3g}")
    for shard in failed:
        print(f"  FAIL {shard['v2_shard']}: {shard['failures'][0]}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
