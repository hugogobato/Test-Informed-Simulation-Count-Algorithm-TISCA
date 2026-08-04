#!/usr/bin/env python3
"""Collect E3 Drive checkpoints and enforce the shard-table contract.

Usage, after copying the per-shard CSVs from the participating Google Drives
into one local directory::

    python collect_shards.py --drive-dir /path/to/E3_drive_exports

The script writes one file per DGP, n, and mode, plus a combined pilot plus
confirmatory file per DGP/n cell, under ``results/E3``.  It fails before writing
outputs if a checkpoint is incomplete, duplicated, assigned to the wrong cell,
or inconsistent with ``shard_table.csv``.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_csv(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        raise RuntimeError(f"empty checkpoint: {path}")
    return fields, rows


def as_int(row, key, path):
    try:
        return int(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(f"invalid {key} in {path}: {row}") from exc


def write_rows(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    here = Path(__file__).resolve().parent
    repo = here.parents[1]
    parser.add_argument("--drive-dir", type=Path, required=True,
                        help="Directory containing copied per-shard Drive CSVs")
    parser.add_argument("--table", type=Path,
                        default=here / "shard_table.csv")
    parser.add_argument("--results-dir", type=Path,
                        default=repo / "results" / "E3")
    args = parser.parse_args()

    table_fields, manifest = read_csv(args.table)
    required = {"dgp", "n", "mode", "seed_start", "seed_end",
                "replication_count", "output_filename"}
    missing_fields = required - set(table_fields)
    if missing_fields:
        raise RuntimeError(f"shard table missing columns: {sorted(missing_fields)}")

    by_cell = defaultdict(list)
    loaded = {}
    for spec in manifest:
        dgp = as_int(spec, "dgp", args.table)
        n = as_int(spec, "n", args.table)
        start = as_int(spec, "seed_start", args.table)
        end = as_int(spec, "seed_end", args.table)
        expected_count = as_int(spec, "replication_count", args.table)
        mode = spec["mode"]
        key = (dgp, n, mode)
        path = args.drive_dir / spec["output_filename"]
        if not path.exists():
            raise RuntimeError(f"missing Drive checkpoint: {path}")
        fields, rows = read_csv(path)
        seeds = [as_int(row, "seed", path) for row in rows]
        expected = set(range(start, end + 1))
        if len(rows) != expected_count:
            raise RuntimeError(
                f"{path.name}: expected {expected_count} rows, found {len(rows)}")
        if len(seeds) != len(set(seeds)):
            raise RuntimeError(f"{path.name}: duplicate seed rows")
        if set(seeds) != expected:
            raise RuntimeError(
                f"{path.name}: seed range mismatch, missing="
                f"{sorted(expected - set(seeds))[:10]}, extra="
                f"{sorted(set(seeds) - expected)[:10]}")
        for row in rows:
            if as_int(row, "dgp", path) != dgp or as_int(row, "n", path) != n:
                raise RuntimeError(f"{path.name}: row metadata disagrees with shard table")
            if row.get("seq_phase") not in (None, "", mode):
                raise RuntimeError(f"{path.name}: seq_phase disagrees with shard table")
        loaded[spec["output_filename"]] = (fields, rows)
        by_cell[key].extend(rows)

    # Validate the table ranges as well as the data files.
    manifest_cells = defaultdict(list)
    for spec in manifest:
        manifest_cells[(int(spec["dgp"]), int(spec["n"]), spec["mode"])].append(spec)

    all_fields = {}
    for key, specs in manifest_cells.items():
        dgp, n, mode = key
        specs.sort(key=lambda r: int(r["seed_start"]))
        for left, right in zip(specs, specs[1:]):
            if int(left["seed_end"]) + 1 != int(right["seed_start"]):
                raise RuntimeError(f"manifest has a seed gap or overlap in {key}")
        rows = by_cell[key]
        seeds = [int(row["seed"]) for row in rows]
        expected = set(range(int(specs[0]["seed_start"]), int(specs[-1]["seed_end"]) + 1))
        if len(rows) != len(expected) or len(seeds) != len(set(seeds)) or set(seeds) != expected:
            raise RuntimeError(f"cell completeness failed for {key}")
        fields = loaded[specs[0]["output_filename"]][0]
        for spec in specs:
            if loaded[spec["output_filename"]][0] != fields:
                raise RuntimeError(f"CSV schemas differ within {key}")
        all_fields[(dgp, n)] = fields
        rows.sort(key=lambda row: int(row["seed"]))
        output = args.results_dir / f"DGP{dgp}_n{n}_{mode}_replications.csv"
        write_rows(output, fields, rows)
        failures = sum(row.get("converged_flag") == "0" for row in rows)
        print(f"{key}: {len(rows)} rows, converged_flag failures={failures}, wrote {output}")

    # Combined cell files are convenient for downstream analysis.  The pilot
    # rows remain explicitly labelled and retain their disjoint emitted seeds.
    for dgp, n in sorted(all_fields):
        combined = []
        fields = all_fields[(dgp, n)]
        for mode in ("pilot", "confirmatory"):
            combined.extend(by_cell.get((dgp, n, mode), []))
        combined.sort(key=lambda row: int(row["seed"]))
        output = args.results_dir / f"DGP{dgp}_n{n}_replications.csv"
        write_rows(output, fields, combined)
        failures = sum(row.get("converged_flag") == "0" for row in combined)
        print(f"(combined) DGP{dgp} n={n}: {len(combined)} rows, "
              f"converged_flag failures={failures}, wrote {output}")

    print("E3 shard collection: PASS")


if __name__ == "__main__":
    main()
