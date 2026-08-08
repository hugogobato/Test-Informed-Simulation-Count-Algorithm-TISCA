#!/usr/bin/env python3
"""Collect E3 shard checkpoints and apply the AMENDMENT 1 pilot/confirmatory split.

Usage, after copying the per-shard CSVs from the participating Google Drives
into one local directory (``notebooks/E3_shards`` already holds them here)::

    python collect_shards.py --drive-dir notebooks/E3_shards

It fails before writing outputs if a checkpoint is incomplete, duplicated,
assigned to the wrong cell, or inconsistent with ``shard_table.csv``.

**The seed blocks this script assembles are not the ones it splits on.**
``run_cell.R`` is a data *generator*: under ``--mode confirmatory`` it draws from
L'Ecuyer master 1 and labels the rows ``0..999``, and under ``--mode pilot`` from
master 2, labelled ``1_000_001..``. That machinery is unchanged and correct, and
nothing here re-runs it.

What ``ANALYSIS_PLAN.md`` AMENDMENT 1 changed is purely a *partition of rows
already collected*: seeds ``0..J0-1`` become the pilot that sizes ``J``, and
seeds ``J0..999`` the confirmatory set that is tested. The two blocks are drawn
from disjoint L'Ecuyer substreams of the same master, so they are independent
samples, which is the only property Algorithm 1 requires of a pilot. The
reserved master-2 block is therefore *not used*: it exists for one cell only and
was produced by a superseded ``run_cell.R`` (see ``CALIBRATION.md`` D1), and its
shard-table rows now carry ``status = superseded`` so this script neither reads
nor requires them.

Outputs under ``results/E3``:

    DGP{d}_n{n}_replications.csv               all 1000 rows, phase-labelled
    DGP{d}_n{n}_amended_pilot_replications.csv seeds 0..J0-1  (sizes J, then discarded)
    DGP{d}_n{n}_amended_confirmatory_replications.csv seeds J0..999 (tested)
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
    parser.add_argument("--pilot-size", type=int, default=100,
                        help="J0: how many of the confirmatory seeds form the "
                             "amended pilot (ANALYSIS_PLAN.md AMENDMENT 1 declares "
                             "100; 50 and 25 are the pre-declared sensitivity rows)")
    args = parser.parse_args()

    table_fields, manifest = read_csv(args.table)
    required = {"dgp", "n", "mode", "seed_start", "seed_end",
                "replication_count", "output_filename"}
    missing_fields = required - set(table_fields)
    if missing_fields:
        raise RuntimeError(f"shard table missing columns: {sorted(missing_fields)}")

    # Superseded rows are skipped outright: their files are absent for three of the
    # four cells and invalid for the fourth, and demanding them made this script
    # unrunnable after AMENDMENT 1. A table without a `status` column is treated as
    # all-active so older checkouts keep working.
    skipped = [s for s in manifest if s.get("status") == "superseded"]
    manifest = [s for s in manifest if s.get("status") != "superseded"]
    for spec in skipped:
        print(f"skip (superseded): {spec['output_filename']}")

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

    # AMENDMENT 1: partition the collected confirmatory block into the pilot that
    # sizes J and the confirmatory set that is tested. This is a row split of data
    # already in hand -- no cell is re-run -- and it is the step the plan calls the
    # "reanalysis".
    J0 = int(args.pilot_size)
    if J0 < 2:
        raise RuntimeError("--pilot-size must be at least 2 (a pilot sd needs df >= 1)")
    for dgp, n in sorted(all_fields):
        combined = []
        fields = all_fields[(dgp, n)]
        for mode in ("pilot", "confirmatory"):
            combined.extend(by_cell.get((dgp, n, mode), []))
        combined.sort(key=lambda row: int(row["seed"]))
        if not combined:
            raise RuntimeError(f"no rows collected for DGP{dgp} n={n}")

        # `seq_phase` records how run_cell.R generated the row; `analysis_phase`
        # records how AMENDMENT 1 uses it. Keeping both is what makes the deviation
        # auditable from the data files alone.
        amended_fields = fields + ["analysis_phase"]
        pilot, conf = [], []
        for row in combined:
            seed = int(row["seed"])
            target = pilot if seed < J0 else conf
            row = dict(row, analysis_phase="pilot" if seed < J0 else "confirmatory")
            target.append(row)
        if len(pilot) != J0:
            raise RuntimeError(
                f"DGP{dgp} n={n}: amended pilot has {len(pilot)} rows, expected {J0}; "
                f"seeds 0..{J0 - 1} are not all present")

        for tag, rows_out in (("", pilot + conf),
                              ("amended_pilot_", pilot),
                              ("amended_confirmatory_", conf)):
            output = args.results_dir / f"DGP{dgp}_n{n}_{tag}replications.csv"
            write_rows(output, amended_fields, rows_out)
        failures = sum(row.get("converged_flag") == "0" for row in combined)
        print(f"(combined) DGP{dgp} n={n}: {len(combined)} rows "
              f"-> pilot {len(pilot)} (seeds 0..{J0 - 1}) + confirmatory {len(conf)} "
              f"(seeds {J0}..{max(int(r['seed']) for r in conf)}), "
              f"converged_flag failures={failures}")

    print(f"E3 shard collection: PASS (J0 = {J0}, pilot rows discarded from inference)")


if __name__ == "__main__":
    main()
