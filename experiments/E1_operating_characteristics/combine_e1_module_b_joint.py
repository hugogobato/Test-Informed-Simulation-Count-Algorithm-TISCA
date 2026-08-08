#!/usr/bin/env python3
"""Combine deterministic joint Module-B shard directories into final E1 outputs.

Each positional input is a directory containing one shard's three CSV files from
``run_e1_module_b_joint.py --shard-index ...``.  The combiner asserts the exact
canonical 660 cell IDs and seeds, rejects every duplicate, and only then writes
the four unsuffixed publication outputs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from run_e1_module_b_joint import (  # noqa: E402
    REQUIRED_MAIN_COLUMNS,
    _atomic_csv,
    _ordered,
    _validate_complete,
    build_grid,
    output_paths,
    write_verdict,
)


def _one(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"expected one {pattern!r} in {directory}, found {matches}")
    return matches[0]


def _read_all(inputs: list[Path], pattern: str) -> tuple[pd.DataFrame, list[Path]]:
    paths = [_one(directory, pattern) for directory in inputs]
    frames = [pd.read_csv(path) for path in paths]
    return pd.concat(frames, ignore_index=True), paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="joint shard output directories")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results" / "E1")
    args = parser.parse_args(argv)
    for directory in args.inputs:
        if not directory.is_dir():
            raise FileNotFoundError(directory)

    results, result_sources = _read_all(
        args.inputs, "module_b_joint_operating_characteristics_shard*.csv"
    )
    covariance, covariance_sources = _read_all(
        args.inputs, "module_b_joint_covariance_shard*.csv"
    )
    manifest, manifest_sources = _read_all(
        args.inputs, "module_b_joint_manifest_shard*.csv"
    )
    for name, frame in (("results", results), ("covariance", covariance), ("manifest", manifest)):
        if frame["cell_id"].astype(str).duplicated().any():
            duplicates = frame.loc[frame["cell_id"].astype(str).duplicated(), "cell_id"].tolist()
            raise RuntimeError(f"duplicate {name} cell IDs: {duplicates[:5]}")

    grid = build_grid()
    expected = {cell["cell_id"] for cell in grid}
    expected_seed = {cell["cell_id"]: int(cell["config"]["seed"]) for cell in grid}
    for frame_name, frame in (("results", results), ("covariance", covariance), ("manifest", manifest)):
        actual_seed = dict(zip(frame["cell_id"].astype(str), frame["seed"].astype(int)))
        bad = [cid for cid, seed in actual_seed.items() if expected_seed.get(cid) != seed]
        if bad:
            raise RuntimeError(f"{frame_name} has non-canonical seeds for {bad[:5]}")

    final = output_paths(args.outdir)
    manifest = manifest.copy()
    manifest["output_paths"] = final.as_json()
    manifest["combined_result_sources"] = ";".join(map(str, result_sources))
    manifest["combined_covariance_sources"] = ";".join(map(str, covariance_sources))
    manifest["combined_manifest_sources"] = ";".join(map(str, manifest_sources))
    results = results.sort_values("cell_id").reset_index(drop=True)
    covariance = covariance.sort_values("cell_id").reset_index(drop=True)
    manifest = manifest.sort_values("cell_id").reset_index(drop=True)
    _validate_complete(results, covariance, manifest, expected)

    _atomic_csv(_ordered(results, REQUIRED_MAIN_COLUMNS), final.operating_characteristics)
    _atomic_csv(covariance, final.covariance)
    _atomic_csv(manifest, final.manifest)
    write_verdict(results, covariance, final.verdict, expected_ids=expected,
                  label="combined two-shard full grid")
    print(f"[PASS] combined {len(results)} unique canonical cells")
    for path in final.__dict__.values():
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
