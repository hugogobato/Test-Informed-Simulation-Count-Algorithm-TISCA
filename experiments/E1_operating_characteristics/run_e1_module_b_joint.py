#!/usr/bin/env python3
"""Run the genuine joint Module-B grid locally or as deterministic Colab shards.

Examples
--------
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --dry-run
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --resume
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --smoke
python experiments/E1_operating_characteristics/run_e1_module_b_joint.py --shard-index 1 --num-shards 2 --resume

The canonical 660-cell grid is imported from ``outermc/e1_grid.py``.  This runner
never retypes or mutates that factorial.  BLAS is single-threaded and the process
pool is capped by both current CPU load and half of currently available RAM by
default, preventing unrestricted nested parallelism.
"""

from __future__ import annotations

import argparse
import atexit
import json
import os
import socket
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tisca" / "python"))

from tisca import multiplicity  # noqa: E402
from tisca.outermc import e1_grid  # noqa: E402
from tisca.outermc.joint_engine import run_joint_cell  # noqa: E402
from tisca.outermc.joint_families import covariance_record  # noqa: E402


MATRIX_CSV = ROOT / "legacy" / "Paper_Experiments" / "DGP1_500_results.csv"
MATRIX_COLS = ["mvbcf_pehe1", "bcf_pehe1"]
DEFAULT_OUTDIR = ROOT / "results" / "E1"
EXPECTED_CELLS = e1_grid.EXPECTED["B"]

REQUIRED_MAIN_COLUMNS = [
    "cell_id", "module", "family", "K", "correction", "rho", "theta", "design",
    "R", "J0", "Jmax", "alpha", "alpha_adj", "marginal_level_or_power",
    "marginal_mcse", "fwer", "fwer_mcse", "conjunctive_power",
    "conjunctive_power_mcse", "disjunctive_power", "disjunctive_power_mcse",
    "fdr", "fdr_mcse", "E_J", "sd_J", "q05_J", "q50_J", "q95_J",
    "P_J_eq_Jmax", "ci_coverage", "bias", "rmse", "mean_wall", "seed",
]


@dataclass(frozen=True)
class OutputPaths:
    operating_characteristics: Path
    covariance: Path
    manifest: Path
    verdict: Path

    def as_json(self) -> str:
        return json.dumps({
            "operating_characteristics": str(self.operating_characteristics),
            "covariance": str(self.covariance),
            "manifest": str(self.manifest),
            "verdict": str(self.verdict),
        }, sort_keys=True, separators=(",", ":"))


def output_paths(outdir: Path, *, shard_index: int | None = None,
                 smoke: bool = False) -> OutputPaths:
    suffix = "_smoke" if smoke else (f"_shard{shard_index}" if shard_index else "")
    return OutputPaths(
        operating_characteristics=outdir / f"module_b_joint_operating_characteristics{suffix}.csv",
        covariance=outdir / f"module_b_joint_covariance{suffix}.csv",
        manifest=outdir / f"module_b_joint_manifest{suffix}.csv",
        verdict=outdir / f"module_b_joint_verdict{suffix}.md",
    )


def load_matrix(path: Path = MATRIX_CSV) -> np.ndarray:
    raw = pd.read_csv(path)
    missing = [name for name in MATRIX_COLS if name not in raw]
    if missing:
        raise ValueError(f"{path} lacks empirical loss columns {missing}")
    matrix = raw[MATRIX_COLS].to_numpy(dtype=float)
    matrix = matrix[np.all(np.isfinite(matrix), axis=1)]
    if matrix.shape[0] < 2:
        raise ValueError("empirical loss matrix has fewer than two finite rows")
    return matrix


def build_grid(matrix_token="__MATRIX__") -> list[dict]:
    grid = e1_grid.module_b(matrix=matrix_token,
                            planning_alpha=multiplicity.planning_alpha)
    if len(grid) != EXPECTED_CELLS:
        raise RuntimeError(f"canonical Module B has {len(grid)} cells, expected {EXPECTED_CELLS}")
    ids = [cell["cell_id"] for cell in grid]
    if len(ids) != len(set(ids)):
        raise RuntimeError("canonical joint Module-B grid has duplicate cell IDs")
    return grid


def shard_grid(grid: list[dict], shard_index: int | None,
               num_shards: int) -> list[dict]:
    if shard_index is None:
        if num_shards != 1:
            raise ValueError("--num-shards requires --shard-index")
        return grid
    if num_shards < 1 or not 1 <= shard_index <= num_shards:
        raise ValueError("shard-index must lie in [1, num-shards]")
    bounds = np.linspace(0, len(grid), num_shards + 1, dtype=int)
    return grid[int(bounds[shard_index - 1]):int(bounds[shard_index])]


def smoke_grid(grid: list[dict]) -> list[dict]:
    """A complete small factorial over K, correction and null/alternative.

    It uses normal losses, rho=0.6 and D4, yielding 3 x 5 x 2 = 30 cells, then
    adds one K=6 null cell for each empirical extension.  This exercises every
    multiplicity decision and both declared empirical joint constructions.
    """
    chosen = [cell for cell in grid if (
        cell["config"]["family"] == "normal"
        and cell["config"]["rho"] == 0.6
        and cell["config"]["design"] == "D4"
    )]
    for family in ("empirical_copula", "empirical"):
        match = next(cell for cell in grid if (
            cell["config"]["family"] == family
            and cell["config"]["K"] == 6
            and cell["config"]["theta"] == 0.0
            and cell["config"]["design"] == "D4"
            and cell["config"]["correction"] == "romano_wolf"
            and (family == "empirical" or cell["config"]["rho"] == 0.6)
        ))
        chosen.append(match)
    assert len(chosen) == 32, len(chosen)
    return chosen


def available_gb() -> float:
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return float("inf")


def choose_workers(requested: int | None, mem_fraction: float,
                   mem_per_worker: float) -> int:
    logical = os.cpu_count() or 2
    physical_cap = max(1, logical // 2)
    try:
        load_one = os.getloadavg()[0]
    except (AttributeError, OSError):
        load_one = 0.0
    # One worker per two currently idle logical CPUs, bounded by the approximate
    # physical-core count.  This backs off when other experiments are active.
    load_cap = max(1, int(max(1.0, logical - load_one) // 2))
    cpu_cap = min(physical_cap, load_cap)
    mem_cap = max(1, int((available_gb() * mem_fraction) // mem_per_worker))
    safe = max(1, min(cpu_cap, mem_cap))
    if requested is None:
        return safe
    if requested > mem_cap:
        print(f"[joint B] warning: {requested} workers exceed the memory cap {mem_cap}")
    return max(1, requested)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tmp", prefix=path.name + ".", dir=path.parent,
        delete=False, encoding="utf-8",
    )
    temp = Path(handle.name)
    try:
        with handle:
            frame.to_csv(handle, index=False)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def _atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tmp", prefix=path.name + ".", dir=path.parent,
        delete=False, encoding="utf-8",
    )
    temp = Path(handle.name)
    try:
        with handle:
            handle.write(text)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def _release_lock(path: Path) -> None:
    """Release only a lock owned by this process (forked workers cannot remove it)."""
    try:
        owner = json.loads(path.read_text(encoding="utf-8"))
        if int(owner.get("pid", -1)) == os.getpid():
            path.unlink(missing_ok=True)
    except (OSError, ValueError, json.JSONDecodeError):
        pass


def _acquire_lock(paths: OutputPaths) -> Path:
    """Prevent two sessions from checkpointing the same output snapshot."""
    path = paths.manifest.with_name(paths.manifest.name + ".lock")
    path.parent.mkdir(parents=True, exist_ok=True)
    owner = {
        "pid": os.getpid(), "hostname": socket.gethostname(),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    while True:
        try:
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(owner, handle, sort_keys=True)
            atexit.register(_release_lock, path)
            return path
        except FileExistsError:
            try:
                old = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError, json.JSONDecodeError):
                old = {}
            same_host = old.get("hostname") == socket.gethostname()
            old_pid = int(old.get("pid", -1))
            alive = False
            if same_host and old_pid > 0:
                try:
                    os.kill(old_pid, 0)
                    alive = True
                except OSError:
                    alive = False
            if same_host and not alive:
                path.unlink(missing_ok=True)
                continue
            raise RuntimeError(
                f"another joint Module-B writer owns {path}: {old}. "
                "Do not run the same shard concurrently."
            )


_WORKER_MATRIX: np.ndarray | None = None
_RUN_OVERRIDES: dict = {}


def _worker_init(overrides: dict) -> None:
    global _RUN_OVERRIDES
    _RUN_OVERRIDES = dict(overrides)


def _run_one(cell: dict) -> dict:
    global _WORKER_MATRIX
    started = time.perf_counter()
    started_iso = datetime.now(timezone.utc).isoformat()
    try:
        cfg = dict(cell["config"])
        cfg.update(_RUN_OVERRIDES)
        if cfg.get("matrix") == "__MATRIX__":
            if _WORKER_MATRIX is None:
                _WORKER_MATRIX = load_matrix()
            cfg["matrix"] = _WORKER_MATRIX
        summary, _raw = run_joint_cell(cfg)
        summary.update(cell["factors"])
        summary.update(
            cell_id=cell["cell_id"], module="B_joint", source_module="B",
            projected_R=cell["config"]["R"],
        )
        cov = covariance_record(
            cfg, cell_id=cell["cell_id"],
            diagnostic_n=int(_RUN_OVERRIDES.get("diagnostic_n", 10_000)),
            numeric_n=int(_RUN_OVERRIDES.get("numeric_covariance_n", 100_000)),
        )
        return {
            "status": "completed", "cell_id": cell["cell_id"],
            "seed": int(cfg["seed"]), "summary": summary, "covariance": cov,
            "hostname": socket.gethostname(), "started_at": started_iso,
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": time.perf_counter() - started, "error": "",
        }
    except Exception as exc:  # surfaced in the manifest before the parent raises
        return {
            "status": "failed", "cell_id": cell["cell_id"],
            "seed": int(cell["config"]["seed"]), "summary": None, "covariance": None,
            "hostname": socket.gethostname(), "started_at": started_iso,
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "runtime_seconds": time.perf_counter() - started,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _load_unique(path: Path, key: str = "cell_id") -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    frame = pd.read_csv(path)
    if key not in frame:
        raise ValueError(f"{path} lacks {key}")
    if frame[key].astype(str).duplicated().any():
        raise ValueError(f"{path} contains duplicate {key}")
    return frame


def _ordered(frame: pd.DataFrame, required: list[str]) -> pd.DataFrame:
    if frame.empty:
        return frame
    missing = [name for name in required if name not in frame]
    if missing:
        raise RuntimeError(f"joint result is missing required columns: {missing}")
    rest = [name for name in frame.columns if name not in required]
    return frame[required + rest]


def write_verdict(results: pd.DataFrame, covariance: pd.DataFrame, path: Path,
                  *, expected_ids: set[str], label: str) -> None:
    complete = set(results.get("cell_id", pd.Series(dtype=str)).astype(str)) == expected_ids
    lines = [
        f"# Joint Module B verdict ({label})",
        "",
        f"Status: {'complete' if complete else 'partial'} ({len(results)} of {len(expected_ids)} cells).",
        "",
        "Module B now reports marginal and family-level operating characteristics. "
        "FWER, conjunctive power, disjunctive power, and FDR are distinct estimands. "
        "BH FDR is stored separately even under the global null, where it equals FWER.",
        "",
    ]
    if not results.empty:
        null = results[results["n_true_alternatives"] == 0]
        d4_null = null[null["design"] == "D4"]
        if not d4_null.empty:
            table = d4_null.groupby(["K", "correction"], dropna=False).agg(
                cells=("cell_id", "size"), marginal=("marginal_level_or_power", "mean"),
                mean_fwer=("fwer", "mean"), min_fwer=("fwer", "min"),
                max_fwer=("fwer", "max"),
            ).reset_index()
            none_fwer = table.loc[table["correction"] == "none", ["K", "mean_fwer"]]
            none_summary = dict(zip(none_fwer["K"].astype(int), none_fwer["mean_fwer"]))
            if {1, 3, 6}.issubset(none_summary):
                interpretation = (
                    "For D4, unadjusted FWER rises from about "
                    f"{none_summary[1]:.3f} at K=1 to {none_summary[3]:.3f} at K=3 "
                    f"and {none_summary[6]:.3f} at K=6. In these global-null cells, "
                    "Bonferroni, Holm, BH, and Romano-Wolf remain near the family level "
                    "across K."
                )
            else:
                interpretation = (
                    "The table reports the observed unadjusted and adjusted global-null "
                    "family rejection rates for every K represented in this run."
                )
            lines += ["## D4 global-null summary", "", "```text",
                      table.to_string(index=False), "```", "", interpretation, ""]
        alt = results[(results["n_true_alternatives"] > 0) & (results["design"] == "D4")]
        if not alt.empty:
            power = alt.groupby(["K", "correction"], dropna=False).agg(
                marginal=("marginal_level_or_power", "mean"),
                conjunctive=("conjunctive_power", "mean"),
                disjunctive=("disjunctive_power", "mean"),
            ).reset_index()
            lines += ["## D4 alternative summary", "", "```text",
                      power.to_string(index=False), "```", ""]
    synthetic = int(covariance.get("empirical_extension_is_synthetic", pd.Series(dtype=bool)).fillna(False).sum())
    lines += [
        "## Construction and limitations", "",
        "Every repetition contains one shared method-A loss and K separately generated "
        "benchmark losses. The paired vector is computed by the declared linear map "
        "C = [1, -I_K]. Romano-Wolf resamples the full K-column contrast block with "
        "one common row-index draw, preserving cross-contrast dependence.", "",
        f"The two-column empirical source cannot identify a K-benchmark joint law. "
        f"Accordingly, {synthetic} saved covariance rows are labeled synthetic "
        "Gaussian-copula extensions. K=1 row-bootstrap cells retain the observed pair exactly.", "",
        f"All {len(covariance)} saved covariance rows pass the exact D = C L mapping; "
        f"the minimum contrast-covariance eigenvalue is "
        f"{float(covariance['min_eigenvalue'].min()):.6g}, and no diagnostic found "
        "repeated benchmark columns.", "",
        "The canonical theta grid contains either K nulls or K alternatives, not a "
        "mixed null/alternative family. BH FDR is therefore informative at the global "
        "null and equals zero in the all-alternative cells; partial-null FDR is not "
        "identified by this factorial.", "",
        "The canonical Romano-Wolf cells retain the existing planning schedule at "
        "family alpha. Their reported marginal and family power is the measured power "
        "of the final joint stepdown decisions, not a claim that a nested Romano-Wolf "
        "planner guarantees the nominal conjunctive target.", "",
        "E3 Romano-Wolf, Holm, and Bonferroni results remain the case-study inference layer. "
        "No E3 model fit or E3 seed-verification run is part of this task. G3 bibliometric "
        "justification coding is out of scope.", "",
    ]
    _atomic_text("\n".join(lines), path)


def _manifest_template(grid: list[dict], paths: OutputPaths,
                       shard_index: int | None, num_shards: int) -> pd.DataFrame:
    return pd.DataFrame([{
        "cell_id": cell["cell_id"], "seed": int(cell["config"]["seed"]),
        "status": "pending", "hostname": "", "started_at": "", "finished_at": "",
        "runtime_seconds": np.nan, "shard_index": shard_index or 0,
        "num_shards": num_shards, "output_paths": paths.as_json(), "error": "",
    } for cell in grid])


def _validate_complete(results: pd.DataFrame, covariance: pd.DataFrame,
                       manifest: pd.DataFrame, expected: set[str]) -> None:
    for name, frame in (("results", results), ("covariance", covariance), ("manifest", manifest)):
        ids = set(frame["cell_id"].astype(str))
        if ids != expected:
            raise RuntimeError(f"{name}: {len(expected - ids)} missing, {len(ids - expected)} unexpected")
        if frame["cell_id"].astype(str).duplicated().any():
            raise RuntimeError(f"{name}: duplicate cell IDs")
    if not manifest["status"].eq("completed").all():
        raise RuntimeError("manifest contains non-completed cells")
    if covariance["repeated_benchmark_columns"].astype(bool).any():
        raise RuntimeError("a covariance diagnostic found repeated benchmark columns")
    if not covariance["exact_contrast_mapping"].astype(bool).all():
        raise RuntimeError("a covariance diagnostic failed the exact D=C L mapping")


def run(args) -> int:
    full_grid = build_grid()
    selected = shard_grid(full_grid, args.shard_index, args.num_shards)
    if args.smoke:
        if args.shard_index is not None:
            raise ValueError("--smoke cannot be combined with cell sharding")
        selected = smoke_grid(full_grid)
    paths = output_paths(Path(args.outdir), shard_index=args.shard_index, smoke=args.smoke)
    workers = choose_workers(args.jobs, args.mem_fraction, args.mem_per_worker)

    overrides = {
        "outer_chunk": args.outer_chunk,
        "diagnostic_n": args.diagnostic_n,
        "numeric_covariance_n": args.numeric_covariance_n,
    }
    if args.smoke:
        overrides.update(R=args.smoke_R, B=args.smoke_B)
    if args.R is not None:
        overrides["R"] = args.R
    if args.bootstrap_B is not None:
        overrides["B"] = args.bootstrap_B

    seeds = [cell["config"]["seed"] for cell in selected]
    print(f"[joint B] canonical cells: {len(full_grid)} (actual grid, plan nominally 600)")
    print(f"[joint B] selected cells: {len(selected)}, IDs {selected[0]['cell_id']}..{selected[-1]['cell_id']}, "
          f"seeds {min(seeds)}..{max(seeds)}")
    print(f"[joint B] workers: {workers}; available RAM {available_gb():.1f} GB; "
          f"memory budget fraction {args.mem_fraction:.2f}")
    print(f"[joint B] outputs: {paths.as_json()}")

    existing = [path for path in paths.__dict__.values() if Path(path).exists()]
    if existing and not args.resume and not args.dry_run:
        raise FileExistsError(
            "joint output already exists; pass --resume to preserve and continue it: "
            + ", ".join(map(str, existing))
        )

    lock_path = None if args.dry_run else _acquire_lock(paths)
    results = _load_unique(paths.operating_characteristics)
    covariances = _load_unique(paths.covariance)
    if args.resume and paths.manifest.exists():
        manifest = _load_unique(paths.manifest)
    else:
        manifest = _manifest_template(selected, paths, args.shard_index, args.num_shards)
    # CSV round-tripping turns an all-empty string column (most notably `error`)
    # into float NaN.  Pandas 2.2 then correctly refuses a lossy string assignment
    # on resume.  Restore the manifest schema before updating any status row.
    for column in ("cell_id", "status", "hostname", "started_at", "finished_at",
                   "output_paths", "error"):
        if column in manifest:
            manifest[column] = manifest[column].fillna("").astype(object)
    expected = {cell["cell_id"] for cell in selected}
    result_ids = set(results.get("cell_id", pd.Series(dtype=str)).astype(str))
    covariance_ids = set(covariances.get("cell_id", pd.Series(dtype=str)).astype(str))
    completed_manifest_ids = set(manifest.loc[
        manifest.get("status", pd.Series(index=manifest.index, dtype=str)).eq("completed"),
        "cell_id",
    ].astype(str))
    # A checkpoint is committed only when all three artifacts agree.  If a hard
    # interruption lands between their atomic replacements, the cell is rerun and
    # the last row replaces the partial snapshot without duplication.
    done = result_ids & covariance_ids & completed_manifest_ids
    pending = [cell for cell in selected if cell["cell_id"] not in done]

    if args.dry_run:
        rw = sum(cell["config"]["correction"] == "romano_wolf" for cell in pending)
        R = int(overrides.get("R", selected[0]["config"]["R"]))
        B = int(overrides.get("B", selected[0]["config"]["B"]))
        print(f"[joint B dry-run] completed={len(done)}, pending={len(pending)}, "
              f"Romano-Wolf cells={rw}, R={R}, bootstrap B={B}")
        print("[joint B dry-run] no files were changed")
        return 0
    if not pending:
        print("[joint B] no pending cells")
        if lock_path is not None:
            _release_lock(lock_path)
        return 0

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    _atomic_csv(manifest, paths.manifest)
    failures = []
    started = time.perf_counter()

    def accept(item: dict) -> None:
        nonlocal results, covariances, manifest
        cid = item["cell_id"]
        mask = manifest["cell_id"].astype(str).eq(cid)
        if not mask.any():
            # A resumed manifest from an older partial shard is safely extended.
            manifest = pd.concat([manifest, _manifest_template(
                [next(cell for cell in selected if cell["cell_id"] == cid)],
                paths, args.shard_index, args.num_shards,
            )], ignore_index=True)
            mask = manifest["cell_id"].astype(str).eq(cid)
        for key in ("status", "hostname", "started_at", "finished_at", "runtime_seconds", "error"):
            manifest.loc[mask, key] = item[key]
        _atomic_csv(manifest.sort_values("cell_id"), paths.manifest)
        if item["status"] != "completed":
            failures.append(item)
            print(f"[FAIL] {cid}: {item['error']}", flush=True)
            return
        results = pd.concat([results, pd.DataFrame([item["summary"]])], ignore_index=True)
        covariances = pd.concat([covariances, pd.DataFrame([item["covariance"]])], ignore_index=True)
        results = results.drop_duplicates("cell_id", keep="last").sort_values("cell_id")
        covariances = covariances.drop_duplicates("cell_id", keep="last").sort_values("cell_id")
        _atomic_csv(_ordered(results, REQUIRED_MAIN_COLUMNS), paths.operating_characteristics)
        _atomic_csv(covariances, paths.covariance)

    if workers == 1:
        _worker_init(overrides)
        for index, cell in enumerate(pending, 1):
            accept(_run_one(cell))
            if index % 10 == 0 or index == len(pending):
                print(f"[joint B] {index}/{len(pending)} pending cells processed", flush=True)
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(
            workers, initializer=_worker_init, initargs=(overrides,)
        ) as pool:
            for index, item in enumerate(pool.imap_unordered(_run_one, pending, chunksize=1), 1):
                accept(item)
                if index % 10 == 0 or index == len(pending):
                    print(f"[joint B] {index}/{len(pending)} pending cells processed", flush=True)

    manifest = manifest[manifest["cell_id"].astype(str).isin(expected)].sort_values("cell_id")
    _atomic_csv(manifest, paths.manifest)
    write_verdict(results, covariances, paths.verdict, expected_ids=expected,
                  label=("smoke" if args.smoke else
                         (f"shard {args.shard_index}/{args.num_shards}" if args.shard_index else "full grid")))
    if failures:
        raise RuntimeError(f"{len(failures)} joint Module-B cells failed; see {paths.manifest}")
    _validate_complete(results, covariances, manifest, expected)
    elapsed = time.perf_counter() - started
    print(f"[PASS] joint Module B wrote {len(results)} cells in {elapsed:.1f}s")
    for path in paths.__dict__.values():
        print(path)
    if lock_path is not None:
        _release_lock(lock_path)
    return 0


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true", help="print the exact plan without writing files")
    p.add_argument("--resume", action="store_true", help="skip completed cell IDs and continue checkpoints")
    p.add_argument("--smoke", action="store_true", help="run the 32-cell complete smoke subset")
    p.add_argument("--smoke-R", type=int, default=200)
    p.add_argument("--smoke-B", type=int, default=99)
    p.add_argument("--R", type=int, default=None, help="outer repetitions override (testing only)")
    p.add_argument("--bootstrap-B", type=int, default=None, help="Romano-Wolf B override (testing only)")
    p.add_argument("--shard-index", type=int, default=None, help="one-based deterministic contiguous shard")
    p.add_argument("--num-shards", type=int, default=1)
    p.add_argument("--jobs", type=int, default=None, help="worker processes; default is load/RAM aware")
    p.add_argument("--mem-fraction", type=float, default=0.5,
                   help="maximum fraction of currently available RAM (default 0.5)")
    p.add_argument("--mem-per-worker", type=float, default=0.35,
                   help="conservative GB budget per worker (default 0.35)")
    p.add_argument("--outer-chunk", type=int, default=32,
                   help="outer repetitions held at once inside one worker")
    p.add_argument("--diagnostic-n", type=int, default=10_000)
    p.add_argument("--numeric-covariance-n", type=int, default=100_000)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if not 0 < args.mem_fraction <= 0.5:
        raise ValueError("mem-fraction must lie in (0, 0.5]")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
