#!/usr/bin/env python3
"""Generate the complete, pre-filled E3 Colab shard set.

The generated notebooks are ready to upload and run from top to bottom.  Each
notebook describes exactly one cell and one contiguous seed range, so there is
no operator configuration cell.  The shard manifest is written to
``experiments/E3_mvbcf_casestudy/shard_table.csv`` and is the source of truth
for collection and completeness checks.

The R library bundle source and SHA256 are generation-time arguments.  The
source may be a direct tarball URL or a public Google Drive folder containing
``tisca_rlib.tar.gz`` and ``tisca_rlib.sha256``.  Before the bundle exists,
``--allow-unbuilt-bundle`` may be used to stage notebooks with a deliberately
invalid sentinel.  Those notebooks fail loudly in their restore cell until a
real bundle is baked in by regenerating the set.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import textwrap
from pathlib import Path


UNBUILT_URL = "UNBUILT_BUNDLE_URL"
UNBUILT_SHA = "UNBUILT_BUNDLE_SHA256"
DEFAULT_BUNDLE_FOLDER_URL = (
    "https://drive.google.com/drive/folders/"
    "1w3quuskj25CBOFCGG0mTRGUHcufPpdb3?usp=sharing"
)
SPEEDUP = 1.62
N500_MINUTES = 6.0
N100_MINUTES = 2.0


def md(cells, source):
    cells.append({"cell_type": "markdown", "metadata": {},
                  "source": source.splitlines(keepends=True)})


def code(cells, source):
    cells.append({"cell_type": "code", "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": source.splitlines(keepends=True)})


def notebook(cells):
    return {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {
            "kernelspec": {"name": "python3", "display_name": "Python 3"},
            "language_info": {"name": "python"},
        },
        "cells": cells,
    }


def _ranges(total, shard_size):
    start = 0
    while start < total:
        end = min(total - 1, start + shard_size - 1)
        yield start, end
        start = end + 1


def build_shards():
    rows = []

    def add_cell(dgp, n, mode, total, shard_size, shard_offset):
        minutes = N500_MINUTES if n == 500 else N100_MINUTES
        if mode == "confirmatory" and n == 100 and total == 1000:
            ranges = [(0, 332), (333, 665), (666, 999)]
        else:
            ranges = _ranges(total, shard_size)
        for local_id, (cli_start, cli_end) in enumerate(
                ranges, start=shard_offset):
            count = cli_end - cli_start + 1
            if mode == "pilot":
                seed_start = cli_start + 1_000_001
                seed_end = cli_end + 1_000_001
                round_name = "Round0"
                seed_label = f"{seed_start}-{seed_end}"
            else:
                seed_start = cli_start
                seed_end = cli_end
                round_name = "Round1"
                seed_label = f"{seed_start:03d}-{seed_end:03d}"
            stem = (f"E3_DGP{dgp}_n{n}_{mode}_shard{local_id:02d}_"
                    f"seeds{seed_label}")
            rows.append({
                "round": round_name,
                "dgp": dgp,
                "n": n,
                "mode": mode,
                "seed_start": seed_start,
                "seed_end": seed_end,
                "cli_seed_start": cli_start,
                "cli_seed_end": cli_end,
                "replication_count": count,
                "projected_hours": f"{count * minutes / SPEEDUP / 60.0:.3f}",
                "notebook_filename": stem + ".ipynb",
                "output_filename": stem + ".csv",
            })

    # Round 0 is one ready-made notebook per pilot cell.  The pilot seed block
    # is emitted by run_cell.R as 1,000,001 ... 1,000,050.
    add_cell(1, 500, "pilot", 50, 50, 1)
    add_cell(2, 500, "pilot", 50, 50, 1)
    add_cell(3, 500, "pilot", 50, 50, 1)
    add_cell(1, 100, "pilot", 50, 50, 1)

    # Round 1 follows the plan budget: ten ~6-hour n=500 shards per cell and
    # three ~6.8-hour n=100 shards for the 1000 confirmatory replications.
    add_cell(1, 500, "confirmatory", 1000, 100, 1)
    add_cell(2, 500, "confirmatory", 1000, 100, 1)
    add_cell(3, 500, "confirmatory", 1000, 100, 1)
    add_cell(1, 100, "confirmatory", 1000, 333, 1)

    # Every cell has an independent contiguous range.  Pilot and confirmatory
    # labels are intentionally disjoint because pilots carry the +1,000,001
    # emitted seed block.
    expected = {"pilot": 50, "confirmatory": 1000}
    for key in {(r["dgp"], r["n"], r["mode"]) for r in rows}:
        cell = sorted((r for r in rows
                       if (r["dgp"], r["n"], r["mode"]) == key),
                      key=lambda r: r["seed_start"])
        assert cell[0]["seed_start"] == (1_000_001 if key[2] == "pilot" else 0)
        for left, right in zip(cell, cell[1:]):
            assert left["seed_end"] + 1 == right["seed_start"], key
        assert sum(r["replication_count"] for r in cell) == expected[key[2]], key
        assert cell[-1]["seed_end"] == (
            1_000_000 + expected[key[2]] if key[2] == "pilot"
            else expected[key[2]] - 1
        )
    assert len(rows) == 37, len(rows)
    assert len({r["notebook_filename"] for r in rows}) == len(rows)
    assert len({r["output_filename"] for r in rows}) == len(rows)
    confirmatory = [r for r in rows if r["mode"] == "confirmatory"]
    for index, row in enumerate(confirmatory):
        row["account_slot"] = f"account{index // 3 + 1:02d}"
        row["session_slot"] = f"session{index % 3 + 1}"
    for index, row in enumerate((r for r in rows if r["mode"] == "pilot"), start=1):
        row["account_slot"] = "round0"
        row["session_slot"] = f"pilot{index}"
    return rows


def _bundle_args(parser):
    parser.add_argument("--bundle-url", default=UNBUILT_URL,
                        help="Direct URL for tisca_rlib.tar.gz.")
    parser.add_argument("--bundle-folder-url", default=None,
                        help=("Public Google Drive folder containing the tarball and sha256 "
                              f"file, e.g. {DEFAULT_BUNDLE_FOLDER_URL}"))
    parser.add_argument("--bundle-sha256", default=UNBUILT_SHA,
                        help="SHA256 of the R library bundle.")
    parser.add_argument("--allow-unbuilt-bundle", action="store_true",
                        help="Stage notebooks with the unbuilt-bundle sentinel.")


def _validated_bundle(args):
    if args.bundle_folder_url:
        if args.bundle_url != UNBUILT_URL:
            raise SystemExit("Use either --bundle-url or --bundle-folder-url, not both")
        if not re.fullmatch(r"https://drive\.google\.com/drive/folders/.+", args.bundle_folder_url):
            raise SystemExit("--bundle-folder-url must be a Google Drive folder URL")
        if not re.fullmatch(r"[0-9a-fA-F]{64}", args.bundle_sha256):
            raise SystemExit("--bundle-sha256 must be a 64-character hexadecimal digest")
        return {"kind": "gdrive_folder", "value": args.bundle_folder_url}, args.bundle_sha256.lower()
    unbuilt = args.bundle_url == UNBUILT_URL or args.bundle_sha256 == UNBUILT_SHA
    if unbuilt and not args.allow_unbuilt_bundle:
        raise SystemExit(
            "The R bundle is not built. Supply --bundle-url or --bundle-folder-url "
            "plus --bundle-sha256, or explicitly pass --allow-unbuilt-bundle to "
            "stage failing notebooks."
        )
    if unbuilt:
        print("WARNING: staging with an unbuilt-bundle sentinel; every notebook "
              "will stop in its bundle assertion cell.")
        return {"kind": "sentinel", "value": UNBUILT_URL}, UNBUILT_SHA
    if not re.fullmatch(r"https?://.+", args.bundle_url):
        raise SystemExit("--bundle-url must be an http(s) URL")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", args.bundle_sha256):
        raise SystemExit("--bundle-sha256 must be a 64-character hexadecimal digest")
    return {"kind": "direct", "value": args.bundle_url}, args.bundle_sha256.lower()


def _shared_cells(bundle_source, bundle_sha, row, repo_root):
    cells = []
    md(cells, f"""\
    # E3 {row['round']} shard

    This notebook is fully pre-filled for **DGP{row['dgp']}**, training
    **n={row['n']}**, mode **{row['mode']}**, and emitted seeds
    **{row['seed_start']}..{row['seed_end']}**. Run cells from top to bottom.
    There is no configuration cell to edit. The CSV checkpoint is written to
    Google Drive after every completed replication.
    """)
    code(cells, textwrap.dedent("""\
        import os, platform, subprocess, time

        def sh(cmd):
            p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return p.stdout.strip()

        print("hostname:", platform.node() or "n/a")
        print("os:", platform.platform())
        print("nproc:", sh("nproc"))
        print("cpu:", sh("grep -m1 -E 'model name' /proc/cpuinfo") or "n/a")
        print(sh("free -g") or "RAM information unavailable")
        print("mc.cores is fixed at 2; model fits are fixed to one thread.")
        """))
    code(cells, textwrap.dedent("""\
        import subprocess
        p = subprocess.run(
            ["bash", "-lc", "apt-get -qq update >/dev/null && "
             "apt-get -qq install -y --no-install-recommends "
             "r-base r-base-dev libcurl4-openssl-dev >/dev/null 2>&1"],
            capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stdout[-1000:])
            print(p.stderr[-2000:])
            raise RuntimeError("R installation failed")
        print(subprocess.check_output(["R", "--version"], text=True).splitlines()[0])
        """))
    code(cells, textwrap.dedent("""\
        from google.colab import drive
        drive.mount("/content/drive")
        print("Drive mounted")
        """))
    if bundle_source["kind"] == "gdrive_folder":
        bundle_restore = textwrap.dedent(f"""\
            import hashlib, os, re, shutil, subprocess, sys
            from pathlib import Path

            BUNDLE_FOLDER_URL = {bundle_source['value']!r}
            BUNDLE_SHA256 = {bundle_sha!r}
            assert re.fullmatch(r"[0-9a-fA-F]{{64}}", BUNDLE_SHA256), \\
                "Bundle SHA256 is missing or malformed; regenerate the notebook."

            BUNDLE_DOWNLOAD_DIR = Path("/content/tisca_bundle_download")
            if BUNDLE_DOWNLOAD_DIR.exists():
                shutil.rmtree(BUNDLE_DOWNLOAD_DIR)
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "-q", "gdown"],
                check=True,
            )
            download = subprocess.run(
                [sys.executable, "-m", "gdown", "--folder", BUNDLE_FOLDER_URL,
                 "--output", str(BUNDLE_DOWNLOAD_DIR), "--remaining-ok"],
                capture_output=True, text=True,
            )
            print(download.stdout[-4000:])
            if download.returncode != 0:
                print(download.stderr[-4000:])
                raise RuntimeError("Google Drive bundle download failed")
            tar_candidates = sorted(BUNDLE_DOWNLOAD_DIR.rglob("tisca_rlib.tar.gz"))
            sha_candidates = sorted(BUNDLE_DOWNLOAD_DIR.rglob("tisca_rlib.sha256"))
            assert len(tar_candidates) == 1, f"expected one tarball, found {{tar_candidates}}"
            assert len(sha_candidates) == 1, f"expected one checksum file, found {{sha_candidates}}"
            published_sha = sha_candidates[0].read_text().split()[0].lower()
            assert published_sha == BUNDLE_SHA256.lower(), \\
                "published tisca_rlib.sha256 differs from the generated notebook"
            download_path = "/content/_dl_tisca_rlib.tar.gz"
            shutil.copy2(tar_candidates[0], download_path)
            with open(download_path, "rb") as f:
                observed_sha = hashlib.sha256(f.read()).hexdigest()
            assert observed_sha == BUNDLE_SHA256.lower(), "R library bundle SHA mismatch"
            if os.path.isdir("/content/tisca_rlib"):
                shutil.rmtree("/content/tisca_rlib")
            subprocess.run(["tar", "xzf", download_path, "-C", "/content"], check=True)
            LIBDIR = "/content/tisca_rlib/rlib"
            assert os.path.isdir(LIBDIR), "bundle did not restore the expected rlib directory"
            print("bundle restored:", LIBDIR, "from", BUNDLE_FOLDER_URL)
        """)
    else:
        bundle_restore = textwrap.dedent(f"""\
            import hashlib, os, re, shutil, subprocess, urllib.request

            BUNDLE_URL = {bundle_source['value']!r}
            BUNDLE_SHA256 = {bundle_sha!r}
            assert BUNDLE_URL not in ("", "{UNBUILT_URL}"), \\
                "Bundle URL is not built; regenerate with a real URL."
            assert re.fullmatch(r"[0-9a-fA-F]{{64}}", BUNDLE_SHA256), \\
                "Bundle SHA256 is missing or malformed; regenerate the notebook."

            DEST = "/content"
            LIBDIR = "/content/tisca_rlib/rlib"
            download_path = "/content/_dl_tisca_rlib.tar.gz"
            urllib.request.urlretrieve(BUNDLE_URL, download_path)
            with open(download_path, "rb") as f:
                observed_sha = hashlib.sha256(f.read()).hexdigest()
            assert observed_sha == BUNDLE_SHA256.lower(), "R library bundle SHA mismatch"
            if os.path.isdir("/content/tisca_rlib"):
                shutil.rmtree("/content/tisca_rlib")
            subprocess.run(["tar", "xzf", download_path, "-C", DEST], check=True)
            assert os.path.isdir(LIBDIR), "bundle did not restore the expected rlib directory"
            print("bundle restored:", LIBDIR)
        """)
    code(cells, bundle_restore)
    code(cells, textwrap.dedent("""\
        import os, urllib.request

        RUNCELL_URL = (
            "https://raw.githubusercontent.com/hugogobato/"
            "Test-Informed-Simulation-Count-Algorithm-TISCA/main/"
            "experiments/E3_mvbcf_casestudy/run_cell.R"
        )
        MVBCF_CPP_URL = (
            "https://raw.githubusercontent.com/Nathan-McJames/MVBCF_Paper/"
            "main/MVBCF_Code.cpp"
        )
        os.makedirs("/content/e3", exist_ok=True)
        urllib.request.urlretrieve(RUNCELL_URL, "/content/e3/run_cell.R")
        # The upstream C++ is downloaded at runtime and is never committed here.
        urllib.request.urlretrieve(MVBCF_CPP_URL, "/content/e3/MVBCF_Code.cpp")
        assert os.path.getsize("/content/e3/run_cell.R") > 1000
        assert os.path.getsize("/content/e3/MVBCF_Code.cpp") > 10000
        print("downloaded run_cell.R and upstream MVBCF_Code.cpp")
        """))
    code(cells, textwrap.dedent("""\
        import os, subprocess

        compile_script = "\\n".join([
            ".libPaths(c('/content/tisca_rlib/rlib', .libPaths()))",
            "if (!requireNamespace('Rcpp', quietly=TRUE) ||",
            "    !requireNamespace('RcppArmadillo', quietly=TRUE) ||",
            "    !requireNamespace('RcppDist', quietly=TRUE)) stop('bundle missing Rcpp dependencies')",
            "library(Rcpp)",
            "sourceCpp('/content/e3/MVBCF_Code.cpp')",
            "stopifnot(is.function(fast_bart))",
            "cat('FAST_BART_OK\\\\n')",
        ])
        with open("/content/e3/compile.R", "w") as f:
            f.write(compile_script)
        p = subprocess.run(["Rscript", "/content/e3/compile.R"],
                           capture_output=True, text=True)
        print(p.stdout[-3000:])
        if p.returncode != 0 or "FAST_BART_OK" not in p.stdout:
            print(p.stderr[-3000:])
            raise RuntimeError("upstream MVBCF C++ compilation failed")
        print("fast_bart() compiled")
        """))
    return cells


def _run_cells(row):
    cells = []
    md(cells, """\
    ## Fixed shard configuration

    The constants below were generated from `shard_table.csv`. They are
    assertions, not operator inputs. A dropped session can be restarted by
    uploading this same notebook, because existing seed rows are skipped.
    """)
    config = textwrap.dedent(f"""\
        import csv, os, subprocess, time

        DGP = {row['dgp']}
        N = {row['n']}
        MODE = {row['mode']!r}
        SHARD_ID = {row['notebook_filename'].split('_shard')[1].split('_')[0]!r}
        CLI_SEED_START = {row['cli_seed_start']}
        CLI_SEED_END = {row['cli_seed_end']}
        EXPECTED_SEED_START = {row['seed_start']}
        EXPECTED_SEED_END = {row['seed_end']}
        MC_CORES = 2
        DRIVE_DIR = "/content/drive/MyDrive/TISCA_E3"
        OUTPUT_CSV = os.path.join(DRIVE_DIR, {row['output_filename']!r})
        GIT_SHA = "not-a-git-build"

        assert DGP in (1, 2, 3)
        assert N in (100, 500)
        assert MODE in ("pilot", "confirmatory")
        assert MC_CORES == 2
        assert CLI_SEED_START == 0 or CLI_SEED_START > 0
        assert CLI_SEED_END >= CLI_SEED_START
        assert EXPECTED_SEED_END - EXPECTED_SEED_START + 1 == CLI_SEED_END - CLI_SEED_START + 1
        os.makedirs(DRIVE_DIR, exist_ok=True)

        def emitted_seed(raw_seed):
            return raw_seed + 1000001 if MODE == "pilot" else raw_seed

        def contiguous_ranges(values):
            values = sorted(values)
            if not values:
                return []
            out = []
            start = previous = values[0]
            for value in values[1:]:
                if value != previous + 1:
                    out.append((start, previous))
                    start = value
                previous = value
            out.append((start, previous))
            return out

        existing = set()
        if os.path.exists(OUTPUT_CSV):
            with open(OUTPUT_CSV, newline="") as f:
                rows = list(csv.DictReader(f))
            existing = {{int(r["seed"]) for r in rows if r.get("seed") not in (None, "")}}
            assert len(existing) == len(rows), \\
                "duplicate checkpoint seeds found; repair the Drive CSV before restarting"
        expected = {{emitted_seed(i) for i in range(CLI_SEED_START, CLI_SEED_END + 1)}}
        assert existing <= expected, "checkpoint contains a seed outside this shard"
        missing_raw = [i for i in range(CLI_SEED_START, CLI_SEED_END + 1)
                       if emitted_seed(i) not in existing]
        print("checkpoint:", len(existing), "rows; missing:", len(missing_raw))
        print("Drive CSV:", OUTPUT_CSV)
        """)
    code(cells, config)
    code(cells, textwrap.dedent("""\
        # Run only missing contiguous ranges. This makes a re-upload idempotent
        # while preserving run_cell.R's shard-offset-invariant RNG streams.
        env = dict(os.environ)
        env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")
        env["TISCA_MVBCF_CPP"] = "/content/e3/MVBCF_Code.cpp"
        env["TISCA_GIT_SHA"] = GIT_SHA
        for raw_start, raw_end in contiguous_ranges(missing_raw):
            cmd = ["Rscript", "/content/e3/run_cell.R", str(DGP), str(N),
                   str(raw_start), str(raw_end), "--out", OUTPUT_CSV,
                   "--cores", str(MC_CORES), "--mode", MODE]
            print("running missing range:", " ".join(cmd))
            t0 = time.time()
            result = subprocess.run(cmd, env=env, capture_output=True, text=True)
            print(result.stdout[-6000:])
            print("range wall-clock: %.1f min" % ((time.time() - t0) / 60.0))
            if result.returncode != 0:
                print(result.stderr[-4000:])
                raise RuntimeError("run_cell.R failed")
        """))
    code(cells, textwrap.dedent("""\
        import csv, os
        with open(OUTPUT_CSV, newline="") as f:
            rows = list(csv.DictReader(f))
        got = [int(r["seed"]) for r in rows]
        expected = list(range(EXPECTED_SEED_START, EXPECTED_SEED_END + 1))
        assert len(got) == len(set(got))
        assert set(got) == set(expected), "shard checkpoint is incomplete"
        failures = sum(r.get("converged_flag") == "0" for r in rows)
        print("completed rows:", len(rows), "of", len(expected))
        print("converged_flag failures:", failures)
        print("checkpoint bytes:", os.path.getsize(OUTPUT_CSV))
        try:
            from google.colab import files
            files.download(OUTPUT_CSV)
            print("Downloaded:", OUTPUT_CSV)
        except Exception as e:
            print("(Not on Colab / download skipped):", e)
        """))
    return cells


def _calibration_cells(row):
    cells = []
    md(cells, """\
    ## P3-T5(e) calibration gate

    This DGP1, n=500 pilot is the first notebook to run. The gate below must
    pass before any Round 1 notebook is started. The other Round 0 pilot
    notebooks are then run for their independent cell pilots.
    """)
    code(cells, textwrap.dedent("""\
        import csv, statistics
        with open(OUTPUT_CSV, newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 50, "the DGP1 n=500 calibration pilot must have 50 rows"
        def values(key):
            return [float(r[key]) for r in rows if r.get(key) not in (None, "")]
        targets = {
            "bcf_pehe1": (9.3, 10.0, "BCF PEHE Y1"),
            "bcf_pehe2": (9.6, 10.3, "BCF PEHE Y2"),
            "bcf_cov951": (0.95, 0.98, "BCF tau 95% coverage Y1"),
            "bcf_cov952": (0.94, 0.98, "BCF tau 95% coverage Y2"),
        }
        passed = True
        for key, (lo, hi, label) in targets.items():
            v = values(key)
            mean = statistics.mean(v)
            se = statistics.stdev(v) / (len(v) ** 0.5)
            ok = lo <= mean <= hi
            passed = passed and ok
            print(f"{label}: mean={{mean:.3f}} se={{se:.3f}} band=[{{lo}},{{hi}}] -> {{'PASS' if ok else 'FAIL'}}")
        print("P3-T5(e) VERDICT:", "PASS" if passed else "FAIL")
        assert passed, "calibration gate failed; diagnose before Round 1"
        """))
    return cells


def _sanity(nb, name):
    for index, cell in enumerate(nb["cells"]):
        source = "".join(cell["source"])
        assert "PASTE_" not in source, f"{name} cell {index} contains a placeholder"
        assert "'''" not in source and '"""' not in source, \
            f"{name} cell {index} contains an unsafe triple quote"


def write_outputs(repo_root, rows, bundle_source, bundle_sha):
    notebook_dir = repo_root / "notebooks" / "E3_shards"
    notebook_dir.mkdir(parents=True, exist_ok=True)
    for old in notebook_dir.glob("E3_*.ipynb"):
        old.unlink()

    row_by_name = {}
    for row in rows:
        cells = _shared_cells(bundle_source, bundle_sha, row, repo_root)
        cells.extend(_run_cells(row))
        if row["mode"] == "pilot" and row["dgp"] == 1 and row["n"] == 500:
            cells.extend(_calibration_cells(row))
        nb = notebook(cells)
        _sanity(nb, row["notebook_filename"])
        path = notebook_dir / row["notebook_filename"]
        path.write_text(json.dumps(nb, indent=1) + "\n")
        row_by_name[row["notebook_filename"]] = nb

    # Preserve the two historical entry points, generated from the exact same
    # code and with real shard parameters baked in.
    first_confirm = next(r for r in rows if r["mode"] == "confirmatory")
    first_pilot = next(r for r in rows
                       if r["mode"] == "pilot" and r["dgp"] == 1 and r["n"] == 500)
    (repo_root / "notebooks" / "E3_mvbcf_shard.ipynb").write_text(
        json.dumps(row_by_name[first_confirm["notebook_filename"]], indent=1) + "\n")
    (repo_root / "notebooks" / "E3_round0_pilot_calibration.ipynb").write_text(
        json.dumps(row_by_name[first_pilot["notebook_filename"]], indent=1) + "\n")

    table_path = repo_root / "experiments" / "E3_mvbcf_casestudy" / "shard_table.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with table_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} shard notebooks to {notebook_dir}")
    print(f"wrote {table_path}")
    print("wrote historical templates from the first confirmatory and calibration shards")


def main():
    parser = argparse.ArgumentParser()
    _bundle_args(parser)
    args = parser.parse_args()
    bundle_url, bundle_sha = _validated_bundle(args)
    repo_root = Path(__file__).resolve().parents[2]
    rows = build_shards()
    write_outputs(repo_root, rows, bundle_url, bundle_sha)


if __name__ == "__main__":
    main()
