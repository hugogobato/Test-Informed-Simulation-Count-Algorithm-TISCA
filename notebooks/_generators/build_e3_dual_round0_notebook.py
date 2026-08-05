#!/usr/bin/env python3
"""Build the paired Round 0 BCF comparison notebook.

The notebook is intentionally separate from the 37-shard E3 generator. It
runs the same DGP1/n=500 pilot seeds through the calibrated stochtree BCF
translation and the original bcf package, producing one paired CSV.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

from build_e3_notebooks import (
    DEFAULT_BUNDLE_FOLDER_URL,
    UNBUILT_SHA,
    UNBUILT_URL,
    _bundle_args,
    _sanity,
    _validated_bundle,
    code,
    md,
    notebook,
)


OUTPUT_NAME = "E3_DGP1_n500_round0_dual_bcf_pilot_seeds1000001-1000050.csv"


def _bundle_restore_cell(bundle_source, bundle_sha):
    if bundle_source["kind"] == "gdrive_folder":
        source = textwrap.dedent(f"""
            import hashlib, os, re, shutil, subprocess, sys
            from pathlib import Path

            BUNDLE_FOLDER_URL = {bundle_source['value']!r}
            BUNDLE_SHA256 = {bundle_sha!r}
            assert re.fullmatch(r"[0-9a-fA-F]{{64}}", BUNDLE_SHA256)
            download_dir = Path("/content/tisca_dual_bundle_download")
            if download_dir.exists():
                shutil.rmtree(download_dir)
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "gdown"], check=True)
            result = subprocess.run(
                [sys.executable, "-m", "gdown", "--folder", BUNDLE_FOLDER_URL,
                 "--output", str(download_dir), "--remaining-ok"],
                capture_output=True, text=True)
            print(result.stdout[-4000:])
            if result.returncode != 0:
                print(result.stderr[-4000:])
                raise RuntimeError("Google Drive bundle download failed")
            tarballs = sorted(download_dir.rglob("tisca_rlib.tar.gz"))
            sidecars = sorted(download_dir.rglob("tisca_rlib.sha256"))
            assert len(tarballs) == 1, f"expected one bundle tarball: {{tarballs}}"
            assert len(sidecars) <= 1, f"expected at most one checksum sidecar: {{sidecars}}"
            if sidecars:
                assert sidecars[0].read_text().split()[0].lower() == BUNDLE_SHA256.lower()
            download_path = "/content/_dl_tisca_dual_rlib.tar.gz"
            shutil.copy2(tarballs[0], download_path)
            with open(download_path, "rb") as f:
                assert hashlib.sha256(f.read()).hexdigest() == BUNDLE_SHA256.lower()
            if os.path.isdir("/content/tisca_rlib"):
                shutil.rmtree("/content/tisca_rlib")
            subprocess.run(["tar", "xzf", download_path, "-C", "/content"], check=True)
            LIBDIR = "/content/tisca_rlib/rlib"
            assert os.path.isdir(LIBDIR)
            print("bundle restored:", LIBDIR)
        """)
    else:
        source = textwrap.dedent(f"""
            import hashlib, os, re, shutil, subprocess, urllib.request

            BUNDLE_URL = {bundle_source['value']!r}
            BUNDLE_SHA256 = {bundle_sha!r}
            assert BUNDLE_URL not in ("", {UNBUILT_URL!r}), "bundle URL is not built"
            assert re.fullmatch(r"[0-9a-fA-F]{{64}}", BUNDLE_SHA256)
            download_path = "/content/_dl_tisca_dual_rlib.tar.gz"
            urllib.request.urlretrieve(BUNDLE_URL, download_path)
            with open(download_path, "rb") as f:
                assert hashlib.sha256(f.read()).hexdigest() == BUNDLE_SHA256.lower()
            if os.path.isdir("/content/tisca_rlib"):
                shutil.rmtree("/content/tisca_rlib")
            subprocess.run(["tar", "xzf", download_path, "-C", "/content"], check=True)
            LIBDIR = "/content/tisca_rlib/rlib"
            assert os.path.isdir(LIBDIR)
            print("bundle restored:", LIBDIR)
        """)
    return source


def build_notebook(bundle_source, bundle_sha):
    cells = []
    md(cells, """
    # E3 Round 0: paired BCF benchmark

    This notebook runs the same 50 DGP1, `n=500` pilot replications through:
    
    1. the calibrated `stochtree::bcf` translation with the paper's 50/20 tree
       structure, standard coding, `num_gfr=0`, 500 burn-in iterations, 500
       posterior draws, and fixed standardized leaf variances. These are
       `1^2/50=0.02` for the prognostic forest and
       `0.375^2/20=0.00703125` for the treatment-effect forest, whose
       forest-level variances are 1 and 0.140625; and
    2. the original `bcf` package using the settings in the paper's
       `GitHub_DGP1.R`.

    Both methods use the same generated data and BART propensity estimate for
    each emitted seed. Their results are written to separate columns in one
    checkpoint CSV. The notebook reports both calibration verdicts without
    stopping merely because one method deviates from the published target.
    """)
    code(cells, textwrap.dedent("""
        import os, platform, subprocess

        def sh(cmd):
            p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return p.stdout.strip()

        print("hostname:", platform.node() or "n/a")
        print("os:", platform.platform())
        print("nproc:", sh("nproc"))
        print("cpu:", sh("grep -m1 -E 'model name' /proc/cpuinfo") or "n/a")
        print(sh("free -g") or "RAM information unavailable")
        """))
    code(cells, textwrap.dedent("""
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
    code(cells, _bundle_restore_cell(bundle_source, bundle_sha))
    code(cells, textwrap.dedent("""
        import os, subprocess, urllib.request

        DUAL_DRIVER_URL = (
            "https://raw.githubusercontent.com/hugogobato/"
            "Test-Informed-Simulation-Count-Algorithm-TISCA/main/"
            "experiments/E3_mvbcf_casestudy/run_cell_round0_dual_bcf.R"
        )
        os.makedirs("/content/e3", exist_ok=True)
        urllib.request.urlretrieve(DUAL_DRIVER_URL, "/content/e3/run_cell_round0_dual_bcf.R")
        assert os.path.getsize("/content/e3/run_cell_round0_dual_bcf.R") > 10000
        with open("/content/e3/run_cell_round0_dual_bcf.R") as f:
            driver_source = f.read()
        required = [
            "num_gfr = 0",
            "sigma2_leaf_init = 1^2 / n_tree_mu",
            "sigma2_leaf_init = 0.375^2 / n_tree_tau",
            "bcf:::predict.bcf",
        ]
        missing = [item for item in required if item not in driver_source]
        assert not missing, f"GitHub main is serving a stale dual driver: {missing}"

        # The shared bundle will contain bcf after the next bundle rebuild. The
        # fallback keeps this notebook runnable with an older bundle as well.
        env = dict(os.environ)
        env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")
        install = subprocess.run(
            ["Rscript", "-e",
             "if (!requireNamespace('bcf', quietly=TRUE)) "
             "install.packages('bcf', lib='/content/tisca_rlib/rlib', "
             "repos='https://cloud.r-project.org')"],
            env=env, capture_output=True, text=True)
        print(install.stdout[-3000:])
        if install.returncode != 0:
            print(install.stderr[-3000:])
            raise RuntimeError("bcf dependency installation/check failed")
        print(subprocess.check_output(
            ["Rscript", "-e",
             ".libPaths(c(Sys.getenv('R_LIBS'), .libPaths())); "
             "cat('bcf=', as.character(packageVersion('bcf')), '\\n'); "
             "cat('stochtree=', as.character(packageVersion('stochtree')), '\\n')"],
            env=env, text=True))
        """))
    code(cells, textwrap.dedent(f"""
        import csv, os, shutil, time

        DGP = 1
        N = 500
        MODE = "pilot"
        CLI_SEED_START = 0
        CLI_SEED_END = 49
        EXPECTED_SEED_START = 1000001
        EXPECTED_SEED_END = 1000050
        MC_CORES = 2
        OUTPUT_DIR = "/content/TISCA_E3"
        OUTPUT_CSV = os.path.join(OUTPUT_DIR, {OUTPUT_NAME!r})
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        def emitted_seed(raw_seed):
            return raw_seed + 1000001

        expected = {{emitted_seed(i) for i in range(CLI_SEED_START, CLI_SEED_END + 1)}}

        def checkpoint_seed(row):
            try:
                return int((row.get("seed") or "").strip())
            except (TypeError, ValueError):
                return None

        if os.path.exists(OUTPUT_CSV):
            with open(OUTPUT_CSV, newline="") as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames
                rows = list(reader)
            if not fieldnames or "seed" not in fieldnames or "converged_flag" not in fieldnames:
                raise RuntimeError("dual checkpoint is missing required columns")
            kept, seen = [], set()
            rejected = {{"malformed": 0, "failed": 0, "duplicate": 0}}
            for row in rows:
                seed = checkpoint_seed(row)
                if seed is None or any(value is None for value in row.values()):
                    rejected["malformed"] += 1
                elif seed not in expected:
                    raise RuntimeError(f"checkpoint contains seed {{seed}} outside this pilot")
                elif row.get("converged_flag") != "1":
                    rejected["failed"] += 1
                elif seed in seen:
                    rejected["duplicate"] += 1
                else:
                    seen.add(seed)
                    kept.append(row)
            if any(rejected.values()):
                backup = f"{{OUTPUT_CSV}}.pre_repair_{{int(time.time())}}.csv"
                shutil.copy2(OUTPUT_CSV, backup)
                temporary = OUTPUT_CSV + ".repair.tmp"
                with open(temporary, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(kept)
                os.replace(temporary, OUTPUT_CSV)
                print("checkpoint repair:", rejected, "backup:", backup)
            existing = seen
        else:
            existing = set()

        lock = OUTPUT_CSV + ".lock"
        if os.path.isdir(lock):
            shutil.rmtree(lock)
        elif os.path.isfile(lock):
            os.remove(lock)
        missing_raw = [i for i in range(CLI_SEED_START, CLI_SEED_END + 1)
                       if emitted_seed(i) not in existing]
        print("checkpoint rows:", len(existing), "missing:", len(missing_raw))
        print("local checkpoint:", OUTPUT_CSV)
        """))
    code(cells, textwrap.dedent("""
        import os, subprocess, time

        def contiguous_ranges(values):
            values = sorted(values)
            if not values:
                return []
            result = []
            start = previous = values[0]
            for value in values[1:]:
                if value != previous + 1:
                    result.append((start, previous))
                    start = value
                previous = value
            result.append((start, previous))
            return result

        env = dict(os.environ)
        env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")
        env["TISCA_GIT_SHA"] = "not-a-git-build"
        for raw_start, raw_end in contiguous_ranges(missing_raw):
            command = ["Rscript", "/content/e3/run_cell_round0_dual_bcf.R",
                       str(DGP), str(N), str(raw_start), str(raw_end),
                       "--out", OUTPUT_CSV, "--cores", str(MC_CORES),
                       "--mode", MODE]
            print("running:", " ".join(command))
            t0 = time.time()
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            print(result.stdout[-8000:])
            print("range wall-clock: %.1f min" % ((time.time() - t0) / 60.0))
            if result.returncode != 0:
                print(result.stderr[-5000:])
                raise RuntimeError("paired BCF driver failed")
        """))
    md(cells, """
    ## Paired Round 0 summary

    The original `bcf` package returns draws by row, whereas stochtree returns
    draws by column. The driver checks both dimensions before calculating
    PEHE, coverage, CRPS, and interval widths. The four target bands below are
    the same DGP1 Table 2 bands used by the existing calibration notebook.
    """)
    code(cells, textwrap.dedent("""
        import csv, math, statistics
        with open(OUTPUT_CSV, newline="") as f:
            rows = list(csv.DictReader(f))
        expected = list(range(EXPECTED_SEED_START, EXPECTED_SEED_END + 1))
        got = [checkpoint_seed(row) for row in rows]
        assert len(rows) == 50, f"expected 50 rows, got {len(rows)}"
        assert len(set(got)) == 50 and sorted(got) == expected
        failures = [row for row in rows if row.get("converged_flag") != "1"]
        assert not failures, [row.get("error_message") for row in failures[:3]]

        targets = {
            "stochtree_bcf1_pehe": (9.3, 10.0, "stochtree BCF PEHE Y1"),
            "stochtree_bcf2_pehe": (9.6, 10.3, "stochtree BCF PEHE Y2"),
            "stochtree_bcf1_cov95": (0.95, 0.98, "stochtree BCF 95% coverage Y1"),
            "stochtree_bcf2_cov95": (0.94, 0.98, "stochtree BCF 95% coverage Y2"),
            "paper_bcf1_pehe": (9.3, 10.0, "paper BCF PEHE Y1"),
            "paper_bcf2_pehe": (9.6, 10.3, "paper BCF PEHE Y2"),
            "paper_bcf1_cov95": (0.95, 0.98, "paper BCF 95% coverage Y1"),
            "paper_bcf2_cov95": (0.94, 0.98, "paper BCF 95% coverage Y2"),
        }
        for key, (lo, hi, label) in targets.items():
            values = [float(row[key]) for row in rows]
            assert all(math.isfinite(value) for value in values)
            mean = statistics.mean(values)
            se = statistics.stdev(values) / (len(values) ** 0.5)
            verdict = "PASS" if lo <= mean <= hi else "FAIL"
            print(f"{{label}}: mean={{mean:.3f}} se={{se:.3f}} band=[{{lo}},{{hi}}] -> {{verdict}}")

        print("\\nPaired Round 0 completed. A FAIL is diagnostic output, not a notebook assertion.")
        print("CSV:", OUTPUT_CSV)
        try:
            from google.colab import files
            files.download(OUTPUT_CSV)
            print("Downloaded:", OUTPUT_CSV)
        except Exception as e:
            print("(Not on Colab / download skipped):", e)
        """))
    return cells


def main():
    parser = argparse.ArgumentParser()
    _bundle_args(parser)
    args = parser.parse_args()
    bundle_source, bundle_sha = _validated_bundle(args)
    repo_root = Path(__file__).resolve().parents[2]
    cells = build_notebook(bundle_source, bundle_sha)
    nb = notebook(cells)
    _sanity(nb, "E3_round0_dual_bcf_pilot.ipynb")
    path = repo_root / "notebooks" / "E3_round0_dual_bcf_pilot.ipynb"
    path.write_text(json.dumps(nb, indent=1) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
