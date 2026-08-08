#!/usr/bin/env python3
"""Generate the two reproducible Colab shards for genuine joint Module B."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _source(text: str) -> list[str]:
    return textwrap.dedent(text).lstrip("\n").splitlines(keepends=True)


def _markdown(cell_id: str, text: str) -> dict:
    return {"id": cell_id, "cell_type": "markdown", "metadata": {}, "source": _source(text)}


def _code(cell_id: str, text: str) -> dict:
    return {
        "id": cell_id, "cell_type": "code", "execution_count": None,
        "metadata": {}, "outputs": [], "source": _source(text),
    }


def notebook(shard: int) -> dict:
    start, end = ((0, 330) if shard == 1 else (330, 660))
    cells = [
        _markdown("joint-intro", f"""
            # E1 joint Module B, shard {shard} of 2

            This notebook runs canonical cells `B_{start:04d}` through
            `B_{end - 1:04d}`. Every outer repetition generates one common
            method-A loss and K distinct benchmark losses, applies the declared
            correction jointly, and checkpoints one completed cell at a time.

            Run the cells in order. The deterministic self-check is a gate: the
            expensive shard does not start unless all joint-path validation tests
            pass. Checkpoints and the per-cell manifest live on Google Drive when
            Drive is available, so restarting this notebook resumes by `cell_id`
            without duplicating completed rows.
        """),
        _code("joint-config", f"""
            import os

            SHARD_INDEX = {shard}
            NUM_SHARDS = 2
            REPO_URL = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
            LOCAL_REPO = "/content/Test-Informed-Simulation-Count-Algorithm-TISCA"
            CLONED_REPO = "/content/TISCA_repo"

            for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                         "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
                os.environ[name] = "1"

            print("Shard:", SHARD_INDEX, "of", NUM_SHARDS)
            print("Expected IDs: B_{start:04d} through B_{end - 1:04d}")
            print("Expected seeds: {400000 + start} through {400000 + end - 1}")
        """),
        _code("joint-environment", """
            import glob
            import json
            import os
            import pathlib
            import socket
            import subprocess
            import sys
            import time

            if os.path.isdir(os.path.join(LOCAL_REPO, "tisca", "python")):
                REPO_ROOT = LOCAL_REPO
            elif os.path.isdir(os.path.join(CLONED_REPO, "tisca", "python")):
                REPO_ROOT = CLONED_REPO
            else:
                subprocess.run(["git", "clone", "--depth", "1", REPO_URL, CLONED_REPO], check=True)
                REPO_ROOT = CLONED_REPO

            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-r",
                            os.path.join(REPO_ROOT, "requirements.txt")], check=True)
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-e", REPO_ROOT], check=True)

            try:
                from google.colab import drive
                drive.mount("/content/drive")
                JOINT_ROOT = "/content/drive/MyDrive/TISCA_E1/module_b_joint"
            except Exception as exc:
                JOINT_ROOT = "/content/TISCA_E1/module_b_joint"
                print("(Drive unavailable; local checkpoint fallback):", exc)

            SHARD_DIR = os.path.join(JOINT_ROOT, f"shard{SHARD_INDEX}")
            COMBINED_DIR = os.path.join(JOINT_ROOT, "combined")
            os.makedirs(SHARD_DIR, exist_ok=True)
            os.makedirs(COMBINED_DIR, exist_ok=True)
            print("Repository:", REPO_ROOT)
            print("Checkpoint directory:", SHARD_DIR)
            print("Host:", socket.gethostname())
        """),
        _markdown("joint-selfcheck-heading", """
            ## Deterministic self-check

            This runs the same eleven joint validation gates used by the project
            test suite, including FWER behavior, BH labeling, covariance PSD,
            shared-control mapping, Romano-Wolf dependence, seed identity, shard
            offsets, and the prohibition on copied benchmark columns.
        """),
        _code("joint-selfcheck", """
            SELFCHECK = os.path.join(
                REPO_ROOT, "experiments", "E1_operating_characteristics",
                "selfcheck_module_b_joint.py",
            )
            subprocess.run([sys.executable, SELFCHECK], cwd=REPO_ROOT, check=True,
                           env=os.environ.copy())
            print("[PASS] pre-grid joint self-check")
        """),
        _markdown("joint-run-heading", """
            ## Run or resume this shard

            The runner imports the 660-cell canonical grid, selects a contiguous
            deterministic half, writes an atomic checkpoint after every cell, and
            stores `cell_id`, seed, status, hostname, runtime, shard metadata, and
            output paths in the manifest. BLAS remains single-threaded and worker
            count is selected from current load and available RAM.
        """),
        _code("joint-run", """
            RUNNER = os.path.join(
                REPO_ROOT, "experiments", "E1_operating_characteristics",
                "run_e1_module_b_joint.py",
            )
            command = [
                sys.executable, RUNNER,
                "--shard-index", str(SHARD_INDEX),
                "--num-shards", str(NUM_SHARDS),
                "--outdir", SHARD_DIR,
                "--resume",
            ]
            print("Running:", " ".join(command))
            subprocess.run(command, cwd=REPO_ROOT, check=True, env=os.environ.copy())

            shard_outputs = sorted(
                glob.glob(os.path.join(SHARD_DIR, "module_b_joint_*_shard*"))
            )
            assert len(shard_outputs) == 4, shard_outputs
            print("[PASS] shard outputs:")
            print("\\n".join(shard_outputs))
        """),
        _markdown("joint-combine-heading", """
            ## Combine when both shards are present

            This cell is safe in either notebook. It combines only when both Drive
            shard directories exist and contain complete outputs. The combiner
            asserts exactly 660 canonical IDs and seeds, with no duplicates.
        """),
        _code("joint-combine", """
            shard1_dir = os.path.join(JOINT_ROOT, "shard1")
            shard2_dir = os.path.join(JOINT_ROOT, "shard2")
            required = [
                os.path.join(shard1_dir, "module_b_joint_operating_characteristics_shard1.csv"),
                os.path.join(shard2_dir, "module_b_joint_operating_characteristics_shard2.csv"),
            ]
            combined_outputs = []
            if all(os.path.exists(path) for path in required):
                COMBINER = os.path.join(
                    REPO_ROOT, "experiments", "E1_operating_characteristics",
                    "combine_e1_module_b_joint.py",
                )
                subprocess.run([
                    sys.executable, COMBINER, shard1_dir, shard2_dir,
                    "--outdir", COMBINED_DIR,
                ], cwd=REPO_ROOT, check=True, env=os.environ.copy())
                combined_outputs = sorted(glob.glob(os.path.join(COMBINED_DIR, "module_b_joint_*")))
                assert len(combined_outputs) == 4, combined_outputs
                print("[PASS] combined 660-cell outputs")
            else:
                print("Other shard not complete yet; combine after it finishes.")
        """),
        _markdown("joint-download-heading", """
            ## Download every generated output

            Drive remains the persistent checkpoint. This final fallback also
            offers every shard output, and every combined output when available,
            through the Colab browser download API.
        """),
        _code("joint-download", """
            def download_output(output_file):
                try:
                    from google.colab import files
                    files.download(output_file)
                    print("Downloaded:", output_file)
                except Exception as e:
                    print("(Not on Colab / download skipped):", e)

            for output_file in shard_outputs + combined_outputs:
                download_output(output_file)
        """),
    ]
    for index, cell in enumerate(cells):
        if cell["cell_type"] == "code":
            compile("".join(cell["source"]), f"joint-shard{shard}-cell{index}", "exec")
    return {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"name": "python3", "display_name": "Python 3"},
            "language_info": {"name": "python"},
            "colab": {"provenance": [], "toc_visible": True},
        },
        "cells": cells,
    }


def main() -> None:
    for shard in (1, 2):
        path = ROOT / f"E1_modB_joint_shard{shard}.ipynb"
        path.write_text(json.dumps(notebook(shard), indent=1) + "\n", encoding="utf-8")
        print("wrote", path)


if __name__ == "__main__":
    main()
