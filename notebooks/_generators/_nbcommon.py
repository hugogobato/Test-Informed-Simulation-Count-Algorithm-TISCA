"""Shared notebook-building helpers for the Phase-3 analysis runners.

``build_e1_modules_nb.py`` and ``build_e3_notebooks.py`` predate this module and
keep their own copies; the four analysis notebooks (P3-T3, P3-T4, P3-T6, P3-T7)
share these instead so their boilerplate exists once.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def md(cells, source):
    cells.append({"cell_type": "markdown", "metadata": {},
                  "source": textwrap.dedent(source).strip("\n").splitlines(keepends=True)})


def code(cells, source):
    cells.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [],
                  "source": textwrap.dedent(source).strip("\n").splitlines(keepends=True)})


def notebook(cells):
    return {"nbformat": 4, "nbformat_minor": 0,
            "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                         "language_info": {"name": "python"}},
            "cells": cells}


def write(name, cells, root=ROOT):
    """Serialise, byte-compile every code cell, and report."""
    nb = notebook(cells)
    out = Path(root) / name
    out.write_text(json.dumps(nb, indent=1) + "\n")
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] == "code":
            compile("".join(cell["source"]), f"{name}:cell-{i}", "exec")
    print("wrote", out, len(nb["cells"]), "cells")
    return out


SETUP = r'''
import os
import sys
import subprocess

import numpy as np
import pandas as pd

REPO_URL = "https://github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA.git"
CANDIDATES = [
    os.path.abspath(os.path.join(os.getcwd(), "..")),          # notebooks/ in a checkout
    os.getcwd(),
    "/content/Test-Informed-Simulation-Count-Algorithm-TISCA",
    "/content/TISCA_repo",
]
REPO_ROOT = next((p for p in CANDIDATES
                  if os.path.isdir(os.path.join(p, "tisca", "python"))), None)
if REPO_ROOT is None:
    subprocess.run(["git", "clone", "--depth", "1", REPO_URL, "/content/TISCA_repo"], check=True)
    REPO_ROOT = "/content/TISCA_repo"
sys.path.insert(0, os.path.join(REPO_ROOT, "tisca", "python"))

RESULTS = os.path.join(REPO_ROOT, "results")
FIGURES = os.path.join(REPO_ROOT, "figures")
os.makedirs(FIGURES, exist_ok=True)
print("repo:", REPO_ROOT)


def download(path):
    """Colab download fallback (standing rule); a no-op off Colab."""
    try:
        from google.colab import files
        files.download(path)
        print("Downloaded:", path)
    except Exception as e:
        print("(Not on Colab / download skipped):", e)
'''

PLOT_STYLE = r'''
import matplotlib
import matplotlib.pyplot as plt

# Print-size legibility: SoutoNeto section 4 asked for larger type in every figure,
# and these are single-column figures in the revised manuscript.
matplotlib.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
    "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
    "axes.spines.right": False, "figure.autolayout": False,
})
# Colour-blind-safe qualitative set, consistent across every Phase-3 figure.
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]
'''
