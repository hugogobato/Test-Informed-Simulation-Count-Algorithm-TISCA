#!/usr/bin/env python3
"""Generate the E3 MVBCF-case-study Colab notebooks (plan P3-T5[d] and P3-T5[e]).

Produces two notebooks under Test-Informed-Simulation-Count-Algorithm-TISCA/
notebooks/:

  E3_mvbcf_shard.ipynb
      The parameterised worker for an E3 confirmatory/pilot shard. Cell 2 is the
      ONLY cell the operator edits: one line per session setting SHARD_ID, DGP,
      N, SEED_START, SEED_END, MC_CORES. Everything else is generic. It:
        - installs R (apt) on the Colab runtime
        - mounts Google Drive (checkpoint target, no shared filesystem)
        - restores the P0-T4 R library bundle from a direct URL + SHA256
        - downloads run_cell.R (this repo) and MVBCF_Code.cpp (the upstream,
          read-only NVBCF repo) and sourceCpp's the cpp
        - runs `Rscript run_cell.R <dgp> <n> <seed_start> <seed_end> --out <drive> \
          --cores <MC_CORES> [--mode pilot|confirmatory]`
        - appends every completed replication to a per-shard CSV on Drive
          (one row per seed), so a killed 8-h session loses at most one row
        - ends with the google.colab.files.download fallback per the user's rule

  E3_round0_pilot_calibration.ipynb
      The Round 0 pilot notebook that runs 50 pilot replications of DGP1/n=500
      first (the only cell the P3-T5[e] calibration gate needs) and then
      evaluates the calibrated stochtree::bcf benchmark against McJames et al.'s
      published DGP1 n=500 Table 2. This gate BLOCKS every confirmatory shard.

Each notebook cell source is stored as a Python list of lines. We build cell
source with plain triple-quoted strings that NEVER contain a nested triple
quote or an unescaped `!` inside an f-string literal (avoid f-strings).
"""
import json, os, textwrap

CELLS_A = []  # shared setup cells for both notebooks
CELLS_B = []  # shard-run cells
CELLS_C = []  # round-0 pilot + calibration cells


def md(cells, src):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)})


def code(cells, src):
    cells.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.splitlines(keepends=True)})


# ---------------------------------------------------------------------------
# Shared setup (both notebooks start identically)
# ---------------------------------------------------------------------------

code(CELLS_A, textwrap.dedent("""\
    import platform, subprocess
    def _sh(cmd):
        return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()
    print("hostname :", platform.node() or "n/a")
    print("os       :", platform.platform())
    print("nproc    :", _sh("nproc"))
    print("cpu      :", _sh("grep -m1 -E 'model name' /proc/cpuinfo") or "no /proc/cpuinfo")
    print("ram:")
    print(_sh("free -g") or "free unavailable")
    """))

code(CELLS_A, textwrap.dedent("""\
    !apt-get -qq update > /dev/null
    !apt-get -qq install -y --no-install-recommends r-base r-base-dev libcurl4-openssl-dev > /dev/null 2>&1
    !R --version | head -1
    print("R install done")
    """))

code(CELLS_A, textwrap.dedent("""\
    from google.colab import drive
    drive.mount("/content/drive")
    print("Drive mounted")
    """))

code(CELLS_A, textwrap.dedent("""\
    import os, subprocess, urllib.request, hashlib

    # Paste the URL + SHA256 printed by the P0-T4 bundle notebook (Cell 7).
    BUNDLE_URL     = "PASTE_BUNDLE_URL"      # direct link to tisca_rlib.tar.gz
    BUNDLE_SHA256  = "PASTE_BUNDLE_SHA256"   # Cell 7 of P0T4_build_rlib_bundle.ipynb
    assert BUNDLE_URL != "PASTE_BUNDLE_URL", "paste BUNDLE_URL"
    assert BUNDLE_SHA256 != "PASTE_BUNDLE_SHA256", "paste BUNDLE_SHA256"

    DEST = "/content"
    LIBDIR = "/content/tisca_rlib/rlib"
    DL = "/content/_dl_tisca_rlib.tar.gz"
    urllib.request.urlretrieve(BUNDLE_URL, DL)
    h = hashlib.sha256(open(DL, "rb").read()).hexdigest()
    assert h == BUNDLE_SHA256, "SHA mismatch"
    if os.path.exists(LIBDIR):
        subprocess.run(["rm", "-rf", "/content/tisca_rlib"])
    subprocess.run(["tar", "xzf", DL, "-C", DEST])
    assert os.path.isdir(LIBDIR)
    print("bundle restored:", LIBDIR)
    print("tisca_rlib.rlib entries:", len(os.listdir(LIBDIR)))
    """))

code(CELLS_A, textwrap.dedent("""\
    import os, subprocess, urllib.request

    # run_cell.R lives in this repo; MVBCF_Code.cpp is the UPSTREAM file (the
    # read-only NVBCF repo). We fetch both at runtime and never copy the .cpp
    # into the TISCA repository (licensing/attribution rule, plan §0.1).
    RUNCELL_URL = "https://raw.githubusercontent.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA/main/experiments/E3_mvbcf_casestudy/run_cell.R"
    MVBCF_CPP_URL = "https://raw.githubusercontent.com/Nathan-McJames/MVBCF_Paper/main/MVBCF_Code.cpp"

    os.makedirs("/content/e3", exist_ok=True)
    urllib.request.urlretrieve(RUNCELL_URL, "/content/e3/run_cell.R")
    urllib.request.urlretrieve(MVBCF_CPP_URL, "/content/e3/MVBCF_Code.cpp")
    assert os.path.getsize("/content/e3/run_cell.R") > 1000
    assert os.path.getsize("/content/e3/MVBCF_Code.cpp") > 10000
    print("run_cell.R :", os.path.getsize("/content/e3/run_cell.R"), "bytes")
    print("MVBCF_Code.cpp :", os.path.getsize("/content/e3/MVBCF_Code.cpp"), "bytes")
    """))

code(CELLS_A, textwrap.dedent("""\
    # Compile the upstream MVBCF C++ once, verify fast_bart() is exposed, and
    # confirm the bundle satisfies Rcpp + RcppArmadillo + RcppDist. Compilation
    # happens once per session (not per replication).
    import subprocess, os
    compile_r = (
        '.libPaths(c("/content/tisca_rlib/rlib", .libPaths()))\\n' +
        'if (!requireNamespace("Rcpp", quietly=TRUE) \\n' +
        '    || !requireNamespace("RcppArmadillo", quietly=TRUE) \\n' +
        '    || !requireNamespace("RcppDist", quietly=TRUE)) \\n' +
        '  stop("bundle missing Rcpp/RcppArmadillo/RcppDist") \\n' +
        'library(Rcpp) \\n' +
        'sourceCpp("/content/e3/MVBCF_Code.cpp") \\n' +
        'stopifnot(is.function(fast_bart)) \\n' +
        'cat("FAST_BART_OK\\\\n") \\n'
    )
    with open("/content/e3/compile.R", "w") as f:
        f.write(compile_r)
    res = subprocess.run(["Rscript", "/content/e3/compile.R"], capture_output=True, text=True)
    print(res.stdout[-3000:])
    if res.returncode != 0 or "FAST_BART_OK" not in res.stdout:
        print("STDERR:", res.stderr[-3000:])
        raise SystemExit("Rcpp compilation of MVBCF_Code.cpp failed")
    print("fast_bart() compiled and available")
    """))

# ---------------------------------------------------------------------------
# E3 confirmatory / general shard: config + run
# ---------------------------------------------------------------------------

md(CELLS_B, textwrap.dedent("""\
    ## E3 worker shard

    **Edit Cell 2 (config) only.** One notebook = one session = one shard of the
    MVBCF case-study re-run. `SEED_START`/`SEED_END` are the small per-master
    replication indices for this shard (confirmatory 0..N-1; pilot 0..J0-1). The
    shard table at the bottom of `REVISION_PLAN.md` §5 tells you which ranges to
    assign. Size each shard for ~6 h of the 8-h budget.

    Output rows are appended to `{DRIVE_DIR}/<cell>_<shardid>.csv` as each
    replication completes; a 8-h kill loses at most one replication.
    """))

code(CELLS_B, textwrap.dedent("""\
    # ================= CONFIG (edit this cell only) =================
    SHARD_ID    = 1        # a short partition label, e.g. 0, 1, 2, ... (just for the filename)
    DGP         = 1        # 1, 2, or 3
    N           = 500      # training sample size: 100 or 500
    SEED_START  = 0        # first replication index for this shard (inclusive)
    SEED_END    = 9        # last replication index for this shard (inclusive)
    MC_CORES    = 2        # Colab is 2 vCPUs; use 2 (plan P0-T3 measured SPEEDUP=1.62)
    MODE        = "confirmatory"  # or "pilot"
    # Directory on Drive where the shard CSV is stored (no shared filesystem;
    # this is the checkpoint target and the source for concat/verification).
    DRIVE_DIR   = "/content/drive/MyDrive/TISCA_E3"
    # ==================================================================
    import os
    assert 1 <= DGP <= 3 and N in (100, 500)
    assert 0 <= SEED_START <= SEED_END
    assert MC_CORES >= 1
    os.makedirs(DRIVE_DIR, exist_ok=True)
    shard_tag = f"DGP{DGP}_n{N}_shard{SHARD_ID}"
    out = os.path.join(DRIVE_DIR, f"{shard_tag}.csv")
    print("shard        :", shard_tag)
    print("mode         :", MODE)
    print("seeds        :", SEED_START, "..", SEED_END)
    print("out csv      :", out)
    print("NOTE: `N` is already handled by run_cell.R via CLI arg; `N` above is "
          "the training size for the DGP, and `n_test` is fixed at 1000 inside the driver.")
    """))

code(CELLS_B, textwrap.dedent("""\
    import os, subprocess, time

    env = dict(os.environ)
    env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")
    env["TISCA_MVBCF_CPP"] = "/content/e3/MVBCF_Code.cpp"
    # R_REPOS directly; the bundle sets none. Keep R_LIBS to the bundle first.
    cmd = [
        "Rscript", "/content/e3/run_cell.R",
        str(DGP), str(N), str(SEED_START), str(SEED_END),
        "--out", out, "--cores", str(MC_CORES), "--mode", MODE,
    ]
    print("running:", " ".join(cmd))
    t0 = time.time()
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    wall = time.time() - t0
    print(res.stdout[-6000:])
    if res.returncode != 0:
        print("STDERR tail:", res.stderr[-4000:])
        raise SystemExit("run_cell failed")
    print("wall-clock: %.1f min" % (wall / 60.0))
    # one-line per-replication timing check (plan re-confirm of SPEEDUP on the
    # real workload): rows already on Drive.
    if os.path.exists(out):
        nrows = sum(1 for _ in open(out)) - 1
        print("rows on Drive so far:", max(nrows, 0))
    """))

code(CELLS_B, textwrap.dedent("""\
    import os
    if not os.path.exists(out):
        raise SystemExit("shard CSV not produced: " + out)
    import csv
    with open(out) as f:
        rows = list(csv.DictReader(f))
    # Emitted `seed` labels differ by mode: confirmatory 0.., pilot 1_000_001..
    emit_start = SEED_START if MODE != "pilot" else 1000000 + SEED_START + 1
    expect = SEED_END - SEED_START + 1
    got = set({int(r["seed"]) for r in rows if r.get("seed")})
    missing = [emit_start + i for i in range(expect) if (emit_start + i) not in got]
    print("expected rows in shard:", expect, "found unique seeds:", len(got))
    print("missing seeds:", missing if missing else "none")
    print("non-converged (converged_flag==0):", sum(1 for r in rows if r.get("converged_flag") == "0"))
    print("shard CSV:", out, os.path.getsize(out), "bytes")
    # Colab download fallback per the standing rule:
    try:
        from google.colab import files
        files.download(out)
        print("Downloaded shard CSV to local machine.")
    except Exception as e:
        print("(not on Colab / download skipped):", e)
    """))

# ---------------------------------------------------------------------------
# Round 0 pilot + P3-T5(e) calibration gate
# ---------------------------------------------------------------------------

md(CELLS_C, textwrap.dedent("""\
    ## Round 0 pilot + P3-T5(e) stochtree::bcf calibration gate

    Run **DGP1 n=500 pilot first** (the only cell the calibration gate needs).
    It is ~5 core-hours. If the calibrated-BCF configuration is wrong, you lose
    one notebook, not four. After an acceptable calibration result, launch the
    remaining pilots (DGP2/DGP3 n=500, DGP1 n=100) in their own shards.

    The gate compares the 50-pilot means of the calibrated `bcf_*` columns to
    McJames et al.'s published DGP1/n=500 Table 2. No confirmatory shard may run
    until all four BCF rows are inside their pass bands (plan P3-T5[e]).
    """))

code(CELLS_C, textwrap.dedent("""\
    import os
    # Round 0 pilot: DGP1, n=500, 50 independent-seed pilot replications.
    DRIVE_DIR = "/content/drive/MyDrive/TISCA_E3"
    DGP, N, MC_CORES, MODE = 1, 500, 2, "pilot"
    SEED_START, SEED_END = 0, 49
    os.makedirs(DRIVE_DIR, exist_ok=True)
    out = os.path.join(DRIVE_DIR, f"DGP1_n500_pilot.csv")
    import subprocess, time
    env = dict(os.environ)
    env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")
    env["TISCA_MVBCF_CPP"] = "/content/e3/MVBCF_Code.cpp"
    cmd = ["Rscript", "/content/e3/run_cell.R", str(DGP), str(N),
           str(SEED_START), str(SEED_END), "--out", out,
           "--cores", str(MC_CORES), "--mode", MODE]
    print("running:", " ".join(cmd))
    t0 = time.time()
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    print("wall-clock: %.1f min" % ((time.time() - t0) / 60.0))
    print(res.stdout[-4000:])
    if res.returncode != 0:
        print("STDERR:", res.stderr[-3000:])
        raise SystemExit("pilot failed")
    """))

code(CELLS_C, textwrap.dedent("""\
    import csv, statistics as st
    with open(out) as f:
        rows = list(csv.DictReader(f))
    n = len(rows)
    if n < 30:
        print("[FAIL] calibration needs >=50 pilot rows; got", n)
    f1 = lambda key: [float(r[key]) for r in rows if r.get(key) not in (None, "")]
    res = {}
    for key in ["bcf_pehe1", "bcf_pehe2", "bcf_cov951", "bcf_cov952"]:
        v = f1(key)
        res[key] = (st.mean(v), st.stdev(v) / (len(v) ** 0.5), len(v)) if v else (float("nan"), float("nan"), 0)
    print("=" * 72)
    print("P3-T5(e) calibration gate vs McJames et al. DGP1 n=500 Table 2")
    print("=" * 72)
    targets = {
        "bcf_pehe1": (9.63, 9.30, 10.00, "BCF PEHE Y1"),
        "bcf_pehe2": (9.96, 9.60, 10.30, "BCF PEHE Y2"),
        "bcf_cov951": (0.97, 0.95, 0.98, "BCF tau 95% cov Y1"),
        "bcf_cov952": (0.96, 0.94, 0.98, "BCF tau 95% cov Y2"),
    }
    ok = True
    for key, (mean, se, nn) in res.items():
        tgt, lo, hi, label = targets[key]
        inside = lo <= mean <= hi
        ok = ok and inside
        print(f"{label:>20}: mean={mean:.3f} se={se:.3f} target={tgt} band=[{lo},{hi}] -> {'OK' if inside else 'OUT OF BAND'}")
    print("-" * 72)
    print("VERDICT:", "PASS - proceed to confirmatory shards" if ok
          else "FAIL - diagnose before launching confirmatory (plan P3-T5[e])")
    """))

code(CELLS_C, textwrap.dedent("""\
    # P3-T5(f): compare fast_bart() vs mvbcf::run_mvbcf() on 10 identical seeds.
    # Placeholder: run both on the same 10 pilot seeds and compare PEHE/coverage.
    print("P3-T5(f) fast_bart-vs-mvbcf equivalence check: "
          "see CALIBRATION.md; record outcome there.")
    """))

# ---------------------------------------------------------------------------
# Assemble the two notebooks
# ---------------------------------------------------------------------------
def build(cells, kernel_py=True):
    return {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                      "language_info": {"name": "python"}},
        "cells": cells,
    }

ROOT = "/home/hugo_souto/Stuff/Submitted_Research/TISCA_paper/Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks"
os.makedirs(ROOT, exist_ok=True)

# E3 shard notebook: shared setup + shard run cells
nb_shard = build(CELLS_A + CELLS_B)
# E3 round-0 pilot + calibration: shared setup + pilot cells
nb_pilot = build(CELLS_A + CELLS_C)

def _sanity(cells, name):
    for i, c in enumerate(cells):
        src = "".join(c["source"])
        assert not any(bad in src for bad in ["'''", '"""']), f"{name} cell {i} has forbidden triple-quote"

_sanity(CELLS_A, "setup")
_sanity(CELLS_B, "shard")
_sanity(CELLS_C, "pilot")

p1 = os.path.join(ROOT, "E3_mvbcf_shard.ipynb")
p2 = os.path.join(ROOT, "E3_round0_pilot_calibration.ipynb")
with open(p1, "w") as f:
    json.dump(nb_shard, f, indent=1)
with open(p2, "w") as f:
    json.dump(nb_pilot, f, indent=1)
print("wrote", p1, len(CELLS_A) + len(CELLS_B), "cells")
print("wrote", p2, len(CELLS_A) + len(CELLS_C), "cells")
