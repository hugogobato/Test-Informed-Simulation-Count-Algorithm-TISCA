#!/usr/bin/env python3
"""Generate the P0-T4 Colab notebook that builds the R library bundle tarball.

Each notebook cell's source is stored as a Python list of lines. We build cell
source via ordinary triple-quoted strings that NEVER contain a nested triple
quote. R scripts that the notebook itself writes are represented as plain
single-line-concatenated strings (no triple quotes inside the cell source),
so the generator never collides with R's own string syntax.

Output: Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks/P0T4_build_rlib_bundle.ipynb
"""
import json, os, textwrap

CELLS = []

def md(src):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)})

def code(src):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.splitlines(keepends=True)})

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
# P0-T4 — Build the R library bundle for the MVBCF case-study re-run

**Task:** `Test-Informed-Simulation-Count-Algorithm-TISCA` revision, task P0-T4.
**Paper:** Souto & Louzada Neto, "Beyond Arbitrary Replications...".
**Purpose:** compile the heavy R packages **once**, archive them as
`tisca_rlib.tar.gz`, and publish that single tarball so that ~33 worker sessions
restore it in ~60 seconds instead of each spending 10-20 minutes compiling from
source (and risking a transient install failure).

**You run this notebook exactly once.** Everything from Cell 3 onward runs on
the Colab runtime and writes the tarball; then you publish it (a Drive share
link or a GitHub release asset on the TISCA repo) and paste the direct URL into
Cell 9 so it can restore-and-verify that URL (the ACCEPTANCE test).

---
## What the bundle must contain (plan P0-T4 / §2)

Authoritative dependency list from `MVBCF-Paper-main/GitHub_DGP1.R` (the
original authors' script), `Paper_Experiments/DGP1.ipynb`, and the `tisca` v2
package plan:

| Package | Source | Reason |
|---|---|---|
| `stochtree` | CRAN | the calibrated BCF benchmark (fast) |
| `dbarts` | CRAN | BART / MVBART |
| `bartCause` | CRAN | treatment-effect BART |
| `skewBART` | GitHub `Seungha-Um/skewBART` | a v1 method |
| `mvbcf` | GitHub `Nathan-McJames/mvbcf` | the MVBCF method column |
| `mvtnorm` | CRAN | DGP sampling |
| `scoringRules` | CRAN | CRPS (the original authors' metric) |
| `matrixStats` | CRAN | per-replication row/col stats |
| `progress` | CRAN | driver progress bars |
| `MCS` | CRAN | cross-check oracle for the v2 MCS layer |

Neither the CRAN nor the GitHub packages are pinned *at install time* — they
resolve to whatever is current when this notebook runs. **Reproducibility comes
from the artefact, not from the install command:** the bundle is built once and
all ~33 worker sessions restore that same tarball, so every shard uses
byte-identical package builds. What has to be recorded is therefore what was
actually resolved, and Cell 6 captures it three ways — `sessionInfo.txt`,
`DEPENDENCIES.csv` (with the `RemoteSha` of each GitHub package) and a real
`renv.lock` — which is what the paper reports (IJDA #14a) and what lets a reader
rebuild the same environment later.
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 1 — Hardware / environment introspection
Record the facts the paper must now report (IJDA #14a): CPU model, `nproc`,
R version, OS. Nothing here gates anything downstream.
"""))

code(textwrap.dedent("""\
import platform, subprocess
def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()
print("hostname :", platform.node() or "n/a")
print("machine  :", platform.machine())
print("os       :", platform.platform())
print("nproc    :", sh("nproc"))
print("cpu      :", sh("grep -m1 -E 'model name' /proc/cpuinfo") or "no /proc/cpuinfo")
print("ram:")
print(sh("free -g") or "free unavailable")
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 2 — Install R
`r-base` from the Debian repos is sufficient. Takes ~1-2 min.
"""))

code(textwrap.dedent("""\
!apt-get -qq update > /dev/null
!apt-get -qq install -y --no-install-recommends r-base r-base-dev gfortran > /dev/null 2>&1
!R --version | head -1
print("R install done")
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 3 — Create the private build library
We install into `/content/rlib` (not the system site-library) so we can `tar`
exactly that directory and so the worker can force `.libPaths()` to use the
bundled versions. We also record the resolved versions in
`/content/rlib_DEPENDENCIES.csv` + `/content/rlib_sessionInfo.txt` for the paper.
"""))

code(textwrap.dedent("""\
import os, subprocess
RLIB = "/content/rlib"
os.makedirs(RLIB, exist_ok=True)

# R settings block (single-line concatenation; no triple-quote nesting needed)
R_LIBPATH = '.libPaths(c("' + RLIB + '", .libPaths()))'
R_OPTIONS = 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 2)'

manifest = (
    R_LIBPATH + "\\n" +
    R_OPTIONS + "\\n" +
    'ip <- installed.packages(lib.loc = "' + RLIB + '")\\n' +
    'd <- data.frame(pkg=ip[,"Package"], version=ip[,"Version"], stringsAsFactors=FALSE)\\n' +
    'write.csv(d[order(d$pkg),], "/content/rlib_DEPENDENCIES.csv", row.names=FALSE)\\n' +
    'sink("/content/rlib_sessionInfo.txt")\\n' +
    'sessionInfo()\\n' +
    'sink()\\n' +
    'cat("MANIFEST_WRITTEN\\\\n")\\n'
)
with open("/content/write_manifest_placeholder.R", "w") as f:
    f.write("# written later after install")
print("RLIB ready:", RLIB, "exists:", os.path.isdir(RLIB))
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 4 — Install the CRAN packages into `/content/rlib`
Order matters: the GitHub packages depend on several of these. We compile with
`Ncpus = 2` (Colab has 2 vCPUs). Takes ~10-15 min.
"""))

code(textwrap.dedent("""\
import subprocess, os, time

RLIB = "/content/rlib"
LB    = '"' + RLIB + '"'                      # R string for the lib path
install = (
    '.libPaths(c(' + LB + ', .libPaths()))\\n' +
    'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 2)\\n' +
    'pkgs <- c("stochtree","dbarts","bartCause","mvtnorm",\\n' +
    '          "scoringRules","matrixStats","progress","MCS")\\n' +
    'ok <- character(0)\\n' +
    'for (p in pkgs) {\\n' +
    '  if (!requireNamespace(p, quietly = TRUE)) {\\n' +
    '    cat("== installing", p, "\\\\n")\\n' +
    '    try(install.packages(p, lib=' + LB + ', quiet = TRUE))\\n' +
    '  }\\n' +
    '  if (requireNamespace(p, quietly = TRUE)) ok <- c(ok, p)\\n' +
    '}\\n' +
    'cat("INSTALLED_OK:", paste(ok, collapse=","), "\\\\n")\\n'
)
with open("/content/install_cran.R", "w") as f:
    f.write(install)
t0 = time.time()
res = subprocess.run(["Rscript", "/content/install_cran.R"], capture_output=True, text=True)
el = time.time() - t0
print(res.stdout[-4000:])
if "INSTALLED_OK:" not in res.stdout:
    print("STDERR (tail):", res.stderr[-3000:])
    raise SystemExit("CRAN install did not report INSTALLED_OK")
print("CRAN install finished in %.1f min" % (el / 60.0))
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 5 — Install the GitHub packages (`skewBART`, `mvbcf`)
Both need `remotes`. They depend on CRAN packages installed above, so they now
compile quickly. Pinned by commit via `remotes::install_github`.
"""))

code(textwrap.dedent("""\
import subprocess, os

RLIB = "/content/rlib"
LB    = '"' + RLIB + '"'
gh = (
    '.libPaths(c(' + LB + ', .libPaths()))\\n' +
    'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 2)\\n' +
    'if (!requireNamespace("remotes", quietly = TRUE))\\n' +
    '  install.packages("remotes", lib=' + LB + ', quiet = TRUE)\\n' +
    'library(remotes)\\n' +
    'for (repo in c("Seungha-Um/skewBART", "Nathan-McJames/mvbcf")) {\\n' +
    '  name <- tail(strsplit(repo, "/")[[1]], 1)\\n' +
    '  if (requireNamespace(name, quietly = TRUE)) { cat("already have", name, "\\\\n"); next }\\n' +
    '  cat("== installing", repo, "\\\\n")\\n' +
    '  install_github(repo, lib=' + LB + ', quiet = TRUE, upgrade = "never", build_vignettes = FALSE)\\n' +
    '}\\n' +
    'cat("GH_DONE\\\\n")\\n'
)
with open("/content/install_gh.R", "w") as f:
    f.write(gh)
res = subprocess.run(["Rscript", "/content/install_gh.R"], capture_output=True,
                     text=True, timeout=3600)
print(res.stdout[-3000:])
if "GH_DONE" not in res.stdout:
    print("STDERR (tail):", res.stderr[-3000:])
    raise SystemExit("GitHub install did not complete")
print("GitHub packages installed")
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 6 — Record the exact environment + versions into the bundle
The paper reports exact package versions (IJDA #14a). We write three manifests
next to `rlib` so they are archived with it:

- `DEPENDENCIES.csv` — package, version, and the `RemoteSha` for the two GitHub
  packages, which is the only thing that identifies *which* `skewBART`/`mvbcf`
  was used;
- `sessionInfo.txt` — the verbatim R session;
- `renv.lock` — a **real** lockfile produced by `renv::snapshot()` over the
  built library. This is the file to commit to the repository root, replacing
  `env/install_R_dependencies.R` as the source of truth. A lockfile written by
  hand without resolved versions is not valid and fails `renv::restore()`,
  which is why the repository carries the installer script until this cell has
  run.
"""))

code(textwrap.dedent("""\
import os, subprocess

RLIB = "/content/rlib"
LB    = '"' + RLIB + '"'
manifest = (
    '.libPaths(c(' + LB + ', .libPaths()))\\n' +
    'options(repos = c(CRAN = "https://cloud.r-project.org"))\\n' +
    'ip <- installed.packages(lib.loc = ' + LB + ')\\n' +
    'pk <- ip[, "Package"]\\n' +
    'sha <- vapply(pk, function(p) {\\n' +
    '  d <- packageDescription(p, lib.loc = ' + LB + ')\\n' +
    '  if (!is.null(d$RemoteSha)) d$RemoteSha else NA_character_ }, character(1))\\n' +
    'd <- data.frame(pkg=pk, version=ip[,"Version"], remote_sha=unname(sha),\\n' +
    '                stringsAsFactors=FALSE)\\n' +
    'write.csv(d[order(d$pkg),], "/content/rlib_DEPENDENCIES.csv", row.names=FALSE)\\n' +
    'sink("/content/rlib_sessionInfo.txt")\\n' +
    'sessionInfo()\\n' +
    'sink()\\n' +
    '# a genuine lockfile, snapshotted from what was actually installed\\n' +
    'if (!requireNamespace("renv", quietly=TRUE))\\n' +
    '  install.packages("renv", lib=' + LB + ', quiet=TRUE)\\n' +
    'ok <- tryCatch({\\n' +
    '  renv::snapshot(library = ' + LB + ', lockfile = "/content/rlib_renv.lock",\\n' +
    '                 type = "all", prompt = FALSE); TRUE\\n' +
    '}, error = function(e) { cat("renv::snapshot failed:", conditionMessage(e), "\\\\n"); FALSE })\\n' +
    'cat("RENV_LOCK_WRITTEN:", ok, "\\\\n")\\n' +
    'cat("MANIFEST_WRITTEN\\\\n")\\n'
)
with open("/content/write_manifest.R", "w") as f:
    f.write(manifest)
res = subprocess.run(["Rscript", "/content/write_manifest.R"], capture_output=True, text=True)
print(res.stdout[-2000:])
if "MANIFEST_WRITTEN" not in res.stdout:
    print("STDERR:", res.stderr[-1500:])
    raise SystemExit("manifest not written")
print("--- DEPENDENCIES.csv ---")
print(open("/content/rlib_DEPENDENCIES.csv").read())
print("--- sessionInfo.txt (head) ---")
print(open("/content/rlib_sessionInfo.txt").read()[:1200])
if os.path.exists("/content/rlib_renv.lock"):
    print("--- renv.lock written; commit this to the repo root ---")
else:
    print("!! renv.lock NOT written: keep env/install_R_dependencies.R as the "
          "source of truth and say so in the paper.")
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 7 — Build the tarball
Archive layout: a single top-level folder `tisca_rlib/` containing `rlib/`,
`DEPENDENCIES.csv`, `sessionInfo.txt`, and `BUNDLE_META.json`. A worker untars
once and points `.libPaths()` at `<dest>/tisca_rlib/rlib`. We write a SHA256 so
workers can verify integrity before loading.
"""))

code(textwrap.dedent("""\
import subprocess, os, hashlib, time, json, platform

RLIB = "/content/rlib"
BUNDLE_DIR = "/content/tisca_rlib"
os.makedirs(BUNDLE_DIR, exist_ok=True)
TAR = "/content/tisca_rlib.tar.gz"

meta = {
    "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "os": platform.platform(),
    "cpu": subprocess.run(["grep","-m1","-E","model name","/proc/cpuinfo"],
                          capture_output=True, text=True).stdout.strip() or "n/a",
    "nproc": subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip(),
    "R": subprocess.run(["Rscript","-e","cat(R.version.string)"],
                        capture_output=True, text=True).stdout.strip(),
}
with open(os.path.join(BUNDLE_DIR,"BUNDLE_META.json"), "w") as f:
    json.dump(meta, f, indent=2)

# move rlib into the bundle dir (fresh copy to avoid clobbering source)
subprocess.run(["cp","-r",RLIB, BUNDLE_DIR+"/rlib"])
subprocess.run(["cp","/content/rlib_DEPENDENCIES.csv", BUNDLE_DIR+"/DEPENDENCIES.csv"])
subprocess.run(["cp","/content/rlib_sessionInfo.txt", BUNDLE_DIR+"/sessionInfo.txt"])
if os.path.exists("/content/rlib_renv.lock"):
    subprocess.run(["cp","/content/rlib_renv.lock", BUNDLE_DIR+"/renv.lock"])

res = subprocess.run(["tar","czf",TAR,"-C","/content","tisca_rlib"],
                     capture_output=True, text=True)
print(res.stderr[-1000:] or "tar OK")
size = os.path.getsize(TAR)
h = hashlib.sha256(open(TAR,"rb").read()).hexdigest()
print("tarball size : %.1f MB" % (size/1e6))
print("SHA256       :", h)
with open("/content/tisca_rlib.sha256","w") as f:
    f.write(h+"\\n")
print("top-level entries of the bundle (first 12):")
for line in subprocess.run(["tar","tzf",TAR], capture_output=True, text=True).stdout.splitlines()[:12]:
    print("  ", line)
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 8 — Publish the tarball

**Pick one.** Recommended Option A (plan P0-T4): *either* a Drive share link,
*or* a GitHub **release asset** on
`github.com/hugogobato/Test-Informed-Simulation-Count-Algorithm-TISCA` (best for
the paper's Data Availability statement). Drop `tisca_rlib.tar.gz` there and
grab the direct-download URL.

Keep the file with `files.download` below, then paste the direct URL into Cell 9.
"""))

code(textwrap.dedent("""\
import os
print("tarball:", "/content/tisca_rlib.tar.gz",
      "%.1f MB" % (os.path.getsize("/content/tisca_rlib.tar.gz")/1e6))

# Colab download fallback (per the user's standing Colab rule):
try:
    from google.colab import files
    files.download("/content/tisca_rlib.tar.gz")
    files.download("/content/tisca_rlib.sha256")
    print("Downloaded tisca_rlib.tar.gz and sha256 to your local machine.")
except Exception as e:
    print("(not on Colab / download skipped):", e)
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Cell 9 — Freshness test (the P0-T4 ACCEPTANCE criterion)

**ACCEPTANCE:** a fresh notebook on an account that did **not** build the bundle
restores it and loads `stochtree`, `dbarts`, `bartCause`, `skewBART` and
`mvtnorm` in under 90 seconds, with zero compilation.

**Run this cell in a FRESH runtime on an account that did not build the bundle.**
That is the whole point of the test, and it means the cell may not read anything
left behind by the build: paste both the `URL` and the `SHA256` printed by
Cell 7 into the two constants below. It downloads, verifies the SHA256, restores
into `/content/tisca_rlib/rlib`, and times `requireNamespace` on the gate
packages.
"""))

code(textwrap.dedent("""\
import subprocess, os, time, urllib.request, hashlib

# Both constants come from Cell 7's output. A fresh account has neither the
# tarball nor the .sha256 file on disk, so they must be pasted, not read.
URL    = ""   # <-- direct download link to tisca_rlib.tar.gz
SHA256 = ""   # <-- the SHA256 Cell 7 printed

if not URL or not SHA256:
    raise SystemExit("Paste both URL and SHA256 from Cell 7 before running the "
                     "acceptance test.")

DEST = "/content"
LIBDIR = "/content/tisca_rlib/rlib"
DL = "/content/_dl_tisca_rlib.tar.gz"

# -- always download: reusing a locally built tarball would not test the thing
#    the acceptance criterion is about (a fresh session restoring from the URL).
t0 = time.time()
urllib.request.urlretrieve(URL, DL)
dl_s = time.time() - t0
print("downloaded in %.1fs (%.1f MB)" % (dl_s, os.path.getsize(DL)/1e6))

# -- integrity --
h = hashlib.sha256(open(DL, "rb").read()).hexdigest()
assert h == SHA256, "SHA mismatch %s vs %s" % (h, SHA256)
print("SHA256 OK")

# -- restore: exactly what a worker shard will run --
t1 = time.time()
if os.path.exists(LIBDIR):
    subprocess.run(["rm","-rf","/content/tisca_rlib"])
subprocess.run(["tar","xzf",DL,"-C",DEST])
assert os.path.isdir(LIBDIR), "restore failed: %s missing" % LIBDIR
print("restored in %.1fs: %s" % (time.time()-t1, LIBDIR))
print("download + restore so far: %.1fs of the 90 s budget" % (time.time()-t0))

check = (
    '.libPaths(c("' + LIBDIR + '", .libPaths()))\\n' +
    'gate <- c("stochtree","dbarts","bartCause","skewBART","mvtnorm",\\n' +
    '          "mvbcf","scoringRules","matrixStats","MCS")\\n' +
    't0 <- Sys.time()\\n' +
    'ok <- sapply(gate, function(p) requireNamespace(p, quietly=TRUE))\\n' +
    'cat("------- gate packages -------\\\\n")\\n' +
    'print(data.frame(pkg=gate, loadable=ok))\\n' +
    'cat("elapsed s:", round(as.numeric(Sys.time()-t0), 2), "\\\\n")\\n' +
    'if (any(!ok)) { cat("MISSING:", paste(gate[!ok], collapse=", "), "\\\\n"); quit(status=2) }\\n'
)
with open("/content/check_bundle.R","w") as f:
    f.write(check)
t0 = time.time()
res = subprocess.run(["Rscript","/content/check_bundle.R"], capture_output=True, text=True)
wall = time.time() - t0
print(res.stdout[-2000:])
if res.returncode != 0:
    print("STDERR:", res.stderr[-1500:])
    print("[FAIL] bundle freshness test failed")
else:
    print("[PASS] gate packages loadable from restored bundle; wall %.1f s" % wall)
"""))

# ---------------------------------------------------------------------------
md(textwrap.dedent("""\
## Acceptance checklist — P0-T4

- [ ] Cell 7 produced `tisca_rlib.tar.gz` + `tisca_rlib.sha256`.
- [ ] You published the tarball to a **Drive share link** or a **GitHub release
      asset** on `Test-Informed-Simulation-Count-Algorithm-TISCA`.
- [ ] Cell 9, run in a **fresh runtime on an account that did not build the
      bundle**, with `URL` and `SHA256` pasted from Cell 7, downloaded and
      restored the bundle and loaded all gate packages in **under 90 s**.
      Running Cell 9 in the build session does not count — it would only prove
      that a local file untars.
- [ ] You recorded the direct `URL`, the explicit `sessionInfo()` (captured in
      the bundle as `sessionInfo.txt`), and any install deviations in
      `experiments/E3_mvbcf_casestudy/CALIBRATION.md`.
- [ ] `renv.lock` from Cell 6 is committed to the repository root, and
      `env/install_R_dependencies.R` now points at it.
- [ ] `DEPENDENCIES.csv` carries a non-`NA` `remote_sha` for `skewBART` and
      `mvbcf` — without it, "which version of the GitHub package" is
      unanswerable and IJDA #14a is not satisfied.

Once this passes, the worker shards (Phase 1+) start with a short
`restore_lib()` helper: `wget <URL> && tar xzf && .libPaths(<dest>/tisca_rlib/rlib)`.
"""))

# ---------------------------------------------------------------------------
nb = {
    "nbformat": 4,
    "nbformat_minor": 0,
    "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                 "language_info": {"name": "python"}},
    "cells": CELLS,
}

out = "/home/hugo_souto/Stuff/Submitted_Research/TISCA_paper/Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks/P0T4_build_rlib_bundle.ipynb"
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    json.dump(nb, f, indent=1)

# sanity: every cell source must not contain a stray triple-double-quote
for i, c in enumerate(CELLS):
    src = "".join(c["source"])
    assert not any(bad in src for bad in ["'''", '"""']), f"cell {i} has forbidden triple-quote"

print("wrote", out, len(CELLS), "cells")
