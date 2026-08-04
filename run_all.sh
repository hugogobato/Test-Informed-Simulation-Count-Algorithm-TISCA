#!/usr/bin/env bash
# TISCA v2 end-to-end pipeline entry point (P0-T1 ACCEPTANCE).
#
# `run_all.sh` must reproduce every number in the paper from raw seeds, from a
# fresh clone and clean environments. Individual experiments are added here as
# they land in Phases 3-4. In its current state it validates the skeleton
# (structure sanity) only; experiment sub-targets are stubbed.
#
# Environment expectation:
#   - Python deps per environment.yml / requirements.txt (for the NumPy
#     outer-MC harness E1, the R-independent part).
#   - R packages from the P0-T4 library bundle on Google Colab, or
#     `Rscript env/install_R_dependencies.R` on a local machine. Once the P0-T4
#     bundle has been built, a real renv.lock is snapshotted from it and
#     committed, and `renv::restore()` becomes the preferred path.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "[run_all] TISCA v2 skeleton sanity check"

# 1) required top-level layout
for d in tisca experiments results figures docs notebooks legacy; do
  [[ -d "$d" ]] || { echo "missing dir: $d"; exit 1; }
done

# 2) the seed/RNG protocol is present
[[ -f docs/seed_rng_protocol.md ]] || { echo "missing docs/seed_rng_protocol.md"; exit 1; }

# 3) Python package imports (E1 outer-MC engine is the R-independent core)
if [[ -f tisca/python/tisca/__init__.py ]]; then
  PY="$(command -v python3 || command -v python || true)"
  if [[ -n "$PY" ]]; then
    PYTHONPATH="$ROOT/tisca/python:${PYTHONPATH:-}" \
      "$PY" -c "import tisca; print('[run_all] tisca python package imports')"
  else
    echo "[run_all] no python interpreter found; skipping package import check"
  fi
fi

# 4) E1 / P2-T3: harness acceptance gate (bivariate normal, rho=0, oracle sigma
#    -> Type I ~ 0.05 and power ~ 0.80, each within 2 MCSE). Suppresses pytest if
#    absent; otherwise runs the focused acceptance test at reduced R for speed.
if command -v python3 >/dev/null 2>&1; then
  echo "[run_all] E1 harness acceptance (P2-T3)"
  PYTHONPATH="$ROOT/tisca/python" \
    python3 experiments/E1_operating_characteristics/run_e1.py --acceptance \
    > results/E1_acceptance.json 2>&1 \
    && echo "[run_all] E1 acceptance PASS (see results/E1_acceptance.json)" \
    || { echo "[run_all] E1 acceptance FAILED"; cat results/E1_acceptance.json; exit 1; }
fi

echo "[run_all] skeleton OK. Experiment pipelines attach here as they land."
