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
#   - R packages from the P0-T4 library bundle on Google Colab, or renv on a
#     local machine (see docs/seed_rng_protocol.md and README Repro section).

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

echo "[run_all] skeleton OK. Experiment pipelines attach here as they land."
