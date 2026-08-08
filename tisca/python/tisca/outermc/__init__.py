"""TISCA v2 outer-Monte-Carlo harness (P2-T3).

A vectorised NumPy engine that repeats the *entire* simulation-design procedure
(pilot -> planning -> confirmatory run -> final adjusted test) ``R`` times and
returns the unconditional operating characteristics. This is the engine behind
experiment E1 (``experiments/E1_operating_characteristics``).

The harness is intentionally pure NumPy/SciPy so it runs on Google Colab without
the R library bundle (see REVISION_PLAN.md P3-T2 and notebooks/E1*.ipynb).

Public surface:
  run_e1(config, ...)            top-level entry point
  make_design(name, config)      design object for one of D1..D6
  families                      supported loss-distribution samplers
  summarize_ocs                collapse raw per-rep results into summary rows
"""

from . import families
from . import joint_families
from .designs import make_design
from .engine import run_e1
from .joint_engine import run_joint_cell
from .oc import summarize_ocs

__all__ = [
    "families", "joint_families", "make_design", "run_e1", "run_joint_cell",
    "summarize_ocs",
]
