"""TISCA v2 - Test-Informed Simulation Count Algorithm (Python reference).

The Python reference implementation is split across the Phase 2 software tasks:

* P2-T1 (this deliverable): ``estimands``, ``planning``, ``inference``,
  ``multiplicity``, ``mcs``, ``procedure``, ``validate``.
* P2-T3 (parallel): ``outermc`` - the outer-MC operating-characteristics engine,
  developed concurrently.

A partial working tree (any phase landing before the others) must never break
``import tisca``, so available submodules are discovered gracefully rather than
imported eagerly. ``from tisca import planning`` (etc.) remains valid for every
module that exists, whether or not it is listed here.
"""

from __future__ import annotations

_SUBMODULES = (
    "estimands",
    "planning",
    "inference",
    "multiplicity",
    "mcs",
    "procedure",
    "validate",
    "outermc",
)

__all__ = []
_imported = []
for _name in _SUBMODULES:
    try:
        __import__(f"{__name__}.{_name}")
        _imported.append(_name)
    except ImportError:
        # The submodule (or one of its transitive dependencies, e.g. the R port
        # or the outer-MC harness) is not available in this checkout; skip it.
        pass

__version__ = "2.0.0"
