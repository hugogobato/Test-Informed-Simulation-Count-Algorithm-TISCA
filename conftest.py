import os
import sys

_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
_PKG_DIR = os.path.join(_REPO_ROOT, "tisca", "python")

if _PKG_DIR not in sys.path:
    sys.path.insert(0, _PKG_DIR)
