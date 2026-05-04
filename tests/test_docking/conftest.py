"""
conftest.py — test_docking package
=====================================
Runs during pytest startup before any test module is imported.
Stubs out module-level imports that test files depend on.
"""

import sys
import types


def _stub(name: str, **attrs) -> types.ModuleType:
    """Insert a lightweight module stub into sys.modules (no-op if present)."""
    if name not in sys.modules:
        mod = types.ModuleType(name)
        for k, v in attrs.items():
            setattr(mod, k, v)
        sys.modules[name] = mod
    return sys.modules[name]


# ---------------------------------------------------------------------------
# tqdm  (imported at module level in np_docking.py)
# ---------------------------------------------------------------------------
_stub("tqdm", tqdm=lambda x, **kw: x)
