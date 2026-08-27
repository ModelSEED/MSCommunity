# -*- coding: utf-8 -*-
"""Batched LP backends. Importing this package auto-registers every
backend whose dependency is installed; missing ones fail silently so
the package import never breaks.

Backends:
    numpy-pdhg   CPU NumPy/SciPy first-order PDHG (always available)
    jax-pdhg     CUDA / Metal / ROCm / TPU / CPU via JAX
    cupy-pdhg    NVIDIA CUDA via CuPy + cuSPARSE
    pdlp         Google OR-Tools PDLP (CPU; high accuracy)

After import, list registered backends with:
    from mscommunity.batched_lp import _SOLVERS
    print(sorted(_SOLVERS))
"""
from __future__ import annotations

# numpy-pdhg has no optional dep beyond scipy and is always loaded.
from . import numpy_pdhg  # noqa: F401

for _mod_name in ("jax_pdhg", "cupy_pdhg", "ortools_pdlp"):
    try:
        __import__(f"mscommunity.backends.{_mod_name}")
    except Exception:
        # The backend module itself records the import failure and only
        # registers if its dependency loaded. Swallow the import error so
        # the package stays usable on systems without the optional dep.
        pass
