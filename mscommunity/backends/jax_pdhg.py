# -*- coding: utf-8 -*-
"""JAX PDHG backend — portable across CUDA, Metal, ROCm, TPU, CPU.

JAX picks the device based on its install: `jax[cuda12]` for NVIDIA,
`jax-metal` for Apple Silicon, `jax[rocm]` for AMD, plain `jax` for CPU.
The same Python code runs unchanged; this module just transfers the
shared `S` and the (B, n) batched bounds / objective vectors to whichever
device JAX selects.

The PDHG inner loop is `jit`-compiled; across LP instances the
per-instance arrays are laid out with the batch axis as the trailing
dimension so column-wise broadcasting maps cleanly onto SpMM, which is
what most accelerators are tuned for.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from mscommunity.batched_lp import BatchedLPSolver, register_batched_solver

from ._pdhg import ArrayBackend, PDHGConfig, run_pdhg_batch

try:
    import jax
    import jax.numpy as jnp
    from jax.experimental.sparse import BCOO
    _JAX_OK = True
    _JAX_ERR = None
except Exception as _exc:  # pragma: no cover — import guard
    jax = None
    jnp = None
    BCOO = None
    _JAX_OK = False
    _JAX_ERR = _exc


def _jax_backend(device=None) -> ArrayBackend:
    if not _JAX_OK:
        raise ImportError(
            f"jax-pdhg backend requires `jax` (install jax[cuda12], "
            f"jax-metal, jax[rocm], or plain jax for CPU). Import failed: {_JAX_ERR}"
        )

    def _to_dev(a):
        arr = jnp.asarray(a)
        return jax.device_put(arr, device) if device is not None else arr

    def _to_host(a):
        return np.asarray(a)

    def _sparse(S):
        # BCOO is JAX's native sparse format; constructed from scipy COO.
        S_coo = S.tocoo()
        return BCOO(
            (jnp.asarray(S_coo.data), jnp.stack([
                jnp.asarray(S_coo.row), jnp.asarray(S_coo.col)
            ], axis=1)),
            shape=S.shape,
        )

    def _spmm(A, X):
        return A @ X

    def _spmm_T(A, Y):
        return A.T @ Y

    return ArrayBackend(
        xp=jnp, to_device=_to_dev, to_host=_to_host,
        sparse_from_scipy=_sparse, spmm=_spmm, spmm_T=_spmm_T,
    )


class JAXPDHGBatchedLPSolver(BatchedLPSolver):
    """First-order PDHG batched LP solver via JAX.

    Device is whatever JAX picks; pass `device="cpu"` / `"gpu"` /
    `"metal"` etc. or a `jax.Device` to pin explicitly. See JAX docs for
    backend-selection details.
    """

    name = "jax-pdhg"

    def __init__(self, max_iters: int = 5000, tol: float = 1e-4,
                 check_every: int = 50, power_iters: int = 30,
                 step_safety: float = 0.9, seed: int = 0, device=None):
        if not _JAX_OK:
            raise ImportError(
                f"jax-pdhg backend requires `jax`. Import failed: {_JAX_ERR}"
            )
        self.config = PDHGConfig(
            max_iters=max_iters, tol=tol, check_every=check_every,
            power_iters=power_iters, step_safety=step_safety, seed=seed,
        )
        if isinstance(device, str):
            device = jax.devices(device)[0]
        self.device = device

    def solve(self, problem, instances):
        return run_pdhg_batch(
            _jax_backend(self.device), problem, instances, self.config
        )


if _JAX_OK:
    register_batched_solver("jax-pdhg", JAXPDHGBatchedLPSolver)
