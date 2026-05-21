# -*- coding: utf-8 -*-
"""CuPy PDHG backend — NVIDIA CUDA only.

Uses `cupy.sparse.csr_matrix` (which wraps cuSPARSE) for the SpMM in the
PDHG inner loop. Same algorithm as the NumPy / JAX backends; the only
difference is array library + memory residency.

Pinned to CUDA: CuPy does not target Metal / ROCm. For cross-vendor
portability use the jax-pdhg backend; for NVIDIA-specific deployments
this is typically faster because cuSPARSE is tightly tuned for it.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from mscommunity.batched_lp import BatchedLPSolver, register_batched_solver

from ._pdhg import ArrayBackend, PDHGConfig, run_pdhg_batch

try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    _CUPY_OK = True
    _CUPY_ERR = None
except Exception as _exc:  # pragma: no cover — import guard
    cp = None
    csp = None
    _CUPY_OK = False
    _CUPY_ERR = _exc


def _cupy_backend(device_id: int = 0) -> ArrayBackend:
    if not _CUPY_OK:
        raise ImportError(
            f"cupy-pdhg backend requires `cupy` (install cupy-cuda12x or "
            f"matching CUDA build). Import failed: {_CUPY_ERR}"
        )

    def _to_dev(a):
        with cp.cuda.Device(device_id):
            return cp.asarray(a)

    def _to_host(a):
        return cp.asnumpy(a)

    def _sparse(S):
        with cp.cuda.Device(device_id):
            return csp.csr_matrix(S.tocsr())

    def _spmm(A, X):
        return A @ X

    def _spmm_T(A, Y):
        return A.T @ Y

    return ArrayBackend(
        xp=cp, to_device=_to_dev, to_host=_to_host,
        sparse_from_scipy=_sparse, spmm=_spmm, spmm_T=_spmm_T,
    )


class CuPyPDHGBatchedLPSolver(BatchedLPSolver):
    """First-order PDHG batched LP solver via CuPy on NVIDIA CUDA."""

    name = "cupy-pdhg"

    def __init__(self, max_iters: int = 5000, tol: float = 1e-4,
                 check_every: int = 50, power_iters: int = 30,
                 step_safety: float = 0.9, seed: int = 0, device: int = 0):
        if not _CUPY_OK:
            raise ImportError(
                f"cupy-pdhg backend requires `cupy`. Import failed: {_CUPY_ERR}"
            )
        self.config = PDHGConfig(
            max_iters=max_iters, tol=tol, check_every=check_every,
            power_iters=power_iters, step_safety=step_safety, seed=seed,
        )
        self.device = int(device)

    def solve(self, problem, instances):
        return run_pdhg_batch(
            _cupy_backend(self.device), problem, instances, self.config
        )


if _CUPY_OK:
    register_batched_solver("cupy-pdhg", CuPyPDHGBatchedLPSolver)
