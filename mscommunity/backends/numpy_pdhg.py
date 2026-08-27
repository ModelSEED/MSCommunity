# -*- coding: utf-8 -*-
"""NumPy / SciPy PDHG backend.

CPU-only. Serves both as a reference for algorithm correctness (the GPU
backends use the same `run_pdhg_batch`, just with different array libs)
and as a fallback when no GPU is available.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from mscommunity.batched_lp import BatchedLPSolver, register_batched_solver

from ._pdhg import ArrayBackend, PDHGConfig, run_pdhg_batch


def _numpy_backend() -> ArrayBackend:
    return ArrayBackend(
        xp=np,
        to_device=lambda a: np.asarray(a),
        to_host=lambda a: np.asarray(a),
        sparse_from_scipy=lambda S: S.tocsr(),
        spmm=lambda A, X: A @ X,
        spmm_T=lambda A, Y: A.T @ Y,
    )


class NumpyPDHGBatchedLPSolver(BatchedLPSolver):
    """First-order PDHG batched LP solver on CPU NumPy.

    Same algorithm as the JAX / CuPy backends — useful as the
    correctness reference and on machines without a GPU.
    """

    name = "numpy-pdhg"

    def __init__(self, max_iters: int = 5000, tol: float = 1e-4,
                 check_every: int = 50, power_iters: int = 30,
                 step_safety: float = 0.9, seed: int = 0):
        self.config = PDHGConfig(
            max_iters=max_iters, tol=tol, check_every=check_every,
            power_iters=power_iters, step_safety=step_safety, seed=seed,
        )

    def solve(self, problem, instances):
        return run_pdhg_batch(_numpy_backend(), problem, instances, self.config)


register_batched_solver("numpy-pdhg", NumpyPDHGBatchedLPSolver)
