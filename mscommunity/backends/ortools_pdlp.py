# -*- coding: utf-8 -*-
"""OR-Tools PDLP backend — serial batch, production-grade accuracy.

Google's PDLP is PDHG plus the heuristics (presolve, restart, adaptive
step size, ratio tests, polishing) that make first-order LP solving
practical at high accuracy. It runs on CPU today via the pip wheel.
GPU-PDLP (`cuPDLP`, `cuPDLP-c`) is a separate effort and ships outside
OR-Tools; if available, swap `_solve_one` to dispatch through it.

This backend solves one LP per instance (no cross-instance batching) but
each solve is much higher accuracy than the PDHG backends and the C++
implementation makes single-LP wall time competitive on real FBA models.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import List, Sequence

import numpy as np
import scipy.sparse as sp

from mscommunity.batched_lp import (
    BatchedLPSolver,
    BatchedSolution,
    register_batched_solver,
)

from ._pdhg import assemble_batch, stack_constraints

try:
    from ortools.pdlp import solvers_pb2
    from ortools.pdlp import solve_log_pb2
    from ortools.pdlp.python import pdlp
    _PDLP_OK = True
    _PDLP_ERR = None
except Exception as _exc:  # pragma: no cover — import guard
    pdlp = None
    solvers_pb2 = None
    solve_log_pb2 = None
    _PDLP_OK = False
    _PDLP_ERR = _exc


def _build_qp(problem, c_row: np.ndarray, lb_row: np.ndarray, ub_row: np.ndarray):
    """Pack one instance into PDLP's QuadraticProgram (LP = QP with no Q)."""
    A_full, full_lb, full_ub, m, e = stack_constraints(problem)
    qp = pdlp.QuadraticProgram()
    # PDLP minimizes c^T x; we already converted to that sign convention in
    # assemble_batch, so flip again here: PDLP's objective_vector is the
    # *minimization* objective.
    qp.objective_vector = (-c_row).astype(np.float64)
    qp.variable_lower_bounds = lb_row.astype(np.float64)
    qp.variable_upper_bounds = ub_row.astype(np.float64)
    A_csc = A_full.tocsc().astype(np.float64)
    qp.constraint_matrix = A_csc
    qp.constraint_lower_bounds = full_lb.astype(np.float64)
    qp.constraint_upper_bounds = full_ub.astype(np.float64)
    return qp


class ORToolsPDLPBatchedLPSolver(BatchedLPSolver):
    """Wrap OR-Tools' PDLP as a batched solver.

    Higher accuracy than the PDHG backends; runs each instance serially
    (or threaded). For GPU PDLP, swap in `cuPDLP` at the `_solve_one`
    seam — the call-site contract is identical.
    """

    name = "pdlp"

    def __init__(self, workers: int = 1,
                 termination_relative_gap: float = 1e-6,
                 termination_iteration_limit: int = 50_000,
                 verbosity_level: int = 0):
        if not _PDLP_OK:
            raise ImportError(
                f"pdlp backend requires `ortools` (`pip install ortools`)."
                f" Import failed: {_PDLP_ERR}"
            )
        self.workers = max(1, int(workers))
        self.params = solvers_pb2.PrimalDualHybridGradientParams()
        # OR-Tools >= 9.x exposes eps_optimal_relative as a scalar float field
        # on termination_criteria (older builds nested it under a message with
        # a `.gap` attribute). Set the scalar directly.
        self.params.termination_criteria.eps_optimal_relative = (
            termination_relative_gap
        )
        self.params.termination_criteria.iteration_limit = (
            termination_iteration_limit
        )
        self.params.verbosity_level = int(verbosity_level)

    def solve(self, problem, instances):
        c_h, lb_h, ub_h, sense_flip = assemble_batch(problem, instances)

        def _one(i):
            qp = _build_qp(problem, c_h[i], lb_h[i], ub_h[i])
            res = pdlp.primal_dual_hybrid_gradient(qp, self.params)
            x = np.asarray(res.primal_solution, dtype=np.float64)
            obj_internal = float(np.dot(c_h[i], x))
            # The termination-reason enum lives in solve_log_pb2 (OR-Tools 9.x).
            status_enum = res.solve_log.termination_reason
            status = (
                "optimal"
                if status_enum == solve_log_pb2.TERMINATION_REASON_OPTIMAL
                else f"pdlp:{int(status_enum)}"
            )
            fluxes = {
                problem.var_names[j]: float(x[j])
                for j in range(problem.n_vars)
            }
            return BatchedSolution(
                id=instances[i].id,
                status=status,
                objective_value=float(sense_flip[i] * obj_internal),
                fluxes=fluxes,
            )

        if self.workers == 1 or len(instances) <= 1:
            return [_one(i) for i in range(len(instances))]
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            return list(pool.map(_one, range(len(instances))))


if _PDLP_OK:
    register_batched_solver("pdlp", ORToolsPDLPBatchedLPSolver)
