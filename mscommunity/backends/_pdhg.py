# -*- coding: utf-8 -*-
"""Shared batched PDHG (primal-dual hybrid gradient) algorithm.

The math is identical across NumPy / CuPy / JAX — only the array module
and the sparse matvec change. Each backend supplies an `ArrayBackend`
with three callables and the rest is library-agnostic.

For each instance i in a batch we solve:

    max  c_i^T x
    s.t.  S x = 0                       (mass balance, shared S)
          A_e x in [a_lo_i, a_hi_i]     (extra row constraints, shared A_e)
          lb_i <= x <= ub_i             (box on net flux)

PDHG updates the augmented saddle point with primal `x` (n, B) and a
stacked dual `y = [y_S; y_E]` (m+e, B). Constant-step variant; we leave
adaptive step sizes / restart heuristics to PDLP-class backends.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp


@dataclass
class PDHGConfig:
    max_iters: int = 5000
    tol: float = 1e-4
    check_every: int = 50
    power_iters: int = 30
    step_safety: float = 0.9
    seed: int = 0


@dataclass
class ArrayBackend:
    """Library-agnostic shim. `xp` is the array namespace; `to_device` and
    `to_host` move tensors; `sparse_from_scipy` returns whatever the
    backend's matmul expects on the left."""

    xp: Any
    to_device: Callable[[np.ndarray], Any]
    to_host: Callable[[Any], np.ndarray]
    sparse_from_scipy: Callable[[sp.spmatrix], Any]
    spmm: Callable[[Any, Any], Any]      # A (m,n)   X (n,B) -> (m,B)
    spmm_T: Callable[[Any, Any], Any]    # A^T (n,m) Y (m,B) -> (n,B)


def assemble_batch(
    problem,  # CommunityProblem (avoid circular import at type-check time)
    instances: Sequence,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build (B, n)-shaped c / lb / ub and per-instance sign trackers.

    Internally we always maximize, so min-sense instances get c flipped
    here and the reported objective is un-flipped on return.
    """
    B = len(instances)
    n = problem.n_vars
    lb = np.tile(problem.lb, (B, 1)).astype(np.float64, copy=True)
    ub = np.tile(problem.ub, (B, 1)).astype(np.float64, copy=True)
    default_sign = 1.0 if problem.sense == "max" else -1.0
    default_c = default_sign * problem.c
    c = np.tile(default_c, (B, 1)).astype(np.float64, copy=True)
    sense_flip = np.ones(B, dtype=np.float64)  # 1 if we report +c^T x, -1 if -c^T x
    for i, inst in enumerate(instances):
        if inst.pfba:
            raise NotImplementedError(
                "pFBA is a two-stage solve; the PDHG backends do not implement"
                " it yet. Use the cpu backend (cobra.pfba) for pfba=True."
            )
        for rxn_id, (lo, hi) in inst.bounds.items():
            j = problem.var_index[rxn_id]
            lb[i, j] = lo
            ub[i, j] = hi
        if inst.objective is not None:
            row = np.zeros(n, dtype=np.float64)
            for rxn_id, coef in inst.objective.items():
                row[problem.var_index[rxn_id]] = float(coef)
            sense = inst.sense or problem.sense
            if sense == "min":
                c[i] = -row
                sense_flip[i] = -1.0
            else:
                c[i] = row
                sense_flip[i] = 1.0
        elif inst.sense == "min" and problem.sense == "max":
            c[i] = -problem.c
            sense_flip[i] = -1.0
        elif inst.sense == "max" and problem.sense == "min":
            c[i] = problem.c
            sense_flip[i] = 1.0
    return c, lb, ub, sense_flip


def stack_constraints(
    problem,
) -> Tuple[sp.csr_matrix, np.ndarray, np.ndarray, int, int]:
    """Stack S (equality, lb=ub=0) and A_extra (boxed) into one matrix.

    Returns A_full (sp.csr), full_lb, full_ub, m (rows of S), e (extra
    rows). The first `m` rows are mass balance with lb=ub=0.
    """
    m = problem.n_mets
    full_lb = np.zeros(m, dtype=np.float64)
    full_ub = np.zeros(m, dtype=np.float64)
    if problem.A_extra is not None:
        A_full = sp.vstack([problem.S, problem.A_extra]).tocsr()
        full_lb = np.concatenate([full_lb, problem.extra_lb])
        full_ub = np.concatenate([full_ub, problem.extra_ub])
        e = problem.n_extra
    else:
        A_full = problem.S.tocsr()
        e = 0
    return A_full, full_lb, full_ub, m, e


def estimate_spectral_norm(
    backend: ArrayBackend,
    A_dev: Any,
    n: int,
    power_iters: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    v_h = rng.standard_normal((n, 1)).astype(np.float64)
    v_h /= np.linalg.norm(v_h)
    v = backend.to_device(v_h)
    for _ in range(power_iters):
        u = backend.spmm(A_dev, v)
        v = backend.spmm_T(A_dev, u)
        norm_v = float(backend.xp.linalg.norm(v))
        if norm_v == 0:
            return 1.0
        v = v / norm_v
    Av = backend.spmm(A_dev, v)
    return float(backend.xp.linalg.norm(Av))


def run_pdhg_batch(
    backend: ArrayBackend,
    problem,
    instances: Sequence,
    config: Optional[PDHGConfig] = None,
):
    """Run batched PDHG on `instances` against `problem`.

    Returns a list[BatchedSolution]; the caller (a backend's `solve()`)
    handles registry + dispatch.
    """
    from mscommunity.batched_lp import BatchedSolution

    cfg = config or PDHGConfig()
    xp = backend.xp
    c_h, lb_h, ub_h, sense_flip = assemble_batch(problem, instances)
    A_full_sp, full_lb_h, full_ub_h, m, e = stack_constraints(problem)
    n = problem.n_vars
    B = len(instances)

    A_dev = backend.sparse_from_scipy(A_full_sp)
    sigma_max = estimate_spectral_norm(
        backend, A_dev, n, cfg.power_iters, cfg.seed
    )
    step = cfg.step_safety / max(sigma_max, 1e-12)
    tau = sigma = step

    c = backend.to_device(c_h.T)        # (n, B)
    lb = backend.to_device(lb_h.T)      # (n, B)
    ub = backend.to_device(ub_h.T)      # (n, B)
    full_lb_d = backend.to_device(full_lb_h.reshape(-1, 1))   # (m+e, 1) broadcast
    full_ub_d = backend.to_device(full_ub_h.reshape(-1, 1))

    x = xp.zeros((n, B), dtype=xp.float64)
    y = xp.zeros((m + e, B), dtype=xp.float64)

    iters_used = 0
    for k in range(cfg.max_iters):
        iters_used = k + 1
        grad = -c + backend.spmm_T(A_dev, y)
        x_new = xp.clip(x - tau * grad, lb, ub)
        xbar = 2.0 * x_new - x
        Ax = backend.spmm(A_dev, xbar)
        # equality rows (b=0) clip [0,0] => Ax kept; inequality rows clip to box
        # PDHG with linear constraint a_lo <= Ax <= a_hi:
        #   y_{k+1} = proj_{[-inf, +inf]}( y_k + sigma * (Ax - b) )  if eq
        #   y_{k+1} = clip( y_k + sigma * Ax - shift, ... )          if ineq
        # Equivalent reformulation: split A by row using a slack r,
        # but simpler closed form for boxed rows:
        #   y_new = y + sigma*Ax, then enforce that y corresponds to
        #   constraint violation outside [a_lo, a_hi]. For boxed:
        #   z = y + sigma*Ax;  if a_lo == a_hi: y_new = z;  else:
        #   y_new = z - sigma * clip(z/sigma, a_lo, a_hi)  (Moreau form
        #   of prox of indicator of box on Ax = u, slack u).
        z = y + sigma * Ax
        u_proj = xp.clip(z / sigma, full_lb_d, full_ub_d)
        y_new = z - sigma * u_proj
        x, y = x_new, y_new

        if (k + 1) % cfg.check_every == 0:
            Ax_cur = backend.spmm(A_dev, x)
            primal_viol = xp.maximum(
                xp.maximum(0.0, full_lb_d - Ax_cur),
                xp.maximum(0.0, Ax_cur - full_ub_d),
            )
            primal_res = xp.linalg.norm(primal_viol, axis=0)
            obj = (c * x).sum(axis=0)
            rel = primal_res / (1.0 + xp.abs(obj))
            if bool((rel < cfg.tol).all()):
                break

    Ax_final = backend.spmm(A_dev, x)
    primal_viol = xp.maximum(
        xp.maximum(0.0, full_lb_d - Ax_final),
        xp.maximum(0.0, Ax_final - full_ub_d),
    )
    primal_res = backend.to_host(xp.linalg.norm(primal_viol, axis=0))
    obj_internal = backend.to_host((c * x).sum(axis=0))
    x_host = backend.to_host(x)   # (n, B)

    out: List[BatchedSolution] = []
    for i, inst in enumerate(instances):
        rel = primal_res[i] / (1.0 + abs(obj_internal[i]))
        status = "optimal" if rel < cfg.tol else "approximate"
        fluxes = {
            problem.var_names[j]: float(x_host[j, i])
            for j in range(n)
        }
        # un-flip objective to user-facing sense
        obj_user = float(sense_flip[i] * obj_internal[i])
        out.append(BatchedSolution(
            id=inst.id,
            status=status,
            objective_value=obj_user,
            fluxes=fluxes,
        ))
    return out
