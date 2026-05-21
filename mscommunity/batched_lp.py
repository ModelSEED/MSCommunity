# -*- coding: utf-8 -*-
"""Batched LP solving for community FBA.

Across many community simulations (one per sample, per medium, per
condition) the stoichiometric matrix S is fixed and only bounds / objective
coefficients vary. This module captures that structure once
(`CommunityProblem.from_model`) and routes many independent LPs through a
pluggable batched-solver backend.

The included `CPUBatchedLPSolver` is a reference backend that re-uses
cobra's per-LP solver (optionally threaded). A future GPU backend can
implement `BatchedLPSolver.solve` against the same contract — consuming the
already-extracted `S`, `lb`, `ub`, `c` plus per-instance patches — without
touching any MSCommunity call site.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from cobra import Model
from cobra.core import Solution


@dataclass
class CommunityProblem:
    """Frozen snapshot of a community model in standard LP form.

    Captured once and reused across many `LPInstance`s that vary only in
    bounds / objective. Backends may consume either the raw matrices
    (`S`, `lb`, `ub`, `c`, plus the boxed extra constraints `A_extra`,
    `extra_lb`, `extra_ub`) or the underlying cobra `model` — both views
    are kept in sync at construction time.

    Constraint extraction caveat: extra row constraints are expressed in
    net-flux form (x = forward - reverse). A constraint with coefficient
    `+K` on `forward_var` and `-K` on `reverse_var` collapses to `+K` on
    `x` (antisymmetric, the cobra default). A constraint with the same
    sign on both halves (e.g. an L1-flux bound) is *non-antisymmetric*
    and cannot be expressed on net flux alone; those rows are dropped and
    `extracted_with_warnings` lists their names. Such constraints need
    split-variable LP form, which is a future enhancement.
    """

    model: Model
    var_names: Tuple[str, ...]
    var_index: Dict[str, int]
    lb: np.ndarray
    ub: np.ndarray
    c: np.ndarray
    sense: str
    S: sp.csr_matrix
    met_names: Tuple[str, ...]
    A_extra: Optional[sp.csr_matrix] = None
    extra_lb: Optional[np.ndarray] = None
    extra_ub: Optional[np.ndarray] = None
    extra_cons_names: Tuple[str, ...] = ()
    extracted_with_warnings: Tuple[str, ...] = ()

    @property
    def n_vars(self) -> int:
        return len(self.var_names)

    @property
    def n_mets(self) -> int:
        return len(self.met_names)

    @property
    def n_extra(self) -> int:
        return len(self.extra_cons_names)

    @classmethod
    def from_model(cls, model: Model) -> "CommunityProblem":
        rxns = list(model.reactions)
        n = len(rxns)
        var_names = tuple(r.id for r in rxns)
        var_index = {name: i for i, name in enumerate(var_names)}
        lb = np.fromiter((r.lower_bound for r in rxns), dtype=float, count=n)
        ub = np.fromiter((r.upper_bound for r in rxns), dtype=float, count=n)

        c = np.array([float(r.objective_coefficient) for r in rxns], dtype=float)
        sense = "min" if model.objective.direction == "min" else "max"

        mets = list(model.metabolites)
        met_names = tuple(m.id for m in mets)
        met_index = {m.id: i for i, m in enumerate(mets)}
        rows: List[int] = []
        cols: List[int] = []
        vals: List[float] = []
        for j, r in enumerate(rxns):
            for met, coef in r.metabolites.items():
                rows.append(met_index[met.id])
                cols.append(j)
                vals.append(float(coef))
        S = sp.csr_matrix((vals, (rows, cols)), shape=(len(mets), n), dtype=float)

        mass_balance_names = set(met_names)
        e_rows: List[int] = []
        e_cols: List[int] = []
        e_vals: List[float] = []
        e_lb: List[float] = []
        e_ub: List[float] = []
        e_names: List[str] = []
        skipped: List[str] = []
        row_i = 0
        for cons in model.constraints:
            if cons.name in mass_balance_names:
                continue
            try:
                coeffs = cons.get_linear_coefficients(model.variables)
            except Exception:
                skipped.append(cons.name)
                continue
            row_entries = []
            antisym = True
            for j, r in enumerate(rxns):
                fwd = float(coeffs.get(r.forward_variable, 0.0))
                rev = float(coeffs.get(r.reverse_variable, 0.0))
                if fwd == 0.0 and rev == 0.0:
                    continue
                if not np.isclose(fwd, -rev):
                    antisym = False
                    break
                row_entries.append((j, fwd))
            if not antisym:
                skipped.append(cons.name)
                continue
            if not row_entries:
                continue
            for j, v in row_entries:
                e_rows.append(row_i)
                e_cols.append(j)
                e_vals.append(v)
            e_lb.append(float(cons.lb) if cons.lb is not None else -np.inf)
            e_ub.append(float(cons.ub) if cons.ub is not None else  np.inf)
            e_names.append(cons.name)
            row_i += 1

        if e_names:
            A_extra = sp.csr_matrix(
                (e_vals, (e_rows, e_cols)), shape=(len(e_names), n), dtype=float,
            )
            extra_lb = np.asarray(e_lb, dtype=float)
            extra_ub = np.asarray(e_ub, dtype=float)
        else:
            A_extra = None
            extra_lb = None
            extra_ub = None

        return cls(
            model=model,
            var_names=var_names,
            var_index=var_index,
            lb=lb,
            ub=ub,
            c=c,
            sense=sense,
            S=S,
            met_names=met_names,
            A_extra=A_extra,
            extra_lb=extra_lb,
            extra_ub=extra_ub,
            extra_cons_names=tuple(e_names),
            extracted_with_warnings=tuple(skipped),
        )


@dataclass
class LPInstance:
    """A per-sample / per-condition variation of the shared community problem."""

    id: str
    bounds: Mapping[str, Tuple[float, float]] = field(default_factory=dict)
    objective: Optional[Mapping[str, float]] = None
    sense: Optional[str] = None
    pfba: bool = False


@dataclass
class BatchedSolution:
    id: str
    status: str
    objective_value: Optional[float]
    fluxes: Dict[str, float]

    def to_cobra_solution(self) -> Solution:
        import pandas as pd

        obj = self.objective_value if self.objective_value is not None else float("nan")
        return Solution(
            objective_value=obj,
            status=self.status,
            fluxes=pd.Series(self.fluxes, dtype=float),
        )


class BatchedLPSolver:
    """Backend contract: solve N LP instances sharing one CommunityProblem."""

    name = "base"

    def solve(
        self,
        problem: CommunityProblem,
        instances: Sequence[LPInstance],
    ) -> List[BatchedSolution]:
        raise NotImplementedError


class CPUBatchedLPSolver(BatchedLPSolver):
    """Reference backend: round-trip each instance through cobra's solver.

    A GPU backend with the same `solve()` signature drops in without
    changing any call site. The structural win (one-shot extraction of S,
    bounds, objective) is already paid; this backend simply doesn't
    exploit batching across instances.
    """

    name = "cpu"

    def __init__(self, workers: int = 1, on_error: str = "record"):
        if on_error not in ("record", "raise"):
            raise ValueError("on_error must be 'record' or 'raise'")
        self.workers = max(1, int(workers))
        self.on_error = on_error

    def solve(
        self,
        problem: CommunityProblem,
        instances: Sequence[LPInstance],
    ) -> List[BatchedSolution]:
        if self.workers == 1 or len(instances) <= 1:
            return [self._solve_one(problem, inst) for inst in instances]
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            return list(pool.map(lambda inst: self._solve_one(problem, inst), instances))

    def _solve_one(
        self, problem: CommunityProblem, inst: LPInstance
    ) -> BatchedSolution:
        model = problem.model
        try:
            with model:
                for rxn_id, (lo, hi) in inst.bounds.items():
                    rxn = model.reactions.get_by_id(rxn_id)
                    rxn.lower_bound = float(lo)
                    rxn.upper_bound = float(hi)
                if inst.objective is not None:
                    expr = sum(
                        float(coef) * model.reactions.get_by_id(rid).flux_expression
                        for rid, coef in inst.objective.items()
                    )
                    direction = inst.sense or problem.sense
                    model.objective = model.problem.Objective(expr, direction=direction)
                elif inst.sense is not None:
                    model.objective.direction = inst.sense
                if inst.pfba:
                    from cobra.flux_analysis import pfba as _pfba
                    sol = _pfba(model)
                else:
                    sol = model.optimize()
        except Exception as exc:
            if self.on_error == "raise":
                raise
            return BatchedSolution(
                id=inst.id,
                status=f"error:{type(exc).__name__}",
                objective_value=None,
                fluxes={},
            )
        return BatchedSolution(
            id=inst.id,
            status=sol.status,
            objective_value=getattr(sol, "objective_value", None),
            fluxes={k: float(v) for k, v in sol.fluxes.items()},
        )


_SOLVERS: Dict[str, type] = {"cpu": CPUBatchedLPSolver}


def register_batched_solver(name: str, cls: type) -> None:
    if not issubclass(cls, BatchedLPSolver):
        raise TypeError("cls must subclass BatchedLPSolver")
    _SOLVERS[name] = cls


def get_batched_solver(name: str = "cpu", **kwargs) -> BatchedLPSolver:
    if name not in _SOLVERS:
        raise KeyError(
            f"Unknown batched LP backend {name!r}. Registered: {sorted(_SOLVERS)}"
        )
    return _SOLVERS[name](**kwargs)


def solve_batch(
    model: Model,
    instances: Sequence[LPInstance],
    backend: str = "cpu",
    problem: Optional[CommunityProblem] = None,
    **backend_kwargs,
) -> List[BatchedSolution]:
    """One-shot helper: snapshot the model once, then run a batch."""
    if problem is None:
        problem = CommunityProblem.from_model(model)
    solver = get_batched_solver(backend, **backend_kwargs)
    return solver.solve(problem, instances)


def media_to_bounds(
    model: Model, medium: Mapping[str, float]
) -> Dict[str, Tuple[float, float]]:
    """Convert a cobra medium dict to a per-exchange bounds patch.

    A cobra medium maps exchange-reaction id → max uptake (a positive
    magnitude); uptake is negative flux on exchanges. The returned patch
    keeps the existing upper bound and sets the lower bound to `-uptake`.
    """
    patch: Dict[str, Tuple[float, float]] = {}
    for ex_id, max_uptake in medium.items():
        rxn = model.reactions.get_by_id(ex_id)
        patch[ex_id] = (-float(max_uptake), float(rxn.upper_bound))
    return patch
