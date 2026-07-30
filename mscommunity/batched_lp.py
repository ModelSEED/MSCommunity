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

import logging
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from cobra import Model
from cobra.core import Solution

logger = logging.getLogger(__name__)


def _preview(names: Sequence[str], limit: int = 3) -> str:
    head = ", ".join(names[:limit])
    return head if len(names) <= limit else f"{head}, ... (+{len(names) - limit} more)"


def _objective_vector(model: Model, rxns: Sequence) -> np.ndarray:
    """Net-flux objective coefficients read off the model's own objective.

    `Reaction.objective_coefficient` only reports a coefficient when the
    objective is antisymmetric over the split variables (`+k` on forward,
    `-k` on reverse) — the form `model.objective = "RXN"` produces.
    MSCommunity builds its objectives over forward variables alone
    (`model.problem.Objective(rxn.forward_variable, ...)`), for which that
    property silently returns 0 for every reaction. Read the objective
    expression directly instead and fold each reaction's forward/reverse
    coefficients into one net-flux coefficient.
    """
    n = len(rxns)
    try:
        obj_coeffs = model.objective.get_linear_coefficients(model.variables)
    except Exception as exc:  # e.g. a quadratic objective
        logger.warning(
            "CommunityProblem.from_model(%s): could not read linear objective "
            "coefficients (%s: %s); falling back to Reaction.objective_coefficient, "
            "which is zero unless the objective is written over net flux.",
            model.id, type(exc).__name__, exc,
        )
        return np.fromiter(
            (float(r.objective_coefficient) for r in rxns), dtype=float, count=n
        )

    c = np.zeros(n, dtype=float)
    approximated: List[str] = []
    for j, r in enumerate(rxns):
        fwd = float(obj_coeffs.get(r.forward_variable, 0.0))
        rev = float(obj_coeffs.get(r.reverse_variable, 0.0))
        if fwd == 0.0 and rev == 0.0:
            continue
        # x = forward - reverse, so an antisymmetric pair (+k, -k) is exactly
        # k*x; a one-sided objective (the MSCommunity form) is k*x over the
        # half-space its bounds allow. Anything else has no exact net-flux
        # coefficient, so take the forward half and say so.
        if fwd != 0.0 and rev != 0.0 and not np.isclose(fwd, -rev):
            approximated.append(r.id)
        c[j] = fwd if fwd != 0.0 else -rev
    if approximated:
        logger.warning(
            "CommunityProblem.from_model(%s): %d objective term(s) load the forward "
            "and reverse halves non-antisymmetrically and have no exact net-flux "
            "coefficient; the forward coefficient was used [%s].",
            model.id, len(approximated), _preview(approximated),
        )
    return c


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
    split-variable LP form, which is a future enhancement. A constraint
    that also loads *auxiliary* optlang variables (the potential /
    delta-G / directionality variables thermodynamic packages add) has no
    net-flux row either — keeping the reaction terms alone would silently
    change the constraint rather than relax it — so those rows are
    dropped as well. Every drop is logged (once per extraction) and named
    in `extracted_with_warnings`; the batched problem is therefore a
    relaxation of the cobra model whenever that tuple is non-empty.
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

        c = _objective_vector(model, rxns)
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

        # forward / reverse split variable -> (reaction column, +1 / -1)
        var_role: Dict[object, Tuple[int, int]] = {}
        for j, r in enumerate(rxns):
            var_role[r.forward_variable] = (j, 1)
            var_role[r.reverse_variable] = (j, -1)

        mass_balance_names = set(met_names)
        e_rows: List[int] = []
        e_cols: List[int] = []
        e_vals: List[float] = []
        e_lb: List[float] = []
        e_ub: List[float] = []
        e_names: List[str] = []
        skipped: List[str] = []
        skipped_unreadable: List[str] = []
        skipped_aux: List[str] = []
        skipped_nonantisym: List[str] = []
        n_candidates = 0
        row_i = 0
        for cons in model.constraints:
            if cons.name in mass_balance_names:
                continue
            n_candidates += 1
            try:
                coeffs = cons.get_linear_coefficients(model.variables)
            except Exception:
                skipped.append(cons.name)
                skipped_unreadable.append(cons.name)
                continue
            # split the row into per-reaction forward / reverse halves; any
            # coefficient on a variable that is not a reaction half is an
            # auxiliary term that net flux cannot carry.
            halves: Dict[int, List[float]] = {}
            has_aux = False
            for var, coef in coeffs.items():
                v = float(coef)
                if v == 0.0:
                    continue
                role = var_role.get(var)
                if role is None:
                    has_aux = True
                    break
                j, sign = role
                pair = halves.setdefault(j, [0.0, 0.0])
                pair[0 if sign > 0 else 1] += v
            if has_aux:
                skipped.append(cons.name)
                skipped_aux.append(cons.name)
                continue
            row_entries = []
            antisym = True
            for j in sorted(halves):
                fwd, rev = halves[j]
                if not np.isclose(fwd, -rev):
                    antisym = False
                    break
                row_entries.append((j, fwd))
            if not antisym:
                skipped.append(cons.name)
                skipped_nonantisym.append(cons.name)
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

        if skipped:
            reasons = []
            if skipped_nonantisym:
                reasons.append(
                    f"{len(skipped_nonantisym)} non-antisymmetric, i.e. same-sign "
                    f"coefficients on the forward/reverse halves - this is what a "
                    f"community-kinetics constraint looks like, and it needs "
                    f"split-variable LP form [{_preview(skipped_nonantisym)}]"
                )
            if skipped_aux:
                reasons.append(
                    f"{len(skipped_aux)} carrying auxiliary (non-reaction) variables "
                    f"such as thermodynamic potential / delta-G / directionality "
                    f"terms, which net flux cannot express "
                    f"[{_preview(skipped_aux)}]"
                )
            if skipped_unreadable:
                reasons.append(
                    f"{len(skipped_unreadable)} whose linear coefficients could not "
                    f"be read (non-linear?) [{_preview(skipped_unreadable)}]"
                )
            logger.warning(
                "CommunityProblem.from_model(%s): dropped %d of %d extra constraint(s); "
                "the batched LP is a RELAXATION of the cobra model and may report a "
                "better objective than model.optimize(). Dropped: %s. Full list in "
                "CommunityProblem.extracted_with_warnings.",
                model.id, len(skipped), n_candidates, "; ".join(reasons),
            )

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

    With `workers > 1` each worker thread solves on its own copy of the
    community model, since a cobra model / optlang solver cannot be
    mutated concurrently.
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
            return [self._solve_one(problem, inst, problem.model) for inst in instances]
        # `with model:` mutates the shared cobra model (and its optlang
        # solver) in place, so workers must not share one. Give each worker
        # thread its own copy, built lazily and one at a time so that at most
        # `workers` copies exist and the copies themselves never race.
        local = threading.local()
        copy_lock = threading.Lock()

        def run(inst: LPInstance) -> BatchedSolution:
            model = getattr(local, "model", None)
            if model is None:
                with copy_lock:
                    model = problem.model.copy()
                local.model = model
            return self._solve_one(problem, inst, model)

        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            return list(pool.map(run, instances))

    def _solve_one(
        self, problem: CommunityProblem, inst: LPInstance, model: Optional[Model] = None
    ) -> BatchedSolution:
        if model is None:
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
    model: Model, medium: Mapping[str, float], close_unlisted: bool = True
) -> Dict[str, Tuple[float, float]]:
    """Convert a cobra medium dict to a per-exchange bounds patch.

    A cobra medium maps exchange-reaction id → max uptake (a positive
    magnitude); uptake is negative flux on an exchange written `met -->`,
    so the patch sets that reaction's lower bound to `-uptake` and keeps
    its upper bound (secretion is untouched). For the rarer exchange
    written `--> met` cobra caps uptake on the *upper* bound instead, and
    so does this function.

    The patch reproduces what `model.medium = medium` does, which means
    it also has to *close* every exchange the medium omits (`close_unlisted`,
    on by default): a medium sweep that only patched the listed exchanges
    would leave the model's baseline uptakes open, so media differing only
    in which exchanges they list would silently give identical solutions.
    Pass `close_unlisted=False` to get a patch that touches only the listed
    reactions.
    """
    patch: Dict[str, Tuple[float, float]] = {}

    def _uptake_bounds(rxn, bound: float) -> Tuple[float, float]:
        # mirrors cobra.Model.medium's set_active_bound
        if rxn.reactants:
            return (-float(bound), float(rxn.upper_bound))
        return (float(rxn.lower_bound), float(bound))

    listed = set()
    for ex_id, max_uptake in medium.items():
        rxn = model.reactions.get_by_id(ex_id)
        listed.add(rxn.id)
        patch[ex_id] = _uptake_bounds(rxn, max_uptake)
    if close_unlisted:
        for rxn in model.exchanges:
            if rxn.id in listed:
                continue
            is_export = rxn.reactants and not rxn.products
            closed = min(0.0, -rxn.lower_bound if is_export else rxn.upper_bound)
            patch[rxn.id] = _uptake_bounds(rxn, closed)
    return patch
