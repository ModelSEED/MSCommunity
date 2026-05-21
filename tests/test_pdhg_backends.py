# -*- coding: utf-8 -*-
"""Tests for the GPU-capable PDHG / PDLP backends.

The NumPy backend is fully exercised (same algorithm as the GPU shims,
just on CPU). The JAX / CuPy / PDLP shims are only tested for their
import-time behavior since their hardware deps aren't installed here.
"""
from __future__ import annotations

import importlib

import numpy as np
import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.batched_lp import (
    BatchedLPSolver,
    CommunityProblem,
    LPInstance,
    get_batched_solver,
    solve_batch,
    _SOLVERS,
)
import mscommunity.backends  # noqa: F401 — triggers registration


def _tiny_model() -> Model:
    m = Model("tiny")
    A = Metabolite("A", compartment="c")
    B = Metabolite("B", compartment="c")
    m.add_metabolites([A, B])
    ex_a = Reaction("EX_A", lower_bound=-10, upper_bound=0)
    ex_a.add_metabolites({A: -1})
    r1 = Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    sink = Reaction("SINK", lower_bound=0, upper_bound=1000)
    sink.add_metabolites({B: -1})
    m.add_reactions([ex_a, r1, sink])
    m.objective = "SINK"
    return m


def test_numpy_pdhg_registered():
    assert "numpy-pdhg" in _SOLVERS


def test_numpy_pdhg_matches_cpu_reference():
    m = _tiny_model()
    caps = [10.0, 5.0, 2.0, 1.0]
    instances = [
        LPInstance(id=f"cap_{c}", bounds={"EX_A": (-c, 0.0)}) for c in caps
    ]
    cpu = solve_batch(m, instances, backend="cpu")
    pdhg = solve_batch(m, instances, backend="numpy-pdhg",
                       max_iters=20_000, tol=1e-5, check_every=100)
    assert [s.status for s in pdhg] == ["optimal"] * len(caps)
    for cpu_s, pd_s in zip(cpu, pdhg):
        assert pd_s.objective_value == pytest.approx(
            cpu_s.objective_value, rel=1e-3, abs=1e-3
        )


def test_numpy_pdhg_objective_patch():
    m = _tiny_model()
    inst = LPInstance(
        id="max_R1", bounds={"EX_A": (-4.0, 0.0)},
        objective={"R1": 1.0}, sense="max",
    )
    [sol] = solve_batch(m, [inst], backend="numpy-pdhg",
                        max_iters=20_000, tol=1e-5)
    assert sol.objective_value == pytest.approx(4.0, rel=1e-3, abs=1e-3)


def test_numpy_pdhg_min_sense_flips_objective_back():
    # The reactions have flux >= 0, so minimizing SINK should give 0.
    m = _tiny_model()
    inst = LPInstance(
        id="min_SINK", bounds={"EX_A": (-4.0, 0.0)},
        objective={"SINK": 1.0}, sense="min",
    )
    [sol] = solve_batch(m, [inst], backend="numpy-pdhg",
                        max_iters=20_000, tol=1e-5)
    assert sol.objective_value == pytest.approx(0.0, abs=1e-3)


def test_problem_extracts_extra_constraints():
    # Force-add a side constraint: SINK <= 3
    m = _tiny_model()
    cons = m.problem.Constraint(
        m.reactions.SINK.flux_expression, name="cap_sink", ub=3.0
    )
    m.add_cons_vars(cons)
    prob = CommunityProblem.from_model(m)
    assert "cap_sink" in prob.extra_cons_names
    assert prob.A_extra is not None
    assert prob.A_extra.shape == (1, 3)
    # PDHG should honor it
    [sol] = solve_batch(m, [LPInstance(id="x", bounds={"EX_A": (-10.0, 0.0)})],
                        backend="numpy-pdhg",
                        max_iters=20_000, tol=1e-5)
    assert sol.objective_value == pytest.approx(3.0, rel=1e-3, abs=1e-3)


def test_pfba_raises_on_pdhg():
    m = _tiny_model()
    inst = LPInstance(id="x", bounds={"EX_A": (-1.0, 0.0)}, pfba=True)
    with pytest.raises(NotImplementedError):
        solve_batch(m, [inst], backend="numpy-pdhg")


def test_optional_backends_clean_failure_when_missing():
    # JAX / CuPy / OR-Tools aren't installed in this env. Their shims
    # must (a) not break package import, (b) raise a clear ImportError
    # at construction time, (c) not be in the registry.
    for name, mod_name in [
        ("jax-pdhg", "mscommunity.backends.jax_pdhg"),
        ("cupy-pdhg", "mscommunity.backends.cupy_pdhg"),
        ("pdlp", "mscommunity.backends.ortools_pdlp"),
    ]:
        mod = importlib.import_module(mod_name)
        # If dep missing, the backend should not be registered
        ok_flag = next(
            (v for k, v in vars(mod).items() if k.endswith("_OK")), None
        )
        if ok_flag is False:
            assert name not in _SOLVERS
            cls = next(
                v for v in vars(mod).values()
                if isinstance(v, type)
                and issubclass(v, BatchedLPSolver)
                and v is not BatchedLPSolver
            )
            with pytest.raises(ImportError):
                cls()
        else:
            # Dep is installed; the backend should be registered
            assert name in _SOLVERS


def test_registry_lists_pdhg_after_import():
    # numpy-pdhg always loads; cpu always loads; pdhg + cpu = >=2 backends
    assert {"cpu", "numpy-pdhg"}.issubset(set(_SOLVERS))
