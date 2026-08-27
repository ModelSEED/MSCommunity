# -*- coding: utf-8 -*-
"""Smoke tests for the batched LP interface.

These exercise the CPU reference backend on a tiny hand-built cobra model
to confirm the shared-S / batched-bounds contract returns the same answers
a per-instance cobra solve would.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.batched_lp import (
    CPUBatchedLPSolver,
    CommunityProblem,
    LPInstance,
    get_batched_solver,
    media_to_bounds,
    register_batched_solver,
    solve_batch,
    BatchedLPSolver,
)


def _tiny_model() -> Model:
    # A -> B -> sink, plus a source EX_A. Max sink flux.
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


def test_problem_extracts_shared_structure():
    m = _tiny_model()
    prob = CommunityProblem.from_model(m)
    assert prob.n_vars == 3
    assert prob.n_mets == 2
    assert set(prob.var_names) == {"EX_A", "R1", "SINK"}
    j_sink = prob.var_index["SINK"]
    assert prob.c[j_sink] == pytest.approx(1.0)
    # mass balance: every column sums into both A and B correctly
    S_dense = prob.S.toarray()
    j_r1 = prob.var_index["R1"]
    i_A = prob.met_names.index("A")
    i_B = prob.met_names.index("B")
    assert S_dense[i_A, j_r1] == -1
    assert S_dense[i_B, j_r1] == 1


def test_batch_matches_per_instance_solve():
    m = _tiny_model()
    # five different EX_A uptake caps
    caps = [10, 5, 2, 1, 0.5]
    instances = [
        LPInstance(id=f"cap_{c}", bounds={"EX_A": (-float(c), 0.0)})
        for c in caps
    ]
    batched = solve_batch(m, instances, backend="cpu")
    assert [b.status for b in batched] == ["optimal"] * len(caps)
    assert [b.objective_value for b in batched] == pytest.approx(caps)
    # ground truth via direct per-instance cobra solve
    for c, sol in zip(caps, batched):
        with m:
            m.reactions.EX_A.lower_bound = -c
            ref = m.optimize()
        assert sol.objective_value == pytest.approx(ref.objective_value)
        for rxn in m.reactions:
            assert sol.fluxes[rxn.id] == pytest.approx(ref.fluxes[rxn.id])
    # state restored after batched run
    assert m.reactions.EX_A.lower_bound == -10


def test_objective_patch_overrides_default():
    m = _tiny_model()
    inst = LPInstance(
        id="max_R1",
        bounds={"EX_A": (-4.0, 0.0)},
        objective={"R1": 1.0},
        sense="max",
    )
    [sol] = solve_batch(m, [inst])
    assert sol.objective_value == pytest.approx(4.0)
    # original model objective is unchanged
    assert "SINK" in str(m.objective.expression)


def test_threaded_backend_matches_sequential():
    m = _tiny_model()
    instances = [
        LPInstance(id=str(i), bounds={"EX_A": (-float(i + 1), 0.0)})
        for i in range(8)
    ]
    seq = solve_batch(m, instances, backend="cpu", workers=1)
    par = solve_batch(m, instances, backend="cpu", workers=4)
    for a, b in zip(seq, par):
        assert a.status == b.status
        assert a.objective_value == pytest.approx(b.objective_value)


def test_error_recorded_not_raised():
    m = _tiny_model()
    bad = LPInstance(id="bad", bounds={"DOES_NOT_EXIST": (0.0, 1.0)})
    [sol] = solve_batch(m, [bad], backend="cpu")
    assert sol.status.startswith("error:")
    assert sol.objective_value is None

    with pytest.raises(KeyError):
        solver = CPUBatchedLPSolver(on_error="raise")
        prob = CommunityProblem.from_model(m)
        solver.solve(prob, [bad])


def test_media_to_bounds_inverts_uptake_sign():
    m = _tiny_model()
    patch = media_to_bounds(m, {"EX_A": 7.5})
    lo, hi = patch["EX_A"]
    assert lo == pytest.approx(-7.5)
    assert hi == pytest.approx(0.0)  # original ub


def test_registry_accepts_custom_backend():
    class Dummy(BatchedLPSolver):
        name = "dummy"

        def solve(self, problem, instances):
            from mscommunity.batched_lp import BatchedSolution

            return [
                BatchedSolution(id=i.id, status="optimal", objective_value=42.0, fluxes={})
                for i in instances
            ]

    register_batched_solver("dummy", Dummy)
    solver = get_batched_solver("dummy")
    out = solver.solve(
        CommunityProblem.from_model(_tiny_model()),
        [LPInstance(id="x")],
    )
    assert out[0].objective_value == 42.0
