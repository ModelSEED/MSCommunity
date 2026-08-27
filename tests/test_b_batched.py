# -*- coding: utf-8 -*-
"""Regression tests for the batched LP extraction / CPU backend (cluster B).

Covers:
  * objective vectors written over forward variables (the MSCommunity form)
  * constraints carrying auxiliary (non-reaction) variables being dropped
    instead of silently emitted without their auxiliary terms
  * a real logging warning whenever rows are dropped
  * the threaded CPU backend not racing on one shared cobra model
  * media_to_bounds agreeing with cobra's ``model.medium`` setter
"""
from __future__ import annotations

import logging

import numpy as np
import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.batched_lp import (
    CommunityProblem,
    LPInstance,
    media_to_bounds,
    solve_batch,
)


def _tiny_model() -> Model:
    # A -> B -> sink, plus a source EX_A. Max sink flux.
    m = Model("tiny_b")
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


# --------------------------------------------------------------------------
# 1. objective extraction
# --------------------------------------------------------------------------

def test_objective_over_forward_variables_is_not_zero():
    """MSCommunity writes objectives over forward variables only."""
    m = _tiny_model()
    sink = m.reactions.SINK
    m.objective = m.problem.Objective(sink.forward_variable, direction="max")
    prob = CommunityProblem.from_model(m)
    assert prob.c[prob.var_index["SINK"]] == pytest.approx(1.0)
    assert prob.c[prob.var_index["R1"]] == pytest.approx(0.0)
    assert prob.sense == "max"


def test_objective_over_net_flux_still_works():
    m = _tiny_model()
    sink = m.reactions.SINK
    for expr in (sink.flux_expression, 1.0 * sink.flux_expression):
        m.objective = m.problem.Objective(expr, direction="max")
        prob = CommunityProblem.from_model(m)
        assert prob.c[prob.var_index["SINK"]] == pytest.approx(1.0)
    # the plain id form (antisymmetric) too
    m.objective = "SINK"
    prob = CommunityProblem.from_model(m)
    assert prob.c[prob.var_index["SINK"]] == pytest.approx(1.0)


def test_multi_term_forward_objective_keeps_weights_and_direction():
    m = _tiny_model()
    sink, r1 = m.reactions.SINK, m.reactions.R1
    m.objective = m.problem.Objective(
        2.0 * sink.forward_variable + 0.5 * r1.forward_variable, direction="min"
    )
    prob = CommunityProblem.from_model(m)
    assert prob.c[prob.var_index["SINK"]] == pytest.approx(2.0)
    assert prob.c[prob.var_index["R1"]] == pytest.approx(0.5)
    assert prob.c[prob.var_index["EX_A"]] == pytest.approx(0.0)
    assert prob.sense == "min"


# --------------------------------------------------------------------------
# 2/3. constraint extraction + warnings
# --------------------------------------------------------------------------

def test_constraint_with_auxiliary_variable_is_dropped_not_mangled():
    m = _tiny_model()
    sink = m.reactions.SINK
    aux = m.problem.Variable("dG_SINK", lb=0, ub=10)
    cons = m.problem.Constraint(
        sink.forward_variable - sink.reverse_variable + aux,
        lb=0, ub=5, name="thermo_SINK",
    )
    m.add_cons_vars([aux, cons])

    prob = CommunityProblem.from_model(m)
    # the row must NOT be emitted with the auxiliary term silently deleted
    assert "thermo_SINK" not in prob.extra_cons_names
    assert prob.A_extra is None
    assert "thermo_SINK" in prob.extracted_with_warnings


def test_pure_reaction_constraint_is_still_retained():
    m = _tiny_model()
    r1, sink = m.reactions.R1, m.reactions.SINK
    cons = m.problem.Constraint(
        r1.flux_expression - 2.0 * sink.flux_expression, lb=-1, ub=1, name="ratio",
    )
    m.add_cons_vars([cons])

    prob = CommunityProblem.from_model(m)
    assert prob.extra_cons_names == ("ratio",)
    assert prob.extracted_with_warnings == ()
    row = prob.A_extra.toarray()[0]
    assert row[prob.var_index["R1"]] == pytest.approx(1.0)
    assert row[prob.var_index["SINK"]] == pytest.approx(-2.0)
    assert prob.extra_lb[0] == pytest.approx(-1.0)
    assert prob.extra_ub[0] == pytest.approx(1.0)


def test_dropped_rows_emit_a_logging_warning(caplog):
    m = _tiny_model()
    sink = m.reactions.SINK
    # a community-kinetics-shaped row: same sign on both halves
    kin = m.problem.Constraint(
        sink.forward_variable + sink.reverse_variable, lb=0, ub=7, name="commkin",
    )
    aux = m.problem.Variable("dG_SINK", lb=0, ub=10)
    thermo = m.problem.Constraint(
        sink.forward_variable - sink.reverse_variable + aux,
        lb=0, ub=5, name="thermo_SINK",
    )
    m.add_cons_vars([kin, aux, thermo])

    with caplog.at_level(logging.WARNING, logger="mscommunity.batched_lp"):
        prob = CommunityProblem.from_model(m)

    assert set(prob.extracted_with_warnings) == {"commkin", "thermo_SINK"}
    warnings = [
        r for r in caplog.records
        if r.name == "mscommunity.batched_lp" and r.levelno >= logging.WARNING
    ]
    # one aggregated message, not one per row
    assert len(warnings) == 1
    msg = warnings[0].getMessage()
    assert "dropped 2" in msg
    assert "commkin" in msg and "thermo_SINK" in msg
    assert "non-antisymmetric" in msg and "auxiliary" in msg


def test_no_warning_when_nothing_is_dropped(caplog):
    m = _tiny_model()
    with caplog.at_level(logging.WARNING, logger="mscommunity.batched_lp"):
        CommunityProblem.from_model(m)
    assert [r for r in caplog.records if r.name == "mscommunity.batched_lp"] == []


# --------------------------------------------------------------------------
# 4. thread safety
# --------------------------------------------------------------------------

def test_threaded_backend_does_not_race_on_shared_model():
    m = _tiny_model()
    n = 24
    instances = [
        LPInstance(id=str(i), bounds={"EX_A": (-float(i + 1), 0.0)}) for i in range(n)
    ]
    expected = [float(i + 1) for i in range(n)]
    for _ in range(3):
        par = solve_batch(m, instances, backend="cpu", workers=8)
        assert [s.status for s in par] == ["optimal"] * n
        assert [s.objective_value for s in par] == pytest.approx(expected)
        # every worker must return the full flux vector for its own instance
        for inst, sol in zip(instances, par):
            assert sol.fluxes["EX_A"] == pytest.approx(inst.bounds["EX_A"][0])
    # the shared model is left exactly as it was
    assert m.reactions.EX_A.bounds == (-10, 0)


def test_threaded_backend_matches_sequential_with_objective_patch():
    m = _tiny_model()
    instances = [
        LPInstance(
            id=str(i),
            bounds={"EX_A": (-float(i + 1), 0.0)},
            objective={"R1": 1.0},
            sense="max",
        )
        for i in range(12)
    ]
    seq = solve_batch(m, instances, backend="cpu", workers=1)
    par = solve_batch(m, instances, backend="cpu", workers=4)
    for a, b in zip(seq, par):
        assert a.id == b.id
        assert a.status == b.status
        assert a.objective_value == pytest.approx(b.objective_value)
    assert "SINK" in str(m.objective.expression)


# --------------------------------------------------------------------------
# 5. media_to_bounds
# --------------------------------------------------------------------------

def _model_with_two_exchanges() -> Model:
    m = _tiny_model()
    ex_b = Reaction("EX_B", lower_bound=-10, upper_bound=0)
    ex_b.add_metabolites({m.metabolites.B: -1})
    m.add_reactions([ex_b])
    return m


def test_media_to_bounds_matches_cobra_medium_setter():
    m = _model_with_two_exchanges()
    medium = {"EX_A": 7.5}
    patch = media_to_bounds(m, medium)

    ref = m.copy()
    ref.medium = medium
    for rxn in ref.reactions:
        lo, hi = patch.get(rxn.id, m.reactions.get_by_id(rxn.id).bounds)
        assert lo == pytest.approx(rxn.lower_bound), rxn.id
        assert hi == pytest.approx(rxn.upper_bound), rxn.id


def test_media_to_bounds_closes_exchanges_the_medium_omits():
    m = _model_with_two_exchanges()
    patch = media_to_bounds(m, {"EX_A": 1.0})
    assert patch["EX_A"] == pytest.approx((-1.0, 0.0))
    # EX_B is open for uptake in the baseline model but absent from the medium
    assert patch["EX_B"][0] == pytest.approx(0.0)
    # opting out reproduces the listed-reactions-only patch
    only = media_to_bounds(m, {"EX_A": 1.0}, close_unlisted=False)
    assert set(only) == {"EX_A"}


def test_media_sweep_distinguishes_media_that_differ_by_an_omission():
    m = _model_with_two_exchanges()
    # B can only be produced from A via R1, or taken up for free through EX_B
    with_b = LPInstance(id="with_b", bounds=media_to_bounds(m, {"EX_A": 1.0, "EX_B": 5.0}))
    without_b = LPInstance(id="without_b", bounds=media_to_bounds(m, {"EX_A": 1.0}))
    sols = solve_batch(m, [with_b, without_b], backend="cpu")
    assert sols[0].objective_value == pytest.approx(6.0)
    assert sols[1].objective_value == pytest.approx(1.0)


def test_media_to_bounds_handles_product_side_exchange():
    m = _tiny_model()
    # an exchange written "--> A": cobra caps uptake on the upper bound
    src = Reaction("SRC_A", lower_bound=0, upper_bound=0)
    src.add_metabolites({m.metabolites.A: 1})
    m.add_reactions([src])
    patch = media_to_bounds(m, {"SRC_A": 3.0}, close_unlisted=False)
    ref = m.copy()
    ref.medium = {"SRC_A": 3.0}
    assert patch["SRC_A"] == pytest.approx(
        (ref.reactions.SRC_A.lower_bound, ref.reactions.SRC_A.upper_bound)
    )


def test_extracted_matrices_stay_consistent():
    m = _tiny_model()
    prob = CommunityProblem.from_model(m)
    assert prob.S.shape == (prob.n_mets, prob.n_vars)
    assert np.asarray(prob.lb).shape == (prob.n_vars,)
    assert np.asarray(prob.c).shape == (prob.n_vars,)
