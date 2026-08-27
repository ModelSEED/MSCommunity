# -*- coding: utf-8 -*-
"""Regression tests for the MSCommunity simulation defects (cluster A).

Everything runs on a synthetic 2-member ModelSEED-style community so the tests
are fast and solver-agnostic:

    memA eats cpd00027 (uptake cap 10)   -> solo growth capacity 10
    memB eats cpd00082 (uptake cap  4)   -> solo growth capacity  4

`bio1` (the community biomass) consumes 0.5 of each member's biomass metabolite,
which is exactly the coupling that made the solo-capacity LPs infeasible when the
`min_comm_growth` floor was installed before they were solved.
"""
from __future__ import annotations

import logging

import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.batched_lp import BatchedLPSolver, BatchedSolution, register_batched_solver
from mscommunity.commhelper import build_from_species_models
from mscommunity.mscommsim import MSCommunity, _pick_qp_backend, _select_biomass_cpd


# --- synthetic ModelSEED-style member models -------------------------------

def _member(mid, carbon, uptake):
    m = Model(mid)
    mets = {}
    for cid, comp in [(carbon, "e0"), (carbon, "c0"), ("cpd00002", "c0"), ("cpd00001", "c0"),
                      ("cpd00008", "c0"), ("cpd00009", "c0"), ("cpd00067", "c0"), ("cpd11416", "c0")]:
        key = f"{cid}_{comp}"
        mets[key] = Metabolite(key, name=key, compartment=comp)
    m.add_metabolites(list(mets.values()))
    ex = Reaction(f"EX_{carbon}_e0", lower_bound=-float(uptake), upper_bound=1000)
    ex.add_metabolites({mets[f"{carbon}_e0"]: -1})
    tr = Reaction("rxn10000_c0", lower_bound=0, upper_bound=1000)
    tr.add_metabolites({mets[f"{carbon}_e0"]: -1, mets[f"{carbon}_c0"]: 1})
    bio = Reaction("bio1", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({mets[f"{carbon}_c0"]: -1, mets["cpd11416_c0"]: 1})
    m.add_reactions([ex, tr, bio])
    m.objective = "bio1"
    return m


def _members():
    """Asymmetric pair: independent carbon sources, solo capacities 10 and 4."""
    return [_member("memA", "cpd00027", 10), _member("memB", "cpd00082", 4)]


def _competing():
    """Both members share one carbon source capped at 10 -> degenerate max-sum face."""
    return [_member("memA", "cpd00027", 10), _member("memB", "cpd00027", 10)]


def _community(members=None):
    return MSCommunity(model=build_from_species_models(members or _members()))


# --- defect 1: solo capacities must be measured without the community floor --

def test_solo_capacities_are_the_true_per_member_maxima():
    comm = _community()
    comm.set_objective(targets=[m.primary_biomass.forward_variable for m in comm.members])
    solo = comm._solo_max_batch()
    assert solo == pytest.approx({"Species0": 10.0, "Species1": 4.0})


def test_regularization_installs_per_member_floors_from_the_true_solo_maxima():
    comm = _community()
    comm.set_objective(targets=[m.primary_biomass.forward_variable for m in comm.members])
    comm.regularization(linear=True)

    cons = {c.name: c for c in comm.util.model.constraints}
    # the floor for the determining solve is still installed (unchanged behavior)
    assert "min_comm_growth" in cons
    assert cons["min_comm_growth"].lb > 0
    # ...and the per-member floors now exist at all, at ratio * TRUE solo capacity.
    # Before the fix min_comm_growth was created first, every solo LP came back
    # infeasible, every capacity was recorded as 0 and no floor was ever built.
    assert "Species0_regularization" in cons and "Species1_regularization" in cons
    ratio = cons["Species0_regularization"].lb / 10.0
    assert ratio == pytest.approx(0.7)
    assert cons["Species1_regularization"].lb == pytest.approx(ratio * 4.0)


def test_regularization_suspends_a_leftover_floor_before_measuring():
    comm = _community()
    comm.set_objective(targets=[m.primary_biomass.forward_variable for m in comm.members])
    # a floor left over from an earlier call must not poison the solo LPs...
    comm.util.create_constraint(comm.util.model.problem.Constraint(
        comm.primary_biomass.flux_expression, name="min_comm_growth", lb=7.5))
    comm.regularization(linear=True, growth_constraint=False)

    cons = {c.name: c for c in comm.util.model.constraints}
    assert cons["Species0_regularization"].lb == pytest.approx(0.7 * 10.0)
    assert cons["Species1_regularization"].lb == pytest.approx(0.7 * 4.0)
    # ...and it must be restored afterwards, since growth_constraint=False means
    # regularization() is not the owner of that constraint.
    assert "min_comm_growth" in cons and cons["min_comm_growth"].lb == pytest.approx(7.5)


def test_regularization_ratio_loop_relaxes_until_feasible():
    # shared carbon capped at 10: both solo maxima are 10, so floors of 0.7*10
    # each are jointly infeasible and the loop must fall back to 0.5.
    comm = _community(_competing())
    comm.set_objective(targets=[m.primary_biomass.forward_variable for m in comm.members])
    comm.regularization(linear=True)
    cons = {c.name: c for c in comm.util.model.constraints}
    assert cons["Species0_regularization"].lb == pytest.approx(5.0)
    assert cons["Species1_regularization"].lb == pytest.approx(5.0)
    assert comm.util.model.slim_optimize() == pytest.approx(10.0)


# --- defect 2: a non-optimal solo LP has no capacity, whatever it reports ----

class _StaleBackend(BatchedLPSolver):
    """Mimics GLPK: a non-optimal solve still reports a stale, positive objective."""

    name = "stale_for_test"

    def __init__(self, workers=1, **kwargs):
        self.workers = workers

    def solve(self, problem, instances):
        return [BatchedSolution(id=i.id, status="infeasible", objective_value=3.75, fluxes={})
                for i in instances]


def test_solo_max_ignores_the_objective_of_a_non_optimal_solve(caplog):
    register_batched_solver("stale_for_test", _StaleBackend)
    comm = _community()
    with caplog.at_level(logging.WARNING, logger="mscommunity.mscommsim"):
        solo = comm._solo_max_batch(backend="stale_for_test")
    assert solo == {"Species0": 0, "Species1": 0}
    assert any("solo-capacity LP" in r.message for r in caplog.records)


# --- defect 3: a sub-optimal solve is loud, and never crashes on the dump ----

def test_suboptimal_solution_is_logged_and_does_not_crash(caplog, monkeypatch, tmp_path):
    comm = _community()
    monkeypatch.chdir(tmp_path)   # the old dump used a bare relative filename
    with comm.util.model:
        comm.util.create_constraint(comm.util.model.problem.Constraint(
            comm.primary_biomass.flux_expression, name="impossible", lb=1e6))
        with caplog.at_level(logging.ERROR, logger="mscommunity.mscommsim"):
            sol = comm.run_fba(None, False)
    assert sol.status != "optimal"
    assert comm.suboptimal_solution is True
    assert any("sub-optimal" in r.message for r in caplog.records)
    # no diagnostics were requested, so nothing was written
    assert not list(tmp_path.iterdir())


def test_print_lp_accepts_a_bare_relative_filename(monkeypatch, tmp_path):
    comm = _community()
    monkeypatch.chdir(tmp_path)
    comm.print_lp("erronous_model.lp")          # used to raise FileNotFoundError('')
    assert (tmp_path / "erronous_model.lp").exists()


# --- defect 4: micom() must not hard-swap the solver -------------------------

def test_micom_keeps_a_qp_capable_incumbent_and_restores_the_solver():
    if _pick_qp_backend() is None:
        pytest.skip("no QP-capable optlang backend installed")
    comm = _community()
    before = comm.util.model.solver.interface.__name__
    sols = comm.micom({"EX_cpd00027_e0": 10, "EX_cpd00082_e0": 4})
    assert len(sols) == 1
    assert comm.util.model.solver.interface.__name__ == before


def test_micom_raises_when_no_qp_backend_is_available(monkeypatch):
    import optlang
    monkeypatch.setattr(optlang, "available_solvers", {}, raising=False)
    comm = _community()
    before = comm.util.model.solver.interface.__name__
    with pytest.raises(RuntimeError, match="QP-capable"):
        comm.micom({"EX_cpd00027_e0": 10})
    assert comm.util.model.solver.interface.__name__ == before


# --- defect 5: the member_models= construction path --------------------------

def test_member_models_construction_end_to_end():
    comm = MSCommunity(member_models=_members())   # used to raise AttributeError
    assert set(comm.abundances) == {"memA", "memB"}
    assert comm.abundances == pytest.approx({"memA": 0.5, "memB": 0.5})
    biomassIDs = {mem.id: mem.biomass_cpd.id for mem in comm.members}
    assert biomassIDs == {"memA": "cpd11416_c1", "memB": "cpd11416_c2"}
    # the resolved metabolites belong to the LIVE (copied) community model
    for mem in comm.members:
        assert comm.util.model.metabolites.get_by_id(mem.biomass_cpd.id) is mem.biomass_cpd
    assert comm.predict_abundances(regularization=True, pfba=False) is not None


def test_select_biomass_cpd_handles_lists_and_scalars():
    met = Metabolite("cpd11416_c1", compartment="c1")
    other = Metabolite("cpd11416_c1", compartment="c1")   # pre-copy twin
    assert _select_biomass_cpd([other], {"cpd11416_c1": met}, {}) is met
    assert _select_biomass_cpd(other, {}, {"cpd11416_c1": met}) is met
    assert _select_biomass_cpd([other], {}, {}) is None


# --- defect 6: __init__ defaults --------------------------------------------

def test_probs_default_is_not_shared_mutable_state():
    import inspect
    assert inspect.signature(MSCommunity.__init__).parameters["probs"].default is None
    a, b = _community(), _community()
    assert a.rxnProbs == {} and b.rxnProbs == {}
    a.rxnProbs["rxn10000_c1"] = 0.5
    assert b.rxnProbs == {}          # the old {} default was shared by every instance


def test_flux_limit_is_still_accepted_even_though_unused():
    comm = MSCommunity(model=build_from_species_models(_members()), flux_limit=42)
    assert comm.abundances == pytest.approx({"Species0": 0.5, "Species1": 0.5})


# --- defect 7: coherent determinism instrumentation --------------------------

def test_determinizer_failure_increments_the_fallback_counter(monkeypatch):
    comm = _community(_competing())
    monkeypatch.setattr(MSCommunity, "_solve_qp_highs", lambda self, hess_vars: None)
    assert comm.pfba_fallback_count == 0
    # unregularized, the shared-carbon max-sum face is flat (any A + B = 10), so the
    # determinizer engages -- and then fails, because the QP backend is stubbed out.
    comm.predict_abundances(regularization=False, pfba=False, determinize=True)
    assert comm.pfba_fell_back is True
    assert comm.pfba_fallback_count == 1        # was stuck at 0 before the fix
    comm.predict_abundances(regularization=False, pfba=False, determinize=True)
    assert comm.pfba_fallback_count == 2


def test_a_converged_row_leaves_the_counter_alone():
    comm = _community()
    comm.predict_abundances(regularization=True, pfba=False, determinize=True)
    assert comm.pfba_fell_back is False
    assert comm.pfba_fallback_count == 0


# --- the hot path still produces sensible relative abundances ----------------

def test_prediction_returns_sensible_relative_abundances():
    comm = _community()
    abund = comm.predict_abundances(regularization=True, pfba=False)
    # memA can pull 10 units of carbon, memB only 4 -> 10/14 vs 4/14
    assert abund["Species0"] == pytest.approx(10 / 14, rel=1e-6)
    assert abund["Species1"] == pytest.approx(4 / 14, rel=1e-6)
    assert sum(abund.values()) == pytest.approx(1.0)
    # the regularization scaffolding is torn down again
    names = {c.name for c in comm.util.model.constraints}
    assert not any("regularization" in n or n == "min_comm_growth" for n in names)


def test_prediction_splits_a_degenerate_face_evenly_under_regularization():
    comm = _community(_competing())
    abund = comm.predict_abundances(regularization=True, pfba=False)
    assert abund["Species0"] == pytest.approx(0.5, abs=1e-6)
    assert abund["Species1"] == pytest.approx(0.5, abs=1e-6)
