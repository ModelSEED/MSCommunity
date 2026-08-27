# -*- coding: utf-8 -*-
"""Regression tests for mscommunity.mskineticsfba.MSKineticsFBA.baseKinFBA."""

from cobra import Model, Metabolite, Reaction
import pytest

from mscommunity.mskineticsfba import MSKineticsFBA

CELL_DRY_G, CELLULAR_L = 1.44e-13, 1e-18
CELL_G_L = CELL_DRY_G / CELLULAR_L
TS_MIN, TOTAL_MIN = 20, 40
TIMESTEP_HR = TS_MIN / 60

# a rate law of 1000*[A_c in M] over an initial 2 mM pool pins R1 at a flux of 2
KINETICS = {"R1": {"src1": {"substituted_rate_law": "1000*A_c", "mets": ["A_c"],
                            "met_id": {"A_c": "A_c"}, "initial_M": {"A_c": 0.002},
                            "metadata": {"Temperature": 30, "pH": 7.2}}}}


def _tiny_model():
    model = Model("tiny")
    A_e = Metabolite("A_e", compartment="e")
    A_c = Metabolite("A_c", compartment="c")
    B_c = Metabolite("B_c", compartment="c")
    B_e = Metabolite("B_e", compartment="e")
    reactions = []
    for rxn_id, stoich, lb, ub in [("EX_A_e", {A_e: -1}, -10, 1000), ("A_t", {A_e: -1, A_c: 1}, 0, 1000),
                                   ("R1", {A_c: -1, B_c: 1}, -5, 5), ("B_t", {B_c: -1, B_e: 1}, 0, 1000),
                                   ("EX_B_e", {B_e: -1}, 0, 1000)]:
        rxn = Reaction(rxn_id, lower_bound=lb, upper_bound=ub)
        rxn.add_metabolites(stoich)
        reactions.append(rxn)
    model.add_reactions(reactions)
    # X_c participates in no reaction, hence its concentration must simply persist across the timesteps
    model.add_metabolites([Metabolite("X_c", compartment="c")])
    model.objective = "EX_B_e"
    return model


def _run(**kwargs):
    kinfba = MSKineticsFBA(_tiny_model())
    conc, fluxes = kinfba.baseKinFBA(kinetics_data=KINETICS, total_min=TOTAL_MIN, ts_min=TS_MIN,
                                     cell_dry_g=CELL_DRY_G, cellular_L=CELLULAR_L,
                                     visualize=False, export=False, **kwargs)
    return kinfba, conc, fluxes


def test_default_call_does_not_raise():
    # < initial_M > defaults to None, which must not be iterated
    kinfba, conc, fluxes = _run()
    assert list(conc.columns) == ["0 min", "20 min", "40 min"]
    assert len(kinfba.sols) == kinfba.parameters["timesteps"] == 2
    assert kinfba.total_min == TOTAL_MIN  # consumed by _visualize


def test_kinetic_flux_reaches_the_linear_program():
    """The rate law must actually constrain the LP.

    The bug was `rxn.lb = rxn.ub = flux`: cobra defines lower_bound/upper_bound,
    so those assignments created stray attributes and the solver never saw the
    rate. Assert the solved flux equals the rate law's value rather than
    inspecting the bounds afterwards, since the pin is scoped to the simulation
    and is deliberately rolled back when baseKinFBA returns.
    """
    kinfba, conc, fluxes = _run()
    for time in ["20 min", "40 min"]:
        assert fluxes.at["R1", time] == pytest.approx(2)
    for sol in kinfba.sols:
        assert sol.fluxes["R1"] == pytest.approx(2)  # 2, not the unconstrained optimum of 5


def test_the_simulation_does_not_mutate_the_callers_model():
    """Pinned bounds and per-timestep constraints must not outlive the call."""
    model = _tiny_model()
    before = model.reactions.get_by_id("R1").bounds
    kinfba = MSKineticsFBA(model)
    kinfba.baseKinFBA(kinetics_data=KINETICS, total_min=TOTAL_MIN, ts_min=TS_MIN,
                      cell_dry_g=CELL_DRY_G, cellular_L=CELLULAR_L, visualize=False, export=False)
    assert model.reactions.get_by_id("R1").bounds == before
    assert not [c for c in model.constraints if c.name.endswith("_conc")]


def test_concentrations_integrate_across_timesteps():
    kinfba, conc, fluxes = _run(initial_M={"X_c": 0.005})
    times = list(conc.columns)
    # an unreacting metabolite retains its concentration instead of being zeroed each timestep
    for time in times:
        assert conc.at["X_c", time] == pytest.approx(5)
    # The medium changes only through the boundary reactions: [X]_t = [X]_(t-1) + v_ex*dt*cell_g_L.
    # Summing over ALL of a metabolite's reactions, as the code used to, evaluates the mass-balance
    # row the solver fixes at zero, so every concentration stayed flat for the whole simulation.
    for previous_time, time in zip(times, times[1:]):
        for met in kinfba.model_util.model.metabolites:
            delta = sum(fluxes.at[rxn.id, time] * TIMESTEP_HR * CELL_G_L
                        for rxn in met.reactions if rxn.boundary)
            assert conc.at[met.id, time] == pytest.approx(conc.at[met.id, previous_time] + delta, abs=1e-9)


def test_the_concentration_update_is_not_identically_zero():
    """Guard against regressing to the mass-balance formulation, which is inert.

    A_e is drawn down through EX_A_e, so its concentration must actually move.
    """
    kinfba, conc, fluxes = _run()
    assert fluxes.at["EX_A_e", "20 min"] != pytest.approx(0)
    assert conc.at["A_e", "20 min"] != pytest.approx(conc.at["A_e", "0 min"])


def test_nonnegativity_constraint_bounds_uptake_by_the_available_amount():
    # A_e must carry a starting concentration for its pool to be tracked at all. The default
    # cell_dry_g/cellular_L is a packed-cell density (1.44e5 g/L) that scales 20 minutes of even a
    # modest flux past any millimolar pool, so use a culture-scale density here instead.
    kinfba = MSKineticsFBA(_tiny_model())
    conc, fluxes = kinfba.baseKinFBA(kinetics_data=KINETICS, total_min=TOTAL_MIN, ts_min=TS_MIN,
                                     cell_dry_g=1e-13, cellular_L=1e-13,  # 1 g/L of biomass
                                     initial_M={"A_e": 0.01}, visualize=False, export=False)
    culture_g_L = 1.0
    model = kinfba.model_util.model
    # The guard is written on the boundary reactions, which are the only ones that move material
    # between the cell and the tracked medium; a bound on the full stoichiometric row would be a
    # bound on the solver's own mass balance and would constrain nothing.
    cons = kinfba.constraints["A_e"][TS_MIN]
    assert cons.ub is None
    ## the floor is -[A_e] at the start of the step, i.e. uptake may not exceed what is present
    assert cons.lb == pytest.approx(-conc.at["A_e", "0 min"])
    boundary = [rxn for rxn in model.metabolites.get_by_id("A_e").reactions if rxn.boundary]
    assert [rxn.id for rxn in boundary] == ["EX_A_e"]
    ## and the solved uptake actually respects it:  v_ex * dt * cell_g_L >= -[A_e]
    uptake = fluxes.at["EX_A_e", "20 min"]
    assert uptake * TIMESTEP_HR * culture_g_L >= -conc.at["A_e", "0 min"] - 1e-6
    ## the pool is drawn down rather than going negative
    assert conc.at["A_e", "20 min"] >= -1e-9


def test_untracked_metabolites_are_not_guarded():
    """A metabolite with no starting concentration must not have its uptake forbidden.

    Bounding uptake by "what is available" for a pool the caller never initialized
    pins the exchange at >= 0, which forbids uptake and makes the medium infeasible
    for the mundane reason that no starting concentration was given.
    """
    kinfba, conc, fluxes = _run()  # nothing in the medium is initialized
    assert "A_e" not in kinfba.constraints
    for sol in kinfba.sols:
        assert sol.status == "optimal"


def test_find_data_match_weights_the_deviations_equally():
    kinfba = MSKineticsFBA(_tiny_model())
    kinfba.kinetics_data = KINETICS
    kinfba.parameters = {"temperature": 25, "pH": 7}
    kinfba.minimum = float("inf")
    match = kinfba._MSKineticsFBA__find_data_match("R1", "src1")
    assert match in ("a", "w")
    # mean(|25-30|/25, |7-7.2|/7) == mean(0.2, 0.0285714)
    assert kinfba.minimum == pytest.approx((5 / 25 + 0.2 / 7) / 2)
