# -*- coding: utf-8 -*-
"""Regression tests for the member-biomass-drain defects.

Three faults lived together in `CommunityMember.__init__` and the two copies of the
community-kinetics constraint:

1. the drain test sat inside an `if "bio" in rxn.id` branch, so the ModelSEED-style
   `SK_cpd11416_c<i>` sink that `build_from_species_models` carries over was never
   recognised and a redundant `DM_cpd11416_c<i>` was built beside it;
2. the classification loop `break`ed on finding the primary biomass, so anything
   ordered after it -- the sink included -- was never examined at all;
3. `MSCommunity.add_commkinetics` and `CommKineticPkg.build_constraint` applied
   different exclusion rules, and neither excluded the drains, so each unit of biomass
   synthesised and discarded bought `kinetic_coef - 1` units of flux budget.
"""
from __future__ import annotations

import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.commhelper import build_from_species_models
from mscommunity.commkineticpkg import member_kinetic_reactions
from mscommunity.mscommsim import MSCommunity


def _member(mid, carbon, uptake):
    """A minimal ModelSEED-style member that already carries its own biomass sink."""
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
    # the sink the real ModelSEED reconstructions ship with, and that the constructor
    # used to walk straight past
    sink = Reaction("SK_cpd11416_c0", lower_bound=0, upper_bound=1000)
    sink.add_metabolites({mets["cpd11416_c0"]: -1})
    m.add_reactions([ex, tr, bio, sink])
    m.objective = "bio1"
    return m


def _community(close_member_drains=False):
    models = [_member("memA", "cpd00027", 10), _member("memB", "cpd00082", 4)]
    comm_model = build_from_species_models(models, abundances={"memA": 0.5, "memB": 0.5})
    return MSCommunity(model=comm_model, ids=["memA", "memB"], kinetic_coeff=1e6,
                       close_member_drains=close_member_drains)


def test_existing_biomass_sink_is_found_instead_of_building_a_second_drain():
    comm = _community()
    for member in comm.members:
        assert member.biomass_drain is not None
        assert member.biomass_drain.id == f"SK_cpd11416_c{member.index}"
    assert not [r.id for r in comm.util.model.reactions if r.id.startswith("DM_cpd11416")]


def test_member_reaction_lists_are_not_truncated_at_the_biomass_reaction():
    """The old `break` cut the scan short, so `self.reactions` -- and the member's
    standalone model -- lost everything ordered after the biomass reaction."""
    comm = _community()
    for member in comm.members:
        ids = {r.id for r in member.reactions}
        assert f"rxn10000_c{member.index}" in ids


def test_the_kinetic_row_excludes_every_biomass_reaction_and_drain():
    comm = _community()
    biomass_ids = {comm.primary_biomass.id} | {m.primary_biomass.id for m in comm.members}
    drain_ids = {m.biomass_drain.id for m in comm.members}
    for member in comm.members:
        counted = {r.id for r in member_kinetic_reactions(comm, member)}
        assert not counted & biomass_ids
        assert not counted & drain_ids
        # ...while the member's own metabolism is still in
        assert f"rxn10000_c{member.index}" in counted


def test_closing_the_drains_pins_member_growth_to_its_abundance():
    # memA can reach its solo capacity of 10 only by discarding biomass: the community
    # biomass reaction needs an equal share from memB, which is capped at 4.
    open_comm = _community(close_member_drains=False)
    open_comm.util.model.objective = "bio2"
    assert open_comm.util.model.slim_optimize() == pytest.approx(10.0)

    shut_a = _community(close_member_drains=True)
    shut_a.util.model.objective = "bio2"
    assert shut_a.util.model.slim_optimize() == pytest.approx(4.0)

    shut = _community(close_member_drains=True)
    for member in shut.members:
        assert member.biomass_drain.bounds == (0, 0)
    shut.util.model.objective = shut.primary_biomass.id
    sol = shut.util.model.optimize()
    assert sol.status == "optimal"
    for member in shut.members:
        assert sol.fluxes[member.primary_biomass.id] == pytest.approx(
            member.abundance * sol.fluxes[shut.primary_biomass.id], abs=1e-6)
