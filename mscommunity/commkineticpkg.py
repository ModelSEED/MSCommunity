# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

logger = logging.getLogger(__name__)


def member_kinetic_reactions(community, member):
    """Yield the reactions whose |flux| counts against ``member``'s kinetic budget.

    A reaction belongs to a member when its compartment index matches the member's,
    which in a ModelSEED community model means that member's cytosol plus the
    transporters reaching the shared ``e0`` pool.  Community-level exchanges sit at
    index 0 and drop out.

    Excluded by identity rather than by compartment:

    * every biomass reaction -- this member's, the other members', and the community
      biomass.  ``FBAHelper.rxn_compartment`` resolves a multi-compartment reaction to
      the first non-extracellular compartment it happens to iterate, so the community
      biomass (c0 + one metabolite per member cytosol) can resolve to an arbitrary
      member and would otherwise land nondeterministically inside somebody's budget.
    * every member's biomass drain, which is a boundary reaction rather than
      metabolism.  Leaving it in let the solver buy (kinetic_coef - 1) units of budget
      for each unit of biomass it synthesised and immediately discarded.

    This is the single definition of membership; both ``CommKineticPkg`` and
    ``MSCommunity.add_commkinetics`` consume it so the two cannot drift apart.
    """
    excluded = {member.primary_biomass.id} if member.primary_biomass is not None else set()
    if getattr(community, "primary_biomass", None) is not None:
        excluded.add(community.primary_biomass.id)
    for other in community.members:
        if other.primary_biomass is not None:  excluded.add(other.primary_biomass.id)
        if other.biomass_drain is not None:  excluded.add(other.biomass_drain.id)
    for rxn in community.util.model.reactions:
        if rxn.id in excluded:  continue
        try:  index = int(FBAHelper.rxn_compartment(rxn)[1:])
        except (ValueError, TypeError, IndexError):  continue
        if index == member.index:  yield rxn


# Base class for FBA packages
class CommKineticPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "community kinetics", {}, {"commKin": "string"})

    def build_package(self, kinetic_coef, community_model, probs=None):
        self.validate_parameters({}, [], {"kinetic_coef": kinetic_coef, "community": community_model})
        cons = {cons.name: cons for cons in self.model.constraints}
        for species in self.parameters["community"].members:
            if species.id+"_commKin" in cons:
                self.model.remove_cons_vars(cons[species.id+"_commKin"])
            self.build_constraint(species, probs)

    def build_constraint(self, species, probs):
        community = self.parameters["community"]
        coef = {species.primary_biomass.forward_variable: -1 * self.parameters["kinetic_coef"],
                species.primary_biomass.reverse_variable: self.parameters["kinetic_coef"]}
        for rxn in member_kinetic_reactions(community, species):
            val = 1 if not isinstance(probs, dict) else probs.get(rxn.id, 1)
            coef[rxn.forward_variable] = coef[rxn.reverse_variable] = val
        return BaseFBAPkg.build_constraint(self, "commKin", None, 0, coef, species.id)
