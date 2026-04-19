# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

# Base class for FBA packages
class CommKineticPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "community kinetics", {}, {"commkin": "string"})

    def build_package(self, kinetic_coef, community_model, probs=None):
        self.validate_parameters({}, [], {"kinetic_coef": kinetic_coef, "community": community_model})
        cons = {cons.name: cons for cons in self.model.constraints}
        for species in self.parameters["community"].members:
            if species.id+"_commKin" in cons:
                self.model.remove_cons_vars(cons[species.id+"_commKin"])
            self.build_constraint(species, probs)

    def build_constraint(self, species, probs):
        coef = {species.primary_biomass.forward_variable: -1 * self.parameters["kinetic_coef"],
                species.primary_biomass.reverse_variable: self.parameters["kinetic_coef"]}
        for rxn in self.model.reactions:
            try:
                rxnIndex = int(FBAHelper.rxn_compartment(rxn)[1:])
            except (ValueError, TypeError, IndexError):
                continue
            if rxnIndex == species.index and rxn != species.primary_biomass:
                val = 1 if not isinstance(probs, dict) else probs.get(rxn.id, 1)
                coef[rxn.forward_variable] = coef[rxn.reverse_variable] = val
        return BaseFBAPkg.build_constraint(self, "commKin", None, 0, coef, species.id)
