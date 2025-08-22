# -*- coding: utf-8 -*-
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.exceptions import ObjectAlreadyDefinedError, FeasibilityError, NoFluxError
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from modelseedpy.core.msatpcorrection import MSATPCorrection
from mscommunity.commhelper import build_from_species_models
from mscommunity.mscommviz import interactions as mscommsim_interactions
from cobra.io import save_matlab_model, write_sbml_model
from itertools import combinations, permutations
from cobra.core.dictlist import DictList
from collections import Counter
from optlang.symbolics import Zero
from cobra.flux_analysis import pfba
from cobra import Reaction, Model
from numpy import array, logspace
from os import makedirs, path
from math import isclose
from pandas import DataFrame
from pprint import pprint
from math import exp
import logging

logger = logging.getLogger(__name__)


class CommunityMember:
    def __init__(self, community, biomass_cpd, ID=None, index=None, abundance=0):
        print(ID, "biomass compound:", biomass_cpd)
        self.community, self.biomass_cpd = community, biomass_cpd
        try:     self.index = int(self.biomass_cpd.compartment[1:])
        except:  self.index = index
        self.abundance = abundance
        if self.biomass_cpd in self.community.primary_biomass.metabolites:
            self.abundance = abs(self.community.primary_biomass.metabolites[self.biomass_cpd])
        if ID is not None:  self.id = ID
        elif "species_name" in self.biomass_cpd.annotation:
            self.id = self.biomass_cpd.annotation["species_name"]
        else:  self.id = f"Species{self.index}"

        logger.info(f"Making atp hydrolysis reaction for species: {self.id}")
        self.model = Model()
        atp_hydrolysis_rxnComp = f"c{self.index}"
        try:
            self.atp_hydrolysis = self.community.util.model.reactions.get_by_id(f"rxn00062_{atp_hydrolysis_rxnComp}")
            print(f"skipping atp hydrolysis rxn00062_{atp_hydrolysis_rxnComp} for {self.id}")
        except:
            atp_rxn = self.community.util.add_atp_hydrolysis(atp_hydrolysis_rxnComp)
            self.atp_hydrolysis = atp_rxn["reaction"]
            self.model.add_reactions([self.atp_hydrolysis])
            print(f"created atp hydrolysis reaction rxn00062_{atp_hydrolysis_rxnComp} for {self.id}")
        self.biomass_drain = self.primary_biomass = None
        self.reactions = []
        for rxn in self.community.util.model.reactions:
            if "bio" in rxn.id:
                mets = {met.id: met for met in rxn.metabolites}
                if self.biomass_cpd.id not in mets:   continue
                met = mets[self.biomass_cpd.id]
                if rxn.metabolites[met] == 1 and len(rxn.metabolites) > 1:  self.primary_biomass = rxn  ;  break
                elif len(rxn.metabolites) == 1 and rxn.metabolites[met] < 0:  self.biomass_drain = rxn
            else:
                rxnComp = FBAHelper.rxn_compartment(rxn)
                if rxnComp is None:  print(f"The reaction {rxn.id} compartment {rxnComp} is undefined.")
                elif rxnComp[1:] == '': print("no compartment", rxn, rxnComp)
                elif int(rxnComp[1:]) == self.index:  self.reactions.append(rxn)

        if self.primary_biomass is None:  print(f"No biomass reaction found for species {self.id}")
        if not self.biomass_drain:
            print(f"Making biomass drain reaction for species: {self.id}")
            self.biomass_drain = Reaction(id=f"DM_{self.biomass_cpd.id}", name=f"DM_{self.biomass_cpd.name}", lower_bound=0, upper_bound=100)
            self.community.util.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({self.biomass_cpd: -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'
        # reactions = self.reactions + [self.primary_biomass, self.biomass_drain]
        # print(Counter([rxn.id for rxn in reactions]))
        # TODO the best way of tracking the models may be to run build_from_species_models inside the MSCommunity class
        ## where the models are available and build available.
        self.model.add_reactions(self.reactions + [self.primary_biomass, self.biomass_drain])
        self.model.add_reactions(self.community.util.exchange_list())
        self.model.medium = self.community.util.model.medium
        self.model.objective = self.primary_biomass.flux_expression

    def disable_species(self):
        for reaction in self.model.reactions:
            reaction_index = FBAHelper.rxn_compartment(reaction)[1:]
            if int(reaction_index) == self.index:  reaction.upper_bound = reaction.lower_bound = 0

    def compute_max_biomass(self):
        if self.primary_biomass is None:  logger.critical("No biomass reaction found for species "+self.id)
        self.community.util.add_objective(self.primary_biomass.flux_expression)
        if self.community.lp_filename:  self.community.print_lp(f"{self.community.lp_filename}_{self.id}_Biomass")
        return self.community.model.optimize()

    def compute_max_atp(self):
        if not self.atp_hydrolysis: logger.critical("No ATP hydrolysis found for species:" + self.id)
        self.community.util.add_objective(Zero, coef={self.atp_hydrolysis.forward_variable: 1})
        if self.community.lp_filename:  self.community.print_lp(f"{self.community.lp_filename}_{self.id}_ATP")
        return self.community.model.optimize()


class MSCommunity:
    def __init__(self, model=None, member_models: list = None, abundances=None, ids=None, kinetic_coeff=750,
                 flux_limit=300, probs={}, climit=None, o2limit=None, lp_filename=None, printing=False):
        assert model is not None or member_models is not None, "Either the community model and the member models must be defined."
        self.lp_filename = lp_filename
        self.printing = printing
        self.gapfillings = {}

        #Define Data attributes as None
        for attr in ['solution', 'biomass_cpd', 'primary_biomass', 'biomass_drain', 'threshold', 'msgapfill', 
                     'element_uptake_limit', 'msdb_path', 'comm_growth', 'threshold']:
            setattr(self, attr, None)
        self.kinCoef = kinetic_coeff
        # defining the models
        if model is None and member_models is not None:
            model, model_tracking = build_from_species_models(member_models, abundances=abundances, printing=printing)
        self.id = model.id
        self.util = MSModelUtil(model, True, None, climit, o2limit)
        self.pkgmgr = MSPackageManager.get_pkg_mgr(self.util.model)
        msid_cobraid_hash = self.util.msid_hash()  # dict of list() of metabolite objects by their msid
        if "cpd11416" not in msid_cobraid_hash:  raise KeyError("Could not find biomass compound for the model.")
        other_biomass_cpds = []
        for self.biomass_cpd in msid_cobraid_hash["cpd11416"]:
            if "c0" in self.biomass_cpd.id:
                for rxn in self.util.model.reactions:
                    if self.biomass_cpd not in rxn.metabolites:  continue
                    # print(self.biomass_cpd, rxn, end=";\t")
                    if rxn.metabolites[self.biomass_cpd] == 1 and len(rxn.metabolites) > 1:
                        if self.primary_biomass:  raise ObjectAlreadyDefinedError(
                            f"The primary biomass {self.primary_biomass} is already defined,"
                            f"hence, the {rxn.id} cannot be defined as the model primary biomass.")
                        if printing:  print('primary biomass defined', rxn.id)
                        self.primary_biomass = rxn
                    elif rxn.metabolites[self.biomass_cpd] < 0 and len(rxn.metabolites) == 1:  self.biomass_drain = rxn
            elif 'c' in self.biomass_cpd.compartment:   other_biomass_cpds.append(self.biomass_cpd)
        
        if ids is None:
            if member_models is not None:   ids = [mem.id for mem in member_models]
            else:  ids = [f"Species{i}" for i in range(len(other_biomass_cpds))]
        if not abundances:
            if member_models is None:
                abundances = {ids[memIndex]: {"biomass_compound": bioCpd, "abundance": 1/len(other_biomass_cpds)}
                              for memIndex, bioCpd in enumerate(other_biomass_cpds)}
            else:
                abundances = {}
                for memID, bioCPD in model.notes["member_biomass_cpds"].items():
                    abundances[memID] = {"abundance": 1/len(other_biomass_cpds)}
                    for met in model.metabolites:
                        if bioCPD.id == met.id:
                            if "biomass_compound" in abundances[memID]:   print("duplicate", bioCPD.id, met.id)
                            abundances[memID].update({"biomass_compound": met})
                            # print(bioCPD, met.id)
                    if "biomass_compound" not in abundances[memID]:   print(f"The {memID} bioCPD was not captured")

        # print()   # this returns the carriage after the tab-ends in the biomass compound printing
        self.members = DictList(CommunityMember(self, info["biomass_compound"], ID, index+1, info["abundance"])
                                for index, (ID, info) in enumerate(abundances.items()))
        # self.members = DictList(
        #     CommunityMember(community=self, biomass_cpd=biomass_cpd, name=ids[memIndex], abundance=abundances[memIndex])
        #     for memIndex, biomass_cpd in enumerate(other_biomass_cpds))
        self.set_abundance(abundances)

        
        # assign the MSCommunity constraints and objective
        self.rxnProbs = probs
        self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff, self, self.rxnProbs)
        # if kinetic_coeff is not None:   self.add_commkinetics(kinetic_coeff, probs)
        

    #Manipulation functions
    def set_abundance(self, abundances):
        #calculate the normalized biomass
        total_abundance = sum(list([content["abundance"] for content in abundances.values()]))
        # map abundances to all species
        for modelID, content in abundances.items():
            if modelID in self.members:  self.members.get_by_id(modelID).abundance = content["abundance"]/total_abundance
        self.abundances = {mem.id: mem.abundance for mem in self.members}
        #remake the primary biomass reaction based on abundances  #TODO what is the purpose of this?
        if self.primary_biomass is None:  logger.critical("Primary biomass reaction not found in community model")
        all_metabolites = {self.primary_biomass.products[0]: 1}
        all_metabolites.update({mem.biomass_cpd: -abundances[mem.id]["abundance"]/total_abundance for mem in self.members})
        self.primary_biomass.add_metabolites(all_metabolites, combine=False)
        self.abundances_set = True

    def set_objective(self, target=None, targets=None, minimize=False):
        targets = targets or [self.util.model.reactions.get_by_id(target or self.primary_biomass.id).flux_expression]
        self.util.model.objective = self.util.model.problem.Objective(sum(targets), direction="max" if not minimize else "min")

    def constrain(self, element_uptake_limit=None, thermo_params=None, msdb_path=None):
        if element_uptake_limit:
            self.element_uptake_limit = element_uptake_limit
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        if thermo_params:
            if msdb_path:
                self.msdb_path = msdb_path
                thermo_params.update({'modelseed_db_path':msdb_path})
                self.pkgmgr.getpkg("FullThermoPkg").build_package(thermo_params)
            else:  self.pkgmgr.getpkg("SimpleThermoPkg").build_package(thermo_params)

    def interactions(self, solution=None, media=None, msdb=None, msdb_path=None, filename=None, figure_format="svg",
                     node_metabolites=True, flux_threshold=1, visualize=True, ignore_mets=None):
        return mscommsim_interactions(self, solution or self.solution, media, flux_threshold, msdb, msdb_path,
                                        visualize, filename, figure_format, node_metabolites, True, ignore_mets)

    def add_commkinetics(self, kinCoef=750, probs={}):  #, abundances):
        self.rxnProbs = probs
        self.kinCoef = kinCoef
        for member in self.members:
            ## remove existing instance of CommKinetics
            consName = f"{member.id}_commKin"
            if consName in self.util.model.constraints:
                print(f"Removing {consName} from {self.util.model.id}")
                self.util.model.remove_cons_vars(self.util.model.constraints[consName])
            ## define the CommKinetics constraint:  kinCoef * bio_f,i > kinCoef * bio_r,i + sum(rxn_i * prob_r) 
            coef = {member.primary_biomass.forward_variable: -kinCoef, member.primary_biomass.reverse_variable: kinCoef}
            for rxn in self.util.model.reactions:
                rxnIndex = int(FBAHelper.rxn_compartment(rxn)[1:])
                if (rxnIndex == member.index and "bio" not in rxn.id):
                    coef[rxn.forward_variable] = coef[rxn.reverse_variable] = self.rxnProbs.get(rxn.id, 1)
            self.util.create_constraint(self.util.model.problem.Constraint(Zero, name=consName, ub=0), coef=coef, printing=True)

    #Utility functions
    def print_lp(self, filename=None):
        filename = filename or self.lp_filename
        makedirs(path.dirname(filename), exist_ok=True)
        with open(filename, 'w') as out:  out.write(str(self.util.model.solver))  ;  out.close()

    def to_sbml(self, export_name):
        makedirs(path.dirname(export_name), exist_ok=True)
        write_sbml_model(self.util.model, export_name)

    #Analysis functions
    def gapfill(self, media = None, target = None, minimize = False, default_gapfill_templates=None, default_gapfill_models=None,
                test_conditions=None, reaction_scores=None, blacklist=None, suffix = None, solver:str="glpk"):
        default_gapfill_templates = default_gapfill_templates or []
        default_gapfill_models = default_gapfill_models or []
        test_conditions, blacklist = test_conditions or [], blacklist or []
        reaction_scores = reaction_scores or {}
        if not target:  target = self.primary_biomass.id
        self.set_objective(target, minimize)
        gfname = FBAHelper.mediaName(media) + "-" + target
        if suffix:  gfname += f"-{suffix}"
        self.gapfillings[gfname] = MSGapfill(self.util.model, default_gapfill_templates, default_gapfill_models,
                                             test_conditions, reaction_scores, blacklist, solver)
        gfresults = self.gapfillings[gfname].run_gapfilling(media, target)
        assert gfresults, f"Gapfilling of {self.util.model.id} in {gfname} towards {target} failed."
        return self.gapfillings[gfname].integrate_gapfill_solution(gfresults)

    def test_individual_species(self, media=None, interacting=True, run_atp=True, run_biomass=True):
        assert run_atp or run_biomass, ValueError("Either the run_atp or run_biomass arguments must be True.")
        # self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if media is not None:  self.util.add_medium(media)
        data = {"Species": [], "Biomass": [], "ATP": []}
        for individual in self.members:
            data["Species"].append(individual.id)
            with self.util.model:
                if not interacting:
                    for other in self.members:
                        if other != individual:  other.disable_species()
                if run_biomass:  data["Biomass"].append(individual.compute_max_biomass())
                if run_atp:  data["ATP"].append(individual.compute_max_atp())
        return DataFrame(data)

    def atp_correction(self, core_template, atp_medias, max_gapfilling=None, gapfilling_delta=0):
        self.atp = MSATPCorrection(self.util.model, core_template, atp_medias, "c0", max_gapfilling, gapfilling_delta)

    def _growth_fraction(self):
        return 0.999 * exp(-0.03 * len(self.members))

    def regularization(self, linear=True):
        self.util.remove_constraint("_regularization")
        # self.util.remove_constraint("min_comm_growth")
        commMax = self.util.model.slim_optimize()
        # self.util.create_constraint(self.util.model.problem.Constraint(
        #     self.primary_biomass.flux_expression, name="min_comm_growth", lb=commMax*self._growth_fraction()), printing=True)
        if linear:
            # TODO create a threshold for the absolute difference between biomasses of the community members
            for threshold in [0.1, 0.5]+list(logspace(-0, 1, 8)):  # growth differences to check
                for mem1, mem2 in permutations(self.members, 2):
                    ## mu_i - mu_j <= threshold   &   mu_j - mu_i <= threshold
                    consName = f"{mem1.id}_{mem2.id}_regularization"
                    # if consName in self.util.model.constraints:
                    #     print(f"Removing {consName} from {self.util.model.id}")
                    #     self.util.model.remove_cons_vars(self.util.model.constraints[consName])
                    coef = {mem1.primary_biomass.forward_variable: 1, mem2.primary_biomass.forward_variable: -1}
                    self.util.create_constraint(self.util.model.problem.Constraint(Zero, name=consName, ub=threshold), coef=coef, printing=True)
                # print()
                sol = self.util.model.optimize()
                if sol.status == "optimal":
                    # if self.printing:
                    print(f"The model {self.util.model.id} is regularized with a max member growth difference of {threshold}")
                    break
        else:
            # TODO create the least-squares method employed in MICOM
            pass
        commCurrent = self.util.model.slim_optimize()
        ## TODO why are these community growths identical?
        print(f"Max community growth ({commMax}) and after regularization ({commCurrent})")
        return threshold

    # TODO evaluate the comparison of this method with MICOM
    def predict_abundances(self, media=None, pfba=True, timeout=60, environName=None, regularization=True, update_abundances=False):
        slimOpt = self.util.model.slim_optimize()
        if isclose(0, slimOpt, abs_tol=1e-3):
            print(f"The model {self.util.model.id} doesn't grow, with a slim_optimize of {slimOpt} in {environName} media")
        # store the original parameters
        ogObj = self.util.model.objective
        ogMedia = self.util.model.medium
        ogTimeout = self.util.model.solver.configuration.timeout
        # simulate the model
        ## maximize the sum of all member biomass reactions)
        self.set_objective(targets=[species.primary_biomass.forward_variable for species in self.members])
        self.util.model.solver.configuration.timeout = timeout
        if regularization:   threshold = self.regularization(linear=True)
        else:   self._remove_regularization_constraints()

        # threshold = self._regularization(linear=True)
        try:    sol = self.run_fba(media, pfba)
        except: sol = self.run_fba(media, pfba=False)
        # if sol.status != "optimal":
        #     raise FeasibilityError(f"The simulation for minimal uptake in {self.util.model.id} was {sol.status}.")
        # reset the model conditions
        self.util.model.solver.configuration.timeout = ogTimeout
        self.util.model.objective = ogObj
        self.util.model.medium = ogMedia
        return self._compute_relative_abundance_from_solution(sol, True, update_abundances)

    def run_fba(self, media=None, pfba=False, fva_reactions=None):
        print("pfba =", pfba)
        if media is not None:  self.util.add_medium(media)
        return self._set_solution(self.util.run_fba(None, pfba, fva_reactions))

    def _compute_relative_abundance_from_solution(self, solution=None, skipNoGrowth=True, update_abundances=False):
        solution = solution or self.solution
        total_growth = sum([solution.fluxes[member.primary_biomass.id] for member in self.members])
        message = f"The total community growth is {total_growth}"
        if self.printing:  print(message)
        if isclose(0, total_growth, abs_tol=1e-3):
            if not skipNoGrowth:   NoFluxError(f"No community growth: {total_growth} in {self.util.model.id}")
            else:    print(message)  ;  return None
        abundances = {member.id: solution.fluxes[member.primary_biomass.id]/total_growth for member in self.members}
        if update_abundances:
            self.set_abundance(abundances)
            if self.printing:  print(f"Updated abundances: {self.abundances}")
        self.comm_growth = sum([solution.fluxes[mem.primary_biomass.id]*mem.abundance for mem in self.members])
        return abundances

    def _set_solution(self, solution):
        if solution.status != "optimal":
            FeasibilityError(f'The solution is sub-optimal, with a(n) {solution} status.')
            self.solution = None
            self.print_lp()
            save_matlab_model(self.util.model, self.util.model.name + ".mat")
        self.solution = solution
        self.member_fluxes, self.memGrowths = {}, {}
        self.comm_growth = 0
        for mem in self.members:
            self.member_fluxes[mem.id] = array([solution.fluxes[rxn.id] for rxn in mem.reactions])
            self.memGrowths[mem.id] = self.solution.fluxes[mem.primary_biomass.id]
            self.comm_growth += self.memGrowths[mem.id] * mem.abundance
        if self.printing:
            print("Member Biomass Fluxes:", self.memGrowths)  # TODO weight the member growths by their abundances
            print("Max member Fluxes:", {memID: growth*self.kinCoef for memID, growth in self.memGrowths.items()})
            print("Total fluxes:", {memID: sum(abs(fluxes)) for memID, fluxes in self.member_fluxes.items()})
        # logger.info(self.util.model.summary())
        return self.solution

    def return_member_models(self):
        ## applicability in disaggregating Filipe's Nitrate reducing community model for the SBI ENIGMA team.
        return [member.model for member in self.members]

    def add_medium(self, media):
        self.util.add_medium(media)
        
