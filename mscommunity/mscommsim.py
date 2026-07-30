# -*- coding: utf-8 -*-
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.exceptions import ObjectAlreadyDefinedError, FeasibilityError, NoFluxError
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from modelseedpy.core.msatpcorrection import MSATPCorrection
from mscommunity.commhelper import build_from_species_models
from mscommunity.commkineticpkg import CommKineticPkg
from mscommunity.mscommviz import interactions as mscommsim_interactions
from mscommunity.batched_lp import (
    BatchedSolution,
    CommunityProblem,
    LPInstance,
    get_batched_solver,
    media_to_bounds,
)
from cobra.io import save_matlab_model, write_sbml_model
from itertools import combinations, permutations
from cobra.core.dictlist import DictList
from collections import Counter
from optlang.symbolics import Zero
from cobra.flux_analysis import pfba
from cobra import Reaction, Model
from numpy import array, logspace, linspace
from os import makedirs, path, environ
from math import isclose, isnan
from icecream import ic
from pandas import DataFrame
from pprint import pprint
from math import exp
import logging

logger = logging.getLogger(__name__)
ic.configureOutput(includeContext=True)

# QP-capable cobra solvers -> optlang.available_solvers capability key. Used by the
# community split-determinizer (min sum mu_i^2) to pick a backend: Gurobi/CPLEX are
# exact; 'hybrid' (HiGHS+OSQP) and 'osqp' are the strictly-open-source path (ADMM,
# ~1e-3 accuracy). GLPK/SciPy reject quadratic objectives outright.
_QP_CAPABLE = {"gurobi": "GUROBI", "cplex": "CPLEX", "hybrid": "OSQP", "osqp": "OSQP"}


def _pick_qp_backend(current=None, prefer=None):
    """Name of an available QP-capable cobra solver: an explicit `prefer`, then the
    `current` solver if it is already QP-capable (avoids a swap), then commercial
    (gurobi/cplex), then strictly open-source (hybrid=HiGHS+OSQP, then osqp). Returns
    None if no QP backend is installed."""
    import optlang
    avail = optlang.available_solvers
    order = ([prefer] if prefer else []) + ([current] if current in _QP_CAPABLE else []) \
            + ["gurobi", "cplex", "hybrid", "osqp"]
    seen = set()
    for s in order:
        if s in seen:
            continue
        seen.add(s)
        if avail.get(_QP_CAPABLE.get(s, ""), False):
            return s
    return None


def _select_biomass_cpd(entry, preferred, fallback):
    """Resolve one member-biomass metabolite from a `member_biomass_cpds` note entry.

    `commhelper.build_from_species_models` stores that note as
    {model_id: LIST of biomass metabolites} (it setdefault(...).append()s every
    cpd11416/"biomass"-named metabolite it renames), and those metabolites belong
    to the pre-copy model. `preferred` maps id -> metabolite of the live community
    model for the member-biomass compounds that feed the primary biomass;
    `fallback` maps id -> any metabolite of the live model. Returns None when the
    entry matches nothing."""
    entries = list(entry) if isinstance(entry, (list, tuple, set)) else [entry]
    hits = [preferred[cpd.id] for cpd in entries if cpd.id in preferred]
    if not hits:  hits = [fallback[cpd.id] for cpd in entries if cpd.id in fallback]
    return hits[0] if hits else None


class CommunityMember:
    def __init__(self, community, biomass_cpd, ID=None, index=None, abundance=0, model=None):
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
        if not model:
            self.model = Model()
        atp_hydrolysis_rxnComp = f"c{self.index}"
        try:
            self.atp_hydrolysis = self.community.util.model.reactions.get_by_id(f"rxn00062_{atp_hydrolysis_rxnComp}")
            print(f"skipping atp hydrolysis rxn00062_{atp_hydrolysis_rxnComp} for {self.id}")
        except:
            atp_rxn = self.community.util.add_atp_hydrolysis(atp_hydrolysis_rxnComp)
            self.atp_hydrolysis = atp_rxn["reaction"]
            if not model:
                self.model.add_reactions([self.atp_hydrolysis.copy()])
            print(f"created atp hydrolysis reaction rxn00062_{atp_hydrolysis_rxnComp} for {self.id}")
        self.biomass_drain = self.primary_biomass = None
        if not model:
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
        self.model.add_reactions([rxn.copy() for rxn in self.reactions + [self.primary_biomass, self.biomass_drain]])
        self.model.add_reactions([rxn.copy() for rxn in self.community.util.exchange_list()])
        self.model.medium = self.community.util.model.medium
        self.model.objective = self.model.reactions.get_by_id(self.primary_biomass.id).flux_expression

    def disable_species(self):
        for reaction in self.community.util.model.reactions:
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
                 flux_limit=300, probs=None, climit=None, o2limit=None, lp_filename=None, printing=False, eleLimits=None, ID=None):
        # `flux_limit` is UNUSED: nothing in the package reads it. It is kept only so
        # existing positional/keyword callers do not break; pass it or not, it has no
        # effect. `probs` was a mutable {} default shared by every instance (mutating
        # one community's rxnProbs would silently edit the default) -> None sentinel.
        assert model is not None or member_models is not None, "Either the community model and the member models must be defined."
        probs = probs if probs is not None else {}
        self.lp_filename = lp_filename
        self.printing = printing
        self.gapfillings = {}

        #Define Data attributes as None
        for attr in ['solution', 'biomass_cpd', 'primary_biomass', 'biomass_drain', 'threshold', 'msgapfill', 
                     'element_uptake_limit', 'msdb_path', 'comm_growth', 'threshold', "memGrowths", "member_fluxes",
                     "suboptimal_solution"]:
            setattr(self, attr, None)
        self.kinCoef = kinetic_coeff
        # Determinism instrumentation: whether the last determining coculture
        # solve had to fall back from pFBA/QP to the plain (solver-degenerate)
        # LP, and how many times that has happened on this community. Callers
        # read these to flag rows whose per-member split is NOT converged.
        self.pfba_fell_back = False
        self.pfba_fallback_count = 0
        # defining the models
        if model is None and member_models is not None:
            model = build_from_species_models(member_models, abundances=abundances, printing=printing)
        self.id = ID or model.id
        # self.modelID_names = model.notes["modelID_names"]
        self.util = MSModelUtil(model, True, None, climit, o2limit)
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
        
        # lookups into the LIVE (possibly copied) community model, since
        # model.notes["member_biomass_cpds"] holds pre-copy metabolite objects
        memberBioCPDs = {cpd.id: cpd for cpd in other_biomass_cpds}
        allModelCPDs = {met.id: met for met in self.util.model.metabolites}
        if ids is None:
            if member_models is not None:   ids = [mem.id for mem in member_models]
            else:  ids = [f"Species{i}" for i in range(len(other_biomass_cpds))]
        # memberIDs_biomass = dict(zip(ids,
        if not abundances:
            if member_models is None:
                abundances = {ids[memIndex]: {"biomass_compound": bioCpd, "abundance": 1/len(other_biomass_cpds)}
                              for memIndex, bioCpd in enumerate(other_biomass_cpds)}
            else:
                abundances = {}
                for memID, bioCPDs in model.notes["member_biomass_cpds"].items():
                    abundances[memID] = {"abundance": 1/len(other_biomass_cpds)}
                    met = _select_biomass_cpd(bioCPDs, memberBioCPDs, allModelCPDs)
                    if met is None:   print(f"The {memID} bioCPD was not captured")
                    else:  abundances[memID].update({"biomass_compound": met})
        elif "abundance" not in list(abundances.values())[0]:
            abundances = {memID:{"abundance": abund,
                                 "biomass_compound": _select_biomass_cpd(
                                     model.notes["member_biomass_cpds"][memID], memberBioCPDs, allModelCPDs)}
                            for memID, abund in abundances.items()}

        # print()   # this returns the carriage after the tab-ends in the biomass compound printing
        self.members = DictList(CommunityMember(self, info["biomass_compound"], ID, index+1, info["abundance"])
                                for index, (ID, info) in enumerate(abundances.items()))
        # self.members = DictList(
        #     CommunityMember(community=self, biomass_cpd=biomass_cpd, name=ids[memIndex], abundance=abundances[memIndex])
        #     for memIndex, biomass_cpd in enumerate(other_biomass_cpds))
        self.set_abundance(abundances)

        # assign the MSCommunity constraints and objective
        self.rxnProbs = probs
        self.pkgmgr = MSPackageManager.get_pkg_mgr(self.util.model)
        kinetic_pkg = CommKineticPkg(self.util.model)
        self.pkgmgr.addpkgobj(kinetic_pkg)
        kinetic_pkg.build_package(kinetic_coeff, self, self.rxnProbs)
        if eleLimits is not None:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(eleLimits)
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

    def set_objective(self, target=None, targets=None, weights=None, minimize=False):
        targets = targets or [self.util.model.reactions.get_by_id(target or self.primary_biomass.id).flux_expression]
        if weights is not None:   targets = [t*w for t, w in zip(targets, weights)]
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
        # a bare relative filename has no dirname; makedirs("") raises FileNotFoundError
        if path.dirname(filename):  makedirs(path.dirname(filename), exist_ok=True)
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
        growth_multiple = 0.999 * exp(-0.03 * len(self.members))
        if self.printing:  print(f"The growth multiple is {growth_multiple}")
        return growth_multiple

    def regularization(self, linear=True, growth_constraint=True, batch_backend="cpu", batch_workers=1):
        self.util.remove_constraint("_regularization")
        # PATCH 5: the per-member solo capacities MUST be measured with the community
        # growth floor suspended. min_comm_growth pins bio1 >= f*max, but bio1 consumes
        # EVERY member's biomass metabolite and nothing else produces them, so a solo
        # LP -- which pins the other members' biomass reactions to (0, 0) -- forces
        # v_bio1 = 0 and contradicts any strictly positive floor: every solo LP came
        # back infeasible and every member's capacity was recorded as garbage/zero, so
        # no per-member floor was ever installed. Measure first (and defensively drop a
        # floor left over from an earlier call, restoring it via the model context),
        # then install the floor for the ratio loop and the determining solve below.
        solo_max = None
        if linear:
            leftoverFloor = [cons for cons in self.util.model.constraints if cons.name == "min_comm_growth"]
            with self.util.model:
                if leftoverFloor:  self.util.model.remove_cons_vars(leftoverFloor)
                solo_max = self._solo_max_batch(backend=batch_backend, workers=batch_workers)
            # the context manager restores any suspended floor on exit
        # PATCH 1: min_comm_growth binds bio1 (community biomass), so its LB
        # must come from the max of bio1 — not the max of whatever objective
        # the caller currently has set (e.g. sum-of-member-biomasses).
        if growth_constraint:
            self.util.remove_constraint("min_comm_growth")
            ogObj_for_bio1 = self.util.model.objective
            self.util.model.objective = self.util.model.problem.Objective(
                self.primary_biomass.flux_expression, direction="max")
            commMax_bio1 = self.util.model.slim_optimize()
            self.util.model.objective = ogObj_for_bio1
            if commMax_bio1 is not None and not isnan(commMax_bio1) and commMax_bio1 > 1e-6:
                self.util.create_constraint(self.util.model.problem.Constraint(
                    self.primary_biomass.flux_expression, name="min_comm_growth",
                    lb=commMax_bio1 * self._growth_fraction()), printing=True)
        if linear:
            # PATCH 4: relative (per-member) regularization. The previous
            # scheme bounded |mu_i - mu_j| <= threshold absolutely, which
            # forced near-equal growth and collapsed real asymmetries
            # (parasitism / commensalism / dominance). Instead, compute each
            # member's solo capacity (what it could reach if the others were
            # forced to zero growth in the same community model) and require
            # each viable member to keep at least `ratio` of its own solo
            # max. Members with zero solo capacity get no constraint, so a
            # surviving member can grow at its full rate while a non-viable
            # partner stays at zero.
            # The N solo-max LPs all share the community S — only bounds and
            # objective differ — so route them through the batched solver. They are
            # solved above, BEFORE min_comm_growth is installed (see PATCH 5).
            # Iterate from tight to loose so we land on the strictest feasible
            for ratio in [0.7, 0.5, 0.3, 0.2, 0.1, 0.05]:
                self.util.remove_constraint("_regularization")
                for mem in self.members:
                    if solo_max[mem.id] > 1e-6:
                        consName = f"{mem.id}_regularization"
                        coef = {mem.primary_biomass.forward_variable: 1}
                        self.util.create_constraint(self.util.model.problem.Constraint(
                            Zero, name=consName, lb=ratio * solo_max[mem.id]), coef=coef, printing=True)
                # BUG FIX (Layer 0 #2): we only need to know whether this ratio is
                # FEASIBLE, not a vertex of the degenerate max-sum face. The old
                # `.optimize()` materialised a full Solution whose per-member flux
                # split is solver-dependent (GLPK vs Gurobi pick different vertices),
                # and in borderline cases that vertex pick could nudge which ratio
                # bucket is accepted. slim_optimize() reads only the objective /
                # status, so the ratio chosen is a pure feasibility decision and the
                # floors handed downstream are solver-independent.
                self.util.model.slim_optimize()
                if self.util.model.solver.status == "optimal":
                    print(f"The model {self.util.model.id} is regularized: each viable member keeps >={ratio*100:.0f}% of its solo capacity")
                    break
        else:
            # TODO create the least-squares method employed in MICOM
            pass
        commCurrent = self.util.model.slim_optimize()
        ic(f"Community growth after regularization ({commCurrent})")
        return None

    def micom(self, media, tradeoff=0.6):
        media = [media] if type(media) == dict else media
        # The second stage minimises sum(mu_i^2), so a QP-capable backend is required.
        # Hard-coding "hybrid" here blew up before anything else ran whenever the
        # incumbent already was QP-capable (a Gurobi model carries LP-method settings
        # that optlang's hybrid interface rejects: "LP Method primal is not valid").
        # Keep the incumbent when it can already do QP, otherwise swap to the best
        # available backend and always swap back.
        ogSolver = self.util.model.solver.interface.__name__.rsplit(".", 1)[-1].replace("_interface", "")
        qpSolver = _pick_qp_backend(current=ogSolver)
        if qpSolver is None:
            raise RuntimeError("micom() minimizes a quadratic objective, but no QP-capable solver "
                               f"is available (any of {sorted(_QP_CAPABLE)}); install gurobi, cplex, "
                               "or osqp (the 'hybrid' HiGHS+OSQP path).")
        swapped = qpSolver != ogSolver
        # the incoming (linear) objective, captured as interface-independent coefficients
        ogObjCoefs = {rxn.id: rxn.objective_coefficient for rxn in self.util.model.reactions
                      if rxn.objective_coefficient}
        ogObjDir = self.util.model.objective.direction
        if swapped:  self.util.model.solver = qpSolver
        try:  solutions = self._micom_inner(media, tradeoff)
        finally:
            # The min-sum(mu^2) objective is still installed here, and a quadratic
            # objective cannot be cloned into a non-QP interface — restoring the solver
            # first would raise "GLPK only supports linear objectives". So put the
            # caller's objective back while the QP backend is still live, then swap.
            ogObjExpr = sum([coef * self.util.model.reactions.get_by_id(rxnID).flux_expression
                             for rxnID, coef in ogObjCoefs.items()]) if ogObjCoefs else Zero
            self.util.model.objective = self.util.model.problem.Objective(ogObjExpr, direction=ogObjDir)
            if swapped:  self.util.model.solver = ogSolver
        return solutions

    def _micom_inner(self, media, tradeoff):
        solutions = []
        for m in media:
            # add the objectives
            if self.abundances is None:
                self.predict_abundances()
            # obj = sum([mem.primary_biomass.forward_variable*mem.abundance for mem in self.members])
            obj = self.primary_biomass.forward_variable
            self.util.model.objective = self.util.model.problem.Objective(obj, direction="max")
            ## maximize the community growth
            sol = self.run_fba(m, pfba=False)
            ### set tradeoff
            biomass = sol.fluxes[self.primary_biomass.id]
            # biomass = sum([sol.fluxes[mem.primary_biomass.id]*mem.abundance for mem in self.members])
            tradeoff_growth = biomass*tradeoff
            ic(biomass, tradeoff)
            self.util.remove_constraint("comm_tradeoff")
            self.util.create_constraint(self.util.model.problem.Constraint(
                self.primary_biomass.forward_variable, name="comm_tradeoff",
                ub=None, lb=tradeoff_growth), printing=True)
            ## minimize the sum of the absolute differences between the member growths
            obj = sum([mem.primary_biomass.forward_variable**2 for mem in self.members])
            self.util.model.objective = self.util.model.problem.Objective(obj, direction="min")
            sol = self.run_fba(m, pfba=False)
            # TODO hard-code a subsequent pFBA optimization that fixes biomasses of all members 
            solutions.append(sol)
        return solutions


    # TODO evaluate the comparison of this method with MICOM
    # TODO implement a community tradeoff and then the objective should be a fraction of the total summed biomass with a minimization of variance
    def predict_abundances(self, media=None, pfba=True, timeout=60,
                           environName=None, regularization=True, update_abundances=False,
                           determinize=False, qp_backend=None):
        print("regularization", regularization)
        # store the original parameters
        ogObj = self.util.model.objective
        ogMedia = self.util.model.medium
        ogTimeout = self.util.model.solver.configuration.timeout
        # apply environment media BEFORE regularization so constraints are calibrated correctly
        if media is not None:
            self.util.add_medium(media)
        slimOpt = self.util.model.slim_optimize()
        # PATCH 3: if kinetic constraints zero out growth on this medium, fall
        # back to a no-kinetics simulation rather than returning None. This
        # rescues low-yield carbons (e.g. Acetate) where Σ|flux| ≤ kinCoef·bio
        # is unsatisfiable at any positive biomass for kinCoef=750. Use the
        # cobra model's context manager so removed constraints auto-restore
        # on exit (so the package's internal tracking stays consistent).
        if isclose(0, slimOpt, abs_tol=1e-3):
            kin_cons = [c for c in self.util.model.constraints if "_commKin" in c.name]
            if kin_cons:
                with self.util.model:
                    self.util.model.remove_cons_vars(kin_cons)
                    slimOpt_no_kin = self.util.model.slim_optimize()
                    if not isclose(0, slimOpt_no_kin, abs_tol=1e-3):
                        print(f"Kinetic constraints disabled for {self.util.model.id} on {environName}: "
                              f"slim_optimize was {slimOpt:.3g}, now {slimOpt_no_kin:.3g}")
                        return self._predict_inner(pfba, timeout, regularization, update_abundances,
                                                    ogObj, ogMedia, ogTimeout, determinize, qp_backend)
                # context exits here, kinetic constraints auto-restored
                print(f"\nThe model {self.util.model.id} doesn't grow even without kinetics on {environName}")
                self.util.model.objective = ogObj
                self.util.model.medium = ogMedia
                self.util.model.solver.configuration.timeout = ogTimeout
                return None
            else:
                print(f"\nThe model {self.util.model.id} doesn't grow, with a slim_optimize of {slimOpt} in {environName} media")
        return self._predict_inner(pfba, timeout, regularization, update_abundances,
                                    ogObj, ogMedia, ogTimeout, determinize, qp_backend)

    def _predict_inner(self, pfba, timeout, regularization, update_abundances,
                       ogObj, ogMedia, ogTimeout, determinize=False, qp_backend=None):
        # maximize the sum of all member biomass reactions
        self.set_objective(targets=[species.primary_biomass.forward_variable for species in self.members])
        self.util.model.solver.configuration.timeout = timeout
        if regularization:   threshold = self.regularization(linear=True)
        else:   self.util.remove_constraint("_regularization")
        # BUG FIX (Layer 0 #1): the old `except: ... pfba=False` silently reverted
        # to the plain, solver-degenerate LP with no trace — masking exactly the
        # hard pairs (and it would swallow a future QP-backend failure too). Catch
        # only real solve errors, log the cause, and record the fallback so callers
        # can flag the row as NOT determinism-converged.
        self.pfba_fell_back = False
        try:
            sol = self.run_fba(None, pfba)
        except Exception as e:
            self._note_fallback(
                "pFBA failed for %s (%s: %s); falling back to plain LP — this "
                "row's per-member split is NOT determinism-converged",
                self.util.model.id, type(e).__name__, e)
            sol = self.run_fba(None, pfba=False)
        # Determinize the per-member split to the unique strictly-convex-QP optimum
        # while the regularization floors are still installed. No-op on the ~98.8%
        # of cells whose split is already unique; on the flat-face ~1.2% it replaces
        # the arbitrary LP vertex with the solver-independent centroid.
        if determinize:
            sol = self._determinize_split(sol, qp_backend)
        self.util.remove_constraint("_regularization")
        self.util.remove_constraint("min_comm_growth")
        self.util.model.solver.configuration.timeout = ogTimeout
        self.util.model.objective = ogObj
        self.util.model.medium = ogMedia
        return self._compute_relative_abundance_from_solution(sol, True, update_abundances)

    def _note_fallback(self, message, *args):
        """Record that THIS row's per-member split is not determinism-converged.

        `pfba_fell_back` is the per-row flag (reset at the top of `_predict_inner`)
        and `pfba_fallback_count` counts ROWS, not events: a row that both loses pFBA
        and fails to determinize is counted once. Previously only `_predict_inner`'s
        pFBA handler touched the counter, so the determinizer's own failure paths were
        invisible and the count under-reported the true non-converged rate."""
        if not self.pfba_fell_back:
            self.pfba_fell_back = True
            self.pfba_fallback_count += 1
        logger.warning(message, *args)

    def _determinize_split(self, sol, qp_backend=None, tol=1e-6, fix_fluxes=True):
        """Refine the per-member coculture-growth split to the UNIQUE strictly-convex
        QP optimum so it is solver-independent. Called from _predict_inner while the
        regularization floors + sum-of-member objective are still installed.

        On the deterministic-construction base ~98.8% of (pair, medium) cells already
        have a unique max-sum-member split; ~1.2% sit on a flat alternate-optima face
        where GLPK and Gurobi pick different vertices (Δ up to ~3.4). A cheap 2-LP
        probe on member[0]'s biomass detects the flat face; only then do we minimise
        sum(mu_i^2) on a QP backend (Gurobi/CPLEX exact, else open-source OSQP/hybrid)
        to select the unique centroid, validate it is feasible on the exact native
        solver, then pin mu* and pFBA on the native solver for a deterministic full
        flux vector. Returns the (possibly new) solution; sets pfba_fell_back if a QP
        backend is needed but unavailable or returns an infeasible point."""
        members = list(self.members)
        if len(members) < 2:
            return sol
        model = self.util.model
        native = model.solver.interface.__name__.rsplit(".", 1)[-1].replace("_interface", "")
        sum_vars = [m.primary_biomass.forward_variable for m in members]
        try:
            C = float(sum(sol.fluxes[m.primary_biomass.id] for m in members))
        except Exception:
            return sol
        if isnan(C) or C <= 1e-9:
            return sol
        # Pin the community total near its max, with a tiny slack BELOW it. Some
        # degenerate cells are a NEAR-flat ridge rather than an exactly-flat face:
        # GLPK and Gurobi stop at points whose totals differ within their optimality
        # tolerance (~1e-5). An exact Σ=C pin collapses that ridge (the probe then
        # reads the split as "unique"); the slack lets the QP traverse the ridge to
        # the centroid. Growth is reduced by at most det_slack (relative ~1e-6).
        det_slack = max(1e-6 * abs(C), 1e-8)
        self.util.remove_constraint("_det_total")
        self.util.create_constraint(model.problem.Constraint(sum(sum_vars), name="_det_total", lb=C - det_slack, ub=None))
        og_obj = model.objective
        try:
            # --- degeneracy probe: range of mu_0 on the fixed-total optimal face ---
            v0 = members[0].primary_biomass.forward_variable
            model.objective = model.problem.Objective(v0, direction="min"); lo = model.slim_optimize()
            model.objective = model.problem.Objective(v0, direction="max"); hi = model.slim_optimize()
            model.objective = og_obj
            if lo is None or hi is None or isnan(lo) or isnan(hi) or (hi - lo) <= tol:
                return sol  # split already unique -> no-op (the common ~98.8% case)

            # --- strictly-convex QP for the unique centroid split ---
            # Solved by HiGHS native QP via matrix extraction (see _solve_qp_highs):
            # it does NOT swap the live model's solver, so it is safe inside the
            # PATCH-3 `with self.util.model:` context (kinetics-disabled cells) and
            # needs no QP-capable optlang backend (GLPK is LP-only). Exact + open-source.
            mu = self._solve_qp_highs(sum_vars)
            if mu is None:
                self._note_fallback("QP determinizer (HiGHS) unavailable or non-optimal "
                                    "for %s; leaving the solver-dependent LP split.",
                                    self.util.model.id)
                return sol

            # --- pin mu* on the native solver (validate feasibility) ---
            def _pin(ptol):
                for m in members:
                    self.util.remove_constraint(f"_det_fix_{m.id}")
                    self.util.create_constraint(model.problem.Constraint(
                        m.primary_biomass.forward_variable, name=f"_det_fix_{m.id}",
                        lb=mu[m.id] - ptol, ub=mu[m.id] + ptol))
                model.slim_optimize()
                return model.solver.status == "optimal"
            if not _pin(1e-6) and not _pin(1e-4):
                self._note_fallback("QP split infeasible on native solver for %s; "
                                    "keeping LP split.", self.util.model.id)
                return sol

            # --- pFBA on native with mu* pinned -> deterministic full flux vector ---
            if fix_fluxes:
                try: final = self.util.run_fba(None, True)
                except Exception: final = model.optimize()
            else:
                final = model.optimize()
            return final
        finally:
            for m in members:
                self.util.remove_constraint(f"_det_fix_{m.id}")
            self.util.remove_constraint("_det_total")
            try: model.objective = og_obj
            except Exception: self.set_objective(targets=sum_vars)

    def _solve_qp_highs(self, hess_vars):
        """min sum(v**2 for v in hess_vars) subject to ALL current model constraints
        and variable bounds, via HiGHS native QP (highspy) by extracting the problem
        matrices — WITHOUT swapping the live model's solver. This is what makes the
        determinizer safe inside the PATCH-3 `with model:` context and independent of
        whether the cobra solver supports quadratics (GLPK does not). Exact and
        strictly open-source. Returns {member_id: mu} or None if highspy is
        unavailable or the QP is not optimal."""
        try:
            import highspy
            import numpy as np
        except Exception:
            return None
        model = self.util.model
        INF = highspy.kHighsInf
        variables = list(model.solver.variables)
        vidx = {v.name: i for i, v in enumerate(variables)}
        n = len(variables)
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("log_to_console", False)
        # Hard wall-clock cap so a pathological/degenerate community QP can never
        # hang the whole scoring run (some KBase-reconstructed communities sent
        # HiGHS into a non-terminating solve). On time-out the status is not
        # kOptimal, so we return None below and fall back to the LP split
        # (flagged via pfba_fell_back) instead of blocking forever.
        h.setOptionValue("time_limit", float(environ.get("QP_TIME_LIMIT", "20")))
        for v in variables:
            h.addVar(-INF if v.lb is None else float(v.lb),
                      INF if v.ub is None else float(v.ub))
        for con in model.solver.constraints:
            items = list(con.get_linear_coefficients(list(con.variables)).items())
            if not items:
                continue
            idx = np.array([vidx[v.name] for v, _ in items], dtype=np.int32)
            val = np.array([float(cf) for _, cf in items], dtype=np.float64)
            lo = -INF if con.lb is None else float(con.lb)
            up =  INF if con.ub is None else float(con.ub)
            h.addRow(lo, up, len(idx), idx, val)
        # objective = sum v^2 = 0.5 x^T Q x with Q_ii = 2 on the member biomass vars
        diag = {vidx[v.name] for v in hess_vars}
        q_start = np.zeros(n + 1, dtype=np.int32)
        q_idx, q_val = [], []
        for col in range(n):
            q_start[col] = len(q_idx)
            if col in diag:
                q_idx.append(col); q_val.append(2.0)
        q_start[n] = len(q_idx)
        try:
            h.passHessian(n, len(q_idx), highspy.HessianFormat.kTriangular,
                          q_start, np.array(q_idx, dtype=np.int32), np.array(q_val, dtype=np.float64))
            h.run()
        except Exception as e:
            logger.warning("HiGHS QP solve failed for %s: %s", self.util.model.id, e)
            return None
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            return None
        cv = h.getSolution().col_value
        return {m.id: float(cv[vidx[m.primary_biomass.forward_variable.name]]) for m in self.members}

    def run_fba(self, media=None, pfba=False, fva_reactions=None):
        # print("pfba =", pfba)
        if media is not None:
            self.util.add_medium(media)
        return self._set_solution(self.util.run_fba(None, pfba, fva_reactions))

    def _comm_growth(self):
        self.comm_growth = 0
        self.member_fluxes, self.memGrowths = {}, {}
        for mem in self.members:
            self.member_fluxes[mem.id] = array([self.solution.fluxes[rxn.id] for rxn in mem.reactions])
            self.memGrowths[mem.id] = self.solution.fluxes[mem.primary_biomass.id]
            self.comm_growth += self.memGrowths[mem.id] * mem.abundance

    def _compute_relative_abundance_from_solution(self, solution=None, skipNoGrowth=True, update_abundances=False):
        if solution is not None:  self._set_solution(solution)
        total_growth = sum([self.solution.fluxes[member.primary_biomass.id] for member in self.members])
        message = f"The total community growth is {total_growth}"
        if self.printing:  ic(message)
        if isclose(0, total_growth, abs_tol=1e-3):
            if not skipNoGrowth:   NoFluxError(f"No community growth: {total_growth} in {self.util.model.id}")
            else:    print(message)  ;  return None
        abundances = {member.id: self.solution.fluxes[member.primary_biomass.id]/total_growth for member in self.members}
        if update_abundances:
            self.set_abundance(abundances)
            if self.printing:  print(f"Updated abundances: {self.abundances}")
        for mem in self.members:
            ic(f"{mem.id} grows {self.solution.fluxes[mem.primary_biomass.id]} with abundance {mem.abundance}")
        return abundances

    def _set_solution(self, solution, dump_diagnostics=None):
        if solution.status != "optimal":
            # This used to be silent AND fatal at once: the FeasibilityError was
            # constructed but never raised (so a sub-optimal solution was consumed as
            # if optimal), while the unconditional diagnostic dump crashed on the bare
            # relative "erronous_model.lp" (print_lp -> makedirs("")). Be loud instead,
            # and only dump when the caller opted in through `lp_filename` (or asked
            # explicitly). Logging rather than raising keeps predict_abundances'
            # documented tolerance of non-growing media intact -- callers detect the
            # bad row through the NaN/zero growth that _compute_relative_abundance_
            # from_solution already screens for.
            self.suboptimal_solution = True
            logger.error("The %s solution is sub-optimal, with a(n) %s status; the fluxes and "
                         "abundances derived from it are NOT trustworthy.",
                         self.util.model.id, solution.status)
            if dump_diagnostics or (dump_diagnostics is None and self.lp_filename):
                dumpDir = path.dirname(self.lp_filename) if self.lp_filename else ""
                try:
                    self.print_lp(path.join(dumpDir, "erronous_model.lp"))
                    save_matlab_model(self.util.model, path.join(
                        dumpDir, f"{self.util.model.name or self.util.model.id}.mat"))
                except Exception as e:
                    logger.warning("Could not dump diagnostics for %s: %s", self.util.model.id, e)
        else:  self.suboptimal_solution = False
        self.solution = solution
        self.exchange_fluxes = {ex.id: self.solution.fluxes[ex.id] for ex in self.util.model.reactions if "EX_" in ex.id}
        self._comm_growth()
        if self.printing:
            ic("Member Biomass Fluxes:", self.memGrowths)  # TODO weight the member growths by their abundances
            ic("Max member Fluxes:", {memID: growth*self.kinCoef for memID, growth in self.memGrowths.items()})
            ic("Total fluxes:", {memID: sum(abs(fluxes)) for memID, fluxes in self.member_fluxes.items()})
        # logger.info(self.util.model.summary())
        return self.solution

    def return_member_models(self):
        ## applicability in disaggregating Filipe's Nitrate reducing community model for the SBI ENIGMA team.
        return [member.model for member in self.members]

    def add_medium(self, media):
        self.util.add_medium(media)

    # --- Batched-LP entry points -----------------------------------------
    # The community S matrix is fixed across samples / conditions; only
    # bounds and objective vary. `extract_problem` snapshots that shared
    # structure once and `solve_batch` routes a list of `LPInstance`s
    # through a pluggable backend (cpu reference today, gpu later) without
    # any call-site changes.

    def extract_problem(self):
        """Snapshot the live community model as a `CommunityProblem`."""
        return CommunityProblem.from_model(self.util.model)

    def solve_batch(self, instances, backend="cpu", problem=None, **backend_kwargs):
        """Solve N independent LPs that share this community's S matrix."""
        if problem is None:
            problem = self.extract_problem()
        solver = get_batched_solver(backend, **backend_kwargs)
        return solver.solve(problem, instances)

    def _solo_max_batch(self, backend="cpu", workers=1):
        """Per-member solo growth max, batched over members.

        Each member's LP zeros out the other members' biomass reactions and
        maximizes its own — same S, different bounds and objective.
        """
        instances = []
        for target in self.members:
            bounds = {}
            for other in self.members:
                if other.id != target.id:
                    bounds[other.primary_biomass.id] = (0.0, 0.0)
            instances.append(LPInstance(
                id=target.id,
                bounds=bounds,
                objective={target.primary_biomass.id: 1.0},
                sense="max",
            ))
        results = self.solve_batch(instances, backend=backend, workers=workers)
        solo = {}
        for r in results:
            # The status is the ONLY trustworthy signal here: on a non-optimal solve
            # GLPK leaves the objective expression evaluated at the STALE incumbent of
            # whatever was solved before (a positive number), while Gurobi returns
            # None -- so reading objective_value alone made the recorded solo capacity
            # solver-dependent garbage. Anything not optimal means "no capacity".
            if r.status != "optimal":
                logger.warning("The solo-capacity LP of %s in %s is %s; recording zero capacity "
                               "(its objective_value %r is not meaningful).",
                               r.id, self.util.model.id, r.status, r.objective_value)
                solo[r.id] = 0
                continue
            v = r.objective_value
            solo[r.id] = v if (v is not None and not isnan(v) and v > 1e-6) else 0
        return solo

    def predict_abundances_batch(self, medias, pfba=False, backend="cpu", workers=1):
        """Predict member abundances across many media in one batched solve.

        Each medium becomes an `LPInstance` whose bounds patch the exchange
        reactions of the shared community model. The objective is the sum
        of member primary biomass fluxes (matches `_predict_inner`).
        """
        problem = self.extract_problem()
        obj_patch = {mem.primary_biomass.id: 1.0 for mem in self.members}
        instances = []
        for i, media in enumerate(medias):
            bounds = media_to_bounds(self.util.model, media)
            instances.append(LPInstance(
                id=f"sample_{i}",
                bounds=bounds,
                objective=obj_patch,
                sense="max",
                pfba=pfba,
            ))
        results = self.solve_batch(instances, backend=backend, problem=problem, workers=workers)
        out = []
        for r in results:
            if r.status != "optimal":
                out.append(None)
                continue
            total = sum(r.fluxes.get(mem.primary_biomass.id, 0.0) for mem in self.members)
            if isclose(0, total, abs_tol=1e-3):
                out.append(None)
                continue
            out.append({
                mem.id: r.fluxes.get(mem.primary_biomass.id, 0.0) / total
                for mem in self.members
            })
        return out

