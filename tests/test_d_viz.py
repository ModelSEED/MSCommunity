"""Regression tests for the visualization/analysis helpers in mscommunity.mscommviz."""
import re
import shutil

import pytest
from cobra import Metabolite, Model, Reaction
from modelseedpy.core.msmodelutl import MSModelUtil
from pandas import DataFrame

from mscommunity import mscommviz


# ------------------------------------------------------------------ fixtures
class _FakeCpd:
    def __init__(self, cid, abbr, name):
        self.id, self.abbr, self.name = cid, abbr, name


class _FakeCompounds:
    def __init__(self, cpds):
        self._cpds = {cpd.id: cpd for cpd in cpds}

    def get_by_id(self, cid):
        return self._cpds[cid]


class _FakeMSDB:
    """The minimal slice of the ModelSEED database that visual_interactions consumes."""
    def __init__(self, cpds):
        self.compounds = _FakeCompounds(cpds)


class _StubMember:
    def __init__(self, ID, index, primary_biomass, biomass_cpd, reactions):
        self.id, self.index = ID, index
        self.primary_biomass, self.biomass_cpd, self.reactions = primary_biomass, biomass_cpd, reactions
        self.abundance = 0.5


class _StubCommunity:
    """Mirrors the MSCommunity API that mscommviz relies upon, sans the heavy construction."""
    def __init__(self, model, members, primary_biomass):
        self.util = MSModelUtil(model)
        self.members, self.primary_biomass = members, primary_biomass
        self.abundances_set, self.printing = True, False
        self.objective_calls = []

    def run_fba(self, media=None, pfba=False, fva_reactions=None):
        return self.util.model.optimize()

    def set_objective(self, target=None, targets=None, weights=None, minimize=False):
        self.objective_calls.append((target, targets, minimize))
        targets = targets or [self.util.model.reactions.get_by_id(
            target or self.primary_biomass.id).flux_expression]
        self.util.model.objective = self.util.model.problem.Objective(
            sum(targets), direction="min" if minimize else "max")


def _toy_community(junk_compartment=False):
    """Two members that share a single environmental substrate."""
    model = Model("toy")
    A_e = Metabolite("A_e0", compartment="e0")
    mets, rxns, memberInfo = [A_e], [], []
    ex = Reaction("EX_A_e0", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A_e: -1})
    rxns.append(ex)
    for index in (1, 2):
        A_c = Metabolite(f"A_c{index}", compartment=f"c{index}")
        bio_cpd = Metabolite(f"cpd11416_c{index}", compartment=f"c{index}")
        transport = Reaction(f"TR_A_c{index}", lower_bound=-1000, upper_bound=1000)
        transport.add_metabolites({A_e: -1, A_c: 1})
        biomass = Reaction(f"bio1_c{index}", lower_bound=0, upper_bound=1000)
        biomass.add_metabolites({A_c: -1, bio_cpd: 1})
        drain = Reaction(f"DM_cpd11416_c{index}", lower_bound=0, upper_bound=1000)
        drain.add_metabolites({bio_cpd: -1})
        mets += [A_c, bio_cpd]  ;  rxns += [transport, biomass, drain]
        memberInfo.append((index, biomass, bio_cpd, [transport, biomass, drain]))
    if junk_compartment:
        ## a reaction whose compartment carries no numeric suffix
        X_c = Metabolite("X_c", compartment="c")
        junk = Reaction("JUNK_c", lower_bound=0, upper_bound=1000)
        junk.add_metabolites({X_c: -1})
        mets.append(X_c)  ;  rxns.append(junk)
    model.add_metabolites(mets)
    model.add_reactions(rxns)
    model.objective = model.reactions.get_by_id("bio1_c1").flux_expression
    members = [_StubMember(f"Species{index}", index, biomass, bio_cpd, memRxns)
               for index, biomass, bio_cpd, memRxns in memberInfo]
    return _StubCommunity(model, members, model.reactions.get_by_id("bio1_c1"))


def _cross_feeding_df():
    """Two cross-fed metabolites whose abbreviations share a three-character prefix."""
    df = DataFrame({"Species1": [3.0, -2.0], "Species2": [-3.0, 2.0], "Environment": [1.0, 1.0]},
                   index=["cpd00027_e0", "cpd00029_e0"])
    df.index.name = "Metabolite/Donor ID"
    return df, _FakeMSDB([_FakeCpd("cpd00027", "glc-D", "D-Glucose"),
                          _FakeCpd("cpd00029", "glcnt", "D-Gluconate")])


def _rank_same_blocks(source):
    """The bodies of every rank=same subgraph in a Graphviz source."""
    blocks, current = [], None
    for line in source.splitlines():
        line = line.strip()
        if line.startswith("subgraph") and line.endswith("{"):  current = []  ;  continue
        if current is None:  continue
        if line == "}":  blocks.append(current)  ;  current = None
        else:  current.append(line)
    return [block for block in blocks if "rank=same" in block]


# ------------------------------------------------------- visual_interactions
@pytest.fixture
def graph_source(tmp_path):
    if shutil.which("dot") is None:  pytest.skip("the Graphviz `dot` executable is unavailable")
    df, msdb = _cross_feeding_df()
    return mscommviz.visual_interactions(df, filename=str(tmp_path / "cross_feeding"),
                                         msdb=msdb, view_figure=False)


def test_metabolite_nodes_are_uniquely_identified(graph_source):
    # both metabolites abbreviate to "glc", yet each must own a distinct node
    metNodes = re.findall(r"^\s*(\S+) \[label=glc .*URL=", graph_source, flags=re.MULTILINE)
    assert len(metNodes) == 2, graph_source
    assert len(set(metNodes)) == 2, f"the two `glc` metabolites collapsed into {metNodes}"
    # neither tooltip may be lost to a collision
    assert graph_source.count('tooltip="cpd00027 ; D-Glucose"') == 1
    assert graph_source.count('tooltip="cpd00029 ; D-Gluconate"') == 1
    # ... and the edges must reference those unique IDs, not the shared abbreviation
    assert " glc " not in graph_source.replace("label=glc ", "")
    for node in metNodes:
        assert re.search(rf"(-> {node} |{node} ->)", graph_source), graph_source


def test_both_member_tiers_are_rank_constrained(graph_source):
    rankedBlocks = _rank_same_blocks(graph_source)
    memberBlocks = [block for block in rankedBlocks
                    if any(line.startswith("S") for line in block)]
    assert len(memberBlocks) == 2, f"expected two member tiers, found {rankedBlocks}"
    tiers = [{re.match(r"(S\d+)", line).group(1) for line in block if line.startswith("S")}
             for block in memberBlocks]
    assert tiers[0] and tiers[1] and not tiers[0] & tiers[1]
    # no member may be declared in the parent graph, which would defeat the rank constraint
    for line in graph_source.splitlines():
        if re.match(r"^\tS\d+ \[label=", line):  pytest.fail(f"member declared outside a tier: {line}")


def test_visual_interactions_does_not_call_display(graph_source):
    # a NameError from the notebook-only `display` would have aborted the fixture
    assert graph_source.startswith("digraph")


# --------------------------------------------------- abundance variability & FBA
def test_abundance_variability_analysis_reports_member_growths():
    comm = _toy_community()
    variability = mscommviz.abundance_variability_analysis(comm, None)
    assert set(variability) == {"Species1", "Species2"}
    for memID, bounds in variability.items():
        assert set(bounds) == {"minVar", "maxVar"}
        assert all(isinstance(float(val), float) for val in bounds.values())
        assert bounds["minVar"] <= bounds["maxVar"]
    assert variability["Species2"]["maxVar"] > 0
    # the objective must be assigned through the (target=, minimize=) signature of MSCommunity
    assert [call[0] for call in comm.objective_calls if call[0] is not None] == [
        "bio1_c1", "bio1_c1", "bio1_c2", "bio1_c2"]


def test_run_fba_enforces_the_minimal_member_growth():
    comm = _toy_community()
    ## Species2 does not contribute to the objective and thus only grows when the floor is enforced
    sol = mscommviz.run_fba(comm, None, minMemGrwoth=1, compute_interactions=False)
    assert sol.fluxes["bio1_c2"] >= 1 - 1e-6
    ## and the floor is idempotent across repeated simulations
    sol = mscommviz.run_fba(comm, None, minMemGrwoth=2, compute_interactions=False)
    assert sol.fluxes["bio1_c2"] >= 2 - 1e-6


def test_minimal_member_growth_floor_does_not_leak_onto_the_model():
    """The floor is scoped to the simulation, not left on the caller's model.

    It is installed inside `with model:` so a later solve — by any other caller,
    with any other objective — is not silently constrained by a floor that some
    earlier run_fba() happened to set.
    """
    comm = _toy_community()
    mscommviz.run_fba(comm, None, minMemGrwoth=2, compute_interactions=False)
    assert not [cons for cons in comm.util.model.constraints
                if cons.name.endswith("minMemGrowth")]


def test_minimal_member_growth_is_opt_in():
    """Default is no floor, because the historical default enforced nothing.

    The old code assigned `member.biomass_cpd.lb`, an attribute cobra's
    Metabolite does not define, so the documented default of 1 never reached the
    LP. Turning it on by default would silently change every existing caller.
    """
    comm = _toy_community()
    sol = mscommviz.run_fba(comm, None, compute_interactions=False)
    assert not [cons for cons in comm.util.model.constraints
                if cons.name.endswith("minMemGrowth")]
    assert sol.fluxes["bio1_c2"] == pytest.approx(0, abs=1e-6)


def test_run_fba_reaches_the_ava_branch():
    comm = _toy_community()
    variability = mscommviz.run_fba(comm, None, ava=True, minMemGrwoth=0)
    assert set(variability) == {"Species1", "Species2"}


def test_run_fba_calls_the_interactions_function(monkeypatch):
    comm = _toy_community()
    monkeypatch.setattr(mscommviz, "visual_interactions", lambda *args, **kwargs: None)
    crossFeeding, exMets = mscommviz.run_fba(comm, None, minMemGrwoth=0)
    assert isinstance(crossFeeding, DataFrame) and isinstance(exMets, DataFrame)


# ------------------------------------------------------- compartment guarding
def test_interactions_tolerates_a_nonnumeric_compartment():
    comm = _toy_community(junk_compartment=True)
    crossFeeding, exMets = mscommviz.interactions(comm, comm.run_fba(), visualize=False)
    assert isinstance(crossFeeding, DataFrame) and isinstance(exMets, DataFrame)


def test_compartment_index_guard():
    model = Model("comps")
    X_c = Metabolite("X_c", compartment="c")
    Y_c3 = Metabolite("Y_c3", compartment="c3")
    junk = Reaction("JUNK_c", lower_bound=0, upper_bound=1000)
    junk.add_metabolites({X_c: -1})
    numbered = Reaction("NUM_c3", lower_bound=0, upper_bound=1000)
    numbered.add_metabolites({Y_c3: -1})
    model.add_metabolites([X_c, Y_c3])
    model.add_reactions([junk, numbered])
    assert mscommviz._compartment_index(junk) is None
    assert mscommviz._compartment_index(numbered) == 3
