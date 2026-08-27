# -*- coding: utf-8 -*-

from __future__ import absolute_import

# import pyximport; pyximport.install(language_level=3)  # improve computational speed

from mscommunity.mscommsim import *
from mscommunity.commkineticpkg import CommKineticPkg
from mscommunity.mscommviz import *
from mscommunity.mskineticsfba import MSKineticsFBA
from mscommunity.batched_lp import (
    BatchedLPSolver,
    BatchedSolution,
    CPUBatchedLPSolver,
    CommunityProblem,
    LPInstance,
    get_batched_solver,
    media_to_bounds,
    register_batched_solver,
    solve_batch,
)
import mscommunity.backends  # auto-registers numpy-pdhg / jax-pdhg / cupy-pdhg / pdlp
