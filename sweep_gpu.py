#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GPU batch-size scaling sweep: characterize jax-pdhg per-LP cost vs batch
width B, to locate the crossover against cobra's ~0.45 s/LP.

Loads the model ONCE, extracts the shared problem ONCE, then times the GPU
batched solver at several B with a fixed iteration budget. Reports per-LP and
per-iteration cost so we can extrapolate a converging-solver crossover.
"""
import sys, os, time
import numpy as np
import cobra
from mscommunity.batched_lp import CommunityProblem, LPInstance, get_batched_solver

ITERS = int(os.environ.get("ITERS", "2000"))
model_path = sys.argv[1]
batch_sizes = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else
                                ["16", "64", "256", "1024", "4096"])]

t0 = time.perf_counter()
model = cobra.io.load_json_model(model_path)
print(f"[load] {model.id} in {time.perf_counter()-t0:.1f}s "
      f"({len(model.reactions)} rxns, {len(model.metabolites)} mets)")
problem = CommunityProblem.from_model(model)
exch = [r.id for r in model.reactions if r.id.startswith("EX_")][:200]
rng = np.random.default_rng(0)

def make_instances(B):
    out = []
    for i in range(B):
        s = float(rng.uniform(0.5, 1.5))
        out.append(LPInstance(id=f"s{i}", bounds={ex: (-10.0*s, 1000.0) for ex in exch},
                              objective={"bio1": 1.0}, sense="max"))
    return out

solver = get_batched_solver("jax-pdhg", max_iters=ITERS, tol=1e-9,
                            check_every=ITERS+1)  # disable early stop: full budget

# warm up (triggers XLA compile for the smallest shape; we still recompile per B)
_ = solver.solve(problem, make_instances(8))

print(f"\n[jax-pdhg GPU scaling] fixed {ITERS} iters, no early stop")
print(f"  {'B':>6} {'total_s':>9} {'ms/LP':>9} {'ms/iter(batch)':>15} {'vs cobra/LP':>12}")
COBRA_PER_LP_MS = 450.0
for B in batch_sizes:
    insts = make_instances(B)
    t0 = time.perf_counter()
    sols = solver.solve(problem, insts)
    dt = time.perf_counter() - t0
    ms_lp = dt / B * 1000
    ms_iter = dt / ITERS * 1000
    ratio = COBRA_PER_LP_MS / ms_lp
    print(f"  {B:>6} {dt:>9.2f} {ms_lp:>9.2f} {ms_iter:>15.3f} {ratio:>11.1f}x")

print(f"\n  (cobra reference: ~{COBRA_PER_LP_MS:.0f} ms/LP, exact, converged)")
print(f"  ratio > 1 means GPU is faster per-LP THAN cobra at this fixed iter budget")
