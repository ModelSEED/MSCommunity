#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Benchmark the batched-LP / GPU acceleration on a real community model.

Usage:
    python bench_gpu.py <model.json> <backend1,backend2,...> [workload]

workload: "solo"  -> per-member solo-max batch (default; the _solo_max_batch
                     workload used by MSCommunity.regularization)
          "media" -> replicate the community objective across N media-like
                     uptake variations (stresses batch width B)

Each backend is one of: cpu, numpy-pdhg, jax-pdhg, cupy-pdhg, pdlp
The 'cpu' backend (cobra/glpk per-LP) is treated as ground truth.
"""
import sys, os, time, json
import numpy as np

MAXITERS = int(os.environ.get("MAXITERS", "20000"))
TOL = float(os.environ.get("TOL", "1e-5"))

t0 = time.perf_counter()
import cobra
from mscommunity.batched_lp import (
    CommunityProblem, LPInstance, get_batched_solver,
)
from mscommunity.batched_lp import _SOLVERS
print(f"[import] mscommunity + cobra in {time.perf_counter()-t0:.2f}s; "
      f"registered backends: {sorted(_SOLVERS)}")

model_path = sys.argv[1]
backends = sys.argv[2].split(",") if len(sys.argv) > 2 else ["cpu", "numpy-pdhg"]
workload = sys.argv[3] if len(sys.argv) > 3 else "solo"

# ---- load model -------------------------------------------------------------
t0 = time.perf_counter()
model = cobra.io.load_json_model(model_path)
print(f"[load] {model.id}: {len(model.reactions)} rxns, "
      f"{len(model.metabolites)} mets, {len(model.compartments)} comps "
      f"in {time.perf_counter()-t0:.2f}s")

# member biomass reactions = bioN with N>=2 (bio1 is community biomass)
import re
member_bios = sorted(
    (r.id for r in model.reactions if re.fullmatch(r"bio\d+", r.id) and r.id != "bio1"),
    key=lambda s: int(s[3:]),
)
print(f"[members] {len(member_bios)} member biomass rxns: {member_bios}")

# ---- extract shared problem once -------------------------------------------
t0 = time.perf_counter()
problem = CommunityProblem.from_model(model)
S = problem.S
print(f"[extract] CommunityProblem: S = {S.shape}, nnz = {S.nnz}, "
      f"extra-cons = {problem.n_extra}, dropped(non-antisym) = "
      f"{len(problem.extracted_with_warnings)} in {time.perf_counter()-t0:.2f}s")
if problem.extracted_with_warnings:
    print(f"          dropped names (first 5): {problem.extracted_with_warnings[:5]}")

# ---- build the batch --------------------------------------------------------
if workload == "solo":
    instances = []
    for tgt in member_bios:
        bounds = {other: (0.0, 0.0) for other in member_bios if other != tgt}
        instances.append(LPInstance(
            id=tgt, bounds=bounds,
            objective={tgt: 1.0}, sense="max",
        ))
else:  # media: vary all exchange uptake caps, community objective (bio1)
    exch = [r.id for r in model.reactions if r.id.startswith("EX_")]
    rng = np.random.default_rng(0)
    instances = []
    n_samples = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    for i in range(n_samples):
        scale = float(rng.uniform(0.5, 1.5))
        bounds = {ex: (-10.0 * scale, 1000.0) for ex in exch[:200]}
        instances.append(LPInstance(
            id=f"sample_{i}", bounds=bounds,
            objective={"bio1": 1.0}, sense="max",
        ))
print(f"[batch] {len(instances)} LP instances ({workload} workload)\n")

# ---- run each backend -------------------------------------------------------
results = {}
for name in backends:
    if name not in _SOLVERS:
        print(f"[{name}] NOT REGISTERED — skipping (dep not installed)")
        continue
    # backend-specific kwargs: give the first-order solvers a generous budget
    kwargs = {}
    if name.endswith("-pdhg"):
        kwargs = dict(max_iters=MAXITERS, tol=TOL, check_every=max(50, MAXITERS // 10))
    try:
        solver = get_batched_solver(name, **kwargs)
        t0 = time.perf_counter()
        sols = solver.solve(problem, instances)
        dt = time.perf_counter() - t0
    except Exception as exc:
        import traceback
        print(f"[{name}] FAILED: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        continue
    obj = {s.id: s.objective_value for s in sols}
    status = {s.id: s.status for s in sols}
    results[name] = (obj, dt, status)
    n_opt = sum(1 for v in status.values() if v == "optimal")
    print(f"[{name}] {dt:.3f}s total | {dt/len(instances)*1000:.1f} ms/LP | "
          f"{n_opt}/{len(instances)} optimal")

# ---- compare to cpu reference ----------------------------------------------
if "cpu" in results:
    ref_obj = results["cpu"][0]
    print(f"\n[compare vs cpu reference]")
    for name, (obj, dt, status) in results.items():
        if name == "cpu":
            continue
        diffs = []
        for k, ref in ref_obj.items():
            got = obj.get(k)
            if ref is None or got is None:
                continue
            denom = max(1.0, abs(ref))
            diffs.append(abs(got - ref) / denom)
        if diffs:
            speedup = results["cpu"][1] / dt
            print(f"  {name:12s}: max rel-err {max(diffs):.2e} | "
                  f"mean {np.mean(diffs):.2e} | speedup vs cpu {speedup:.2f}x")

# print a few per-member values for sanity
if workload == "solo" and "cpu" in results:
    print(f"\n[solo-max values per member]")
    print(f"  {'member':8s} " + " ".join(f"{n:>14s}" for n in results))
    for m in member_bios:
        row = " ".join(
            f"{(results[n][0].get(m) if results[n][0].get(m) is not None else float('nan')):14.4f}"
            for n in results
        )
        print(f"  {m:8s} {row}")
