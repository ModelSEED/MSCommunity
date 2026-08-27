# -*- coding: utf-8 -*-
"""Regression tests for the JAX PDHG backend's precision contract.

`_pdhg.run_pdhg_batch` asks for float64 everywhere and the shared
docstring claims the math is identical across NumPy / CuPy / JAX. JAX
defaults to 32-bit and silently truncates float64 requests, so the JAX
backend used to run the whole kernel in single precision. These tests
pin the fixed behavior: double precision by default, a loud failure or
warning when it is unavailable, and no process-global config mutation.
"""
from __future__ import annotations

import inspect
import os
import warnings

import numpy as np
import pytest
from cobra import Metabolite, Model, Reaction

from mscommunity.backends._pdhg import ArrayBackend, _check_double_precision
from mscommunity.batched_lp import LPInstance, solve_batch, _SOLVERS
import mscommunity.backends  # noqa: F401 — triggers registration

jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")


def _jax_can_run() -> bool:
    """True once JAX can actually execute an op, falling back to CPU.

    A CUDA plugin whose driver does not match makes *every* JAX op raise,
    including `jax.devices("cpu")`. The precision contract under test is
    device-independent, so retry on CPU before giving up. JAX_PLATFORMS is
    only read by JAX, and only when it has no working backend yet.
    """
    for _ in range(2):
        try:
            float(jnp.zeros(1).sum())
            return True
        except Exception:
            os.environ["JAX_PLATFORMS"] = "cpu"
    return False


pytestmark = pytest.mark.skipif(
    not _jax_can_run(), reason="no working JAX backend on this machine"
)


def _chain_model(n: int = 60) -> Model:
    """Linear metabolic chain: EX_M0 -> M0 -> ... -> M(n-1) -> SINK."""
    m = Model("chain")
    mets = [Metabolite(f"M{i}", compartment="c") for i in range(n)]
    m.add_metabolites(mets)
    rxns = []
    ex = Reaction("EX_M0", lower_bound=-10, upper_bound=0)
    ex.add_metabolites({mets[0]: -1})
    rxns.append(ex)
    for i in range(n - 1):
        r = Reaction(f"R{i}", lower_bound=0, upper_bound=1000)
        r.add_metabolites({mets[i]: -1, mets[i + 1]: 1})
        rxns.append(r)
    sk = Reaction("SINK", lower_bound=0, upper_bound=1000)
    sk.add_metabolites({mets[-1]: -1})
    rxns.append(sk)
    m.add_reactions(rxns)
    m.objective = "SINK"
    return m


# The contract under test is precision, not convergence, and the float32 /
# float64 gap is already several orders of magnitude wide after a few hundred
# iterations. 4000 iterations cost ~30 s across this file — a >10x increase on
# the whole suite — for no additional discrimination.
_KW = dict(max_iters=300, tol=1e-12, check_every=1000)


def _instances():
    return [
        LPInstance(id=f"cap_{c}", bounds={"EX_M0": (-c, 0.0)})
        for c in (10.0, 5.0, 2.0, 1.0)
    ]


@pytest.mark.skipif("jax-pdhg" not in _SOLVERS, reason="jax backend not registered")
def test_jax_pdhg_matches_numpy_to_double_precision():
    m = _chain_model()
    insts = _instances()
    ref = solve_batch(m, insts, backend="numpy-pdhg", **_KW)
    jx = solve_batch(m, insts, backend="jax-pdhg", **_KW)
    for a, b in zip(ref, jx):
        # float32 agreement was ~1e-6; float64 agreement is ~1e-15.
        assert abs(a.objective_value - b.objective_value) < 1e-10
        for rxn_id, flux in a.fluxes.items():
            assert abs(flux - b.fluxes[rxn_id]) < 1e-10


@pytest.mark.skipif("jax-pdhg" not in _SOLVERS, reason="jax backend not registered")
def test_jax_pdhg_emits_no_float64_truncation_warning():
    m = _chain_model(12)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solve_batch(m, _instances(), backend="jax-pdhg", **_KW)
    truncated = [w for w in caught if "truncated to dtype float32" in str(w.message)]
    assert not truncated, [str(w.message) for w in truncated]


@pytest.mark.skipif("jax-pdhg" not in _SOLVERS, reason="jax backend not registered")
def test_jax_x64_scope_does_not_leak_into_host_process():
    # The scope must be entered/exited around the solve only: a host
    # application's JAX dtype defaults have to survive untouched.
    before = jax.dtypes.canonicalize_dtype(np.float64)
    solve_batch(_chain_model(12), _instances(), backend="jax-pdhg", **_KW)
    assert jax.dtypes.canonicalize_dtype(np.float64) == before


@pytest.mark.skipif("jax-pdhg" not in _SOLVERS, reason="jax backend not registered")
def test_jax_x64_opt_out_warns_loudly():
    m = _chain_model(12)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sols = solve_batch(m, _instances(), backend="jax-pdhg", x64=False, **_KW)
    assert len(sols) == 4
    if jax.dtypes.canonicalize_dtype(np.float64) != np.float64:
        # Only meaningful when the ambient JAX default really is 32-bit:
        # opting out must still say, unmistakably, that it is single.
        msgs = [str(w.message) for w in caught if w.category is RuntimeWarning]
        assert any("float32" in msg for msg in msgs), msgs


def test_check_double_precision_raises_on_silent_downgrade():
    backend = ArrayBackend(
        xp=np, to_device=lambda a: a, to_host=np.asarray,
        sparse_from_scipy=lambda s: s,
        spmm=lambda a, x: a @ x, spmm_T=lambda a, y: a.T @ y,
        name="fake",
    )
    with pytest.raises(RuntimeError, match="float64"):
        _check_double_precision(backend, np.zeros(3, dtype=np.float32))
    # float64 passes silently
    _check_double_precision(backend, np.zeros(3, dtype=np.float64))


def test_check_double_precision_warns_when_backend_opted_in():
    backend = ArrayBackend(
        xp=np, to_device=lambda a: a, to_host=np.asarray,
        sparse_from_scipy=lambda s: s,
        spmm=lambda a, x: a @ x, spmm_T=lambda a, y: a.T @ y,
        allow_float32=True, name="fake",
    )
    with pytest.warns(RuntimeWarning, match="float32"):
        _check_double_precision(backend, np.zeros(3, dtype=np.float32))


def test_jax_backend_does_not_jit_compile_the_inner_loop():
    """The docstring used to claim the loop was `jit`-compiled; it is not.

    Assert the behavior rather than the prose, so rewording the docstring does
    not break the test and actually jitting the loop later does.
    """
    from mscommunity.backends import jax_pdhg

    source = inspect.getsource(jax_pdhg)
    assert "jax.jit" not in source and "@jit" not in source
