# -*- coding: utf-8 -*-
"""JAX PDHG backend — portable across CUDA, Metal, ROCm, TPU, CPU.

JAX picks the device based on its install: `jax[cuda12]` for NVIDIA,
`jax-metal` for Apple Silicon, `jax[rocm]` for AMD, plain `jax` for CPU.
The same Python code runs unchanged; this module just transfers the
shared `S` and the (B, n) batched bounds / objective vectors to whichever
device JAX selects.

The PDHG inner loop is *not* `jit`-compiled: it is the plain Python loop
in `_pdhg.run_pdhg_batch`, so every iteration dispatches JAX ops eagerly
and the periodic convergence poll forces a host sync. Across LP instances
the per-instance arrays are laid out with the batch axis as the trailing
dimension so column-wise broadcasting maps cleanly onto SpMM, which is
what most accelerators are tuned for.

Precision: the shared kernel is written in float64, but JAX defaults to
32-bit and silently truncates any float64 request. This backend therefore
runs the solve inside `jax.enable_x64()`, a scoped (thread-local) context
rather than the process-global `jax.config.update("jax_enable_x64", True)`
-- so a host application that embeds this package keeps whatever default
precision it chose, and only our kernel sees 64-bit dtypes. If 64-bit is
unavailable (e.g. `jax-metal`, which has no float64 support), the solve
raises rather than quietly running in single precision; pass `x64=False`
to accept float32 and get a loud warning instead.
"""
from __future__ import annotations

import contextlib

import numpy as np
import scipy.sparse as sp

from mscommunity.batched_lp import BatchedLPSolver, register_batched_solver

from ._pdhg import ArrayBackend, PDHGConfig, run_pdhg_batch

try:
    import jax
    import jax.numpy as jnp
    from jax.experimental.sparse import BCOO
    _JAX_OK = True
    _JAX_ERR = None
except Exception as _exc:  # pragma: no cover — import guard
    jax = None
    jnp = None
    BCOO = None
    _JAX_OK = False
    _JAX_ERR = _exc


def _x64_scope():
    """Thread-local 64-bit scope, or a no-op if this JAX is too old.

    `jax.enable_x64()` sets the flag only for the calling thread and only
    for the duration of the `with` block, which is what makes it safe to
    use inside a library: unlike `jax.config.update("jax_enable_x64", True)`
    it does not change the dtype defaults of the host application. Arrays
    must be *created* inside the scope to be 64-bit, so the whole solve
    runs in it.
    """
    if hasattr(jax, "enable_x64"):
        return jax.enable_x64()
    return contextlib.nullcontext()


def _x64_active() -> bool:
    """True if float64 survives canonicalization in the current scope."""
    try:
        return jax.dtypes.canonicalize_dtype(np.float64) == np.float64
    except Exception:  # pragma: no cover — defensive
        return False


def _jax_backend(device=None, allow_float32: bool = False) -> ArrayBackend:
    if not _JAX_OK:
        raise ImportError(
            f"jax-pdhg backend requires `jax` (install jax[cuda12], "
            f"jax-metal, jax[rocm], or plain jax for CPU). Import failed: {_JAX_ERR}"
        )

    def _to_dev(a):
        arr = jnp.asarray(a)
        return jax.device_put(arr, device) if device is not None else arr

    def _to_host(a):
        return np.asarray(a)

    def _sparse(S):
        # BCOO is JAX's native sparse format; constructed from scipy COO.
        S_coo = S.tocoo()
        return BCOO(
            (jnp.asarray(S_coo.data), jnp.stack([
                jnp.asarray(S_coo.row), jnp.asarray(S_coo.col)
            ], axis=1)),
            shape=S.shape,
        )

    def _spmm(A, X):
        return A @ X

    def _spmm_T(A, Y):
        return A.T @ Y

    return ArrayBackend(
        xp=jnp, to_device=_to_dev, to_host=_to_host,
        sparse_from_scipy=_sparse, spmm=_spmm, spmm_T=_spmm_T,
        allow_float32=allow_float32, name="jax-pdhg",
    )


class JAXPDHGBatchedLPSolver(BatchedLPSolver):
    """First-order PDHG batched LP solver via JAX.

    Device is whatever JAX picks; pass `device="cpu"` / `"gpu"` /
    `"metal"` etc. or a `jax.Device` to pin explicitly. See JAX docs for
    backend-selection details.

    `x64=True` (the default) runs the solve inside `jax.enable_x64()` so
    the kernel really is double precision, matching the NumPy and CuPy
    backends. The scope is thread-local and exited before `solve()`
    returns, so it does not alter the enclosing application's JAX dtype
    defaults. Set `x64=False` only if you want single precision (or are on
    a device such as Metal that has no float64); results then agree with
    the NumPy backend to ~1e-6 instead of ~1e-15 and a `RuntimeWarning`
    says so on every solve.
    """

    name = "jax-pdhg"

    def __init__(self, max_iters: int = 5000, tol: float = 1e-4,
                 check_every: int = 50, power_iters: int = 30,
                 step_safety: float = 0.9, seed: int = 0, device=None,
                 x64: bool = True):
        if not _JAX_OK:
            raise ImportError(
                f"jax-pdhg backend requires `jax`. Import failed: {_JAX_ERR}"
            )
        self.config = PDHGConfig(
            max_iters=max_iters, tol=tol, check_every=check_every,
            power_iters=power_iters, step_safety=step_safety, seed=seed,
        )
        if isinstance(device, str):
            device = jax.devices(device)[0]
        self.device = device
        self.x64 = bool(x64)

    def solve(self, problem, instances):
        if not self.x64:
            # This does NOT force single precision, it only declines to request
            # double: on a host that already enabled x64 globally the kernel
            # still runs in float64. It suppresses the downgrade error, nothing
            # more.
            return run_pdhg_batch(
                _jax_backend(self.device, allow_float32=True),
                problem, instances, self.config,
            )
        with _x64_scope():
            if not _x64_active():
                raise RuntimeError(
                    "jax-pdhg could not enable 64-bit mode, so the kernel "
                    "would silently run in float32 while the NumPy and CuPy "
                    "backends run in float64. Upgrade JAX (`jax.enable_x64` "
                    "is required) or use a device with float64 support "
                    "(`jax-metal` has none). Pass x64=False to stop requesting "
                    "double precision and accept whatever JAX is configured for."
                )
            # Suppress nothing: any residual truncation is still caught by
            # the float64 guard inside run_pdhg_batch.
            return run_pdhg_batch(
                _jax_backend(self.device), problem, instances, self.config
            )


if _JAX_OK:
    register_batched_solver("jax-pdhg", JAXPDHGBatchedLPSolver)
