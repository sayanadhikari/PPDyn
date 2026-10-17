"""External field grids: reading PINC / ppdyn files and solving for the
field around conducting objects (ported and extended from the debris branch).
"""
from __future__ import annotations

import re
import time

import h5py
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numba import njit, prange

from . import geometry, kernels


class FieldGrid:
    """Cell-node grid with E stored as (nx, ny, nz, 3)."""

    def __init__(self, origin, spacing, E, units="SI", phi=None, mask=None):
        self.origin = np.asarray(origin, float)
        self.spacing = np.asarray(spacing, float)
        self.E = np.ascontiguousarray(E, dtype=np.float64)
        self.units = units
        self.phi = phi
        self.mask = mask

    @property
    def shape(self):
        return self.E.shape[:3]

    def axes(self):
        return [self.origin[d] + self.spacing[d] * np.arange(self.shape[d]) for d in range(3)]

    def save(self, path, attrs=None):
        with h5py.File(path, "w") as f:
            f.attrs["format"] = "ppdyn-field"
            f.attrs["units"] = self.units
            f.attrs["origin"] = self.origin
            f.attrs["spacing"] = self.spacing
            for k, v in (attrs or {}).items():
                f.attrs[k] = v
            f.create_dataset("E", data=self.E, compression="gzip", compression_opts=4)
            for name, arr in (("x", 0), ("y", 1), ("z", 2)):
                f.create_dataset(name, data=self.axes()[arr])
            if self.phi is not None:
                f.create_dataset("phi", data=self.phi, compression="gzip", compression_opts=4)
            if self.mask is not None:
                f.create_dataset("mask", data=self.mask.astype(np.uint8), compression="gzip", compression_opts=4)


def _numeric_key(k):
    m = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", k)
    return float(m[-1]) if m else -1.0


def load(path, fmt="auto", step=-1, center=None, center_grid=True) -> FieldGrid:
    with h5py.File(path, "r") as f:
        if fmt == "auto":
            fmt = "ppdyn" if f.attrs.get("format", "") == "ppdyn-field" else "pinc"
        if fmt == "ppdyn":
            return FieldGrid(f.attrs["origin"], f.attrs["spacing"], f["E"][:], str(f.attrs.get("units", "SI")),
                             f["phi"][:] if "phi" in f else None)
        # PINC: one dataset per stored time step, each (nx, ny, nz, 3)
        keys = sorted((k for k in f.keys() if isinstance(f[k], h5py.Dataset)), key=_numeric_key)
        if not keys:
            raise ValueError(f"{path}: no field datasets found")
        ds = f[keys[step]]
        E = ds[:]
        if E.ndim != 4 or E.shape[-1] != 3:
            raise ValueError(f"{path}:{keys[step]} has shape {E.shape}, expected (nx, ny, nz, 3)")
        dx = float(np.ravel(f.attrs.get("Axis denormalization factor", [1.0]))[0])
        qf = f.attrs.get("Quantity denormalization factor", None)
        if qf is not None:
            E = E * float(np.ravel(qf)[0])
    n = np.array(E.shape[:3])
    spacing = np.full(3, dx)
    c = np.zeros(3) if center is None else np.asarray(center, float)
    origin = c - 0.5 * (n - 1) * spacing if center_grid else np.zeros(3)
    return FieldGrid(origin, spacing, E, "SI")


# ----------------------------------------------------------------------
# field solver
# ----------------------------------------------------------------------
def solve(cfg: dict, prims: list[dict], log=print) -> FieldGrid:
    """Solve (lap - 1/lambda^2) phi = 0 with phi = V_k on object k and phi = 0
    on the domain boundary. Everything in the input unit system."""
    fs = cfg["field_solver"]
    n = np.array([int(v) for v in fs["grid"]])
    if n.size != 3 or (n < 3).any():
        raise ValueError("field_solver.grid needs three integers >= 3")
    box = np.asarray(cfg["domain"]["box"], float)
    if cfg["domain"]["shape"] == "sphere":
        box = np.full(3, 2 * cfg["domain"]["radius"])
    c = np.asarray(cfg["domain"]["center"], float)
    h = box / (n - 1)
    origin = c - box / 2
    axes = [origin[d] + h[d] * np.arange(n[d]) for d in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)

    kinds, pars, _, _ = geometry.kernel_arrays(prims)
    phi_fixed = np.full(pts.shape[0], np.nan)
    if len(prims):
        d, who = _sdf_grid_kernel(pts, kinds, pars)
        inside = d <= 0.5 * h.min()
        pot = np.array([p["potential"] for p in prims])
        phi_fixed[inside] = pot[who[inside]]
    else:
        inside = np.zeros(pts.shape[0], bool)

    idx = np.arange(pts.shape[0]).reshape(tuple(n))
    boundary = np.zeros(tuple(n), bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    boundary = boundary.ravel()
    phi_fixed[boundary & ~inside] = 0.0
    known = ~np.isnan(phi_fixed)
    free = np.flatnonzero(~known)
    log(f"field solver: grid {tuple(n)}, {inside.sum()} object nodes, {free.size} unknowns")

    lam = fs["debye_length"]
    inv_h2 = 1.0 / h ** 2
    diag = 2 * inv_h2.sum() + (1.0 / lam ** 2 if lam > 0 else 0.0)
    pos_of = -np.ones(pts.shape[0], np.int64)
    pos_of[free] = np.arange(free.size)
    rows, cols, vals = [np.arange(free.size)], [np.arange(free.size)], [np.full(free.size, diag)]
    rhs = np.zeros(free.size)
    ijk = np.array(np.unravel_index(free, tuple(n))).T
    for d in range(3):
        for sgn in (-1, 1):
            nb = ijk.copy()
            nb[:, d] += sgn
            nb_lin = np.ravel_multi_index(nb.T, tuple(n))   # free nodes are never on the boundary
            nf = pos_of[nb_lin]
            isfree = nf >= 0
            rows.append(np.arange(free.size)[isfree])
            cols.append(nf[isfree])
            vals.append(np.full(isfree.sum(), -inv_h2[d]))
            rhs[~isfree] += inv_h2[d] * phi_fixed[nb_lin[~isfree]]
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(free.size, free.size))
    t0 = time.time()
    x = _spd_solve(A, rhs, fs["tolerance"], log)
    log(f"field solver: solved in {time.time() - t0:.2f} s")
    phi = phi_fixed.copy()
    phi[free] = x
    phi = phi.reshape(tuple(n))
    gx, gy, gz = np.gradient(phi, *h)
    E = -np.stack([gx, gy, gz], -1)
    E[inside.reshape(tuple(n))] = 0.0
    return FieldGrid(origin, h, E, cfg["units"]["system"], phi, inside.reshape(tuple(n)))


@njit(parallel=True, cache=True)
def _sdf_grid_kernel(pts, kinds, pars):
    n = pts.shape[0]
    d = np.empty(n)
    who = np.empty(n, np.int64)
    for i in prange(n):
        d[i], who[i] = kernels.sdf_union(kinds, pars, pts[i, 0], pts[i, 1], pts[i, 2])
    return d, who


def _spd_solve(A, b, tol, log):
    try:
        import pyamg
        ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        res = []
        x = ml.solve(b, tol=tol, accel="cg", residuals=res, maxiter=500)
        log(f"field solver: AMG-CG, {len(res) - 1} iterations")
        return x
    except ImportError:
        pass
    M = sp.diags(1.0 / A.diagonal())
    it = [0]

    def cb(_):
        it[0] += 1
    x, info = spla.cg(A, b, rtol=tol, maxiter=20000, M=M, callback=cb)
    if info != 0:
        log(f"field solver: CG did not fully converge (info={info})")
    log(f"field solver: Jacobi-CG, {it[0]} iterations (install pyamg for speed)")
    return x
