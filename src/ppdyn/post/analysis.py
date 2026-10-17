"""Analysis of PPDyn output files."""
from __future__ import annotations

import json
import math

import h5py
import numpy as np
from numba import njit, prange

from ..units import KB, Scales


class Run:
    """Lazy reader for a PPDyn HDF5 file.

    >>> r = Run("output/ocp.h5")
    >>> r.times, r.position(-1), r.history["temperature"]
    """

    def __init__(self, path):
        self.path = path
        self.f = h5py.File(path, "r")
        a = self.f.attrs
        self.n = int(a["N"])
        self.lo = np.asarray(a["box_lo"])
        self.hi = np.asarray(a["box_hi"])
        self.L = self.hi - self.lo
        self.periodic = np.asarray(a["periodic"]).astype(bool)
        self.kappa = float(a["kappa"])
        self.dt = float(a["dt"])
        u = json.loads(a["units"])
        self.scales = Scales(u["length"], u["charge"], u["mass"], u["si"])
        self.input_text = a.get("input", "")
        self.species_names = json.loads(a.get("species_names", '["dust"]'))
        self.species = self.f["particles/species"][:]
        self.radius = self.f["particles/radius"][:]
        self.mass = self.f["particles/mass"][:]
        self.charge0 = self.f["particles/charge0"][:]
        s = self.f["snapshots"]
        self.times = s["time"][:]
        self.steps = s["step"][:]
        self.quantities = [k for k in s.keys() if k not in ("time", "step")]
        h = self.f["history"]
        cols = list(h.attrs["columns"])
        data = h["data"][:]
        self.history = {c: data[:, i] for i, c in enumerate(cols)}
        self.final = {k: v[:] for k, v in self.f["final"].items()} if "final" in self.f else {}

    def __len__(self):
        return len(self.times)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self.f.close()

    def _get(self, name, k):
        if name not in self.f["snapshots"]:
            raise KeyError(f"'{name}' was not written; add it to output.quantities")
        return self.f["snapshots"][name][k].astype(np.float64) if name != "state" else self.f["snapshots"][name][k]

    def position(self, k=-1):
        return self._get("position", k)

    def velocity(self, k=-1):
        return self._get("velocity", k)

    def charge(self, k=-1):
        return self._get("charge", k) if "charge" in self.quantities else self.charge0

    def state(self, k=-1):
        if "state" in self.quantities:
            return self._get("state", k)
        return np.zeros(self.n, np.uint8)

    def active(self, k=-1):
        return self.state(k) == 0

    # ---------------------------------------------------------------- unit helpers
    def si(self, values, dim):
        """Convert reduced values to SI (identity for reduced-unit runs)."""
        return values * self.scales.factor(dim)

    def temperature_kelvin(self, T_reduced):
        return T_reduced * self.scales.temperature


# ======================================================================
@njit(parallel=True, cache=True)
def _rdf_hist(pos, L, periodic, rmax, nbins):
    """Pair-distance histogram; each chunk owns a private row (race-free)."""
    n = pos.shape[0]
    nch = 64
    rows = np.zeros((nch, nbins), np.int64)
    dr = rmax / nbins
    r2max = rmax * rmax
    for c in prange(nch):
        for i in range(c, n, nch):
            for j in range(i + 1, n):
                d2 = 0.0
                for d in range(3):
                    x = pos[i, d] - pos[j, d]
                    if periodic[d]:
                        x -= L[d] * math.floor(x / L[d] + 0.5)
                    d2 += x * x
                if d2 < r2max:
                    b = int(math.sqrt(d2) / dr)
                    if b < nbins:
                        rows[c, b] += 1
    return rows.sum(axis=0)


def rdf(pos, L, periodic, rmax=None, nbins=200, select=None, max_particles=20000, rng=None):
    """Pair correlation g(r). Non-periodic axes are not edge corrected, so use
    rmax well below the cloud size for confined systems."""
    pos = np.asarray(pos, float)
    sel = np.ones(len(pos), bool) if select is None else np.asarray(select, bool)
    idx = np.flatnonzero(sel)
    if len(idx) > max_particles:
        rng = rng or np.random.default_rng(0)
        idx = np.sort(rng.choice(idx, max_particles, replace=False))
    p = np.ascontiguousarray(pos[idx])
    L = np.asarray(L, float)
    per = np.asarray(periodic, np.bool_)
    if rmax is None:
        rmax = 0.5 * L[per].min() if per.any() else 0.25 * float(np.ptp(p, axis=0).min())
    hist = _rdf_hist(p, L, per, float(rmax), int(nbins))
    edges = np.linspace(0, rmax, nbins + 1)
    r = 0.5 * (edges[1:] + edges[:-1])
    n = len(p)
    ext = np.where(per, L, np.ptp(p, axis=0))
    flat = (~per) & (ext < 0.05 * rmax)
    if flat.sum() == 1:
        # quasi-2D layer (monolayer crystal, deposit on a wall): areal normalisation
        shell = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
        area = float(np.prod(ext[~flat]))
        g = 2.0 * hist / (n * (n / area) * shell)
    else:
        shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        vol = float(np.prod(np.maximum(ext, 1e-300)))
        g = 2.0 * hist / (n * (n / vol) * shell)
    return r, g


def unwrap(positions, L, periodic):
    """Unwrap a (nt, N, 3) trajectory across periodic boundaries."""
    x = np.array(positions, float, copy=True)
    for d in range(3):
        if periodic[d]:
            dx = np.diff(x[:, :, d], axis=0)
            jumps = -L[d] * np.round(dx / L[d])
            x[1:, :, d] += np.cumsum(jumps, axis=0)
    return x


def msd(run: Run, select=None, max_frames=400):
    k = np.unique(np.linspace(0, len(run) - 1, min(len(run), max_frames)).astype(int))
    traj = np.stack([run.position(i) for i in k])
    if select is not None:
        traj = traj[:, select]
    x = unwrap(traj, run.L, run.periodic)
    d = x - x[0]
    return run.times[k] - run.times[k[0]], (d ** 2).sum(-1).mean(-1)


def diffusion_coefficient(t, m, fit_from=0.5):
    """Long-time diffusion from MSD = 6 D t (fit on the last part)."""
    sel = t >= fit_from * t[-1]
    if sel.sum() < 2:
        return float("nan")
    return float(np.polyfit(t[sel], m[sel], 1)[0] / 6.0)


def velocity_distribution(v, mass=None, bins=80):
    """Histogram of one velocity component plus the matching Maxwellian."""
    comp = v.ravel()
    T = float(np.mean((mass[:, None] if mass is not None else 1.0) * v ** 2))
    m_mean = float(np.mean(mass)) if mass is not None else 1.0
    hist, edges = np.histogram(comp, bins=bins, density=True)
    c = 0.5 * (edges[1:] + edges[:-1])
    s2 = T / m_mean
    maxw = np.exp(-c ** 2 / (2 * s2)) / np.sqrt(2 * np.pi * s2) if s2 > 0 else np.zeros_like(c)
    return c, hist, maxw, T


def density_profile(pos, axis=2, bins=60, lo=None, hi=None, select=None):
    p = pos if select is None else pos[select]
    x = p[:, axis]
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    h, e = np.histogram(x, bins=bins, range=(lo, hi))
    return 0.5 * (e[1:] + e[:-1]), h


def radial_profile(pos, center=(0, 0, 0), bins=60, rmax=None, select=None, cylindrical=False):
    p = (pos if select is None else pos[select]) - np.asarray(center)
    r = np.linalg.norm(p[:, :2] if cylindrical else p, axis=1)
    rmax = r.max() if rmax is None else rmax
    h, e = np.histogram(r, bins=bins, range=(0, rmax))
    if cylindrical:
        vol = np.pi * (e[1:] ** 2 - e[:-1] ** 2)
    else:
        vol = 4 / 3 * np.pi * (e[1:] ** 3 - e[:-1] ** 3)
    return 0.5 * (e[1:] + e[:-1]), h / vol


def coupling(run: Run, k=-1):
    """Effective coupling Gamma = <Q>^2 / (a T) in reduced units (a = 1)."""
    T = run.history["temperature"]
    t = run.history["time"]
    Tk = float(np.interp(run.times[k], t, T))
    q = run.charge(k)
    act = run.active(k)
    qm = float(np.mean(np.abs(q[act]))) if act.any() else 0.0
    return qm ** 2 / Tk if Tk > 0 else math.inf


def fall_times(run: Run):
    """Times at which absorbed particles reached a wall (dust-fall)."""
    t = run.final.get("absorb_time")
    if t is None:
        return np.array([])
    return t[t >= 0]
