"""Species sampling: sizes, masses, charges, initial positions and velocities.

Everything here works in the *input* unit system; the simulation converts
the resulting arrays to reduced units afterwards.
"""
from __future__ import annotations

import math
import os

import numpy as np

from .units import EPS0, KB


def _redraw_positive(rng, draw, n, minimum):
    x = draw(n)
    bad = x <= minimum
    tries = 0
    while bad.any():
        x[bad] = draw(int(bad.sum()))
        bad = x <= minimum
        tries += 1
        if tries > 1000:
            raise ValueError("size distribution keeps producing values below radius.min; check mean/std")
    return x


def sample_size(spec: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    """Sample the size variable (radius, or relative mass for variable = mass)."""
    d = spec["distribution"]
    mn = spec["min"]
    if d == "fixed":
        return np.full(n, spec["mean"], dtype=float)
    if d == "normal":
        return _redraw_positive(rng, lambda k: rng.normal(spec["mean"], spec["std"], k), n, mn)
    if d == "lognormal":
        mu = math.log(spec["mean"]) if spec["mean"] > 0 else 0.0
        return rng.lognormal(mu, spec["std"], n)
    if d == "uniform":
        return rng.uniform(spec["low"], spec["high"], n)
    vals = spec["values"]
    w = np.asarray(spec["weights"], float)
    if len(vals) < 2:
        raise ValueError(f"radius.distribution = {d!r} needs radius.values with two entries")
    w = w[: len(vals)] / w[: len(vals)].sum()
    pick = rng.choice(len(vals), size=n, p=w)
    if d == "binary":
        return np.asarray(vals, float)[pick]
    stds = spec["stds"] or [0.0] * len(vals)
    out = np.empty(n)
    for k in range(len(vals)):
        sel = pick == k
        out[sel] = _redraw_positive(rng, lambda c, k=k: rng.normal(vals[k], stds[k], c), int(sel.sum()), mn)
    return out


def species_properties(sp: dict, rng, si: bool):
    """Return radius, mass, charge arrays for one species (input units)."""
    n = sp["count"]
    rs = sp["radius"]
    x = sample_size(rs, n, rng)
    if rs["variable"] == "mass":
        # legacy: x is mass in units of species.mass; radius ~ x^(1/3)
        rel = np.cbrt(np.maximum(x, 0.0))
        rref = sp["reference_radius"]
        radius = rel * rref
        mass = sp["mass"] * x
        if sp["density"] > 0 and rref > 0:
            mass = sp["density"] * 4.0 / 3.0 * math.pi * radius ** 3
    else:
        radius = x
        rref = sp["reference_radius"] if sp["reference_radius"] > 0 else (radius.mean() if n else 0.0)
        rel = radius / rref if rref > 0 else np.ones(n)
        if sp["density"] > 0:
            mass = sp["density"] * 4.0 / 3.0 * math.pi * radius ** 3
        else:
            mass = sp["mass"] * rel ** 3
    ch = sp["charge"]
    if ch["model"] == "potential":
        charge = 4.0 * math.pi * EPS0 * radius * ch["potential"]
    else:
        charge = ch["value"] * rel ** ch["exponent"] if ch["exponent"] != 0 else np.full(n, ch["value"])
    if np.any(mass <= 0):
        raise ValueError(f"species '{sp['name']}': non-positive masses; set mass, density or a positive radius")
    return radius, mass, charge, rel


def _lattice(kind, n, size, center, rng, jitter):
    basis = {"sc": [(0, 0, 0)],
             "bcc": [(0, 0, 0), (0.5, 0.5, 0.5)],
             "fcc": [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]}[kind]
    vol = np.prod(size)
    cell = (vol * len(basis) / max(n, 1)) ** (1 / 3)
    m = np.maximum(np.ceil(np.asarray(size) / cell).astype(int), 1)
    idx = np.stack(np.meshgrid(*[np.arange(k) for k in m], indexing="ij"), -1).reshape(-1, 3)
    pts = (idx[:, None, :] + np.asarray(basis)[None]).reshape(-1, 3) * (np.asarray(size) / m)
    pts = pts[:n] - np.asarray(size) / 2 + center
    if len(pts) < n:
        raise ValueError("lattice too small for the requested count")
    if jitter > 0:
        pts += rng.normal(0, jitter, pts.shape)
    return pts


def _read_positions(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npy":
        return np.load(path)
    if ext in (".h5", ".hdf5"):
        import h5py
        with h5py.File(path, "r") as f:
            if "snapshots/position" in f:
                return f["snapshots/position"][-1][:]
            return f["position"][:]
    return np.loadtxt(path, delimiter="," if ext == ".csv" else None)


def place(sp: dict, n: int, domain: dict, rng, inside_fn=None, base_dir="."):
    pl = sp["placement"]
    center = np.asarray(pl["center"] if pl["center"] is not None else domain["center"], float)
    size = np.asarray(pl["size"] if pl["size"] is not None else domain["box"], float)
    region = pl["region"]
    if region == "file":
        p = _read_positions(os.path.join(base_dir, pl["file"]))
        if len(p) < n:
            raise ValueError(f"{pl['file']} holds {len(p)} positions, species needs {n}")
        return np.asarray(p[:n], float)
    if region == "lattice":
        return _lattice(pl["lattice"], n, size, center, rng, pl["jitter"])

    def draw(k):
        if region == "box":
            return center + (rng.random((k, 3)) - 0.5) * size
        if region in ("sphere", "shell"):
            rmax = pl["r_max"] if pl["r_max"] > 0 else (domain["radius"] or size.min() / 2)
            rmin = pl["r_min"] if region == "shell" else 0.0
            u = rng.normal(size=(k, 3))
            u /= np.linalg.norm(u, axis=1)[:, None]
            r = (rmin ** 3 + rng.random(k) * (rmax ** 3 - rmin ** 3)) ** (1 / 3)
            return center + u * r[:, None]
        if region == "plane":
            p = center + (rng.random((k, 3)) - 0.5) * size
            ax = "xyz".index(pl["axis"])
            p[:, ax] = pl["position"]
            return p
        raise ValueError(region)

    pts = draw(n)
    if inside_fn is not None and pl["avoid_objects"] and n:
        for _ in range(200):
            bad = inside_fn(pts)
            if not bad.any():
                break
            pts[bad] = draw(int(bad.sum()))
    return pts


def velocities(sp: dict, mass: np.ndarray, rng, kT_unit: float):
    """kT_unit converts velocity.temperature into the energy unit used with mass."""
    n = len(mass)
    vs = sp["velocity"]
    if vs["distribution"] == "zero" or n == 0:
        v = np.zeros((n, 3))
    elif vs["distribution"] == "uniform":
        v = (rng.random((n, 3)) - 0.5) * np.asarray(vs["vmax"])
    else:
        kT = vs["temperature"] * kT_unit
        v = rng.normal(size=(n, 3)) * np.sqrt(np.maximum(kT, 0.0) / mass)[:, None]
    if vs["remove_net_momentum"] and n > 1 and vs["distribution"] != "zero":
        v -= (mass[:, None] * v).sum(0) / mass.sum()
    v += np.asarray(vs["drift"])
    return v


def epstein_gamma(radius_m, mass_kg, pressure, Tn_K, mass_amu, delta):
    """Epstein neutral drag rate [1/s] for grains of radius r in a gas at pressure p."""
    from .units import AMU
    mn = mass_amu * AMU
    vth = math.sqrt(8.0 * KB * Tn_K / (math.pi * mn))
    n_n = pressure / (KB * Tn_K)
    return delta * 4.0 / 3.0 * math.pi * radius_m ** 2 * n_n * mn * vth / mass_kg
