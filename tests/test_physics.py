"""Physics and numerics tests. Run with `pytest -q`."""
import math

import h5py
import numpy as np
import pytest

from ppdyn import kernels, loads
from ppdyn.simulation import Simulation, auto_cutoff


def make(text, **kw):
    return Simulation(loads("[run]\nprogress = false\n" + text), quiet=True, **kw)


BASE = """
[domain]
box = [{L}, {L}, {L}]
[[species]]
count = {n}
[species.velocity]
temperature = {T}
[interaction]
kappa = {kappa}
method = "{method}"
[time]
dt = {dt}
steps = {steps}
"""


def box_for(n):
    return (4 / 3 * math.pi * n) ** (1 / 3)


# ----------------------------------------------------------------------
def test_cells_match_direct():
    n = 3000
    L = box_for(n)
    a = make(BASE.format(L=L, n=n, T=0.01, kappa=1.5, method="cells", dt=0.01, steps=1))
    b = make(BASE.format(L=L, n=n, T=0.01, kappa=1.5, method="direct", dt=0.01, steps=1))
    b.pos[:] = a.pos
    b.rc = a.rc
    a.compute_forces()
    b.compute_forces()
    assert np.allclose(a.force, b.force, rtol=1e-10, atol=1e-12)
    assert np.isclose(a.pe.sum(), b.pe.sum(), rtol=1e-12)


def test_cells_match_direct_open_slab():
    """Non-periodic thin slab with one cell along z."""
    txt = """
[domain]
box = [40, 40, 3]
[boundary]
x = "periodic"
y = "periodic"
z = "reflecting"
[[species]]
count = 2500
[interaction]
kappa = 0.8
method = "{m}"
[time]
steps = 1
"""
    a = make(txt.format(m="cells"))
    b = make(txt.format(m="direct"))
    b.pos[:] = a.pos
    b.rc = a.rc
    a.compute_forces(); b.compute_forces()
    assert np.allclose(a.force, b.force, rtol=1e-10, atol=1e-12)


def test_pair_force_analytic():
    s = make("""
[domain]
box = [40, 40, 40]
[boundary]
x = "reflecting"
y = "reflecting"
z = "reflecting"
[[species]]
count = 2
[interaction]
kappa = 0.7
method = "direct"
[time]
steps = 1
""")
    s.pos[:] = [[0, 0, 0], [1.3, 0, 0]]
    s.compute_forces()
    r, k = 1.3, 0.7
    f = (1 + k * r) * math.exp(-k * r) / r ** 2
    assert np.isclose(s.force[1, 0], f) and np.isclose(s.force[0, 0], -f)
    assert np.isclose(s.pe.sum(), math.exp(-k * r) / r - math.exp(-k * s.rc) / s.rc)


def test_energy_conservation_second_order():
    """NVE drift shrinks ~dt^2 (velocity Verlet)."""
    n = 500
    L = box_for(n)
    drifts = []
    for dt in (0.04, 0.02):
        s = make(BASE.format(L=L, n=n, T=0.005, kappa=1.0, method="direct", dt=dt, steps=1))
        # equilibrate-ish start: put particles on a jittered lattice to avoid close pairs
        g = np.stack(np.meshgrid(*[np.arange(8)] * 3, indexing="ij"), -1).reshape(-1, 3)[:n]
        s.pos[:] = (g + 0.5) * (L / 8) - L / 2 + np.random.default_rng(1).normal(0, 0.05, (n, 3))
        s.compute_forces()
        E = []
        for _ in range(int(round(4.0 / dt))):
            s.advance()
            kx, ky, kz, _ = kernels.kinetic(s.vel, s.m, s.state)
            E.append(kx + ky + kz + s.pe.sum())
        E = np.asarray(E)
        drifts.append(np.ptp(E) / abs(E.mean()))
    print("relative energy fluctuation", drifts)
    assert drifts[0] < 1e-3
    assert 2.5 < drifts[0] / drifts[1] < 6.0


def test_langevin_reaches_target_temperature():
    s = make("""
[domain]
box = [30, 30, 30]
[[species]]
count = 4000
mass = 2.0
[interaction]
model = "none"
[neutrals]
drag = "constant"
gamma = 2.0
temperature = 0.3
[time]
dt = 0.05
steps = 1
""")
    s.compute_forces()
    T = []
    for i in range(400):
        s.advance()
        if i > 100:
            T.append(s.diagnostics()["temperature"])
    assert abs(np.mean(T) - 0.3) / 0.3 < 0.02


def test_free_fall_and_absorbing_floor():
    s = make("""
[domain]
box = [10, 10, 10]
[boundary]
x = "periodic"
y = "periodic"
z = ["absorbing", "reflecting"]
[[species]]
count = 1
charge.value = 0.5
[species.velocity]
distribution = "zero"
[species.placement]
region = "plane"
axis = "z"
position = 4.0
[interaction]
model = "none"
[forces.gravity]
g = [0, 0, -2.0]
[time]
dt = 0.001
t_end = 5
""")
    s.compute_forces()
    while s.state[0] == 0:
        s.advance()
    t_exact = math.sqrt(2 * 9.0 / 2.0)
    assert abs(s.absorb_t[0] - t_exact) < 2e-3
    assert s.q[0] == 0.0 and np.isclose(s.absorb_q[0], 0.5)
    assert s.state[0] == kernels.STUCK


def test_sphere_reflection_keeps_particles_inside():
    s = make("""
[domain]
shape = "sphere"
radius = 5.0
[[species]]
count = 500
[species.placement]
region = "sphere"
r_max = 4.9
[species.velocity]
temperature = 4.0
[interaction]
model = "none"
[time]
dt = 0.05
steps = 400
""")
    s.compute_forces()
    k0 = s.diagnostics()["kinetic"]
    for _ in range(400):
        s.advance()
    assert (np.linalg.norm(s.pos, axis=1) <= 5.0 + 1e-12).all()
    assert np.isclose(s.diagnostics()["kinetic"], k0)          # elastic walls


def test_object_reflect_and_absorb():
    txt = """
[domain]
box = [20, 20, 20]
[boundary]
x = "open"
y = "open"
z = "open"
[[objects]]
type = "{t}"
center = [0, 0, 0]
size = [4, 4, 4]
radius = 2.0
length = 4.0
action = "{a}"
[[species]]
count = 400
[species.placement]
region = "box"
center = [0, 0, 6]
size = [3, 3, 0.1]
[species.velocity]
distribution = "zero"
drift = [0, 0, -1.0]
[interaction]
model = "none"
[time]
dt = 0.01
t_end = 8
"""
    for t in ("box", "sphere", "cylinder"):
        s = make(txt.format(t=t, a="reflect"))
        s.compute_forces()
        inside_any = False
        for _ in range(800):
            s.advance()
            d = np.array([kernels.sdf_union(s.obj_kinds, s.obj_pars, *p)[0] for p in s.pos[s.state == 0]])
            inside_any |= (d < -1e-9).any()
        assert not inside_any
        s = make(txt.format(t=t, a="absorb"))
        s.compute_forces()
        for _ in range(800):
            s.advance()
        assert (s.state == kernels.STUCK).sum() > 300
        assert (s.hit_obj[s.state == kernels.STUCK] == 0).all()


def test_head_on_collision_exchanges_velocities():
    s = make("""
[domain]
box = [20, 20, 20]
[boundary]
x = "reflecting"
y = "reflecting"
z = "reflecting"
[[species]]
count = 2
[species.radius]
mean = 0.5
[interaction]
model = "none"
[collisions]
enabled = true
[time]
dt = 0.01
steps = 1
""")
    s.pos[:] = [[-2, 0, 0], [2, 0, 0]]
    s.vel[:] = [[1, 0, 0], [-1, 0, 0]]
    s.compute_forces()
    for _ in range(400):
        s.advance()
    assert np.allclose(s.vel, [[-1, 0, 0], [1, 0, 0]])
    assert np.allclose((s.m[:, None] * s.vel).sum(0), 0)


def test_trilinear_exact_for_linear_field():
    rng = np.random.default_rng(0)
    nx, ny, nz = 9, 7, 8
    o, d = np.array([-2.0, -1.5, -1.0]), np.array([0.5, 0.5, 0.25])
    X, Y, Z = np.meshgrid(*(o[k] + d[k] * np.arange(n) for k, n in enumerate((nx, ny, nz))), indexing="ij")
    E = np.stack([1 + 2 * X - Y, 3 * Z, X + Y + Z], -1)
    pts = rng.uniform(o, o + d * (np.array([nx, ny, nz]) - 1), (200, 3))
    for p in pts:
        e = kernels._trilinear(p[0], p[1], p[2], o, d, E, 0)
        exact = [1 + 2 * p[0] - p[1], 3 * p[2], p.sum()]
        assert np.allclose(e, exact)
    assert np.allclose(kernels._trilinear(10.0, 0.0, 0.0, o, d, E, 0), 0)


def test_auto_cutoff():
    for k in (0.5, 1.0, 2.0):
        rc = auto_cutoff(k, 1e-4)
        f = lambda r: (1 + k * r) * math.exp(-k * r) / r ** 2
        assert np.isclose(f(rc) / f(1.0), 1e-4, rtol=1e-6)
    assert auto_cutoff(0.0, 1e-4) == math.inf


def test_si_matches_reduced():
    """The same physical system in SI and reduced units gives identical dynamics."""
    a, Q, m = 3e-4, 1e-15, 1e-13
    from ppdyn.units import Scales, KB
    sc = Scales(a, Q, m, True)
    T_K = 0.02 * sc.temperature
    si = make(f"""
[units]
system = "SI"
length = {a}
charge = {Q}
mass = {m}
[domain]
box = [{10 * a}, {10 * a}, {10 * a}]
[[species]]
count = 200
mass = {m}
charge.value = {Q}
charge.exponent = 0.0
[species.velocity]
temperature = {T_K}
[interaction]
debye_length = {a / 1.3}
method = "direct"
[time]
dt = {0.02 * sc.time}
steps = 50
""")
    red = make("""
[domain]
box = [10, 10, 10]
[[species]]
count = 200
[species.velocity]
temperature = 0.02
[interaction]
kappa = 1.3
method = "direct"
[time]
dt = 0.02
steps = 50
""")
    assert np.allclose(si.pos, red.pos) and np.allclose(si.vel, red.vel)
    assert np.isclose(si.kappa, 1.3) and np.isclose(si.dt, 0.02)
    for s in (si, red):
        s.compute_forces()
        for _ in range(50):
            s.advance()
    assert np.allclose(si.pos, red.pos, atol=1e-10)


def test_epstein_rate():
    from ppdyn.particles import epstein_gamma
    # 3.5 um melamine grain, 30 Pa argon: about 50 1/s
    r = 3.5e-6
    m = 1510 * 4 / 3 * math.pi * r ** 3
    g = epstein_gamma(r, m, 30.0, 300.0, 39.95, 1.44)
    assert 40 < g < 60
