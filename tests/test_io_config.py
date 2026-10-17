"""Input handling, I/O, restart, field solver, legacy conversion and CLI."""
import math
import os
import subprocess
import sys

import h5py
import numpy as np
import pytest

import ppdyn
from ppdyn import ConfigError, fields, geometry, loads
from ppdyn.cli import main as cli
from ppdyn.legacy import convert

HERE = os.path.dirname(__file__)
EXAMPLES = os.path.join(os.path.dirname(ppdyn.__file__), "examples")


def test_unknown_key_suggestion():
    with pytest.raises(ConfigError, match="Did you mean 'kappa'"):
        loads("[[species]]\ncount = 10\n[interaction]\nkapa = 1.0\n")


def test_unknown_section_and_choice():
    with pytest.raises(ConfigError, match="Unknown section"):
        loads("[[species]]\ncount=1\n[interation]\n")
    with pytest.raises(ConfigError, match="not one of"):
        loads("[[species]]\ncount=1\n[boundary]\nx = \"periodc\"\n")


def test_booleans_are_real_booleans():
    """Legacy bug: bool('False') was True. TOML booleans cannot do that."""
    cfg = loads("[[species]]\ncount=1\n[output]\nxdmf = false\n")
    assert cfg["output"]["xdmf"] is False
    with pytest.raises(ConfigError):
        loads("[[species]]\ncount=1\n[output]\nxdmf = \"False\"\n")


def test_species_required():
    with pytest.raises(ConfigError, match="species"):
        loads("[time]\ndt = 0.1\n")


def _small(tmp_path, extra="", steps=40):
    txt = f"""
[run]
name = "t"
progress = false
[domain]
box = [12, 12, 12]
[[species]]
count = 300
[species.velocity]
temperature = 0.05
[neutrals]
drag = "constant"
gamma = 0.3
temperature = 0.05
[time]
dt = 0.02
steps = {steps}
[output]
dir = '{tmp_path / 'out'}'
every = 10
history_every = 5
checkpoint_every = 20
quantities = ["position", "velocity", "charge", "state"]
{extra}
"""
    return txt


def test_output_file_layout(tmp_path):
    sim, path = ppdyn.run(loads(_small(tmp_path)), progress=False)
    with h5py.File(path) as f:
        assert f["snapshots/position"].shape == (5, 300, 3)
        assert f["history/data"].shape[0] == 9
        assert "input" in f.attrs
    assert os.path.exists(path.replace(".h5", ".xdmf"))
    hist = np.genfromtxt(os.path.join(tmp_path, "out", "history.csv"), delimiter=",", names=True)
    assert len(hist) == 9


def test_restart_is_bitwise(tmp_path):
    """40 steps straight == 20 steps + restart from checkpoint + 20 steps (incl. Langevin noise)."""
    full, _ = ppdyn.run(loads(_small(tmp_path / "a", steps=40)), progress=False)
    half, _ = ppdyn.run(loads(_small(tmp_path / "b", steps=20)), progress=False)
    ck = tmp_path / "b" / "out" / "checkpoint.h5"
    cfg = loads(_small(tmp_path / "b", steps=40))
    cfg["run"]["restart"] = str(ck)
    rest, _ = ppdyn.run(cfg, progress=False)
    assert rest.step == 40
    assert np.array_equal(full.pos, rest.pos) and np.array_equal(full.vel, rest.vel)


def test_post_processing(tmp_path):
    from ppdyn.post import analysis, plots
    _, path = ppdyn.run(loads(_small(tmp_path)), progress=False)
    with analysis.Run(path) as r:
        rr, g = analysis.rdf(r.position(-1), r.L, r.periodic)
        assert np.isclose(g[len(g) // 2:].mean(), 1.0, atol=0.15)
        idx = plots.report(r, str(tmp_path / "rep"))
    assert os.path.exists(idx)
    for f in ("history.png", "snapshot.png", "rdf.png", "vdf.png", "viewer.html"):
        assert os.path.exists(tmp_path / "rep" / f), f


def test_rdf_ideal_gas():
    from ppdyn.post.analysis import rdf
    rng = np.random.default_rng(0)
    L = np.array([10.0, 10.0, 10.0])
    pos = rng.uniform(-5, 5, (4000, 3))
    r, g = rdf(pos, L, [True] * 3, rmax=4.0, nbins=20)
    assert np.allclose(g[2:], 1.0, atol=0.08)


def test_field_solver_screened_sphere(tmp_path):
    """Conducting sphere at V: phi = V R/r exp(-(r-R)/lambda)."""
    cfg = loads("""
[[species]]
count = 1
[domain]
box = [8, 8, 8]
[[objects]]
type = "sphere"
radius = 1.0
potential = 2.0
[field_solver]
grid = [81, 81, 81]
debye_length = 1.0
tolerance = 1e-10
""")
    g = fields.solve(cfg, geometry.expand(cfg["objects"]), log=lambda *a: None)
    x = g.axes()[0]
    i0 = len(x) // 2
    phi = g.phi[:, i0, i0]
    for r in (1.5, 2.0, 2.5):
        k = np.argmin(abs(x - r))
        exact = 2.0 / x[k] * math.exp(-(x[k] - 1.0))
        # finite box (phi = 0 at |x| = 4) and staircase sphere: a few percent
        assert abs(phi[k] - exact) / exact < 0.06, (r, phi[k], exact)


def test_pinc_reader(tmp_path):
    p = tmp_path / "E.grid.h5"
    E = np.zeros((5, 6, 7, 3))
    E[..., 2] = 3.0
    with h5py.File(p, "w") as f:
        f.attrs["Axis denormalization factor"] = [0.01]
        f.attrs["Quantity denormalization factor"] = [2.0]
        f["/n=10.0"] = E * 0.5
        f["/n=100.0"] = E
    g = fields.load(str(p))
    assert g.E.shape == (5, 6, 7, 3) and np.allclose(g.E[..., 2], 6.0)
    assert np.allclose(g.spacing, 0.01) and np.allclose(g.origin, -0.5 * np.array([4, 5, 6]) * 0.01)


LEGACY = {
    "main": """[simbox]
Lx = 15.0
Ly = 15.0
Lz = 20.0
[particles]
N = 200
Vxmax = 1.0
Vymax = 1.0
Vzmax = 1.0
Temp = 0.1
dist = True
mean = 1
stdDev = 0.2
maxwell_load = True
[screening]
k = 0.5
[gravity]
g_0 = 0.0
[cutoff radius]
rc = 1e-5
[boundary]
btype = periodic
[time]
tmax = 1.0
dt = 0.01
[diagnostics]
dumpPeriod = 5
dumpData = True
vtkData = False
realTime = False
[directory]
dataDir = data
[options]
parallelMode = True
""",
    "dust_fall": """[simbox]
Lx = 15.0
Ly = 15.0
Lz = 20.0
[particles]
N = 200
Vxmax = 1.0
Vymax = 1.0
Vzmax = 1.0
Temp = 0.1
dist = True
mean = 1
stdDev = 0.2
[screening]
k = 0.0
[gravity]
g_0= 0.0063
[cutoff radius]
rc = 1e-5
[boundary]
btype = mixed
[time]
tmax = 1.0
dt = 0.010
[diagnostics]
dumpPeriod = 5
vtkData = True
[directory]
dataDir = data
""",
    "dust_void": """[simbox]
Lx = 100
Ly = 100
Lz = 1
[particles]
N = 300
Vxmax = 1.0
Vymax = 1.0
Vzmax = 1.0
Temp = 0.1
dist = False
mean = 1
stdDev = 0.2
f0 = 1e-13
f_flow = 2.5e-14
[constant]
KB = 1.38e-23
[screening]
k = 0.5
lambda_c = 20
Gamma = 5
a = 3e-4
[dust]
Td = 9283.0
md = 1.7e-13
[neutral]
nu = 1.6e-4
Tn = 3481.2
[gravity]
g_0= 0.0
[cutoff radius]
rc = 1e-5
[boundary]
btype = periodic
[time]
tmax = 1.0
dt = 0.001
[diagnostics]
dumpPeriod = 5
[directory]
dataDir = data
""",
    "dust_exp": """scope = default
[simbox]
Lx  = 1.85
Ly  = 1.85
Lz  = 1.85
r_min = 1.2
r_max = 1.5
[particles]
N     = 300
Temp  = 0
dist = binary
stdDev = 2.2e-05, 2.9e-05
mass = 4.29e-14
a = 0.000426, 0.000438
Q = 2e-9
[plasma]
E = 0
k_fel = 3.16e-02
k_fion = 3.16e-02
Ffield = staticfields/F_digitized.npz
[neutrals]
k_drag = 6e-07
[screening]
k = 24
[gravity]
g_0= 0.0
[boundary]
btype = periodic
geometry = radial
[time]
tmax  = 1.5e1
dt    = 0.5e-2
[diagnostics]
dumpPeriod  = 5
vtkData     = False
[directory]
path = data
""",
}


@pytest.mark.parametrize("name", sorted(LEGACY))
def test_legacy_conversion_validates(name):
    text, notes = convert(LEGACY[name], name)
    cfg = loads(text)
    assert cfg["species"][0]["count"] in (200, 300)
    if name == "dust_fall":
        assert cfg["boundary"]["z"] == ["absorbing", "reflecting"]
        assert cfg["species"][0]["placement"]["region"] == "plane"
        assert cfg["forces"]["harmonic"]["omega2"][2] == pytest.approx(0.0063)
    if name == "dust_void":
        assert cfg["forces"]["void"] and cfg["neutrals"]["drag"] == "constant"
    if name == "dust_exp":
        assert cfg["domain"]["shape"] == "sphere" and len(cfg["forces"]["radial_table"]) == 2
    if name == "main":
        assert cfg["output"]["vtk"] is False      # legacy bool('False') bug fixed
        assert cfg["domain"]["box"] == [30.0, 30.0, 40.0]


def test_legacy_main_runs(tmp_path):
    text, _ = convert(LEGACY["main"], "legacy_main")
    p = tmp_path / "legacy.toml"
    p.write_text(text.replace('dir = "data"', f"dir = '{tmp_path / 'out'}'"))
    assert cli(["run", str(p), "-q"]) == 0


@pytest.mark.parametrize("example", ["yukawa_crystal", "dust_fall", "dust_void", "dust_cloud", "debris_cubesat"])
def test_examples_run_briefly(tmp_path, example):
    assert cli(["init", example, str(tmp_path / example)]) == 0
    inp = str(tmp_path / example / "input.toml")
    extra = ["--set", "time.steps=3", "--set", "output.every=1", "--set", "field_solver.grid=[24, 24, 24]"]
    if example in ("dust_cloud", "dust_void"):
        extra += ["--set", "species.0.count=1500"]
    assert cli(["run", inp, "-q"] + extra) == 0


def test_docs_cover_schema():
    from ppdyn.docgen import markdown
    md = markdown()
    for key in ("[[species]]", "[forces.void]", "[[forces.radial_table]]", "[field]", "[[objects]]", "debye_length"):
        assert key in md


def test_cli_entry_point():
    out = subprocess.run([sys.executable, "-m", "ppdyn", "--version"], capture_output=True, text=True)
    assert out.returncode == 0 and "PPDyn" in out.stdout
