"""3D rendering (needs the optional PyVista dependency and an offscreen OpenGL context)."""
import os

import numpy as np
import pytest

pv = pytest.importorskip("pyvista")
PIL = pytest.importorskip("PIL")

import ppdyn  # noqa: E402
from ppdyn import loads  # noqa: E402
from ppdyn.cli import main as cli  # noqa: E402
from ppdyn.post import Run, render  # noqa: E402

pytestmark = pytest.mark.render


def _run(tmp_path, extra=""):
    txt = f"""
[run]
name = "small"
progress = false
[domain]
box = [12, 12, 12]
[boundary]
z = ["absorbing", "reflecting"]
[[species]]
count = 200
[species.radius]
distribution = "binary"
values = [1.0, 1.3]
[species.velocity]
temperature = 0.05
[forces.gravity]
g = [0.0, 0.0, -0.5]
[[objects]]
type = "sphere"
center = [0.0, 0.0, -2.0]
radius = 1.5
action = "absorb"
[time]
dt = 0.02
steps = 60
[output]
dir = '{tmp_path / 'out'}'
every = 10
quantities = ["position", "velocity", "charge", "state"]
{extra}
"""
    sim, path = ppdyn.run(loads(txt), progress=False)
    return path


def _not_blank(path):
    from PIL import Image
    a = np.asarray(Image.open(path).convert("RGB"), float)
    return a.std() > 5.0


def test_scene_reads_geometry_back_from_the_file(tmp_path):
    with Run(_run(tmp_path)) as r:
        sc = render.Scene(r)
        assert sc.name == "small"
        assert (2, 0) in sc.absorbing
        assert len(sc.prims) == 1 and sc.prims[0]["type"] == "sphere"
        assert sc.time_text(-1).startswith("t = ")


@pytest.mark.parametrize("color", ["speed", "radius", "displacement", "outward", "state", "charge", "z"])
def test_colour_values(tmp_path, color):
    with Run(_run(tmp_path)) as r:
        v, kind, label, unit = render.colour_values(render.Scene(r), -1, color)
        assert v.shape == (r.n,) and np.isfinite(v).all()
        assert kind in ("continuous", "categorical")
        if color == "radius":
            assert kind == "categorical" and len(label) == 2       # binary sizes get a legend


def test_snapshot_compare_and_movie(tmp_path):
    with Run(_run(tmp_path)) as r:
        p = render.snapshot(r, str(tmp_path / "a.png"), color="radius", cut="quadrant", trails=3, size=(480, 320),
                            title="test", subtitle="{time}")
        assert _not_blank(p)
        c = render.compare(r, str(tmp_path / "c.png"), frames=(0, -1), color="displacement", size=(320, 300))
        from PIL import Image
        assert Image.open(c).size[0] > 600
        g = render.movie(r, str(tmp_path / "m.gif"), every=3, orbit=30, size=(320, 200), width=200)
        assert Image.open(g).n_frames >= 3


def test_unknown_option_is_rejected(tmp_path):
    with Run(_run(tmp_path)) as r:
        with pytest.raises(TypeError, match="unknown render option"):
            render.snapshot(r, str(tmp_path / "x.png"), colour="speed")


def test_cli_render(tmp_path):
    h5 = _run(tmp_path)
    out = str(tmp_path / "cli.png")
    assert cli(["render", h5, "-o", out, "--size", "400x300", "--color", "state", "--azimuth", "10"]) == 0
    assert os.path.exists(out) and _not_blank(out)
    gif = str(tmp_path / "cli.gif")
    assert cli(["render", h5, "--movie", gif, "--every", "3", "--width", "160"]) == 0
    assert os.path.exists(gif)


def test_every_bundled_example_has_a_preset():
    examples = os.path.join(os.path.dirname(ppdyn.__file__), "examples")
    names = {d for d in os.listdir(examples) if os.path.isdir(os.path.join(examples, d)) and not d.startswith("_")}
    assert names <= set(render.PRESETS)
    for p in render.PRESETS.values():
        assert set(p["look"]) <= set(render.DEFAULTS)
