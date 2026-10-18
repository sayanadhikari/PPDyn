"""Renders the README gallery from the five example runs.

    pip install pyvista scipy
    python render_all.py --prepare RUNDIR     # init + run the examples (about 8 min on 2 cores)
    python render_all.py RUNDIR [scene ...]   # render (about 8 min)

Scenes: cloud, crystal, void, fall, debris (default: all). Images and GIFs are
written to docs/media next to this script.
"""
import os, subprocess, sys, tempfile
import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *  # noqa
from ppdyn.post import Run

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TMP = tempfile.mkdtemp(prefix="ppdyn_render_")


# denser output where the renders draw trajectories
PREPARE = {"yukawa_crystal": [], "dust_void": [], "dust_cloud": [],
           "dust_fall": ["--set", "output.every=5"], "debris_cubesat": ["--set", "output.every=10"]}


def prepare(rundir):
    os.makedirs(rundir, exist_ok=True)
    for name, extra in PREPARE.items():
        d = os.path.join(rundir, name)
        if not os.path.isdir(d):
            subprocess.run(["ppdyn", "init", name, d], check=True)
        print("running", name, flush=True)
        subprocess.run(["ppdyn", "run", "input.toml", "-q", *extra], cwd=d, check=True)


def run(rundir, name, file=None):
    return Run(os.path.join(rundir, name, "output", f"{file or name}.h5"))


def min_image(d, L, periodic):
    d = d.copy()
    for a in range(3):
        if periodic[a]:
            d[:, a] -= L[a] * np.round(d[:, a] / L[a])
    return d


# ====================================================================== crystal
def crystal(rundir):
    r = run(rundir, "yukawa_crystal")
    p0 = r.position(0)
    c = (r.lo + r.hi) / 2

    def disp(k):  # displacement from the lattice site, through the periodic images
        steps = np.unique(np.linspace(0, k, 25).astype(int))
        tot = np.zeros_like(p0); prev = p0
        for s in steps[1:]:
            cur = r.position(s); tot += min_image(cur - prev, r.L, r.periodic); prev = cur
        return np.linalg.norm(tot, axis=1)

    def scene(k, path, size=(1100, 1000), az=0.0, e=118.0):
        k = k % len(r)
        p = r.position(k) - c
        d = disp(k) if k else np.zeros(r.n)
        pl = plotter(size)
        grains(pl, p, 0.44, d, cmap=HEAT, clim=(0, 4.0), res=16)
        h = r.L / 2
        pl.add_mesh(pv.Box((-h[0], h[0], -h[1], h[1], -h[2], h[2])).extract_all_edges(), color=SLATE, opacity=0.35, line_width=1.2)
        lights(pl, key=(0.7, -1.2, 1.4))
        a = np.radians(-38 + az)
        pl.camera_position = [(e * np.cos(a), e * np.sin(a), 0.62 * e), (0, 0, -1.8), (0, 0, 1)]
        pl.camera.view_angle = 19.5
        finish(pl, path, ssao_radius=1.2)

    a, b = f"{TMP}/c0.png", f"{TMP}/c1.png"
    scene(0, a); scene(-1, b)
    label(a, "bcc crystal, Γ = 400", "t = 0, every grain on its lattice site")
    label(b, "heated to Γ = 100", f"t = {r.times[-1]:.0f} ω₀⁻¹, the lattice has melted")
    out = f"{OUT}/crystal_melting.png"
    hstack([a, b], out, gap=4)
    colorbar(out, HEAT, "0", "4 a", "displacement from lattice site")
    # animation: melting while the camera turns
    frames = []
    ks = list(range(0, len(r), 2))
    for i, k in enumerate(ks):
        f = f"{TMP}/cf{i:03d}.png"
        scene(k, f, size=(1000, 640), az=i * 2.0, e=92.0)
        label(f, None, f"t = {r.times[k]:5.1f} ω₀⁻¹", corner=(28, 24))
        frames.append(f)
    gif(frames + [frames[-1]] * 8, f"{OUT}/crystal_melting.gif", fps=10, width=500, colors=112)
    print(out)


# ====================================================================== void
ORDER = LinearSegmentedColormap.from_list("order", [ORANGE, "#f7a073", CREAM, "#7fd3e8", CYAN, "#127fa3"])


def psi6(p):
    """|psi_6| of each grain from its Delaunay neighbours in the plane."""
    from scipy.spatial import Delaunay
    tri = Delaunay(p[:, :2])
    ind, nb = tri.vertex_neighbor_vertices
    out = np.zeros(len(p))
    for i in range(len(p)):
        j = nb[ind[i]:ind[i + 1]]
        d = p[j, :2] - p[i, :2]
        out[i] = abs(np.exp(6j * np.arctan2(d[:, 1], d[:, 0])).mean())
    return out


def void(rundir):
    r = run(rundir, "dust_void")
    p0 = r.position(0)
    rhat = p0[:, :2] / np.maximum(np.hypot(p0[:, 0], p0[:, 1]), 1e-9)[:, None]
    cum = {0: np.zeros_like(p0)}
    tot = np.zeros_like(p0); prev = p0
    for k in range(1, len(r)):
        cur = r.position(k); tot = tot + min_image(cur - prev, r.L, r.periodic); prev = cur; cum[k] = tot.copy()
    vmax = 6.0

    def scene(k, path, size=(1800, 1000)):
        k = k % len(r)
        p = r.position(k)
        push = (cum[k][:, :2] * rhat).sum(1)          # outward radial displacement
        pl = plotter(size)
        grains(pl, p, 0.62, push, cmap=HEAT, clim=(0, vmax), res=14)
        floor = pv.Plane(center=(0, 0, -4.2), i_size=r.L[0] * 1.6, j_size=r.L[1] * 1.6, i_resolution=48, j_resolution=48)
        pl.add_mesh(floor.extract_all_edges(), color="#2b3a55", opacity=0.55, line_width=1)
        lights(pl, key=(0.5, -1.0, 1.6))
        e = 165.0; a = np.radians(-62)
        pl.camera_position = [(e * np.cos(a), e * np.sin(a), 0.55 * e), (0, 6, -6), (0, 0, 1)]
        pl.camera.view_angle = 34
        finish(pl, path, ssao_radius=1.0)

    out = f"{OUT}/dust_void.png"
    scene(-1, out)
    label(out, "Void opening in a dust monolayer",
          "5000 grains, κ = 0.5, an ion-drag force pushes grains out of the centre\ncolour = how far each grain has been pushed outward")
    colorbar(out, HEAT, "0", f"{vmax:.0f} a", "outward displacement")
    frames = []
    for i, k in enumerate(range(len(r))):
        f = f"{TMP}/vf{i:03d}.png"
        scene(k, f, size=(1000, 640))
        label(f, None, f"t = {r.times[k]:5.1f} ω₀⁻¹", corner=(26, 22))
        frames.append(f)
    gif(frames + [frames[-1]] * 10, f"{OUT}/dust_void.gif", fps=10, width=500, colors=96)
    print(out)


# ====================================================================== cloud (hero)
def cloud(rundir):
    r = run(rundir, "dust_cloud")
    p = r.position(-1); big = r.radius > np.median(r.radius)
    keep = ~((p[:, 0] > 0) & (p[:, 1] < 0) & (p[:, 2] > -8))
    R = r.L[0] / 2

    def scene(path, size, az=0.0, cut=True):
        kk = keep if cut else np.ones(r.n, bool)
        pl = plotter(size)
        grains(pl, p[kk], 0.42, big[kk].astype(float), cmap=[CYAN, ORANGE], clim=(0, 1), res=14)
        pl.add_mesh(pv.Sphere(radius=R, theta_resolution=120, phi_resolution=80), color="#9fb7d9", opacity=0.045, smooth_shading=True)
        lights(pl, key=(0.6, -1.0, 1.4))
        e = np.hypot(62, 66); a = np.radians(-46.8 + az)
        pl.camera_position = [(e * np.cos(a), e * np.sin(a), 44), (0, 0, -2), (0, 0, 1)]
        pl.camera.view_angle = 30
        finish(pl, path)

    out = f"{OUT}/hero_dust_cloud.png"
    scene(out, (2000, 1000))
    label(out, "Binary dust cloud in microgravity",
          "20 000 grains, cut open to show the void.\nSmall grains line the void,\nlarger grains settle outside.")
    im = Image.open(out); d = ImageDraw.Draw(im); f = font(22, "Regular")
    for y, col, t in ((880, CYAN, "3.45 µm"), (918, ORANGE, "3.55 µm")):
        d.ellipse((1760, y - 9, 1778, y + 9), fill=col); d.text((1792, y), t, font=f, fill="#c9d4e3", anchor="lm")
    im.save(out)
    print(out)


# ====================================================================== dust fall
def fall(rundir):
    r = run(rundir, "dust_fall")
    rad = r.radius * r.scales.length * 1e6           # µm
    disp = 0.10 + 0.034 * rad
    h = r.L / 2; c = (r.lo + r.hi) / 2
    tail = 16

    def scene(k, path, size=(1400, 1000), az=0.0):
        k = k % len(r)
        pl = plotter(size)
        p = r.position(k) - c
        alive = r.active(k)
        if alive.any():
            grains(pl, p[alive], disp[alive], rad[alive], cmap=HEAT, clim=(2.5, 6.0), res=12)
            # streaks over the last frames
            ks = list(range(max(0, k - tail), k + 1))
            if len(ks) > 1:
                P = np.stack([r.position(j) - c for j in ks], 1)[alive]
                ok = (np.abs(np.diff(P, axis=1)).max(axis=(1, 2)) < 0.4 * r.L.min())   # no wrap-around jumps
                P = P[ok]; rr_ = rad[alive][ok]
                n, m = P.shape[:2]
                lines = np.hstack([np.full((n, 1), m), np.arange(n * m).reshape(n, m)]).ravel()
                poly = pv.PolyData(P.reshape(-1, 3), lines=lines)
                poly["s"] = np.repeat(rr_, m)
                pl.add_mesh(poly, scalars="s", cmap=HEAT, clim=(2.5, 6.0), opacity=0.6, line_width=2.5, show_scalar_bar=False)
        landed = ~alive
        if landed.any():
            q = p[landed].copy(); q[:, 2] = -h[2] + disp[landed]
            grains(pl, q, disp[landed], rad[landed], cmap=HEAT, clim=(2.5, 6.0), res=10)
        plate = pv.Box((-h[0] * 1.08, h[0] * 1.08, -h[1] * 1.08, h[1] * 1.08, -h[2] - 0.6, -h[2]))
        pl.add_mesh(plate, color="#8a96a8", ambient=0.2, diffuse=0.7, specular=0.8, specular_power=60, smooth_shading=False)
        pl.add_mesh(pv.Box((-h[0], h[0], -h[1], h[1], -h[2], h[2])).extract_all_edges(), color=SLATE, opacity=0.35, line_width=1.2)
        lights(pl, key=(0.6, -1.1, 1.5))
        e = 66.0; a = np.radians(-55 + az)
        pl.camera_position = [(e * np.cos(a), e * np.sin(a), 0.42 * e), (0, 0, -1.5), (0, 0, 1)]
        pl.camera.view_angle = 30
        finish(pl, path, ssao_radius=0.8)

    tsc = r.scales.time * 1e3
    kmid = int(np.argmin(abs(r.times * tsc - 44.0)))
    out = f"{OUT}/dust_fall.png"
    scene(kmid, out)
    label(out, "Charged grains falling onto an electrode",
          f"2000 melamine grains in 10 Pa argon, t = {r.times[kmid] * tsc:.0f} ms\nlarge grains fall faster and land first")
    colorbar(out, HEAT, "2.5 µm", "6 µm", "grain radius")
    frames = []
    for i, k in enumerate(range(0, len(r), 2)):
        f = f"{TMP}/ff{i:03d}.png"
        scene(k, f, size=(1000, 640))
        label(f, None, f"t = {r.times[k] * tsc:4.1f} ms", corner=(26, 22))
        frames.append(f)
    gif(frames + [frames[-1]] * 10, f"{OUT}/dust_fall.gif", fps=14, width=500, colors=112)
    print(out)


# ====================================================================== CubeSat debris
def cubesat_meshes(input_path, scale):
    from ppdyn import config, geometry
    cfg = config.load(input_path)
    out = []
    for k, pr in enumerate(geometry.expand(cfg["objects"])):
        cc = np.asarray(pr["center"]) / scale
        if pr["type"] == "box":
            sz = np.asarray(pr["size"]) / scale / 2
            out.append(("body" if k == 0 else "wing", pv.Box((cc[0] - sz[0], cc[0] + sz[0], cc[1] - sz[1], cc[1] + sz[1], cc[2] - sz[2], cc[2] + sz[2]))))
        elif pr["type"] == "cylinder":
            ax = np.eye(3)["xyz".index(pr["axis"])]
            out.append(("antenna", pv.Cylinder(center=cc, direction=ax, radius=pr["radius"] / scale * 2.5, height=pr["length"] / scale, resolution=24)))
        else:
            out.append(("body", pv.Sphere(center=cc, radius=pr["radius"] / scale)))
    return out


def add_cubesat(pl, meshes):
    for kind, m in meshes:
        if kind == "body":
            pl.add_mesh(m, color="#b8bfc9", ambient=0.2, diffuse=0.75, specular=0.9, specular_power=70)
            b = m.bounds
            # solar cells on the four long faces: thin dark tiles slightly proud of the surface
            for face in range(4):
                for j in range(3):
                    z0 = b[4] + (b[5] - b[4]) * (j + 0.08) / 3; z1 = b[4] + (b[5] - b[4]) * (j + 0.92) / 3
                    e = 0.01
                    x0, x1 = b[0] + 0.12 * (b[1] - b[0]), b[1] - 0.12 * (b[1] - b[0])
                    y0, y1 = b[2] + 0.12 * (b[3] - b[2]), b[3] - 0.12 * (b[3] - b[2])
                    tiles = [(x0, x1, b[2] - e, b[2]), (x0, x1, b[3], b[3] + e), (b[0] - e, b[0], y0, y1), (b[1], b[1] + e, y0, y1)]
                    t = tiles[face]
                    pl.add_mesh(pv.Box((t[0], t[1], t[2], t[3], z0, z1)), color="#1b2c55", ambient=0.15, specular=1.0, specular_power=90)
        elif kind == "wing":
            pl.add_mesh(m, color="#16264d", ambient=0.12, diffuse=0.7, specular=1.0, specular_power=100)
            b = m.bounds
            g = pv.Plane(center=((b[0] + b[1]) / 2, (b[2] + b[3]) / 2, b[5] + 0.01), i_size=b[1] - b[0], j_size=b[3] - b[2],
                         i_resolution=8, j_resolution=4)
            pl.add_mesh(g.extract_all_edges(), color="#8fa6c9", opacity=0.55, line_width=1.2)
            g2 = g.translate((0, 0, -(b[5] - b[4]) - 0.02))
            pl.add_mesh(g2.extract_all_edges(), color="#8fa6c9", opacity=0.55, line_width=1.2)
        else:
            pl.add_mesh(m, color="#dfe4ea", specular=1.0, specular_power=80)


def debris(rundir):
    name = "debris_cubesat"
    r = run(rundir, name)
    L0 = r.scales.length
    meshes = cubesat_meshes(os.path.join(rundir, name, "input.toml"), L0)
    import h5py
    with h5py.File(os.path.join(rundir, name, "field.h5")) as f:
        x, y, z, phi = f["x"][:] / L0, f["y"][:] / L0, f["z"][:] / L0, f["phi"][:]
    grid = pv.RectilinearGrid(x, y, z); grid["phi"] = phi.ravel(order="F")
    sheath = [grid.contour([v], scalars="phi").smooth(n_iter=60) for v in (-0.35, -1.5)]
    rad = r.radius * L0 * 1e6
    disp = 0.03 + 0.008 * rad
    tail = 14
    hit_t = r.final["absorb_time"]

    def scene(k, path, size=(1800, 1000), az=0.0, e=34.0):
        k = k % len(r)
        pl = plotter(size)
        add_cubesat(pl, meshes)
        pl.add_mesh(sheath[0], color=CYAN, opacity=0.07, smooth_shading=True, specular=0.4)
        pl.add_mesh(sheath[1], color=CYAN, opacity=0.12, smooth_shading=True, specular=0.4)
        p = r.position(k); alive = r.active(k)
        if alive.any():
            grains(pl, p[alive], disp[alive], color="#e9eef6", res=10)
            ks = list(range(max(0, k - tail), k + 1))
            if len(ks) > 1:
                P = np.stack([r.position(j) for j in ks], 1)[alive]
                n, m = P.shape[:2]
                lines = np.hstack([np.full((n, 1), m), np.arange(n * m).reshape(n, m)]).ravel()
                poly = pv.PolyData(P.reshape(-1, 3), lines=lines)
                poly["a"] = np.tile(np.linspace(0, 1, m), n)
                pl.add_mesh(poly, scalars="a", cmap=["#0d1626", "#9fd9ea"], opacity=0.5, line_width=1.6, show_scalar_bar=False)
        dep = (~alive) & (r.final["hit_object"] >= 0)
        if dep.any():
            grains(pl, p[dep], disp[dep], color=ORANGE, res=10)
        lights(pl, key=(0.9, -1.0, 1.1), dist=80)
        a = np.radians(-58 + az)
        pl.camera_position = [(e * np.cos(a), e * np.sin(a), 0.48 * e), (0, 0, -2.2), (0, 0, 1)]
        pl.camera.view_angle = 32
        finish(pl, path, ssao_radius=0.4)

    tsc = r.scales.time
    out = f"{OUT}/debris_cubesat.png"
    kk = int(np.argmin(abs(r.times * tsc - 6.0)))
    scene(kk, out)
    label(out, "Charged micro-debris on a 3U CubeSat",
          f"3000 aluminium grains at +2 V drift onto a −10 V spacecraft\nt = {r.times[kk] * tsc:.0f} s, translucent shells: sheath potential −1.5 V and −0.35 V")
    im = Image.open(out); d = ImageDraw.Draw(im); f = font(22, "Regular")
    for yy, col, t in ((900, "#e9eef6", "drifting"), (938, ORANGE, "deposited")):
        d.ellipse((1620, yy - 9, 1638, yy + 9), fill=col); d.text((1652, yy), t, font=f, fill="#c9d4e3", anchor="lm")
    im.save(out)
    frames = []
    for i, k in enumerate(range(0, len(r), 5)):
        f = f"{TMP}/df{i:03d}.png"
        scene(k, f, size=(1000, 640), az=i * 0.8, e=44.0)
        label(f, None, f"t = {r.times[k] * tsc:4.1f} s", corner=(26, 22))
        frames.append(f)
    gif(frames + [frames[-1]] * 10, f"{OUT}/debris_cubesat.gif", fps=12, width=500, colors=112)
    print(out)


SCENES = {"cloud": cloud, "crystal": crystal, "void": void, "fall": fall, "debris": debris}


if __name__ == "__main__":
    if sys.argv[1] == "--prepare":
        prepare(sys.argv[2]); sys.exit()
    rundir = sys.argv[1]
    for s in (sys.argv[2:] or SCENES):
        SCENES[s](rundir)
    stills = {"cloud": "hero_dust_cloud", "crystal": "crystal_melting", "void": "dust_void",
              "fall": "dust_fall", "debris": "debris_cubesat"}
    for s in (sys.argv[2:] or SCENES):
        compress_png(f"{OUT}/{stills[s]}.png")
