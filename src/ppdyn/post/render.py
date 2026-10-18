"""Publication-quality 3D renders and animations of PPDyn runs (PyVista / VTK).

Everything needed to draw a run is read back from its HDF5 file: particles,
the stored input (domain, absorbing walls, objects such as a CubeSat) and, if
present, the field file next to the input (sheath isosurfaces).

    from ppdyn.post import Run, render
    with Run("output/dust_cloud.h5") as r:
        render.snapshot(r, "cloud.png", color="radius", cut="quadrant")
        render.movie(r, "cloud.gif", orbit=90)
        render.gallery(r, "media/")          # the README images of a bundled example

or from the shell: ``ppdyn render output/dust_cloud.h5 [--movie cloud.gif] [--gallery media/]``.
Needs ``pip install "ppdyn[render]"``. On a headless machine VTK renders
offscreen through EGL or OSMesa; on Linux CI use a virtual display (xvfb).
"""
from __future__ import annotations

import os
import tempfile

import numpy as np

from .analysis import Run

NAVY, CYAN, ORANGE, SLATE, CREAM = "#0f1b2d", "#19b5d8", "#f2662f", "#7d8ea6", "#f4efe6"
BG_BOTTOM, BG_TOP = "#070b12", "#1d2a40"
HEAT_COLORS = ["#127fa3", CYAN, "#bfe9f2", CREAM, "#f7a073", ORANGE, "#b8361a"]
COLORS = ("speed", "radius", "species", "displacement", "outward", "charge", "state", "z", "none")
STATE_COLORS = {0: "#e9eef6", 1: ORANGE}          # drifting, deposited


def _pv():
    try:
        import pyvista as pv
    except ImportError as exc:  # pragma: no cover
        raise ImportError('3D rendering needs PyVista: pip install "ppdyn[render]"') from exc
    pv.OFF_SCREEN = True
    return pv


def _heat():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("ppdyn_heat", HEAT_COLORS)


# ====================================================================== scene description
class Scene:
    """Geometry of a run in reduced units: box, sphere, absorbing walls, objects, field."""

    def __init__(self, run: Run, field_file=None):
        from .. import config, geometry
        self.run = run
        self.lo, self.hi, self.L = run.lo, run.hi, run.L
        self.center = (run.lo + run.hi) / 2
        self.size = float(run.L.max())
        try:
            cfg = config.loads(run.input_text) if run.input_text else None
        except Exception:  # pragma: no cover - inputs from other versions
            cfg = None
        self.cfg = cfg
        self.name = cfg["run"]["name"] if cfg else os.path.splitext(os.path.basename(run.path))[0]
        L0 = run.scales.length
        self.si = bool(cfg and cfg["units"]["system"] == "SI")
        s = L0 if self.si else 1.0
        self.sphere = None
        if cfg and cfg["domain"]["shape"] == "sphere":
            self.sphere = (np.asarray(cfg["domain"]["center"], float) / s, cfg["domain"]["radius"] / s)
        self.absorbing = []                        # (axis, side) of absorbing walls
        if cfg:
            for a, ax in enumerate("xyz"):
                for side in (0, 1):
                    if cfg["boundary"][ax][side] == "absorbing":
                        self.absorbing.append((a, side))
        self.prims = []
        if cfg and cfg.get("objects"):
            for k, p in enumerate(geometry.expand(cfg["objects"])):
                q = dict(p)
                q["center"] = np.asarray(p["center"], float) / s
                if p["type"] == "box":
                    q["size"] = np.asarray(p["size"], float) / s
                else:
                    q["radius"] = (p.get("radius") or 0.0) / s
                    q["length"] = (p.get("length") or 0.0) / s
                q["role"] = ("body" if k == 0 else "wing" if p["type"] == "box" else "antenna") \
                    if cfg["objects"][p["parent"]]["type"] == "cubesat" else "object"
                self.prims.append(q)
        self.field = None
        f = field_file or (cfg.get("field") or {}).get("file") if cfg else field_file
        if f:
            for base in (os.path.dirname(os.path.dirname(os.path.abspath(run.path))),
                         os.path.dirname(os.path.abspath(run.path)), os.getcwd()):
                p = f if os.path.isabs(f) else os.path.join(base, f)
                if os.path.exists(p):
                    self.field = p
                    break

    # ------------------------------------------------------------------ units
    def time_text(self, k):
        t = self.run.times[k % len(self.run)]
        if not self.si:
            return f"t = {t:.0f} ω₀⁻¹" if t >= 10 else f"t = {t:.1f} ω₀⁻¹"
        ts = t * self.run.scales.time
        return f"t = {ts * 1e3:.0f} ms" if ts < 1 else f"t = {ts:.0f} s"

    def radius_um(self):
        return self.run.radius * self.run.scales.length * 1e6 if self.si else self.run.radius

    def frame_at(self, time):
        """Snapshot index closest to `time` (seconds for SI runs, ω₀⁻¹ otherwise)."""
        t = self.run.times * (self.run.scales.time if self.si else 1.0)
        return int(np.argmin(abs(t - time)))


# ====================================================================== particle colours
def _min_image(d, L, periodic):
    d = d.copy()
    for a in range(3):
        if periodic[a]:
            d[:, a] -= L[a] * np.round(d[:, a] / L[a])
    return d


class _Displacement:
    """Unwrapped displacement since frame 0, accumulated through periodic images."""

    def __init__(self, run):
        self.run, self.cache, self.tot, self.k = run, {0: np.zeros((run.n, 3))}, np.zeros((run.n, 3)), 0
        self.prev = run.position(0)

    def __call__(self, k):
        k %= len(self.run)
        if k in self.cache:
            return self.cache[k]
        if k < self.k:                       # restart from the beginning
            self.tot, self.k, self.prev = np.zeros((self.run.n, 3)), 0, self.run.position(0)
        while self.k < k:
            self.k += 1
            cur = self.run.position(self.k)
            self.tot = self.tot + _min_image(cur - self.prev, self.run.L, self.run.periodic)
            self.prev = cur
            self.cache[self.k] = self.tot.copy()
        return self.cache[k]


def colour_values(scene: Scene, k: int, color: str, disp=None):
    """Returns (values, kind, label, unit) with kind 'continuous' or 'categorical'."""
    r = scene.run
    k %= len(r)
    if color == "speed":
        v = r.velocity(k)
        sp = np.linalg.norm(v, axis=1)
        if scene.si:
            sp = sp * r.scales.length / r.scales.time
        return sp, "continuous", "speed", "m/s" if scene.si else "a ω₀"
    if color == "radius":
        ru = scene.radius_um()
        uniq = np.unique(np.round(ru, 4))
        if len(uniq) <= 4:
            idx = np.searchsorted(uniq, np.round(ru, 4)).astype(float)
            labels = [f"{u:.2f} µm" if scene.si else f"r = {u:g}" for u in uniq]
            return idx, "categorical", labels, None
        return ru, "continuous", "grain radius", "µm" if scene.si else "a"
    if color == "species":
        return r.species.astype(float), "categorical", list(r.species_names), None
    if color in ("displacement", "outward"):
        d = (disp or _Displacement(r))(k)
        if color == "displacement":
            return np.linalg.norm(d, axis=1), "continuous", "displacement from start", "a"
        p0 = r.position(0)[:, :2] - scene.center[:2]
        rh = p0 / np.maximum(np.linalg.norm(p0, axis=1), 1e-12)[:, None]
        return (d[:, :2] * rh).sum(1), "continuous", "outward displacement", "a"
    if color == "charge":
        q = r.charge(k)
        return (q * r.scales.charge / 1.602176634e-19 if scene.si else q), "continuous", "charge", "e" if scene.si else "Q₀"
    if color == "state":
        return (r.state(k) != 0).astype(float), "categorical", ["drifting", "deposited"], None
    if color == "z":
        return r.position(k)[:, 2], "continuous", "height", "a"
    return np.zeros(r.n), "categorical", [""], None


# ====================================================================== drawing primitives
def _plotter(size):
    pv = _pv()
    pl = pv.Plotter(off_screen=True, window_size=tuple(int(v) for v in size), lighting="none")
    pl.set_background(BG_BOTTOM, top=BG_TOP)
    return pl


def _lights(pl, center, dist, key=(1.0, -1.6, 1.2), key_int=1.25):
    pv = _pv()
    c = np.asarray(center, float)
    k = np.asarray(key, float); k /= np.linalg.norm(k)
    fill = np.array([-k[0], k[1] * 0.4, -0.3]); fill /= np.linalg.norm(fill)
    rim = np.array([-k[0] * 0.3, -k[1], 0.9]); rim /= np.linalg.norm(rim)
    for d, i, col in ((k, key_int, "#fff4e6"), (fill, 0.45, "#a9c8ff"), (rim, 0.6, "#cfe3ff")):
        pl.add_light(pv.Light(position=tuple(c + dist * d), focal_point=tuple(c), intensity=i, color=col,
                              light_type="scene light"))
    pl.add_light(pv.Light(intensity=0.18, light_type="headlight"))


SHADE = dict(smooth_shading=True, ambient=0.16, diffuse=0.9, specular=0.5, specular_power=48)


def _spheres(pl, pos, radius, scalars=None, cmap=None, clim=None, color=None, res=14, opacity=1.0):
    pv = _pv()
    if len(pos) == 0:
        return
    pts = pv.PolyData(np.ascontiguousarray(pos, float))
    pts["r"] = np.broadcast_to(np.asarray(radius, float), (len(pos),)).copy()
    sph = pv.Sphere(radius=1.0, theta_resolution=res, phi_resolution=max(8, res * 2 // 3))
    if scalars is not None:
        pts["s"] = np.asarray(scalars, float)
    g = pts.glyph(geom=sph, scale="r", orient=False)
    if scalars is not None:
        pl.add_mesh(g, scalars="s", cmap=cmap, clim=clim, show_scalar_bar=False, opacity=opacity, **SHADE)
    else:
        pl.add_mesh(g, color=color, opacity=opacity, **SHADE)


def _lines(pl, P, scalars, cmap, clim, opacity, width):
    pv = _pv()
    n, m = P.shape[:2]
    if n == 0 or m < 2:
        return
    cells = np.hstack([np.full((n, 1), m), np.arange(n * m).reshape(n, m)]).ravel()
    poly = pv.PolyData(P.reshape(-1, 3), lines=cells)
    poly["s"] = scalars
    pl.add_mesh(poly, scalars="s", cmap=cmap, clim=clim, opacity=opacity, line_width=width, show_scalar_bar=False)


def _box(lo, hi):
    return _pv().Box((lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]))


def _draw_objects(pl, scene: Scene):
    pv = _pv()
    for p in scene.prims:
        c = p["center"]
        if p["type"] == "box":
            h = p["size"] / 2
            m = _box(c - h, c + h)
        elif p["type"] == "sphere":
            m = pv.Sphere(center=c, radius=p["radius"], theta_resolution=48, phi_resolution=32)
        else:
            ax = np.eye(3)["xyz".index(p["axis"])]
            rr = p["radius"] * (2.5 if p["role"] == "antenna" else 1.0)     # thin masts drawn visible
            m = pv.Cylinder(center=c, direction=ax, radius=rr, height=p["length"], resolution=32)
        if p["role"] == "body" and p["type"] == "box":
            pl.add_mesh(m, color="#b8bfc9", ambient=0.2, diffuse=0.75, specular=0.9, specular_power=70)
            _solar_cells(pl, m.bounds)
        elif p["role"] == "wing":
            pl.add_mesh(m, color="#16264d", ambient=0.12, diffuse=0.7, specular=1.0, specular_power=100)
            b = m.bounds
            for z in (b[5] + 0.01, b[4] - 0.01):
                g = pv.Plane(center=((b[0] + b[1]) / 2, (b[2] + b[3]) / 2, z), i_size=b[1] - b[0], j_size=b[3] - b[2],
                             i_resolution=8, j_resolution=4)
                pl.add_mesh(g.extract_all_edges(), color="#8fa6c9", opacity=0.55, line_width=1.2)
        elif p["role"] == "antenna":
            pl.add_mesh(m, color="#dfe4ea", specular=1.0, specular_power=80)
        else:
            pl.add_mesh(m, color="#aab3bf", ambient=0.2, diffuse=0.75, specular=0.8, specular_power=60)


def _solar_cells(pl, b):
    """Dark solar-cell tiles on the four long faces of a CubeSat body."""
    e = 0.01
    x0, x1 = b[0] + 0.12 * (b[1] - b[0]), b[1] - 0.12 * (b[1] - b[0])
    y0, y1 = b[2] + 0.12 * (b[3] - b[2]), b[3] - 0.12 * (b[3] - b[2])
    for j in range(3):
        z0 = b[4] + (b[5] - b[4]) * (j + 0.08) / 3
        z1 = b[4] + (b[5] - b[4]) * (j + 0.92) / 3
        for t in ((x0, x1, b[2] - e, b[2]), (x0, x1, b[3], b[3] + e), (b[0] - e, b[0], y0, y1), (b[1], b[1] + e, y0, y1)):
            pl.add_mesh(_box((t[0], t[2], z0), (t[1], t[3], z1)), color="#1b2c55", ambient=0.15, specular=1.0,
                        specular_power=90)


def _sheath_meshes(scene: Scene, levels):
    import h5py
    pv = _pv()
    L0 = scene.run.scales.length if scene.si else 1.0
    with h5py.File(scene.field) as f:
        x, y, z, phi = f["x"][:] / L0, f["y"][:] / L0, f["z"][:] / L0, f["phi"][:]
    grid = pv.RectilinearGrid(x, y, z)
    grid["phi"] = phi.ravel(order="F")
    out = []
    for v in levels:
        m = grid.contour([v], scalars="phi")
        if m.n_points:
            out.append(m.smooth(n_iter=60))
    return out


# ====================================================================== 2D overlays (PIL)
def _font(size, weight="Medium"):
    from PIL import ImageFont
    for path in (f"/usr/share/fonts/opentype/inter/Inter-{weight}.otf",):
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    from matplotlib import font_manager
    w = {"SemiBold": "semibold", "Regular": "normal", "Medium": "medium"}.get(weight, "normal")
    return ImageFont.truetype(font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight=w)), size)


def _overlay(path, title=None, sub=None, corner=(36, 30), legend=None, cbar=None, scale=1.0):
    from PIL import Image, ImageDraw
    im = Image.open(path).convert("RGB")
    d = ImageDraw.Draw(im)
    W, H = im.size
    s = scale
    x, y = int(corner[0] * s), int(corner[1] * s)
    if title:
        d.text((x, y), title, font=_font(int(34 * s), "SemiBold"), fill="#eef3fa")
        y += int(46 * s)
    if sub:
        for line in sub.split("\n"):
            d.text((x, y), line, font=_font(int(22 * s), "Regular"), fill="#9fb0c8")
            y += int(31 * s)
    f = _font(int(21 * s), "Regular")
    if legend:
        x0, y0 = W - int(250 * s), H - int(40 * s) - int(38 * s) * (len(legend) - 1)
        for i, (col, text) in enumerate(legend):
            yy = y0 + int(38 * s) * i
            r_ = int(9 * s)
            d.ellipse((x0 - r_, yy - r_, x0 + r_, yy + r_), fill=col)
            d.text((x0 + int(22 * s), yy), text, font=f, fill="#c9d4e3", anchor="lm")
    if cbar:
        cmap, lo, hi, text = cbar
        length = int(360 * s)
        x0, y0 = W - length - int(50 * s), H - int(70 * s)
        bar = (np.array([cmap(t)[:3] for t in np.linspace(0, 1, length)]) * 255).astype(np.uint8)
        im.paste(Image.fromarray(np.repeat(bar[None], max(4, int(10 * s)), 0)), (x0, y0))
        d.text((x0, y0 - int(28 * s)), text, font=f, fill="#c9d4e3")
        d.text((x0, y0 + int(16 * s)), lo, font=f, fill="#9fb0c8")
        d.text((x0 + length, y0 + int(16 * s)), hi, font=f, fill="#9fb0c8", anchor="ra")
    im.save(path)


def compress_png(path, colors=256):
    """Palette PNG with dithering: about 3x smaller, visually identical for these renders."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    im.quantize(colors=colors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG).save(path, optimize=True)


def _fmt(v):
    return f"{v:.3g}"


# ====================================================================== core frame renderer
DEFAULTS = dict(color="speed", cmap=None, clim=None, radius="auto", cut=None, trails=0, trail_opacity=0.5,
                azimuth=-50.0, elevation=28.0, distance=None, view_angle=30.0, focus=(0.0, 0.0, 0.0),
                box=None, floor=False, objects=True, sheath=None, size=(1600, 1000), res=14,
                title=None, subtitle=None, legend=True, colorbar=True, ssao=True, ssao_radius=None,
                key_light=(1.0, -1.6, 1.2))


class Renderer:
    """Holds per-run state (displacements, sheath meshes) so movies do not redo work."""

    def __init__(self, run: Run, field_file=None, **opts):
        unknown = set(opts) - set(DEFAULTS)
        if unknown:
            raise TypeError(f"unknown render option(s): {', '.join(sorted(unknown))}")
        self.scene = Scene(run, field_file)
        self.o = {**DEFAULTS, **opts}
        if self.o["color"] not in COLORS:
            raise ValueError(f"color must be one of {', '.join(COLORS)}")
        self.disp = _Displacement(run)
        self._sheath = None
        self._radius = self._display_radius()
        self._clim = self.o["clim"]
        if self._clim is None:
            vals = [colour_values(self.scene, k, self.o["color"], self.disp)[0] for k in (0, -1)]
            kind = colour_values(self.scene, 0, self.o["color"], self.disp)[1]
            if kind == "continuous":
                allv = np.concatenate(vals)
                self._clim = (float(np.percentile(allv, 1)), float(np.percentile(allv, 99)))

    def _display_radius(self):
        sc, r = self.scene, self.scene.run
        rad = self.o["radius"]
        if isinstance(rad, (int, float)):
            return np.full(r.n, float(rad))
        if isinstance(rad, (tuple, list)):                       # (base, per µm or per a of grain radius)
            return rad[0] + rad[1] * sc.radius_um()
        # auto: a third of the median nearest-neighbour distance, scaled mildly with grain size
        from scipy.spatial import cKDTree
        p = r.position(-1)[r.active(-1)] if r.active(-1).any() else r.position(0)
        nn = np.median(cKDTree(p).query(p, k=2)[0][:, 1]) if len(p) > 1 else sc.size / 20
        base = float(np.clip(0.33 * nn, 0.002 * sc.size, 0.02 * sc.size))
        ru = r.radius / max(r.radius.mean(), 1e-300)
        return base * np.clip(ru, 0.5, 2.0) ** 0.5

    def camera(self, azimuth_offset=0.0, distance=None):
        o, sc = self.o, self.scene
        a = np.radians(o["azimuth"] + azimuth_offset)
        e = np.radians(o["elevation"])
        dist = distance or o["distance"]
        if dist is None:                           # frame the whole domain
            rb = sc.sphere[1] if sc.sphere is not None else 0.5 * float(np.linalg.norm(sc.L))
            dist = 1.08 * rb / (sc.size * np.sin(np.radians(o["view_angle"]) / 2))
        d = dist * sc.size
        f = sc.center + np.asarray(o["focus"], float) * sc.size
        pos = f + d * np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
        return [tuple(pos), tuple(f), (0, 0, 1)]

    def frame(self, k, path, size=None, azimuth_offset=0.0, distance=None, title=None, subtitle=None,
              overlays=True, scale=None):
        pv = _pv()
        o, sc, r = self.o, self.scene, self.scene.run
        k %= len(r)
        size = size or o["size"]
        pl = _plotter(size)
        vals, kind, label, unit = colour_values(sc, k, o["color"], self.disp)
        if kind == "continuous":
            cmap = o["cmap"] or _heat()
        else:
            cats = label if isinstance(label, list) else [label]
            pal = o["cmap"] if isinstance(o["cmap"], (list, tuple)) else (
                [STATE_COLORS[0], STATE_COLORS[1]] if o["color"] == "state" else [CYAN, ORANGE, "#9be15d", "#c084fc"])
            cmap = list(pal[:max(len(cats), 1)]) if len(cats) > 1 else [pal[0]]
        clim = self._clim if kind == "continuous" else (0, max(len(cmap) - 1, 1))
        state = r.state(k)
        pos = r.position(k).copy()
        rad = self._radius
        show = state != 2
        # stuck grains: keep them inside the domain so they rest on walls
        stuck = state == 1
        for a in range(3):
            if not r.periodic[a]:
                pos[stuck, a] = np.clip(pos[stuck, a], r.lo[a] + rad[stuck], r.hi[a] - rad[stuck])
        if o["cut"] == "quadrant":
            rel = pos - sc.center
            az = np.degrees(np.arctan2(rel[:, 1], rel[:, 0]))
            dz = (az - (o["azimuth"] + azimuth_offset) + 180) % 360 - 180
            show &= ~((np.abs(dz) < 45) & (rel[:, 2] > -0.15 * sc.size))
        elif o["cut"] == "half":
            rel = pos - sc.center
            a = np.radians(o["azimuth"] + azimuth_offset)
            show &= (rel[:, 0] * np.cos(a) + rel[:, 1] * np.sin(a)) < 0
        sv = vals if not (kind == "categorical" and len(cmap) == 1) else None
        _spheres(pl, pos[show], rad[show], sv[show] if sv is not None else None, cmap=cmap if sv is not None else None,
                 clim=clim, color=cmap[0] if sv is None else None, res=o["res"])
        # trails
        nt = int(o["trails"])
        if nt > 0 and k > 0:
            ks = list(range(max(0, k - nt), k + 1))
            alive = (state == 0) & show
            P = np.stack([r.position(j) for j in ks], 1)[alive]
            ok = np.abs(np.diff(P, axis=1)).max(axis=(1, 2)) < 0.4 * r.L.min()        # drop periodic wraps
            P = P[ok]
            m = P.shape[1]
            if kind == "continuous":
                _lines(pl, P, np.repeat(vals[alive][ok], m), cmap, clim, o["trail_opacity"], 2.5)
            else:
                _lines(pl, P, np.tile(np.linspace(0, 1, m), len(P)), ["#0d1626", "#9fd9ea"], (0, 1),
                       o["trail_opacity"], 1.6)
        # domain decorations
        box = o["box"] if o["box"] is not None else (sc.sphere is None and not sc.prims)
        h = r.L / 2
        c = sc.center
        if box:
            pl.add_mesh(_box(c - h, c + h).extract_all_edges(), color=SLATE, opacity=0.35, line_width=1.2)
        if sc.sphere is not None:
            pl.add_mesh(pv.Sphere(center=sc.sphere[0], radius=sc.sphere[1], theta_resolution=120, phi_resolution=80),
                        color="#9fb7d9", opacity=0.045, smooth_shading=True)
        for a, side in sc.absorbing:
            lo, hi = c - h * 1.08, c + h * 1.08
            t = 0.03 * r.L[a]
            if side == 0:
                lo[a], hi[a] = r.lo[a] - t, r.lo[a]
            else:
                lo[a], hi[a] = r.hi[a], r.hi[a] + t
            pl.add_mesh(_box(lo, hi), color="#8a96a8", ambient=0.2, diffuse=0.7, specular=0.8, specular_power=60)
        if o["floor"]:
            fl = pv.Plane(center=(c[0], c[1], r.lo[2] - 0.2 * r.L[2]), i_size=r.L[0] * 1.6, j_size=r.L[1] * 1.6,
                          i_resolution=48, j_resolution=48)
            pl.add_mesh(fl.extract_all_edges(), color="#2b3a55", opacity=0.55, line_width=1)
        if o["objects"] and sc.prims:
            _draw_objects(pl, sc)
        if o["sheath"] and sc.field:
            if self._sheath is None:
                self._sheath = _sheath_meshes(sc, o["sheath"])
            for i, m in enumerate(self._sheath):
                pl.add_mesh(m, color=CYAN, opacity=0.07 + 0.05 * i, smooth_shading=True, specular=0.4)
        _lights(pl, sc.center, 4 * sc.size, key=o["key_light"])
        pl.camera_position = self.camera(azimuth_offset, distance)
        pl.camera.view_angle = o["view_angle"]
        if o["ssao"]:
            sr = o["ssao_radius"] or 3.0 * float(np.median(rad))
            pl.enable_ssao(radius=sr, bias=0.03, kernel_size=128, blur=True)
        pl.enable_anti_aliasing("ssaa")
        pl.screenshot(path)
        pl.close()
        if not overlays:
            return path
        s = scale or min(1.0, size[0] / 1300)
        s = max(s, 0.55)
        t = title if title is not None else o["title"]
        st = subtitle if subtitle is not None else o["subtitle"]
        fill = {"time": sc.time_text(k), "n": r.n, "name": sc.name}
        legend = cbar = None
        if kind == "categorical" and o["legend"] and isinstance(label, list) and len(label) > 1:
            legend = list(zip(cmap, label))
        if kind == "continuous" and o["colorbar"]:
            u = f" {unit}" if unit else ""
            cbar = (cmap, _fmt(clim[0]), _fmt(clim[1]) + u, label)
        _overlay(path, t.format(**fill) if t else None, st.format(**fill) if st else None, legend=legend, cbar=cbar,
                 scale=s)
        return path


# ====================================================================== public API
def snapshot(run: Run, path, frame=-1, time=None, field_file=None, compress=True, **opts):
    """Render one snapshot to PNG. `time` (s for SI runs, ω₀⁻¹ otherwise) overrides `frame`."""
    R = Renderer(run, field_file, **opts)
    k = R.scene.frame_at(time) if time is not None else frame
    R.frame(k, path)
    if compress:
        compress_png(path)
    return path


def compare(run: Run, path, frames=(0, -1), titles=None, subtitles=None, field_file=None, compress=True, **opts):
    """Side-by-side panels of the same view at several frames, with one shared colour bar."""
    from PIL import Image
    opts.setdefault("size", (1100, 1000))
    R = Renderer(run, field_file, **{**opts, "colorbar": False, "legend": False})
    tmp = tempfile.mkdtemp(prefix="ppdyn_render_")
    panels = []
    for i, k in enumerate(frames):
        p = os.path.join(tmp, f"p{i}.png")
        R.frame(k, p, title=(titles or [None] * len(frames))[i] or R.scene.time_text(k),
                subtitle=(subtitles or [None] * len(frames))[i] or "", scale=1.0)
        panels.append(Image.open(p).convert("RGB"))
    W = sum(i.width for i in panels) + 4 * (len(panels) - 1)
    out = Image.new("RGB", (W, panels[0].height), BG_BOTTOM)
    x = 0
    for im in panels:
        out.paste(im, (x, 0)); x += im.width + 4
    out.save(path)
    vals, kind, label, unit = colour_values(R.scene, frames[-1], R.o["color"], R.disp)
    if kind == "continuous" and opts.get("colorbar", True):
        _overlay(path, cbar=(R.o["cmap"] or _heat(), _fmt(R._clim[0]), _fmt(R._clim[1]) + (f" {unit}" if unit else ""),
                             label), scale=1.0)
    if compress:
        compress_png(path)
    return path


def movie(run: Run, path, frames=None, every=1, orbit=0.0, fps=12, width=500, colors=112, size=(1000, 640),
          distance=None, hold=10, field_file=None, **opts):
    """Animated GIF (or MP4 when `path` ends in .mp4 and ffmpeg is available); the camera may orbit."""
    from PIL import Image
    opts.setdefault("title", None)
    R = Renderer(run, field_file, **{**opts, "colorbar": opts.get("colorbar", False), "legend": opts.get("legend", False)})
    ks = list(frames) if frames is not None else list(range(0, len(run), max(1, int(every))))
    tmp = tempfile.mkdtemp(prefix="ppdyn_movie_")
    files = []
    for i, k in enumerate(ks):
        f = os.path.join(tmp, f"f{i:04d}.png")
        R.frame(k, f, size=size, azimuth_offset=orbit * i / max(len(ks) - 1, 1), distance=distance,
                title="", subtitle="{time}", scale=1.0)
        files.append(f)
    files += [files[-1]] * hold
    if path.endswith(".mp4"):
        import imageio.v2 as imageio
        with imageio.get_writer(path, fps=fps, codec="libx264", quality=8) as w:
            for f in files:
                w.append_data(np.asarray(Image.open(f).convert("RGB")))
        return path
    ims = []
    for f in files:
        im = Image.open(f).convert("RGB")
        if width:
            im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
        ims.append(im.quantize(colors=colors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
    ims[0].save(path, save_all=True, append_images=ims[1:], duration=int(1000 / fps), loop=0, optimize=True)
    return path


# ====================================================================== presets for the bundled examples
PRESETS = {
    "dust_cloud": dict(
        look=dict(color="radius", radius=0.42, cut="quadrant", key_light=(0.6, -1.0, 1.4), ssao_radius=1.5, azimuth=-46.8, elevation=25.9, distance=1.854,
                  focus=(0, 0, -0.0368), view_angle=30, size=(2000, 1000),
                  title="Binary dust cloud in microgravity",
                  subtitle="20 000 grains, cut open to show the void.\nSmall grains line the void,\nlarger grains settle outside."),
        still="hero_dust_cloud.png"),
    "yukawa_crystal": dict(
        look=dict(color="displacement", radius=0.44, clim=(0, 4), key_light=(0.7, -1.2, 1.4), ssao_radius=1.2, azimuth=-38, elevation=31.8, distance=5.70,
                  focus=(0, 0, -0.074), view_angle=19.5, res=16, ssao=True),
        compare=dict(file="crystal_melting.png", frames=(0, -1),
                     titles=["bcc crystal, Γ = 400", "heated to Γ = 100"],
                     subtitles=["t = 0, every grain on its lattice site", "t = 150 ω₀⁻¹, the lattice has melted"]),
        movie=dict(file="crystal_melting.gif", every=2, orbit=60, distance=4.44, fps=10, hold=8)),
    "dust_void": dict(
        look=dict(color="outward", radius=0.62, clim=(0, 6), key_light=(0.5, -1.0, 1.6), ssao_radius=1.0, floor=True, box=False, azimuth=-62, elevation=28.8,
                  distance=1.569, focus=(0, 0.05, -0.05), view_angle=34, size=(1800, 1000),
                  title="Void opening in a dust monolayer",
                  subtitle="5000 grains, κ = 0.5, an ion-drag force pushes grains out of the centre\n"
                           "colour = how far each grain has been pushed outward"),
        still="dust_void.png", movie=dict(file="dust_void.gif", every=1, fps=10, colors=96)),
    "dust_fall": dict(
        look=dict(color="radius", radius=(0.10, 0.034), clim=(2.5, 6.0), key_light=(0.6, -1.1, 1.5), ssao_radius=0.8, trails=16, trail_opacity=0.6, res=12,
                  azimuth=-55, elevation=22.8, distance=3.524, focus=(0, 0, -0.0739), view_angle=30, size=(1400, 1000),
                  title="Charged grains falling onto an electrode",
                  subtitle="2000 melamine grains in 10 Pa argon, {time}\nlarge grains fall faster and land first"),
        still="dust_fall.png", still_time=0.044, movie=dict(file="dust_fall.gif", every=2, fps=14)),
    "debris_cubesat": dict(
        look=dict(color="state", radius=(0.03, 0.008), key_light=(0.9, -1.0, 1.1), ssao_radius=0.4, trails=14, res=10, sheath=(-0.35, -1.5), azimuth=-58,
                  elevation=25.6, distance=1.622, focus=(0, 0, -0.0946), view_angle=32, size=(1800, 1000),
                  title="Charged micro-debris on a 3U CubeSat",
                  subtitle="3000 aluminium grains at +2 V drift onto a −10 V spacecraft\n"
                           "{time}, translucent shells: sheath potential −1.5 V and −0.35 V"),
        still="debris_cubesat.png", still_time=6.0,
        movie=dict(file="debris_cubesat.gif", every=5, orbit=48, distance=2.099, fps=12)),
}


def preset_for(run: Run):
    name = Scene(run).name
    return PRESETS.get(name), name


def gallery(run: Run, outdir, which=("still", "compare", "movie"), field_file=None):
    """Render the README images of a bundled example (chosen by the run name)."""
    p, name = preset_for(run)
    if p is None:
        raise ValueError(f"no gallery preset for run '{name}'; presets: {', '.join(PRESETS)}")
    os.makedirs(outdir, exist_ok=True)
    look = p["look"]
    out = []
    if "still" in which and "still" in p:
        out.append(snapshot(run, os.path.join(outdir, p["still"]), time=p.get("still_time"), field_file=field_file, **look))
    if "compare" in which and "compare" in p:
        c = p["compare"]
        out.append(compare(run, os.path.join(outdir, c["file"]), c["frames"], c["titles"], c["subtitles"],
                           field_file=field_file, **look))
    if "movie" in which and "movie" in p:
        m = dict(p["movie"])
        f = m.pop("file")
        args = {**{k: v for k, v in look.items() if k not in ("size", "title", "subtitle")}, **m}
        out.append(movie(run, os.path.join(outdir, f), field_file=field_file, **args))
    return out
