"""Command line interface: ``ppdyn <command> ...``"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

from . import __version__

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "examples")


def _examples():
    return sorted(d for d in os.listdir(EXAMPLES_DIR) if os.path.isdir(os.path.join(EXAMPLES_DIR, d)))


def _example_doc(name):
    with open(os.path.join(EXAMPLES_DIR, name, "input.toml")) as fh:
        first = fh.readline().strip()
    return first.lstrip("# ").strip()


def _parse_value(text):
    from .config import tomllib
    try:
        return tomllib.loads(f"v = {text}")["v"]
    except Exception:
        return text


def _apply_overrides(raw: dict, sets):
    for s in sets or []:
        if "=" not in s:
            raise SystemExit(f"--set expects key=value, got {s!r}")
        path, val = s.split("=", 1)
        keys = path.strip().split(".")
        node = raw
        for k in keys[:-1]:
            if isinstance(node, list):
                node = node[int(k)]
            else:
                if k == "species" and k not in node:
                    node[k] = [{}]
                node = node.setdefault(k, {})
        if isinstance(node, list):
            node = node[0]
        node[keys[-1]] = _parse_value(val.strip())
    return raw


def _load(path, sets=None):
    from . import config
    with open(path, "rb") as fh:
        raw = config.tomllib.load(fh)
    raw = _apply_overrides(raw, sets)
    cfg = config.complete(raw)
    cfg["_source"] = os.path.abspath(path)
    with open(path, encoding="utf-8") as fh:
        cfg["_text"] = fh.read() + ("".join(f"\n# --set {s}" for s in (sets or [])))
    return cfg


# ----------------------------------------------------------------------
def cmd_run(a):
    from .simulation import run
    cfg = _load(a.input, a.set)
    if a.threads:
        cfg["run"]["threads"] = a.threads
    if a.output:
        cfg["output"]["dir"] = a.output
    else:
        base = os.path.dirname(os.path.abspath(a.input))
        if not os.path.isabs(cfg["output"]["dir"]):
            cfg["output"]["dir"] = os.path.join(base, cfg["output"]["dir"])
    if a.restart:
        cfg["run"]["restart"] = os.path.abspath(a.restart)
    sim, path = run(cfg, progress=not a.quiet and cfg["run"]["progress"])
    print(f"output: {path}")
    if a.report:
        from .post.analysis import Run
        from .post.plots import report
        with Run(path) as r:
            idx = report(r, os.path.join(os.path.dirname(path), "report"))
        print(f"report: {idx}")
    return 0


def cmd_check(a):
    from .simulation import Simulation
    cfg = _load(a.input, a.set)
    sim = Simulation(cfg, quiet=True)
    print(sim.summary())
    print("input OK")
    return 0


def cmd_init(a):
    if not a.example:
        print("available examples:")
        for e in _examples():
            print(f"  {e:<16} {_example_doc(e)}")
        print("\nusage: ppdyn init <example> [directory]")
        return 0
    if a.example not in _examples():
        raise SystemExit(f"unknown example {a.example!r}; choose from {', '.join(_examples())}")
    dest = a.directory or a.example
    if os.path.exists(dest) and os.listdir(dest) and not a.force:
        raise SystemExit(f"{dest} exists and is not empty (use --force)")
    shutil.copytree(os.path.join(EXAMPLES_DIR, a.example), dest, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("output", "__pycache__"))
    print(f"created {dest}/  ->  cd {dest} && ppdyn run input.toml")
    return 0


def cmd_field(a):
    from . import fields, geometry
    cfg = _load(a.input, a.set)
    prims = geometry.expand(cfg["objects"])
    if not prims:
        raise SystemExit("ppdyn field needs at least one [[objects]] table")
    g = fields.solve(cfg, prims)
    base = os.path.dirname(os.path.abspath(a.input))
    out = a.output or os.path.join(base, cfg["field_solver"]["output"])
    g.save(out, {"debye_length": cfg["field_solver"]["debye_length"]})
    print(f"field written to {out}  (|E|max = {abs(g.E).max():.4g})")
    return 0


def cmd_convert(a):
    from .legacy import convert
    with open(a.legacy) as fh:
        text, notes = convert(fh.read(), name=os.path.splitext(os.path.basename(a.legacy))[0])
    out = a.output or os.path.splitext(a.legacy)[0] + ".toml"
    with open(out, "w") as fh:
        fh.write(text)
    print(f"wrote {out}")
    for n in notes:
        print(f"  note: {n}")
    return 0


def cmd_plot(a):
    from .post import plots
    from .post.analysis import Run
    with Run(a.file) as r:
        out = a.output or os.path.join(os.path.dirname(os.path.abspath(a.file)), "figures")
        os.makedirs(out, exist_ok=True)
        k = a.frame
        what = a.what
        jobs = {
            "history": lambda: plots.history(r, os.path.join(out, "history.png"), a.si),
            "snapshot": lambda: plots.snapshot(r, os.path.join(out, f"snapshot_{k}.png"), k, a.color, a.si),
            "rdf": lambda: plots.rdf_plot(r, os.path.join(out, "rdf.png"), k, a.si),
            "vdf": lambda: plots.vdf_plot(r, os.path.join(out, "vdf.png"), k, a.si),
            "msd": lambda: plots.msd_plot(r, os.path.join(out, "msd.png"), a.si),
            "profiles": lambda: plots.profiles_plot(r, os.path.join(out, "profiles.png"), k, a.si),
            "absorption": lambda: plots.falltime_plot(r, os.path.join(out, "absorption.png"), a.si),
        }
        if what == "report":
            print(plots.report(r, out, a.si))
            return 0
        names = list(jobs) if what == "all" else [what]
        for n in names:
            p = jobs[n]()
            if p:
                print(p)
    return 0


def cmd_movie(a):
    from .post import plots
    from .post.analysis import Run
    with Run(a.file) as r:
        out = a.output or os.path.splitext(a.file)[0] + (".mp4" if a.view else ".mp4")
        print(plots.movie(r, out, a.color, a.view, a.fps, a.frames, a.si))
    return 0


def cmd_render(a):
    from .post.analysis import Run
    try:
        from .post import render
    except ImportError as exc:  # pragma: no cover
        print(f"error: {exc}", file=sys.stderr)
        return 2
    with Run(a.file) as r:
        if a.gallery:
            for p in render.gallery(r, a.gallery, field_file=a.field):
                print(p)
            return 0
        opts = {}
        preset, name = render.preset_for(r)
        if preset and a.preset != "none":
            opts.update({k: v for k, v in preset["look"].items() if k not in ("title", "subtitle")})
            if a.titles:
                opts.update({k: preset["look"][k] for k in ("title", "subtitle") if k in preset["look"]})
        for key in ("color", "cut", "trails", "azimuth", "elevation", "distance", "view_angle"):
            v = getattr(a, key)
            if v is not None:
                opts[key] = v
        if a.radius is not None:
            opts["radius"] = a.radius
        if a.clim:
            opts["clim"] = tuple(a.clim)
        if a.sheath:
            opts["sheath"] = tuple(a.sheath)
        if a.title is not None:
            opts["title"] = a.title
        if a.subtitle is not None:
            opts["subtitle"] = a.subtitle.replace("\\n", "\n")
        if a.size:
            w, h = (int(v) for v in a.size.lower().split("x"))
            opts["size"] = (w, h)
        if a.movie:
            opts.pop("title", None); opts.pop("subtitle", None)
            still_size = opts.pop("size", None)
            size = still_size if a.size else (1000, 640)
            print(render.movie(r, a.movie, every=a.every, orbit=a.orbit, fps=a.fps, width=a.width, size=size,
                               field_file=a.field, **opts))
            return 0
        out = a.output or os.path.splitext(a.file)[0] + "_render.png"
        print(render.snapshot(r, out, frame=a.frame, time=a.time, field_file=a.field, **opts))
    return 0


def cmd_view(a):
    from .post.analysis import Run
    from .post.viewer import viewer_html
    with Run(a.file) as r:
        out = a.output or os.path.splitext(a.file)[0] + "_viewer.html"
        print(viewer_html(r, out, a.color, a.points, a.frames, a.si))
    return 0


def cmd_export(a):
    import numpy as np
    from . import io
    from .post.analysis import Run
    if a.format == "xdmf":
        print(io.write_xdmf(a.file))
        return 0
    with Run(a.file) as r:
        out = a.output or os.path.join(os.path.dirname(os.path.abspath(a.file)), a.format)
        os.makedirs(out, exist_ok=True)
        frames = range(len(r)) if a.frame is None else [a.frame]
        for k in frames:
            p = r.position(k)
            data = {"charge": r.charge(k), "radius": r.radius, "species": r.species.astype(float),
                    "state": r.state(k).astype(float)}
            if "velocity" in r.quantities:
                data["velocity"] = r.velocity(k)
            if a.format == "vtk":
                io.write_vtk(os.path.join(out, f"particles_{k:06d}.vtk"), p, data)
            else:
                cols = [p] + [v if v.ndim == 2 else v[:, None] for v in data.values()]
                hdr = "x,y,z," + ",".join(
                    (f"{n}x,{n}y,{n}z" if v.ndim == 2 else n) for n, v in data.items())
                np.savetxt(os.path.join(out, f"particles_{k:06d}.csv"), np.hstack(cols), delimiter=",",
                           header=hdr, comments="")
        print(out)
    return 0


def cmd_info(a):
    import numpy as np
    from .post.analysis import Run, coupling
    with Run(a.file) as r:
        h = r.history
        print(f"{a.file}: {r.n} particles, {len(r)} snapshots, t = {r.times[0]:.4g} .. {r.times[-1]:.4g}")
        print(f"  units: {r.scales.describe()}")
        print(f"  box: {np.round(r.lo, 4).tolist()} .. {np.round(r.hi, 4).tolist()}, periodic {r.periodic.tolist()}")
        print(f"  kappa = {r.kappa:.4g}, final T = {h['temperature'][-1]:.5g}, Gamma_eff = {coupling(r):.4g}")
        e = h["total"]
        print(f"  energy change (max-min)/|E0| = {(e.max() - e.min()) / max(abs(e[0]), 1e-300):.3g}"
              + ("" if h["n_stuck"][-1] or h["n_removed"][-1] else "  (conserved only without drag, thermostat or walls)"))
        if h["n_stuck"][-1] or h["n_removed"][-1]:
            print(f"  stuck {int(h['n_stuck'][-1])}, removed {int(h['n_removed'][-1])}, "
                  f"collected charge {h['collected_charge'][-1]:.5g}")
        if "final" in r.f and "wall_time" in r.f["final"].attrs:
            print(f"  wall time {r.f['final'].attrs['wall_time']:.2f} s")
    return 0


def cmd_watch(a):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import animation
    path = a.history if a.history.endswith(".csv") else os.path.join(a.history, "history.csv")
    fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

    def update(_):
        try:
            d = np.genfromtxt(path, delimiter=",", names=True)
        except Exception:
            return
        if d.size < 2:
            return
        for x in ax:
            x.clear()
        ax[0].plot(d["time"], d["kinetic"], label="kinetic")
        ax[0].plot(d["time"], d["potential"], label="potential")
        ax[0].plot(d["time"], d["total"], label="total")
        ax[0].legend()
        ax[1].plot(d["time"], d["temperature"])
        ax[1].set_ylabel("temperature")
        ax[1].set_xlabel("time")
        ax[0].set_title(f"step {int(d['step'][-1])}")

    _ani = animation.FuncAnimation(fig, update, interval=1000, cache_frame_data=False)
    plt.show()
    return 0


def cmd_docs(a):
    from .docgen import markdown
    text = markdown()
    if a.output:
        with open(a.output, "w") as fh:
            fh.write(text)
        print(a.output)
    else:
        print(text)
    return 0


def cmd_bench(a):
    import time
    import numba
    from .config import loads
    from .simulation import Simulation
    print(f"PPDyn {__version__}, numba {numba.__version__}, {numba.get_num_threads()} threads")
    for n in a.sizes:
        L = (4 / 3 * 3.141592653589793 * n) ** (1 / 3)
        cfg = loads(f"""
[run]
progress = false
[domain]
box = [{L}, {L}, {L}]
[[species]]
count = {n}
[species.velocity]
temperature = 0.01
[interaction]
kappa = {a.kappa}
method = "{a.method}"
[time]
dt = 0.01
steps = {a.steps}
""")
        sim = Simulation(cfg, quiet=True)
        sim.compute_forces()
        sim.advance()  # compile
        t0 = time.perf_counter()
        for _ in range(a.steps):
            sim.advance()
        dt = time.perf_counter() - t0
        print(f"  N = {n:>8}: {dt / a.steps * 1e3:9.3f} ms/step, {n * a.steps / dt / 1e6:7.2f} M particle-steps/s ({sim.method})")
    return 0


# ----------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(prog="ppdyn", description="PPDyn: parallel molecular dynamics for dusty plasmas")
    p.add_argument("--version", action="version", version=f"PPDyn {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("run", help="run a simulation")
    s.add_argument("input")
    s.add_argument("-t", "--threads", type=int, default=0)
    s.add_argument("-o", "--output", help="output directory (overrides output.dir)")
    s.add_argument("--set", action="append", metavar="KEY=VALUE", help="override an input value, e.g. --set time.dt=0.005")
    s.add_argument("--restart", help="checkpoint file to continue from")
    s.add_argument("--report", action="store_true", help="write figures and an HTML report when done")
    s.add_argument("-q", "--quiet", action="store_true")
    s.set_defaults(fn=cmd_run)

    s = sub.add_parser("check", help="validate an input file and print the run summary")
    s.add_argument("input")
    s.add_argument("--set", action="append", metavar="KEY=VALUE")
    s.set_defaults(fn=cmd_check)

    s = sub.add_parser("init", help="copy a ready-to-run example")
    s.add_argument("example", nargs="?")
    s.add_argument("directory", nargs="?")
    s.add_argument("--force", action="store_true")
    s.set_defaults(fn=cmd_init)

    s = sub.add_parser("field", help="solve the electrostatic field around [[objects]] and write a field file")
    s.add_argument("input")
    s.add_argument("-o", "--output")
    s.add_argument("--set", action="append", metavar="KEY=VALUE")
    s.set_defaults(fn=cmd_field)

    s = sub.add_parser("convert", help="convert a legacy input.ini to TOML")
    s.add_argument("legacy")
    s.add_argument("-o", "--output")
    s.set_defaults(fn=cmd_convert)

    s = sub.add_parser("plot", help="figures from an output file")
    s.add_argument("file")
    s.add_argument("what", nargs="?", default="all",
                   choices=["all", "report", "history", "snapshot", "rdf", "vdf", "msd", "profiles", "absorption"])
    s.add_argument("-o", "--output")
    s.add_argument("-k", "--frame", type=int, default=-1)
    s.add_argument("-c", "--color", default="speed", choices=["speed", "charge", "radius", "species", "z"])
    s.add_argument("--si", action="store_true", help="label axes in SI units (SI runs only)")
    s.set_defaults(fn=cmd_plot)

    s = sub.add_parser("movie", help="render an animation (mp4 with ffmpeg, otherwise gif)")
    s.add_argument("file")
    s.add_argument("-o", "--output")
    s.add_argument("-c", "--color", default="speed", choices=["speed", "charge", "radius", "species", "z"])
    s.add_argument("--view", default="3d", choices=["3d", "xy", "xz", "yz"])
    s.add_argument("--fps", type=int, default=20)
    s.add_argument("--frames", type=int, default=300)
    s.add_argument("--si", action="store_true")
    s.set_defaults(fn=cmd_movie)

    s = sub.add_parser("render", help="publication-quality 3D render (PNG) or animation (GIF/MP4) with PyVista",
                       description="Ray-traced-looking 3D renders of a run: shaded grains, trails, walls, objects "
                                   "and sheath isosurfaces. Runs of the bundled examples use their gallery look "
                                   "unless --preset none; every option below overrides it.")
    s.add_argument("file")
    s.add_argument("-o", "--output", help="PNG file (default <file>_render.png)")
    s.add_argument("--movie", metavar="OUT", help="write an animation instead (.gif, or .mp4 with ffmpeg)")
    s.add_argument("--gallery", metavar="DIR", help="write the README images of a bundled example into DIR")
    s.add_argument("-k", "--frame", type=int, default=-1)
    s.add_argument("-t", "--time", type=float, help="snapshot nearest to this time (s for SI runs, 1/omega0 otherwise)")
    s.add_argument("-c", "--color", choices=["speed", "radius", "species", "displacement", "outward", "charge", "state", "z", "none"])
    s.add_argument("--clim", type=float, nargs=2, metavar=("LO", "HI"))
    s.add_argument("--radius", type=float, help="display radius of the grains in a (default: automatic)")
    s.add_argument("--cut", choices=["quadrant", "half"], help="cut the particle cloud open towards the camera")
    s.add_argument("--trails", type=int, help="draw trajectories over the last N snapshots")
    s.add_argument("--sheath", type=float, nargs="+", metavar="V", help="potential isosurfaces (V) from the field file")
    s.add_argument("--field", help="field file (default: [field] file of the run's input)")
    s.add_argument("--azimuth", type=float); s.add_argument("--elevation", type=float)
    s.add_argument("--distance", type=float, help="camera distance in units of the largest box side")
    s.add_argument("--view-angle", dest="view_angle", type=float)
    s.add_argument("--size", help="image size WxH in pixels, e.g. 1600x1000")
    s.add_argument("--title"); s.add_argument("--subtitle", help="{time}, {n} and {name} are filled in")
    s.add_argument("--titles", action="store_true", help="use the preset's title and subtitle")
    s.add_argument("--preset", default="auto", choices=["auto", "none"])
    s.add_argument("--every", type=int, default=1, help="movie: use every N-th snapshot")
    s.add_argument("--orbit", type=float, default=0.0, help="movie: camera turns this many degrees")
    s.add_argument("--fps", type=int, default=12)
    s.add_argument("--width", type=int, default=500, help="movie: GIF width in pixels")
    s.set_defaults(fn=cmd_render)

    s = sub.add_parser("view", help="interactive 3D HTML viewer with a time slider")
    s.add_argument("file")
    s.add_argument("-o", "--output")
    s.add_argument("-c", "--color", default="speed", choices=["speed", "charge", "radius", "species", "z"])
    s.add_argument("--points", type=int, default=20000)
    s.add_argument("--frames", type=int, default=60)
    s.add_argument("--si", action="store_true")
    s.set_defaults(fn=cmd_view)

    s = sub.add_parser("export", help="export snapshots to VTK, CSV or (re)write the XDMF index")
    s.add_argument("file")
    s.add_argument("format", choices=["vtk", "csv", "xdmf"])
    s.add_argument("-o", "--output")
    s.add_argument("-k", "--frame", type=int)
    s.set_defaults(fn=cmd_export)

    s = sub.add_parser("info", help="summarise an output file")
    s.add_argument("file")
    s.set_defaults(fn=cmd_info)

    s = sub.add_parser("watch", help="live energy/temperature plot of a running simulation")
    s.add_argument("history", help="history.csv or the output directory")
    s.set_defaults(fn=cmd_watch)

    s = sub.add_parser("docs", help="print the input reference (Markdown)")
    s.add_argument("-o", "--output")
    s.set_defaults(fn=cmd_docs)

    s = sub.add_parser("bench", help="quick performance benchmark")
    s.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 100000])
    s.add_argument("--steps", type=int, default=20)
    s.add_argument("--kappa", type=float, default=1.0)
    s.add_argument("--method", default="auto", choices=["auto", "cells", "direct"])
    s.set_defaults(fn=cmd_bench)
    return p


def main(argv=None):
    from .config import ConfigError
    args = build_parser().parse_args(argv)
    try:
        return args.fn(args)
    except ConfigError as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
