"""Figures, movies and a one-page HTML report."""
from __future__ import annotations

import html
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import animation  # noqa: E402

from . import analysis  # noqa: E402
from .analysis import Run  # noqa: E402

INK = "#1f2933"
C = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c", "#0891b2"]
STYLE = {
    "figure.dpi": 110, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.alpha": 0.25, "axes.prop_cycle": matplotlib.cycler(color=C),
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "legend.frameon": False, "figure.constrained_layout.use": True,
}


class Units:
    """Axis labels and conversion factors for reduced or SI display."""

    def __init__(self, run: Run, si: bool):
        self.si = si and run.scales.si
        sc = run.scales
        if self.si:
            self.x, self.xl = sc.length * 1e3, "mm"
            self.t, self.tl = sc.time, "s"
            self.T, self.Tl = sc.temperature, "K"
            self.v, self.vl = sc.length * sc.omega * 1e3, "mm/s"
            self.e, self.el = sc.energy, "J"
            self.q, self.ql = sc.charge / 1.602176634e-19, "e"
        else:
            self.x, self.xl = 1.0, "a"
            self.t, self.tl = 1.0, r"$\omega_0^{-1}$"
            self.T, self.Tl = 1.0, r"$Q_0^2/(4\pi\epsilon_0 a)$"
            self.v, self.vl = 1.0, r"$a\omega_0$"
            self.e, self.el = 1.0, r"$Q_0^2/(4\pi\epsilon_0 a)$"
            self.q, self.ql = 1.0, r"$Q_0$"


def _save(fig, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def color_values(run: Run, k: int, color: str):
    if color == "speed":
        if "velocity" in run.quantities:
            return np.linalg.norm(run.velocity(k), axis=1), "speed"
        color = "species"
    if color == "charge":
        return run.charge(k), "charge"
    if color == "radius":
        return run.radius, "radius"
    if color == "z":
        return run.position(k)[:, 2], "z"
    return run.species.astype(float), "species"


# ----------------------------------------------------------------------
def history(run: Run, path, si=False):
    u = Units(run, si)
    h = run.history
    t = h["time"] * u.t
    absorbing = h["n_stuck"].max() > 0 or h["n_removed"].max() > 0
    n = 3 if absorbing else 2
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(n, 1, figsize=(7, 2.4 * n), sharex=True)
        ax[0].plot(t, h["kinetic"] * u.e, label="kinetic")
        ax[0].plot(t, h["potential"] * u.e, label="potential")
        ax[0].plot(t, h["total"] * u.e, label="total", color=INK, lw=1.2)
        ax[0].set_ylabel(f"energy [{u.el}]")
        ax[0].legend(ncol=3, loc="best")
        for c, lab in (("tx", "$T_x$"), ("ty", "$T_y$"), ("tz", "$T_z$")):
            ax[1].plot(t, h[c] * u.T, label=lab, lw=1)
        ax[1].plot(t, h["temperature"] * u.T, color=INK, lw=1.4, label="T")
        ax[1].set_ylabel(f"temperature [{u.Tl}]")
        ax[1].set_yscale("log" if np.all(h["temperature"][1:] > 0) else "linear")
        ax[1].legend(ncol=4)
        if absorbing:
            ax[2].plot(t, h["n_active"], label="active")
            ax[2].plot(t, h["n_stuck"], label="stuck")
            ax[2].plot(t, h["n_removed"], label="removed")
            ax[2].set_ylabel("particles")
            ax2 = ax[2].twinx()
            ax2.plot(t, h["collected_charge"] * u.q, color=C[3], ls="--", label="collected charge")
            ax2.set_ylabel(f"collected charge [{u.ql}]", color=C[3])
            ax2.grid(False)
            ax[2].legend(loc="center left")
        ax[-1].set_xlabel(f"time [{u.tl}]")
        return _save(fig, path)


def snapshot(run: Run, path, k=-1, color="speed", si=False, max_points=200_000):
    u = Units(run, si)
    p = run.position(k) * u.x
    c, clab = color_values(run, k, color)
    st = run.state(k)
    keep = st != 2
    idx = np.flatnonzero(keep)
    if len(idx) > max_points:
        idx = np.sort(np.random.default_rng(0).choice(idx, max_points, replace=False))
    p, c = p[idx], c[idx]
    size = np.clip(6.0 * (run.radius[idx] / run.radius.max()) ** 2, 0.5, 12) if run.radius.max() > 0 else 2.0
    if len(p) > 20000:
        size = np.asarray(size) * 0.3
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(13.5, 4.4))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.08)
        ax3 = fig.add_subplot(gs[0], projection="3d")
        sc = ax3.scatter(p[:, 0], p[:, 1], p[:, 2], c=c, s=size, cmap="viridis", depthshade=True, linewidths=0)
        ax3.set_xlabel(f"x [{u.xl}]"); ax3.set_ylabel(f"y [{u.xl}]")
        ax3.set_zlabel(f"z [{u.xl}]", labelpad=2)
        lo, hi = run.lo * u.x, run.hi * u.x
        ax3.set_xlim(lo[0], hi[0]); ax3.set_ylim(lo[1], hi[1]); ax3.set_zlim(lo[2], hi[2])
        span = np.maximum(hi - lo, 1e-12)
        ax3.set_box_aspect(tuple(np.maximum(span, 0.15 * span.max())))   # keep flat boxes readable
        for j, (a, b) in enumerate(((0, 1), (0, 2))):
            ax = fig.add_subplot(gs[j + 1])
            ax.scatter(p[:, a], p[:, b], c=c, s=size, cmap="viridis", linewidths=0)
            ax.set_xlabel(f"{'xyz'[a]} [{u.xl}]"); ax.set_ylabel(f"{'xyz'[b]} [{u.xl}]")
            ratio = span[a] / span[b]
            ax.set_aspect("equal" if 0.125 < ratio < 8 else "auto")
            ax.set_xlim(lo[a], hi[a]); ax.set_ylim(lo[b], hi[b])
        cb = fig.colorbar(sc, ax=fig.axes, shrink=0.8, pad=0.02)
        cb.set_label(clab)
        fig.suptitle(f"t = {run.times[k] * u.t:.4g} {u.tl}", color=INK)
        return _save(fig, path)


def rdf_plot(run: Run, path, k=-1, si=False):
    u = Units(run, si)
    sel = run.state(k) != 2
    if sel.sum() < 10:
        return None
    r, g = analysis.rdf(run.position(k), run.L, run.periodic, select=sel)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6, 3.6))
        ax.plot(r * u.x, g, color=C[0])
        ax.axhline(1, color=INK, lw=0.8, ls=":")
        ax.set_xlabel(f"r [{u.xl}]"); ax.set_ylabel("g(r)")
        ax.set_title(f"pair correlation, t = {run.times[k] * u.t:.4g} {u.tl}")
        return _save(fig, path)


def vdf_plot(run: Run, path, k=-1, si=False):
    u = Units(run, si)
    act = run.active(k)
    if act.sum() < 10:
        return None
    v = run.velocity(k)[act]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6, 3.6))
        for d in range(3):
            c, h, mw, T = analysis.velocity_distribution(v[:, d:d + 1], run.mass[act])
            ax.semilogy(c * u.v, h / u.v, label=f"$v_{'xyz'[d]}$", lw=1.2)
        c, h, mw, T = analysis.velocity_distribution(v, run.mass[act])
        ax.semilogy(c * u.v, mw / u.v, color=INK, ls="--", lw=1, label="Maxwellian")
        ax.set_ylim(max(mw.max() * 1e-5, 1e-12) / u.v, None)
        ax.set_xlabel(f"v [{u.vl}]"); ax.set_ylabel("f(v)")
        ax.legend()
        return _save(fig, path)


def msd_plot(run: Run, path, si=False):
    u = Units(run, si)
    t, m = analysis.msd(run)
    D = analysis.diffusion_coefficient(t, m)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6, 3.6))
        ax.loglog(t[1:] * u.t, m[1:] * u.x ** 2, color=C[0])
        ax.set_xlabel(f"t [{u.tl}]"); ax.set_ylabel(f"MSD [{u.xl}$^2$]")
        ax.set_title(f"mean squared displacement, D = {D * u.x ** 2 / u.t:.3g} {u.xl}$^2$/{u.tl}")
        return _save(fig, path)


def profiles_plot(run: Run, path, k=-1, si=False):
    u = Units(run, si)
    p = run.position(k)
    sel = run.state(k) != 2
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
        for d in range(3):
            x, h = analysis.density_profile(p, d, lo=run.lo[d], hi=run.hi[d], select=sel)
            ax[0].plot(x * u.x, h, label="xyz"[d])
        ax[0].set_xlabel(f"position [{u.xl}]"); ax[0].set_ylabel("count per bin"); ax[0].legend()
        c = 0.5 * (run.lo + run.hi)
        r, n = analysis.radial_profile(p, c, select=sel, rmax=0.5 * run.L.min())
        ax[1].plot(r * u.x, n / u.x ** 3, color=C[0])
        ax[1].set_xlabel(f"distance from centre [{u.xl}]"); ax[1].set_ylabel(f"density [{u.xl}$^{{-3}}$]")
        return _save(fig, path)


def falltime_plot(run: Run, path, si=False):
    u = Units(run, si)
    ft = analysis.fall_times(run)
    if len(ft) == 0:
        return None
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
        ax[0].hist(ft * u.t, bins=40, color=C[0])
        ax[0].set_xlabel(f"absorption time [{u.tl}]"); ax[0].set_ylabel("particles")
        hit = run.final["absorb_time"] >= 0
        if run.radius.max() > 0:
            ax[1].scatter(run.radius[hit] * u.x, ft * u.t, s=3, color=C[1])
            ax[1].set_xlabel(f"radius [{u.xl}]"); ax[1].set_ylabel(f"absorption time [{u.tl}]")
        else:
            ax[1].scatter(np.abs(run.charge0[hit]) * u.q, ft * u.t, s=3, color=C[1])
            ax[1].set_xlabel(f"initial |charge| [{u.ql}]"); ax[1].set_ylabel(f"absorption time [{u.tl}]")
        return _save(fig, path)


# ----------------------------------------------------------------------
def movie(run: Run, path, color="speed", view="3d", fps=20, max_frames=300, si=False, max_points=50_000):
    """Animate the run. view: 3d, xy, xz, yz. Writes mp4 (ffmpeg) or gif."""
    u = Units(run, si)
    frames = np.unique(np.linspace(0, len(run) - 1, min(len(run), max_frames)).astype(int))
    rng = np.random.default_rng(0)
    idx = np.arange(run.n) if run.n <= max_points else np.sort(rng.choice(run.n, max_points, replace=False))
    lo, hi = run.lo * u.x, run.hi * u.x
    vals = [color_values(run, k, color)[0][idx] for k in (frames[0], frames[-1])]
    vmin, vmax = min(v.min() for v in vals), max(v.max() for v in vals)
    s = 2.0 if len(idx) < 5000 else 0.6
    with plt.rc_context(STYLE):
        if view == "3d":
            fig = plt.figure(figsize=(6.4, 6))
            ax = fig.add_subplot(projection="3d")
            ax.set_box_aspect(tuple(np.maximum(hi - lo, 1e-12)))
            ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
            p = run.position(frames[0])[idx] * u.x
            sc = ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=s, c=vals[0], cmap="viridis", vmin=vmin, vmax=vmax, linewidths=0)
        else:
            a, b = ("xyz".index(view[0]), "xyz".index(view[1]))
            fig, ax = plt.subplots(figsize=(6.4, 6))
            ax.set_xlim(lo[a], hi[a]); ax.set_ylim(lo[b], hi[b]); ax.set_aspect("equal")
            ax.set_xlabel(f"{view[0]} [{u.xl}]"); ax.set_ylabel(f"{view[1]} [{u.xl}]")
            p = run.position(frames[0])[idx] * u.x
            sc = ax.scatter(p[:, a], p[:, b], s=s, c=vals[0], cmap="viridis", vmin=vmin, vmax=vmax, linewidths=0)
        title = ax.set_title("")

        def update(fi):
            k = frames[fi]
            p = run.position(k)[idx] * u.x
            c = color_values(run, k, color)[0][idx]
            if view == "3d":
                sc._offsets3d = (p[:, 0], p[:, 1], p[:, 2])
            else:
                sc.set_offsets(p[:, [a, b]])
            sc.set_array(c)
            title.set_text(f"t = {run.times[k] * u.t:.4g} {u.tl}")
            return sc, title

        ani = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        if path.endswith(".mp4") and animation.writers.is_available("ffmpeg"):
            ani.save(path, writer=animation.FFMpegWriter(fps=fps, bitrate=4000))
        else:
            if path.endswith(".mp4"):
                path = path[:-4] + ".gif"
            ani.save(path, writer=animation.PillowWriter(fps=fps))
        plt.close(fig)
    return path


# ----------------------------------------------------------------------
def report(run: Run, outdir, si=None, viewer=True):
    """All standard figures plus an index.html that ties them together.
    SI runs are labelled in SI units unless si=False."""
    si = run.scales.si if si is None else si
    os.makedirs(outdir, exist_ok=True)
    figs = []

    def add(title, fn, *a, **k):
        try:
            p = fn(*a, **k)
            if p:
                figs.append((title, os.path.basename(p)))
        except Exception as exc:  # keep the report going
            figs.append((title + f" (failed: {exc})", None))

    j = lambda name: os.path.join(outdir, name)
    add("Energy and temperature history", history, run, j("history.png"), si)
    add("Final snapshot", snapshot, run, j("snapshot.png"), -1, "speed", si)
    if len(run) > 1:
        add("Initial snapshot", snapshot, run, j("snapshot_initial.png"), 0, "speed", si)
    add("Pair correlation", rdf_plot, run, j("rdf.png"), -1, si)
    if "velocity" in run.quantities:
        # last snapshot that still has moving particles
        kv = max([k for k in range(len(run)) if run.active(k).sum() >= 10] or [0])
        add("Velocity distribution", vdf_plot, run, j("vdf.png"), kv, si)
    if len(run) > 3:
        add("Mean squared displacement", msd_plot, run, j("msd.png"), si)
    add("Density profiles", profiles_plot, run, j("profiles.png"), -1, si)
    add("Absorption statistics", falltime_plot, run, j("absorption.png"), si)
    viewer_file = None
    if viewer:
        try:
            from .viewer import viewer_html
            viewer_file = os.path.basename(viewer_html(run, j("viewer.html")))
        except ImportError:
            pass
    h = run.history
    rows = [("particles", run.n), ("snapshots", len(run)), ("final time", f"{run.times[-1]:.5g}"),
            ("kappa", f"{run.kappa:.4g}"), ("final temperature", f"{h['temperature'][-1]:.5g}"),
            ("coupling Gamma (final)", f"{analysis.coupling(run):.4g}"),
            ("energy drift (rel.)", f"{(h['total'].max() - h['total'].min()) / max(abs(h['total'][0]), 1e-300):.3g}"),
            ("units", run.scales.describe() if run.scales.si else "reduced")]
    body = ["<!doctype html><meta charset='utf-8'><title>PPDyn report</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 16px;color:#1f2933}"
            "img{max-width:100%;border:1px solid #e5e7eb;border-radius:6px}td{padding:2px 12px}"
            "h2{margin-top:1.6em;font-size:1.1em}</style>",
            f"<h1>PPDyn report: {html.escape(os.path.basename(run.path))}</h1><table>"]
    body += [f"<tr><td>{html.escape(str(a))}</td><td>{html.escape(str(b))}</td></tr>" for a, b in rows]
    body.append("</table>")
    if viewer_file:
        body.append(f"<p><a href='{viewer_file}'>Open the interactive 3D viewer</a></p>")
    for title, fn in figs:
        body.append(f"<h2>{html.escape(title)}</h2>")
        if fn:
            body.append(f"<img src='{fn}'>")
        elif "failed" not in title:
            body.pop()
    body.append(f"<h2>Input</h2><pre>{html.escape(run.input_text)}</pre>")
    with open(j("index.html"), "w") as fh:
        fh.write("\n".join(body))
    return j("index.html")
