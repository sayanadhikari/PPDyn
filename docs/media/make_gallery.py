"""Build docs/media/gallery.png from the five example runs.

Usage: python make_gallery.py <dir containing the example folders>
(each run with `ppdyn run <example>/input.toml`).
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ppdyn.post import Run, rdf  # noqa: E402
from ppdyn.post.plots import STYLE, C, INK  # noqa: E402

root = sys.argv[1] if len(sys.argv) > 1 else "."
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gallery.png")
path = lambda e: os.path.join(root, e, "output", f"{e}.h5")

with plt.rc_context(STYLE):
    fig = plt.figure(figsize=(16, 3.9))
    gs = fig.add_gridspec(1, 5, wspace=0.28)

    # 1 crystal: g(r) liquid vs crystal
    ax = fig.add_subplot(gs[0])
    with Run(path("yukawa_crystal")) as r:
        for k, lab, col in ((0, "bcc crystal", C[0]), (-1, "melted", C[1])):
            x, g = rdf(r.position(k), r.L, r.periodic, rmax=6, nbins=150)
            ax.plot(x, g, color=col, label=f"{lab}, t = {r.times[k]:.0f}")
    ax.set_ylim(0, 6)
    ax.set_xlabel("r / a"); ax.set_ylabel("g(r)")
    ax.legend(fontsize=8, loc="upper right", frameon=True, framealpha=0.95, edgecolor="none")
    ax.set_title("Melting a Yukawa crystal", fontsize=10)

    # 2 dust fall: absorption time vs radius
    ax = fig.add_subplot(gs[1])
    with Run(path("dust_fall")) as r:
        t = r.final["absorb_time"]
        hit = t >= 0
        sc = r.scales
        ax.scatter(r.radius[hit] * sc.length * 1e6, t[hit] * sc.time * 1e3, s=2, color=C[1])
    ax.set_xlabel("grain radius [µm]"); ax.set_ylabel("landing time [ms]")
    ax.set_title("Dust fall onto an electrode", fontsize=10)

    # 3 dust void: top view
    ax = fig.add_subplot(gs[2])
    with Run(path("dust_void")) as r:
        p = r.position(-1)
        ax.scatter(p[:, 0], p[:, 1], s=0.4, color=INK, linewidths=0)
        ax.set_xlim(r.lo[0], r.hi[0]); ax.set_ylim(r.lo[1], r.hi[1])
    ax.set_aspect("equal"); ax.set_xlabel("x / a"); ax.set_ylabel("y / a")
    ax.set_title("Void in a dust monolayer", fontsize=10)

    # 4 dust cloud: radial distribution per grain size
    ax = fig.add_subplot(gs[3])
    with Run(path("dust_cloud")) as r:
        R = np.linalg.norm(r.position(-1), axis=1) * r.scales.length * 1e3
        big = r.radius > np.median(r.radius)
        bins = np.linspace(R.min() - 0.2, R.max() + 0.2, 50)
        ax.hist(R[~big], bins, color=C[0], alpha=0.8, label="3.45 µm")
        ax.hist(R[big], bins, color=C[1], alpha=0.8, label="3.55 µm")
    ax.set_xlabel("distance from centre [mm]"); ax.set_ylabel("grains")
    ax.legend(fontsize=8)
    ax.set_title("Size segregation around a void", fontsize=10)

    # 5 debris on a CubeSat
    ax = fig.add_subplot(gs[4])
    with Run(path("debris_cubesat")) as r:
        p = r.position(-1) * r.scales.length * 100
        st = r.state(-1)
        ax.scatter(p[st == 0, 0], p[st == 0, 2], s=1, color=C[0], linewidths=0, label="drifting")
        ax.scatter(p[st == 1, 0], p[st == 1, 2], s=1, color=C[1], linewidths=0, label="deposited")
        lo, hi = r.lo * r.scales.length * 100, r.hi * r.scales.length * 100
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[2], hi[2])
    ax.set_aspect("equal"); ax.set_xlabel("x [cm]"); ax.set_ylabel("z [cm]")
    ax.legend(fontsize=8, markerscale=6, loc="upper left")
    ax.set_title("Charged debris on a 3U CubeSat", fontsize=10)

    fig.savefig(out, dpi=130, bbox_inches="tight")
print(out)
