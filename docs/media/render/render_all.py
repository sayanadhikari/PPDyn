"""Rebuilds the README gallery with the PPDyn command-line tools only.

    pip install "ppdyn[render]"
    python docs/media/render/render_all.py RUNDIR            # init + run the five examples, then render

Every image and animation is written by ``ppdyn render <run>.h5 --gallery docs/media``;
the look of each example (colours, camera, cut-away, trails, sheath) is the preset
of the same name in ``ppdyn.post.render.PRESETS``. Takes about 8 min to run the
examples and 8 min to render on two cores. Pass --render-only to skip the runs.
"""
import os
import subprocess
import sys

MEDIA = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# denser output where the renders draw trajectories
EXAMPLES = {"dust_cloud": [], "yukawa_crystal": [], "dust_void": [],
            "dust_fall": ["--set", "output.every=5"], "debris_cubesat": ["--set", "output.every=10"]}


def sh(*cmd, cwd=None):
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def main(rundir, render_only=False):
    os.makedirs(rundir, exist_ok=True)
    for name, extra in EXAMPLES.items():
        d = os.path.join(rundir, name)
        if not render_only:
            if not os.path.isdir(d):
                sh("ppdyn", "init", name, d)
            sh("ppdyn", "run", "input.toml", "-q", *extra, cwd=d)
        sh("ppdyn", "render", os.path.join(d, "output", f"{name}.h5"), "--gallery", MEDIA)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0] if args else "gallery_runs", "--render-only" in sys.argv)
