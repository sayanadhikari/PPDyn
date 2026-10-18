"""Generates the PPDyn logo (SVG + PNG): a hexagonal plasma crystal of charged
grains with Debye-screening halos, one grain on a curved trajectory, and the
wordmark. Text is converted to outlines, so the SVGs need no fonts.

    pip install cairosvg fonttools
    python make_logo.py
"""
import math
import os

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY, CYAN, ORANGE, SLATE = "#0f1b2d", "#19b5d8", "#f2662f", "#7d8ea6"
FONT_DIR = "/usr/share/fonts/opentype/inter"
_fonts = {}


# ---------------------------------------------------------------- text as outlines
def _font(weight):
    if weight not in _fonts:
        _fonts[weight] = TTFont(os.path.join(FONT_DIR, f"Inter-{weight}.otf"))
    return _fonts[weight]


def text_paths(runs, x, y, size, weight="ExtraBold", tracking=0.0):
    """runs = [(text, colour), ...]; returns (svg, advance width)."""
    f = _font(weight)
    gs, cmap, hmtx = f.getGlyphSet(), f.getBestCmap(), f["hmtx"]
    s = size / f["head"].unitsPerEm
    out, pen_x = [], x
    for text, col in runs:
        pen = SVGPathPen(gs)
        for ch in text:
            g = cmap[ord(ch)]
            gs[g].draw(TransformPen(pen, (s, 0, 0, -s, pen_x, y)))
            pen_x += hmtx[g][0] * s + tracking * size
        out.append(f'<path d="{pen.getCommands()}" fill="{col}"/>')
    return "".join(out), pen_x - x - tracking * size


# ---------------------------------------------------------------- emblem
def emblem(cx, cy, R, dark):
    ink = "#e8eef6" if dark else NAVY
    o = []
    uid = "d" if dark else "l"
    o.append(f'<defs>'
             f'<radialGradient id="halo{uid}"><stop offset="0" stop-color="{CYAN}" stop-opacity="0.34"/>'
             f'<stop offset="0.45" stop-color="{CYAN}" stop-opacity="0.12"/>'
             f'<stop offset="1" stop-color="{CYAN}" stop-opacity="0"/></radialGradient>'
             f'<radialGradient id="glow{uid}"><stop offset="0" stop-color="{ORANGE}" stop-opacity="0.45"/>'
             f'<stop offset="1" stop-color="{ORANGE}" stop-opacity="0"/></radialGradient>'
             f'<linearGradient id="trail{uid}" x1="0" y1="0" x2="1" y2="0">'
             f'<stop offset="0" stop-color="{ORANGE}" stop-opacity="0"/>'
             f'<stop offset="1" stop-color="{ORANGE}" stop-opacity="1"/></linearGradient>'
             f'<clipPath id="hex{uid}"><path d="{hexagon(cx, cy, R * 0.98, round_r=R * 0.12)}"/></clipPath></defs>')
    # hexagonal frame with a faint fill
    o.append(f'<path d="{hexagon(cx, cy, R, round_r=R * 0.12)}" fill="{CYAN}" fill-opacity="{0.08 if dark else 0.05}" '
             f'stroke="{ink}" stroke-width="{R * 0.045:.2f}" stroke-linejoin="round"/>')
    a = R * 0.42                                   # lattice constant
    sites = [(cx, cy)] + [(cx + a * math.cos(math.radians(60 * k + 30)), cy + a * math.sin(math.radians(60 * k + 30)))
                          for k in range(6)]
    g = [f'<g clip-path="url(#hex{uid})">']
    # second-shell grains, partly outside: the crystal continues beyond the frame
    for k in range(6):
        ang = math.radians(60 * k)
        x, y = cx + a * math.sqrt(3) * math.cos(ang), cy + a * math.sqrt(3) * math.sin(ang)
        g.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{a * 0.62:.2f}" fill="url(#halo{uid})"/>')
        g.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{a * 0.13:.2f}" fill="{SLATE}" fill-opacity="0.55"/>')
    # Debye halos and bonds of the first shell
    for x, y in sites:
        g.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{a * 0.62:.2f}" fill="url(#halo{uid})"/>')
    for k in range(6):
        x1, y1 = sites[1 + k]; x2, y2 = sites[1 + (k + 1) % 6]
        for (p, q) in (((cx, cy), (x1, y1)), ((x1, y1), (x2, y2))):
            g.append(f'<line x1="{p[0]:.2f}" y1="{p[1]:.2f}" x2="{q[0]:.2f}" y2="{q[1]:.2f}" stroke="{ink}" '
                     f'stroke-opacity="0.22" stroke-width="{R * 0.012:.2f}"/>')
    g.append('</g>')
    o += g
    # one grain leaves its site: curved trajectory ending outside the lattice, upper right
    x0, y0 = sites[4]                                # lower-left site (empty now: a vacancy)
    o.append(f'<circle cx="{x0:.2f}" cy="{y0:.2f}" r="{a * 0.17:.2f}" fill="none" stroke="{ink}" stroke-opacity="0.45" '
             f'stroke-width="{R * 0.012:.2f}" stroke-dasharray="{R * 0.03:.1f} {R * 0.03:.1f}"/>')
    xe, ye = cx + R * 0.86, cy - R * 0.70
    c1 = (x0 + R * 0.05, y0 - R * 0.62)
    c2 = (cx + R * 0.30, cy - R * 0.78)

    def bez(t):
        return ((1 - t) ** 3 * x0 + 3 * (1 - t) ** 2 * t * c1[0] + 3 * (1 - t) * t * t * c2[0] + t ** 3 * xe,
                (1 - t) ** 3 * y0 + 3 * (1 - t) ** 2 * t * c1[1] + 3 * (1 - t) * t * t * c2[1] + t ** 3 * ye)
    n = 60
    for i in range(n):                              # tapered, colour-graded trail
        t0, t1 = i / n, (i + 1) / n
        (ax, ay), (bx, by) = bez(t0), bez(t1)
        col = mix(CYAN, ORANGE, t1 ** 0.8)
        o.append(f'<line x1="{ax:.2f}" y1="{ay:.2f}" x2="{bx:.2f}" y2="{by:.2f}" stroke="{col}" '
                 f'stroke-width="{R * (0.012 + 0.030 * t1):.2f}" stroke-linecap="round" stroke-opacity="{0.35 + 0.65 * t1:.2f}"/>')
    # grains of the first shell (bright, on top of the bonds)
    for i, (x, y) in enumerate(sites):
        if i == 4:
            continue
        o.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{a * 0.20:.2f}" fill="{CYAN}"/>')
        o.append(f'<circle cx="{x - a * 0.06:.2f}" cy="{y - a * 0.06:.2f}" r="{a * 0.07:.2f}" fill="#ffffff" fill-opacity="0.45"/>')
    # the moving grain, with a glow
    o.append(f'<circle cx="{xe:.2f}" cy="{ye:.2f}" r="{a * 0.55:.2f}" fill="url(#glow{uid})"/>')
    o.append(f'<circle cx="{xe:.2f}" cy="{ye:.2f}" r="{a * 0.23:.2f}" fill="{ORANGE}"/>')
    o.append(f'<circle cx="{xe - a * 0.07:.2f}" cy="{ye - a * 0.07:.2f}" r="{a * 0.08:.2f}" fill="#ffffff" fill-opacity="0.5"/>')
    return "".join(o)


def mix(a, b, t):
    ca = [int(a[i:i + 2], 16) for i in (1, 3, 5)]
    cb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(ca, cb))


def hexagon(cx, cy, R, round_r=0.0):
    """Pointy-top hexagon path with rounded corners."""
    pts = [(cx + R * math.cos(math.radians(60 * k - 90)), cy + R * math.sin(math.radians(60 * k - 90))) for k in range(6)]
    if round_r <= 0:
        return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
    d = []
    for k in range(6):
        p0, p1, p2 = pts[k - 1], pts[k], pts[(k + 1) % 6]
        u = [(p0[i] - p1[i]) for i in range(2)]; v = [(p2[i] - p1[i]) for i in range(2)]
        lu, lv = math.hypot(*u), math.hypot(*v)
        a = (p1[0] + u[0] / lu * round_r, p1[1] + u[1] / lu * round_r)
        b = (p1[0] + v[0] / lv * round_r, p1[1] + v[1] / lv * round_r)
        d.append(("M" if k == 0 else "L") + f"{a[0]:.2f},{a[1]:.2f} Q{p1[0]:.2f},{p1[1]:.2f} {b[0]:.2f},{b[1]:.2f}")
    return " ".join(d) + " Z"


# ---------------------------------------------------------------- assembly
def svg(w, h, body, bg=None):
    rect = f'<rect width="{w}" height="{h}" rx="{min(w, h) * 0.08:.0f}" fill="{bg}"/>' if bg else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
            f'{rect}{body}</svg>')


def wordmark(x, y, size, dark):
    ink = "#f4f7fb" if dark else NAVY
    sub = "#9fb0c6" if dark else "#51627a"
    big, w = text_paths([("PP", ink), ("Dyn", ORANGE)], x, y, size, "ExtraBold", tracking=-0.01)
    tag, tw = text_paths([("MOLECULAR DYNAMICS FOR DUSTY PLASMAS AND CHARGED DEBRIS", sub)],
                         x + size * 0.03, y + size * 0.40, size * 0.118, "Medium", tracking=0.16)
    return big + tag, max(w, tw + size * 0.03)


def main():
    files = {}
    # light: for white pages; dark: on a navy card; on_dark: transparent, for dark-mode pages
    for dark, name, bg in ((False, "ppdyn_logo", None), (True, "ppdyn_logo_dark", NAVY), (True, "ppdyn_logo_on_dark", None)):
        wm, w = wordmark(330, 212, 168, dark)
        W = int(330 + w + 60)
        files[name + ".svg"] = (svg(W, 360, emblem(170, 180, 140, dark) + wm, bg), 1400)
    files["ppdyn_icon.svg"] = (svg(512, 512, emblem(256, 262, 205, True), NAVY), 512)
    files["ppdyn_icon_light.svg"] = (svg(512, 512, emblem(256, 262, 205, False)), 512)
    import cairosvg
    for name, (text, width) in files.items():
        p = os.path.join(HERE, name)
        with open(p, "w") as fh:
            fh.write(text)
        cairosvg.svg2png(bytestring=text.encode(), write_to=p.replace(".svg", ".png"), output_width=width)
        print("wrote", name)


if __name__ == "__main__":
    main()
