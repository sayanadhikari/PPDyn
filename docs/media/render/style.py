"""Shared look for the PPDyn README renders (PyVista / VTK, offscreen)."""
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True
NAVY, CYAN, ORANGE, SLATE, CREAM = "#0f1b2d", "#19b5d8", "#f2662f", "#7d8ea6", "#f4efe6"
BG_BOTTOM, BG_TOP = "#070b12", "#1d2a40"


def plotter(size=(1600, 1000)):
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    pl.set_background(BG_BOTTOM, top=BG_TOP)
    return pl


def lights(pl, key=(1.0, -1.6, 1.2), center=(0, 0, 0), dist=200.0, key_int=1.25):
    c = np.asarray(center, float)
    k = np.asarray(key, float); k /= np.linalg.norm(k)
    fill = np.array([-k[0], k[1] * 0.4, -0.3]); fill /= np.linalg.norm(fill)
    rim = np.array([-k[0] * 0.3, -k[1], 0.9]); rim /= np.linalg.norm(rim)
    for d, i, col in ((k, key_int, "#fff4e6"), (fill, 0.45, "#a9c8ff"), (rim, 0.6, "#cfe3ff")):
        pl.add_light(pv.Light(position=tuple(c + dist * d), focal_point=tuple(c), intensity=i,
                              color=col, light_type="scene light"))
    pl.add_light(pv.Light(intensity=0.18, light_type="headlight"))


def grains(pl, pos, radius, scalars=None, cmap=None, clim=None, color=None, res=14,
           roughness=0.32, metallic=0.05, opacity=1.0, pbr=False):
    pts = pv.PolyData(np.ascontiguousarray(pos, float))
    r = np.broadcast_to(np.asarray(radius, float), (len(pos),)).copy()
    pts["r"] = r
    if scalars is not None:
        pts["s"] = np.asarray(scalars, float)
    sph = pv.Sphere(radius=1.0, theta_resolution=res, phi_resolution=max(8, res * 2 // 3))
    g = pts.glyph(geom=sph, scale="r", orient=False)
    if pbr:
        kw = dict(pbr=True, metallic=metallic, roughness=roughness, smooth_shading=True, opacity=opacity)
    else:
        kw = dict(smooth_shading=True, opacity=opacity, ambient=0.16, diffuse=0.9, specular=0.5,
                  specular_power=48)
    if scalars is not None:
        pl.add_mesh(g, scalars="s", cmap=cmap, clim=clim, show_scalar_bar=False, **kw)
    else:
        pl.add_mesh(g, color=color, **kw)
    return g


def finish(pl, path, ssao=True, ssao_radius=1.5, aa=True):
    if ssao:
        pl.enable_ssao(radius=ssao_radius, bias=0.03, kernel_size=128, blur=True)
    if aa:
        pl.enable_anti_aliasing("ssaa")
    pl.screenshot(path)
    pl.close()


# ---------------------------------------------------------------- 2D finishing (labels, panels)
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

HEAT = LinearSegmentedColormap.from_list("ppdyn_heat", ["#127fa3", CYAN, "#bfe9f2", CREAM, "#f7a073", ORANGE, "#b8361a"])
FONT_DIR = "/usr/share/fonts/opentype/inter"


def font(size, weight="Medium"):
    try:
        return ImageFont.truetype(f"{FONT_DIR}/Inter-{weight}.otf", size)
    except OSError:
        return ImageFont.load_default()


def label(path, title=None, sub=None, corner=(36, 30)):
    """Write a small title and subtitle onto a rendered image (in place)."""
    im = Image.open(path).convert("RGB")
    d = ImageDraw.Draw(im)
    x, y = corner
    if title:
        d.text((x, y), title, font=font(34, "SemiBold"), fill="#eef3fa")
        y += 46
    if sub:
        for line in sub.split("\n"):
            d.text((x, y), line, font=font(22, "Regular"), fill="#9fb0c8")
            y += 31
    im.save(path)


def colorbar(path, cmap, lo, hi, text, pos=None, length=360):
    """Slim horizontal colour bar drawn onto an image (in place)."""
    im = Image.open(path).convert("RGB")
    W, H = im.size
    x0, y0 = pos if pos else (W - length - 50, H - 70)
    bar = (np.array([cmap(t)[:3] for t in np.linspace(0, 1, length)]) * 255).astype(np.uint8)
    im.paste(Image.fromarray(np.repeat(bar[None], 10, 0)), (x0, y0))
    d = ImageDraw.Draw(im)
    f = font(19, "Regular")
    d.text((x0, y0 - 28), text, font=f, fill="#c9d4e3")
    d.text((x0, y0 + 16), lo, font=f, fill="#9fb0c8")
    d.text((x0 + length, y0 + 16), hi, font=f, fill="#9fb0c8", anchor="ra")
    im.save(path)


def hstack(paths, out, gap=0):
    ims = [Image.open(p).convert("RGB") for p in paths]
    W = sum(i.width for i in ims) + gap * (len(ims) - 1)
    c = Image.new("RGB", (W, max(i.height for i in ims)), BG_BOTTOM)
    x = 0
    for i in ims:
        c.paste(i, (x, 0)); x += i.width + gap
    c.save(out)


def gif(frames, out, fps=15, width=None, colors=192):
    ims = []
    for f in frames:
        im = Image.open(f).convert("RGB")
        if width:
            im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
        ims.append(im.quantize(colors=colors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
    ims[0].save(out, save_all=True, append_images=ims[1:], duration=int(1000 / fps), loop=0, optimize=True)


def compress_png(path, colors=256):
    """256-colour palette PNG with dithering: about 3x smaller, visually identical for these renders."""
    im = Image.open(path).convert("RGB")
    im.quantize(colors=colors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG).save(path, optimize=True)
