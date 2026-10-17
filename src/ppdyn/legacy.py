"""Convert legacy PPDyn input.ini files (main, dust-fall, dust-void, dust-exp,
debris_dynamics branches) to the TOML input format.

Legacy runs used a domain [-L, L] per axis, so box lengths are doubled. The
conversion keeps the legacy physics where it was intentional and records every
judgement call as a comment in the generated file.
"""
from __future__ import annotations

import configparser
import math


def _num(s, default=None):
    if s is None:
        return default
    s = str(s).split(";")[0].split("#")[0].strip().strip("'\"")
    try:
        return float(s)
    except ValueError:
        return default


def _bool(s, default=False):
    if s is None:
        return default
    s = str(s).split(";")[0].split("#")[0].strip().lower()
    return s in ("true", "1", "yes", "on")


def _str(s, default=""):
    if s is None:
        return default
    return str(s).split(";")[0].split("#")[0].strip().strip("'\"")


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int,)) and not isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return repr(v) if (abs(v) < 1e-4 or abs(v) >= 1e6) and v != 0 else f"{v:.10g}"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    return '"' + str(v).replace('"', '\\"') + '"'


class _Out:
    def __init__(self):
        self.lines = []

    def section(self, name, comment=None):
        self.lines.append("")
        if comment:
            for c in comment.split("\n"):
                self.lines.append(f"# {c}")
        self.lines.append(f"[{name}]")

    def kv(self, k, v, comment=None):
        self.lines.append(f"{k} = {_fmt(v)}" + (f"   # {comment}" if comment else ""))

    def note(self, text):
        for c in text.split("\n"):
            self.lines.append(f"# {c}")

    def text(self):
        return "\n".join(self.lines).strip() + "\n"


def convert(ini_text: str, name="legacy") -> tuple[str, list[str]]:
    cp = configparser.ConfigParser(inline_comment_prefixes=(";", "#"), strict=False)
    cp.optionxform = str
    # some legacy files put keys before any section
    cp.read_string("[__top__]\n" + ini_text)
    g = lambda sec, key, d=None: cp.get(sec, key, fallback=d) if cp.has_section(sec) else d
    notes = []
    o = _Out()
    o.note(f"Converted from a legacy PPDyn input.ini by `ppdyn convert`.\nReview every line marked NOTE.")

    Lx, Ly, Lz = (_num(g("simbox", k), 10.0) for k in ("Lx", "Ly", "Lz"))
    N = int(_num(g("particles", "N"), 1000))
    btype = _str(g("boundary", "btype"), "periodic")
    geometry = _str(g("boundary", "geometry"), "box")
    kappa = _num(g("screening", "k"), 0.0)
    g0 = _num(g("gravity", "g_0"), 0.0)
    rc = _num(g("cutoff radius", "rc"), 0.0)
    tmax = _num(g("time", "tmax"), 10.0)
    dt = _num(g("time", "dt"), 0.01)
    every = int(_num(g("diagnostics", "dumpPeriod"), 10))
    vtk = _bool(g("diagnostics", "vtkData"))
    datadir = _str(g("directory", "dataDir"), "") or _str(g("directory", "path"), "data")
    is_void = cp.has_option("particles", "f0") if cp.has_section("particles") else False
    is_exp = cp.has_section("plasma")
    is_debris = cp.has_option("options", "PIC_data") if cp.has_section("options") else False

    o.section("run")
    o.kv("name", name)
    o.kv("seed", 99999999)

    if is_debris:
        o.section("units", "NOTE: the debris_dynamics branch normalised PINC SI fields with Qref/mref.\nRunning in SI keeps that physics explicit.")
        o.kv("system", "SI")
        notes.append("debris: units switched to SI; set box, dt and species values in SI before running")

    o.section("domain", "Legacy domain was [-L, L] per axis, so box = 2 L.")
    if geometry == "radial":
        o.kv("shape", "sphere")
        o.kv("radius", Lx)
    else:
        o.kv("box", [2 * Lx, 2 * Ly, 2 * Lz])

    o.section("boundary")
    if geometry == "radial":
        o.kv("sphere", "reflecting", "NOTE: legacy radial boundary mirrored particles through the origin")
        notes.append("dust-exp radial boundary mapped to a reflecting sphere")
    elif btype == "periodic":
        o.kv("x", "periodic"); o.kv("y", "periodic"); o.kv("z", "periodic")
    elif btype == "reflecting":
        o.kv("x", "reflecting"); o.kv("y", "reflecting")
        o.kv("z", "reflecting", "NOTE: legacy only reflected at the top; floor made reflecting")
    elif btype == "mixed":
        o.kv("x", "periodic"); o.kv("y", "periodic")
        o.kv("z", ["absorbing", "reflecting"], "dust-fall: absorbing ground, reflecting top")
        o.kv("on_absorb", "stick")
        o.kv("neutralize", True)

    o.section("time")
    o.kv("dt", dt)
    o.kv("t_end", tmax)

    # ------------------------------------------------------------ species
    o.lines.append("")
    o.lines.append("[[species]]")
    o.kv("name", "dust")
    o.kv("count", N)
    dist = _str(g("particles", "dist"), "False")
    if is_exp:
        mass = _num(g("particles", "mass"), 1.0)
        Q = _num(g("particles", "Q"), 1.0)
        a_vals = [float(x) for x in _str(g("particles", "a"), "1").split(",")]
        stds = [float(x) for x in _str(g("particles", "stdDev"), "0").split(",")]
        a_ref = sum(a_vals) / len(a_vals)
        o.kv("mass", mass, "mass at the reference radius (legacy m ~ a^3)")
        o.kv("reference_radius", a_ref, "legacy scaled Q, m with a / mean(a)")
        o.lines.append("[species.radius]")
        if dist == "binary":
            o.kv("distribution", "binary"); o.kv("values", a_vals)
        elif dist == "mean2":
            o.kv("distribution", "bimodal"); o.kv("values", a_vals); o.kv("stds", stds)
        elif dist == "uniform":
            o.kv("distribution", "uniform"); o.kv("low", a_vals[0]); o.kv("high", a_vals[-1])
        else:
            o.kv("distribution", "normal"); o.kv("mean", a_vals[0]); o.kv("std", stds[0])
        o.lines.append("[species.charge]")
        o.kv("value", Q); o.kv("exponent", 1.0)
        rmin, rmax = _num(g("simbox", "r_min"), 0.0), _num(g("simbox", "r_max"), Lx)
        o.lines.append("[species.placement]")
        o.kv("region", "shell"); o.kv("r_min", rmin); o.kv("r_max", rmax)
    elif _bool(dist) or dist in ("gaussian",):
        o.kv("mass", 1.0)
        o.lines.append("[species.radius]")
        o.note("legacy dist = True: m ~ N(mean, stdDev), Q = m^(2/3)")
        o.kv("distribution", "normal")
        o.kv("variable", "mass")
        o.kv("mean", _num(g("particles", "mean"), 1.0))
        o.kv("std", _num(g("particles", "stdDev"), 0.2))
        o.kv("min", 1e-6)
        o.lines.append("[species.charge]")
        o.kv("value", 1.0); o.kv("exponent", 2.0)
    elif dist == "uniform_radius":
        o.kv("density", _num(g("particles", "density"), 2000.0))
        o.lines.append("[species.radius]")
        o.kv("distribution", "uniform")
        o.kv("low", _num(g("particles", "r_min"), 1e-3)); o.kv("high", _num(g("particles", "r_max"), 1e-2))
        o.lines.append("[species.charge]")
        o.kv("value", 1.0, "NOTE: set the reference charge in coulomb (legacy Qref)")
    if btype == "mixed" and not is_exp:
        o.lines.append("[species.placement]")
        o.kv("region", "plane"); o.kv("axis", "z"); o.kv("position", Lz, "dust-fall: released from the top")
    temp = _num(g("particles", "Temp"), 0.0)
    o.lines.append("[species.velocity]")
    if is_void:
        o.kv("distribution", "uniform")
        o.kv("vmax", [_num(g("particles", k), 1.0) for k in ("Vxmax", "Vymax", "Vzmax")])
    elif _bool(g("particles", "maxwell_load"), True):
        o.kv("distribution", "maxwell")
        o.kv("temperature", temp ** 2, "legacy sampled v ~ N(0, Temp), i.e. T = Temp^2 for m = 1")
    else:
        o.kv("distribution", "uniform")
        o.kv("vmax", [_num(g("particles", k), 1.0) for k in ("Vxmax", "Vymax", "Vzmax")])

    # ------------------------------------------------------------ interaction
    o.section("interaction")
    o.kv("model", "yukawa" if kappa > 0 else "coulomb")
    o.kv("kappa", kappa)
    if kappa == 0:
        notes.append("kappa = 0 (pure Coulomb): uses a minimum-image O(N^2) sum like the legacy code")

    # ------------------------------------------------------------ forces
    forces = []
    if g0 != 0 and btype in ("reflecting", "mixed") and not is_debris:
        forces.append(("harmonic", [("omega2", [0.0, 0.0, g0]), ("origin", [0.0, 0.0, -Lz])],
                       "NOTE: legacy 'gravity' was acceleration -(z + Lz) g_0 (a spring anchored at the floor).\n"
                       "For uniform gravity use [forces.gravity] g = [0, 0, -g] instead."))
        notes.append("g_0 mapped to the legacy spring-like vertical force (forces.harmonic)")
    if g0 != 0 and is_debris:
        forces.append(("gravity", [("g", [0.0, 0.0, -g0])], None))
    if is_void:
        a = _num(g("screening", "a"), 1.0)
        KB = _num(g("constant", "KB"), 1.38e-23)
        Td = _num(g("dust", "Td"), 1.0)
        Gam = _num(g("screening", "Gamma"), 1.0)
        ff = a / (KB * Td * Gam)
        f0 = _num(g("particles", "f0"), 0.0)
        fflow = _num(g("particles", "f_flow"), 0.0)
        lam = _num(g("screening", "lambda_c"), 1.0)
        forces.append(("void", [("amplitude", f0 * ff), ("decay_length", lam), ("scale_with", "none")],
                       "dust-void: F = f0 a/(kB Td Gamma) exp(-r/lambda_c) r_hat.\n"
                       "NOTE: the legacy time dependence sin(w0 dt) was constant (it used dt, not t);\n"
                       "set modulation and omega here for a genuine time-modulated void."))
        forces.append(("constant", [("F", [0.0, fflow * ff, 0.0])], "dust-void ion flow force"))
    if is_exp:
        ff = _str(g("plasma", "Ffield"), "")
        if ff:
            forces.append(("radial_table", [("file", ff), ("r_key", "arr_0"), ("f_key", "arr_1"),
                                            ("factor", _num(g("plasma", "k_fel"), 1.0)), ("scale_with", "radius")],
                           "dust-exp electric force profile (force on a reference grain, ~ a)"))
            forces.append(("radial_table", [("file", ff), ("r_key", "arr_0"), ("f_key", "arr_2"),
                                            ("factor", _num(g("plasma", "k_fion"), 1.0)), ("scale_with", "radius2")],
                           "dust-exp ion drag profile (~ a^2)"))
    for nm, kvs, com in forces:
        if nm == "radial_table":
            o.lines.append("")
            if com:
                o.note(com)
            o.lines.append("[[forces.radial_table]]")
        else:
            o.section(f"forces.{nm}", com)
        for k, v in kvs:
            o.kv(k, v)

    # ------------------------------------------------------------ neutrals
    if is_void:
        nu = _num(g("neutral", "nu"), 0.0)
        Tn = _num(g("neutral", "Tn"), 0.0)
        Td = _num(g("dust", "Td"), 1.0)
        Gam = _num(g("screening", "Gamma"), 1.0)
        o.section("neutrals", "dust-void neutral drag with a proper Langevin thermal force.\n"
                              "NOTE: the legacy random kick was a constant (no random numbers); fixed here.")
        o.kv("drag", "constant"); o.kv("gamma", nu)
        o.kv("temperature", Tn / (Td * Gam), "Tn / (Td Gamma) in reduced units")
        notes.append("dust-void Langevin force fixed (legacy kick was deterministic)")
    if is_exp:
        kd = _num(g("neutrals", "k_drag"), 0.0)
        if kd > 0:
            o.section("neutrals", "NOTE: legacy drag force was -k_drag a^2 v. Choose an equivalent rate gamma = k_drag a^2 / m,\nor switch to SI units and drag = \"epstein\" with the gas pressure.")
            mass = _num(g("particles", "mass"), 1.0)
            a_vals = [float(x) for x in _str(g("particles", "a"), "1").split(",")]
            a_ref = sum(a_vals) / len(a_vals)
            o.kv("drag", "constant"); o.kv("gamma", kd * a_ref ** 2 / mass)
            notes.append("dust-exp drag mapped to a constant rate at the reference radius")

    # ------------------------------------------------------------ collisions
    if rc > 0 and btype != "periodic" and not is_void and not is_exp:
        o.section("collisions", "Legacy reflecting runs had elastic collisions for r < 2 rc.\nNOTE: give the species a radius (rc) to keep them.")
        o.kv("enabled", False)

    # ------------------------------------------------------------ field
    if is_debris and _bool(g("options", "PIC_data")):
        pic = _str(g("directory", "picDir"), "PIC_data")
        o.section("field", "debris_dynamics: PINC electric field")
        o.kv("file", pic.rstrip("/") + "/E.grid.h5")
        o.kv("format", "pinc")
        if _bool(g("options", "object_data")):
            o.note("NOTE: object.grid.h5 voxel masks are replaced by analytic objects, e.g.\n[[objects]]\ntype = \"cubesat\"")
            notes.append("debris object mask: describe the object with [[objects]] (e.g. type = \"cubesat\")")

    o.section("output")
    o.kv("dir", datadir or "output")
    o.kv("every", every)
    o.kv("vtk", vtk)
    return o.text(), notes
