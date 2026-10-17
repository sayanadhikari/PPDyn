"""Solid objects as signed-distance primitives, including a CubeSat preset."""
from __future__ import annotations

import numpy as np

from .kernels import OBJ_BOX, OBJ_CYL, OBJ_SPHERE

# CubeSat envelopes (m), Cal Poly CDS rev 14.1 (same table as the debris branch)
CUBESAT_SIZES = {"1U": (0.1000, 0.1000, 0.1135), "2U": (0.1000, 0.1000, 0.2270),
                 "3U": (0.1000, 0.1000, 0.3405), "6U": (0.1000, 0.2263, 0.3660)}
WING_THICKNESS = 0.002
ANTENNA_RADIUS = 0.002


def expand(objects: list[dict]) -> list[dict]:
    """Turn composite objects (cubesat) into primitives. Each primitive keeps
    the parent index in ``parent`` so hits can be attributed."""
    prims = []
    for k, ob in enumerate(objects):
        base = {"action": ob["action"], "restitution": ob["restitution"],
                "potential": ob["potential"], "parent": k}
        c = np.asarray(ob["center"], float)
        if ob["type"] != "cubesat":
            prims.append({**base, "type": ob["type"], "center": c, "size": np.asarray(ob["size"], float),
                          "radius": ob["radius"], "length": ob["length"], "axis": ob["axis"]})
            continue
        s = ob["scale"]
        wx, wy, wz = (v * s for v in CUBESAT_SIZES[ob["cubesat_size"]])
        t, L = WING_THICKNESS * s, ob["wing_length"] * s
        prims.append({**base, "type": "box", "center": c, "size": np.array([wx, wy, wz])})
        if ob["wings"] == "body_mounted":
            for sg in (-1, 1):
                prims.append({**base, "type": "box", "center": c + [sg * (wx / 2 + t / 2), 0, 0],
                              "size": np.array([t, wy, wz])})
        elif ob["wings"] == "deployed":
            for sg in (-1, 1):
                prims.append({**base, "type": "box", "center": c + [sg * (wx / 2 + L / 2), 0, 0],
                              "size": np.array([L, wy, t])})
        if ob["antenna"]:
            prims.append({**base, "type": "cylinder", "center": c + [0, 0, -wz / 2 - L / 2],
                          "radius": ANTENNA_RADIUS * s, "length": L, "axis": "z"})
    return prims


def kernel_arrays(prims: list[dict]):
    """Arrays consumed by kernels.apply_objects / sdf_union."""
    n = len(prims)
    kinds = np.zeros(n, np.int64)
    pars = np.zeros((max(n, 1), 7))
    action = np.zeros(n, np.int64)
    rest = np.ones(n)
    for k, p in enumerate(prims):
        pars[k, :3] = p["center"]
        if p["type"] == "box":
            kinds[k] = OBJ_BOX
            pars[k, 3:6] = np.asarray(p["size"]) / 2
        elif p["type"] == "sphere":
            kinds[k] = OBJ_SPHERE
            pars[k, 3] = p["radius"]
        else:
            kinds[k] = OBJ_CYL
            pars[k, 3] = p["radius"]
            pars[k, 4] = p["length"] / 2
            pars[k, 6] = "xyz".index(p["axis"])
        action[k] = 0 if p["action"] == "reflect" else 1
        rest[k] = p["restitution"]
    return kinds, pars[:n], action, rest


def characteristic_size(prims) -> float:
    s = []
    for p in prims:
        if p["type"] == "box":
            s.append(float(np.min(p["size"])))
        else:
            s.append(2 * float(p["radius"]))
    return min(s) if s else 1.0
