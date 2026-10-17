"""Output: one HDF5 file per run, a live CSV history, ParaView XDMF/VTK, checkpoints."""
from __future__ import annotations

import json
import os

import h5py
import numpy as np

HISTORY_COLUMNS = ["step", "time", "kinetic", "kx", "ky", "kz", "potential", "total",
                   "temperature", "tx", "ty", "tz", "n_active", "n_stuck", "n_removed",
                   "collected_charge", "collisions"]

PER_PARTICLE = ("position", "velocity", "charge", "state", "force")


class Writer:
    def __init__(self, cfg, sim):
        out = cfg["output"]
        self.dir = out["dir"]
        os.makedirs(self.dir, exist_ok=True)
        self.name = cfg["run"]["name"]
        self.path = os.path.join(self.dir, f"{self.name}.h5")
        self.hist_path = os.path.join(self.dir, "history.csv")
        self.quantities = [q for q in out["quantities"] if q in PER_PARTICLE]
        unknown = set(out["quantities"]) - set(PER_PARTICLE)
        if unknown:
            raise ValueError(f"output.quantities: unknown entries {sorted(unknown)}; choose from {PER_PARTICLE}")
        if "position" not in self.quantities:
            self.quantities.insert(0, "position")
        self.dtype = np.float32 if out["precision"] == "float32" else np.float64
        comp = out["compression"]
        self.comp = None if comp == "none" else comp
        self.xdmf = out["xdmf"]
        self.vtk = out["vtk"]
        resume = bool(cfg["run"]["restart"]) and os.path.exists(self.path)
        self.f = h5py.File(self.path, "a" if resume else "w")
        if not resume:
            self._init_file(cfg, sim)
        self.hist = open(self.hist_path, "a" if resume else "w", buffering=1)
        if not resume:
            self.hist.write(",".join(HISTORY_COLUMNS) + "\n")
        self.n_vtk = len(self.f["snapshots/time"])

    def _init_file(self, cfg, sim):
        f = self.f
        f.attrs["code"] = "PPDyn"
        from . import __version__
        f.attrs["version"] = __version__
        f.attrs["input"] = cfg.get("_text", "")
        f.attrs["units"] = json.dumps(sim.scales.as_dict())
        f.attrs["unit_system_input"] = cfg["units"]["system"]
        f.attrs["N"] = sim.n
        f.attrs["dt"] = sim.dt
        f.attrs["box_lo"] = sim.lo
        f.attrs["box_hi"] = sim.hi
        f.attrs["periodic"] = sim.periodic
        f.attrs["kappa"] = sim.kappa
        f.attrs["species_names"] = json.dumps([s["name"] for s in cfg["species"]])
        p = f.create_group("particles")
        p["species"] = sim.species
        p["radius"] = sim.rad
        p["mass"] = sim.m
        p["charge0"] = sim.q.copy()
        s = f.create_group("snapshots")
        s.create_dataset("time", (0,), maxshape=(None,), dtype="f8", chunks=(256,))
        s.create_dataset("step", (0,), maxshape=(None,), dtype="i8", chunks=(256,))
        n = sim.n
        chunk_n = max(1, min(n, 1 << 16))
        for qn in self.quantities:
            if qn in ("charge",):
                shape, dt = (0, n), self.dtype
            elif qn == "state":
                shape, dt = (0, n), np.uint8
            else:
                shape, dt = (0, n, 3), self.dtype
            s.create_dataset(qn, shape, maxshape=(None,) + shape[1:], dtype=dt,
                             chunks=(1,) + (chunk_n,) + shape[2:], compression=self.comp)
        h = f.create_group("history")
        h.attrs["columns"] = HISTORY_COLUMNS
        h.create_dataset("data", (0, len(HISTORY_COLUMNS)), maxshape=(None, len(HISTORY_COLUMNS)),
                         dtype="f8", chunks=(512, len(HISTORY_COLUMNS)))

    # ------------------------------------------------------------------
    def snapshot(self, sim):
        s = self.f["snapshots"]
        k = s["time"].shape[0]
        for name in ("time", "step"):
            s[name].resize((k + 1,))
        s["time"][k] = sim.t
        s["step"][k] = sim.step
        data = {"position": sim.pos, "velocity": sim.vel, "charge": sim.q, "state": sim.state,
                "force": sim.force}
        for qn in self.quantities:
            ds = s[qn]
            ds.resize((k + 1,) + ds.shape[1:])
            ds[k] = data[qn].astype(ds.dtype, copy=False)
        if self.vtk:
            write_vtk(os.path.join(self.dir, "vtk", f"{self.name}_{k:06d}.vtk"), sim.pos,
                      {"velocity": sim.vel, "charge": sim.q, "radius": sim.rad,
                       "state": sim.state.astype(np.float64), "species": sim.species.astype(np.float64)})
        self.f.flush()

    def history(self, row: dict):
        vals = [row.get(c, np.nan) for c in HISTORY_COLUMNS]
        ds = self.f["history/data"]
        k = ds.shape[0]
        ds.resize((k + 1, ds.shape[1]))
        ds[k] = vals
        self.hist.write(",".join(f"{v:.10g}" for v in vals) + "\n")

    def finalize(self, sim):
        g = self.f.require_group("final")
        for name, arr in (("absorb_time", sim.absorb_t), ("absorb_charge", sim.absorb_q),
                          ("hit_object", sim.hit_obj), ("state", sim.state)):
            if name in g:
                del g[name]
            g[name] = arr
        g.attrs["time"] = sim.t
        g.attrs["step"] = sim.step
        g.attrs["wall_time"] = sim.wall_time
        self.f.flush()
        if self.xdmf:
            write_xdmf(self.path)

    def close(self):
        self.hist.close()
        self.f.close()


# ----------------------------------------------------------------------
def write_xdmf(h5path: str, xdmf_path: str | None = None):
    """Temporal collection of point clouds referencing the HDF5 directly."""
    xdmf_path = xdmf_path or os.path.splitext(h5path)[0] + ".xdmf"
    with h5py.File(h5path, "r") as f:
        s = f["snapshots"]
        times = s["time"][:]
        n = f.attrs["N"]
        qs = [k for k in ("velocity", "charge", "state", "force") if k in s]
        prec = 4 if s["position"].dtype == np.float32 else 8
    h5 = os.path.basename(h5path)
    lines = ['<?xml version="1.0" ?>', '<Xdmf Version="3.0">', "<Domain>",
             '<Grid Name="particles" GridType="Collection" CollectionType="Temporal">']

    def slab(name, k, comps):
        dims = f"{n} {comps}" if comps > 1 else f"{n}"
        full = f"{len(times)} {n} {comps}" if comps > 1 else f"{len(times)} {n}"
        start = f"{k} 0 0" if comps > 1 else f"{k} 0"
        stride = "1 1 1" if comps > 1 else "1 1"
        count = f"1 {n} {comps}" if comps > 1 else f"1 {n}"
        r = len(start.split())
        p = 1 if name == "state" else prec
        nt = "UChar" if name == "state" else "Float"
        return (f'<DataItem ItemType="HyperSlab" Dimensions="{dims}">'
                f'<DataItem Dimensions="3 {r}" Format="XML">{start} {stride} {count}</DataItem>'
                f'<DataItem Dimensions="{full}" NumberType="{nt}" Precision="{p}" Format="HDF">'
                f'{h5}:/snapshots/{name}</DataItem></DataItem>')

    for k, t in enumerate(times):
        lines.append(f'<Grid Name="t{k}" GridType="Uniform"><Time Value="{t:.10g}"/>')
        lines.append(f'<Topology TopologyType="Polyvertex" NumberOfElements="{n}"/>')
        lines.append(f'<Geometry GeometryType="XYZ">{slab("position", k, 3)}</Geometry>')
        for q in qs:
            comps = 1 if q in ("charge", "state") else 3
            at = "Scalar" if comps == 1 else "Vector"
            lines.append(f'<Attribute Name="{q}" AttributeType="{at}" Center="Node">{slab(q, k, comps)}</Attribute>')
        for q, path in (("radius", "particles/radius"), ("species", "particles/species")):
            lines.append(f'<Attribute Name="{q}" AttributeType="Scalar" Center="Node">'
                         f'<DataItem Dimensions="{n}" NumberType="Float" Precision="8" Format="HDF">{h5}:/{path}</DataItem></Attribute>')
        lines.append("</Grid>")
    lines += ["</Grid>", "</Domain>", "</Xdmf>"]
    with open(xdmf_path, "w") as fh:
        fh.write("\n".join(lines))
    return xdmf_path


def write_vtk(path: str, pos, point_data: dict):
    """Legacy binary VTK polydata (vertices) without external dependencies."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    n = len(pos)
    with open(path, "wb") as fh:
        fh.write(b"# vtk DataFile Version 3.0\nPPDyn particles\nBINARY\nDATASET POLYDATA\n")
        fh.write(f"POINTS {n} double\n".encode())
        fh.write(np.asarray(pos, ">f8").tobytes())
        fh.write(f"\nVERTICES {n} {2 * n}\n".encode())
        fh.write(np.column_stack([np.ones(n, ">i4"), np.arange(n, dtype=">i4")]).astype(">i4").tobytes())
        fh.write(f"\nPOINT_DATA {n}\n".encode())
        for name, arr in point_data.items():
            arr = np.asarray(arr)
            if arr.ndim == 2:
                fh.write(f"VECTORS {name} double\n".encode())
            else:
                fh.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n".encode())
            fh.write(arr.astype(">f8").tobytes())
            fh.write(b"\n")


# ----------------------------------------------------------------------
CHECKPOINT_FIELDS = ("pos", "vel", "q", "m", "rad", "rr", "state", "species", "absorb_t", "absorb_q",
                     "hit_obj", "gamma", "Tn")


def save_checkpoint(path, sim):
    tmp = path + ".tmp"
    with h5py.File(tmp, "w") as f:
        for k in CHECKPOINT_FIELDS:
            f[k] = getattr(sim, k)
        f.attrs["t"] = sim.t
        f.attrs["step"] = sim.step
        f.attrs["collected_charge"] = sim.collected
        f.attrs["rng"] = json.dumps(sim.rng.bit_generator.state)
        f.attrs["scales"] = json.dumps(sim.scales.as_dict())
    os.replace(tmp, path)


def load_checkpoint(path):
    with h5py.File(path, "r") as f:
        data = {k: f[k][:] for k in CHECKPOINT_FIELDS if k in f}
        data["t"] = float(f.attrs["t"])
        data["step"] = int(f.attrs["step"])
        data["collected_charge"] = float(f.attrs["collected_charge"])
        data["rng"] = json.loads(f.attrs["rng"])
        data["scales"] = json.loads(f.attrs["scales"])
    return data
