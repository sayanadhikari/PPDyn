"""Simulation setup and the BAOAB time loop."""
from __future__ import annotations

import copy
import math
import os
import shutil
import sys
import time

import numpy as np

from . import fields, geometry, io, kernels, particles
from .schema import SCHEMA
from .units import KB, Scales, convert_tree

SCALE_CODES = {"none": 0, "charge": 1, "mass": 2, "radius": 3, "radius2": 4}
BC_CODES = {"periodic": kernels.BC_PERIODIC, "reflecting": kernels.BC_REFLECT,
            "absorbing": kernels.BC_ABSORB, "open": kernels.BC_OPEN}


def auto_cutoff(kappa: float, tol: float) -> float:
    """Distance where the Yukawa force falls to tol times its value at r = a."""
    if kappa <= 0:
        return math.inf
    f1 = (1 + kappa) * math.exp(-kappa)
    lo, hi = 1.0, 1.0
    while (1 + kappa * hi) * math.exp(-kappa * hi) / hi ** 2 > tol * f1:
        hi *= 2
        if hi > 1e8:
            return math.inf
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if (1 + kappa * mid) * math.exp(-kappa * mid) / mid ** 2 > tol * f1:
            lo = mid
        else:
            hi = mid
    return hi


class Simulation:
    def __init__(self, cfg: dict, log=print, quiet=False):
        self.cfg_input = cfg
        self.log = (lambda *a, **k: None) if quiet else log
        self.base = os.path.dirname(cfg.get("_source", "")) if cfg.get("_source", "<string>") != "<string>" else "."
        cfg = copy.deepcopy(cfg)
        self.si = cfg["units"]["system"] == "SI"
        self.rng = np.random.default_rng(cfg["run"]["seed"])
        threads = cfg["run"]["threads"]
        if threads > 0:
            import numba
            numba.set_num_threads(min(threads, numba.config.NUMBA_NUM_THREADS))

        # -------------------------------------------------- objects (input units)
        self.prims = geometry.expand(cfg["objects"])
        kinds, pars, _, _ = geometry.kernel_arrays(self.prims)
        inside_in = (lambda p: kernels.inside_objects(np.ascontiguousarray(p, float), kinds, pars)) if self.prims else None

        # -------------------------------------------------- species (input units)
        R, M, Q, REL, P, V, S, G = [], [], [], [], [], [], [], []
        for k, sp in enumerate(cfg["species"]):
            r, m, q, rel = particles.species_properties(sp, self.rng, self.si)
            p = particles.place(sp, sp["count"], cfg["domain"], self.rng, inside_in, self.base)
            v = particles.velocities(sp, m, self.rng, KB if self.si else 1.0)
            R.append(r); M.append(m); Q.append(q); REL.append(rel); P.append(p); V.append(v)
            S.append(np.full(sp["count"], k, np.int32))
            G.append(self._drag_rate_input(cfg, r, m))
        cat = lambda xs, shape=(0,): np.concatenate(xs) if xs else np.zeros(shape)
        rad, m, q, rel = cat(R), cat(M), cat(Q), cat(REL)
        pos, vel, gam_in = cat(P, (0, 3)), cat(V, (0, 3)), cat(G)
        self.species = cat(S).astype(np.int32)
        self.n = len(m)
        if self.n == 0:
            raise ValueError("no particles: every species has count = 0")

        # -------------------------------------------------- scales
        self.scales = self._make_scales(cfg, m, q)
        sc = self.scales
        for name, table in SCHEMA.items():
            if name in ("species", "objects"):
                continue
            if table.array:
                for it in cfg[name]:
                    convert_tree(it, table, sc)
            else:
                convert_tree(cfg[name], table, sc)
        self.cfg = cfg
        fL = sc.factor("L")
        self.pos = np.ascontiguousarray(pos / fL)
        self.vel = np.ascontiguousarray(vel / sc.factor("V"))
        self.q = np.ascontiguousarray(q / sc.factor("Q"))
        self.m = np.ascontiguousarray(m / sc.factor("M"))
        self.rad = np.ascontiguousarray(rad / fL)
        self.rr = np.ascontiguousarray(rel, dtype=float)
        for p_ in self.prims:
            p_["center"] = np.asarray(p_["center"]) / fL
            if "size" in p_:
                p_["size"] = np.asarray(p_["size"]) / fL
            for key in ("radius", "length"):
                if key in p_ and p_[key] is not None:
                    p_[key] = p_[key] / fL
        self.obj_kinds, self.obj_pars, self.obj_action, self.obj_rest = geometry.kernel_arrays(self.prims)
        self.obj_h = 1e-4 * geometry.characteristic_size(self.prims) if self.prims else 1e-6

        # -------------------------------------------------- state arrays
        self.state = np.zeros(self.n, np.int8)
        self.absorb_t = -np.ones(self.n)
        self.absorb_q = np.zeros(self.n)
        self.hit_obj = -np.ones(self.n, np.int32)
        self.force = np.zeros((self.n, 3))
        self.pe = np.zeros(self.n)
        self.dv = np.zeros((self.n, 3))
        self.t = 0.0
        self.step = 0
        self.collected = 0.0
        self.ncoll = 0
        self.wall_time = 0.0

        self._setup_domain()
        self._setup_interaction()
        self._setup_drag(gam_in)
        self._setup_forces()
        dt = cfg["time"]["dt"]
        self.dt = dt
        self.nsteps = cfg["time"]["steps"] or int(math.ceil(cfg["time"]["t_end"] / dt - 1e-9))
        self.collisions = cfg["collisions"]["enabled"]
        if self.collisions and self.rad.max() <= 0:
            raise ValueError("collisions.enabled needs particles with radius > 0")
        if cfg["run"]["restart"]:
            self._restart(os.path.join(self.base, cfg["run"]["restart"]))

    # ------------------------------------------------------------------ setup helpers
    def _drag_rate_input(self, cfg, r, m):
        nt = cfg["neutrals"]
        if nt["drag"] == "constant":
            return np.full(len(m), nt["gamma"])
        if nt["drag"] == "epstein":
            Tn = nt["temperature"] if nt["temperature"] > 0 else 300.0
            return particles.epstein_gamma(r, m, nt["pressure"], Tn, nt["mass_amu"], nt["delta"])
        return np.zeros(len(m))

    def _make_scales(self, cfg, m, q):
        if not self.si:
            return Scales()
        u = cfg["units"]
        a = u["length"]
        if a <= 0:
            dom = cfg["domain"]
            vol = (4 / 3 * math.pi * dom["radius"] ** 3) if dom["shape"] == "sphere" else float(np.prod(dom["box"]))
            a = (3 * vol / (4 * math.pi * len(m))) ** (1 / 3)
        Q0 = u["charge"] if u["charge"] > 0 else float(np.mean(np.abs(q)))
        if Q0 <= 0:
            Q0 = 1.602176634e-19
        m0 = u["mass"] if u["mass"] > 0 else float(np.mean(m))
        return Scales(length=a, charge=Q0, mass=m0, si=True)

    def _setup_domain(self):
        d = self.cfg["domain"]
        c = np.asarray(d["center"], float)
        bnd = self.cfg["boundary"]
        if d["shape"] == "sphere":
            R = d["radius"]
            self.sphere = (c, R, BC_CODES[bnd["sphere"]])
            self.lo, self.hi = c - R, c + R
            self.bc = np.full((3, 2), kernels.BC_OPEN, np.int64)
            self.periodic = np.zeros(3, np.int64)
        else:
            self.sphere = None
            L = np.asarray(d["box"], float)
            self.lo, self.hi = c - L / 2, c + L / 2
            self.bc = np.array([[BC_CODES[s] for s in bnd[ax]] for ax in "xyz"], np.int64)
            self.periodic = (self.bc[:, 0] == kernels.BC_PERIODIC).astype(np.int64)
        self.L = self.hi - self.lo
        self.remove = bnd["on_absorb"] == "remove"
        self.neutralize = bnd["neutralize"]
        self.wall_rest = bnd["restitution"]

    def _setup_interaction(self):
        it = self.cfg["interaction"]
        self.pair = it["model"] != "none"
        if it["model"] == "coulomb":
            kappa = 0.0
        elif it["debye_length"] > 0:
            kappa = 1.0 / it["debye_length"]
        else:
            kappa = it["kappa"]
        self.kappa = kappa
        qmax = np.abs(self.q).max() or 1.0
        rc = it["cutoff"] if it["cutoff"] > 0 else auto_cutoff(kappa, it["tolerance"])
        half = np.where(self.periodic == 1, self.L / 2, np.inf)
        if rc > half.min():
            if it["method"] == "cells" or (it["cutoff"] > 0):
                self.log(f"note: cutoff {rc:.3g} exceeds half the periodic box; clipped to {half.min():.3g}")
            rc = min(rc, half.min())
        self.rc = rc
        method = it["method"]
        if method == "auto":
            method = "cells" if np.isfinite(rc) else "direct"
        self.method = method
        if self.pair and method == "direct" and self.n > 20000:
            self.log(f"warning: O(N^2) direct sum with N = {self.n}; increase kappa, set a cutoff or use a larger box for the cell list")
        if self.pair and kappa == 0 and self.periodic.any():
            self.log("warning: unscreened Coulomb with periodic boundaries uses a minimum-image sum (no Ewald); prefer yukawa")
        # cells of size >= rc/reach; reach 2 halves the volume searched per particle
        self.reach = 2
        rcell = (rc if np.isfinite(rc) else self.L.max()) / self.reach
        self.nc = np.minimum(np.maximum((self.L / rcell).astype(np.int64), 1), 1024)
        rmax = self.rad.max() if self.n else 0.0
        self.nc_coll = np.minimum(np.maximum((self.L / max(2 * rmax, 1e-12)).astype(np.int64), 1), 512)
        self.qmax = qmax

    def _setup_drag(self, gam_in):
        nt = self.cfg["neutrals"]
        th = self.cfg["thermostat"]
        self.gamma = gam_in * self.scales.factor("T") if self.si else gam_in.copy()  # rate * time unit
        self.Tn = np.full(self.n, nt["temperature"] if nt["noise"] else 0.0)
        self.th_type = th["type"]
        self.th_until = th["until"]
        self._langevin_on = None
        self._update_ou()

    def _update_ou(self):
        th = self.cfg["thermostat"]
        on = self.th_type == "langevin" and (self.th_until < 0 or self.t < self.th_until)
        if on == self._langevin_on:
            return
        self._langevin_on = on
        g = self.gamma.copy()
        T = self.Tn * self.gamma
        if on:
            g += th["gamma"]
            T += th["temperature"] * th["gamma"]
        Teff = np.divide(T, g, out=np.zeros_like(g), where=g > 0)
        c1 = np.exp(-g * self.cfg["time"]["dt"])
        self.c1 = c1
        self.c2 = np.sqrt(np.maximum((1 - c1 ** 2) * Teff / self.m, 0.0))
        self.use_ou = bool((g > 0).any())

    def _setup_forces(self):
        f = self.cfg["forces"]
        sc = self.scales
        ep = np.zeros(kernels.EP_SIZE)
        ef = np.zeros(kernels.EF_SIZE, np.int64)
        K = kernels
        if f["gravity"]:
            ef[K.EF_G] = 1; ep[K.EP_G:K.EP_G + 3] = f["gravity"]["g"]
        if f["efield"]:
            ef[K.EF_E] = 1; ep[K.EP_E:K.EP_E + 3] = f["efield"]["E"]
        if f["linear_efield"]:
            le = f["linear_efield"]
            ef[K.EF_LE] = 1
            ep[K.EP_LE0:K.EP_LE0 + 3] = le["E0"]; ep[K.EP_LEG:K.EP_LEG + 3] = le["gradient"]
            ep[K.EP_LEO:K.EP_LEO + 3] = le["origin"]
        if f["harmonic"]:
            ef[K.EF_H] = 1
            ep[K.EP_HW:K.EP_HW + 3] = f["harmonic"]["omega2"]; ep[K.EP_HO:K.EP_HO + 3] = f["harmonic"]["origin"]
        if f["constant"]:
            ef[K.EF_C] = 1; ef[K.EF_CS] = SCALE_CODES[f["constant"]["scale_with"]]
            ep[K.EP_CF:K.EP_CF + 3] = f["constant"]["F"]
        if f["void"]:
            v = f["void"]
            ef[K.EF_V] = 1; ef[K.EF_VS] = SCALE_CODES[v["scale_with"]]; ef[K.EF_VCYL] = int(v["cylindrical"])
            ep[K.EP_VA] = v["amplitude"]; ep[K.EP_VL] = v["decay_length"]
            ep[K.EP_VC:K.EP_VC + 3] = v["center"]; ep[K.EP_VM] = v["modulation"]; ep[K.EP_VW] = v["omega"]
        # radial tables
        rs, fs, off, cs, ss = [], [], [0], [], []
        for tab in f["radial_table"]:
            r, F = self._read_table(tab)
            if tab["kind"] == "efield":
                F = F / sc.factor("E")
                code = SCALE_CODES["charge"]
            else:
                F = F / sc.factor("F")
                code = SCALE_CODES[tab["scale_with"]]
            rs.append(r); fs.append(F); off.append(off[-1] + len(r)); cs.append(tab["center"]); ss.append(code)
        self.tab = (np.concatenate(rs) if rs else np.zeros(1), np.concatenate(fs) if fs else np.zeros(1),
                    np.asarray(off, np.int64), np.asarray(cs, float).reshape(-1, 3),
                    np.asarray(ss, np.int64))
        # grid field
        fld = self.cfg["field"]
        if fld:
            path = os.path.join(self.base, fld["file"])
            if not os.path.exists(path) and fld["auto_solve"] and self.cfg_input["objects"]:
                self.log(f"field: {fld['file']} not found, solving it from [[objects]] and [field_solver]")
                gs = fields.solve(self.cfg_input, geometry.expand(self.cfg_input["objects"]), log=self.log)
                gs.save(path, {"debye_length": self.cfg_input["field_solver"]["debye_length"]})
            g = fields.load(path, fld["format"], fld["step"],
                            np.asarray(self.cfg_input["domain"]["center"]), fld["center_grid"])
            if g.units == "SI" and not self.si:
                raise ValueError(f"{fld['file']} is in SI units; set units.system = \"SI\" to use it")
            if g.units == "reduced" and self.si:
                raise ValueError(f"{fld['file']} is in reduced units but the run uses SI input")
            fL = sc.factor("L")
            self.grid = (g.origin / fL, g.spacing / fL, np.ascontiguousarray(g.E * fld["scale"] / sc.factor("E")))
            ef[K.EF_GRID] = 1
            ef[K.EF_GOUT] = {"zero": 0, "clamp": 1, "periodic": 2}[fld["outside"]]
            self.log(f"field: {os.path.basename(path)} grid {g.shape}, |E|max = {np.abs(g.E).max():.4g} ({g.units})")
        else:
            self.grid = (np.zeros(3), np.ones(3), np.zeros((2, 2, 2, 3)))
        self.ep, self.ef = ep, ef

    def _read_table(self, tab):
        path = os.path.join(self.base, tab["file"])
        if path.endswith(".npz"):
            d = np.load(path)
            r, F = np.asarray(d[tab["r_key"]], float), np.asarray(d[tab["f_key"]], float)
        else:
            arr = np.loadtxt(path, delimiter="," if path.endswith(".csv") else None, comments="#")
            r, F = arr[:, 0], arr[:, 1]
        r = r * tab["r_unit"]
        F = F * tab["factor"]
        keep = r >= 0          # symmetric (legacy dust-exp) tables: use the r >= 0 half
        r, F = r[keep], F[keep]
        o = np.argsort(r)
        return np.ascontiguousarray(r[o] / self.scales.factor("L")), np.ascontiguousarray(F[o])

    def _restart(self, path):
        d = io.load_checkpoint(path)
        if len(d["pos"]) != self.n:
            raise ValueError(f"checkpoint has {len(d['pos'])} particles, input has {self.n}")
        for k in io.CHECKPOINT_FIELDS:
            if k in d:
                setattr(self, k, np.ascontiguousarray(d[k]))
        self.t, self.step, self.collected = d["t"], d["step"], d["collected_charge"]
        self.rng.bit_generator.state = d["rng"]
        self._langevin_on = None
        self._update_ou()
        self.log(f"restarted from {path} at t = {self.t:.6g} (step {self.step})")

    # ------------------------------------------------------------------ physics
    def compute_forces(self):
        if self.pair:
            if self.method == "cells":
                start, order = kernels.build_cells(self.pos, self.state, self.lo, self.L, self.nc)
                kernels.pair_cells(self.pos, self.q, self.state, self.force, self.pe, self.kappa, self.rc,
                                   self.lo, self.L, self.periodic, self.nc, self.reach, start, order)
            else:
                kernels.pair_direct(self.pos, self.q, self.state, self.force, self.pe, self.kappa, self.rc,
                                    self.L, self.periodic)
        else:
            self.force[:] = 0.0
            self.pe[:] = 0.0
        tr, tf, toff, tc, ts = self.tab
        go, gd, g = self.grid
        kernels.external_forces(self.pos, self.q, self.m, self.rr, self.state, self.force, self.t,
                                self.ep, self.ef, tr, tf, toff, tc, ts, go, gd, g)

    def apply_walls(self):
        if self.sphere is not None:
            c, R, code = self.sphere
            kernels.apply_sphere(self.pos, self.vel, self.q, self.state, self.absorb_t, self.absorb_q, self.t,
                                 c, R, code, self.wall_rest, self.remove, self.neutralize)
        else:
            kernels.apply_box(self.pos, self.vel, self.q, self.state, self.absorb_t, self.absorb_q, self.t,
                              self.lo, self.hi, self.bc, self.wall_rest, self.remove, self.neutralize)
        if len(self.prims):
            kernels.apply_objects(self.pos, self.vel, self.q, self.state, self.absorb_t, self.absorb_q,
                                  self.hit_obj, self.t, self.obj_kinds, self.obj_pars, self.obj_action,
                                  self.obj_rest, self.obj_h, self.remove, self.neutralize)

    def advance(self):
        dt = self.dt
        K = kernels
        K.kick(self.vel, self.force, self.m, self.state, 0.5 * dt)
        K.drift(self.pos, self.vel, self.state, 0.5 * dt)
        if self.use_ou:
            xi = self.rng.standard_normal((self.n, 3))
            K.ou_step(self.vel, self.state, self.c1, self.c2, xi)
        K.drift(self.pos, self.vel, self.state, 0.5 * dt)
        self.t += dt
        self.step += 1
        self.apply_walls()
        self.compute_forces()
        K.kick(self.vel, self.force, self.m, self.state, 0.5 * dt)
        if self.collisions:
            start, order = K.build_cells(self.pos, self.state, self.lo, self.L, self.nc_coll)
            self.ncoll += K.collide(self.pos, self.vel, self.m, self.rad, self.state,
                                    self.cfg["collisions"]["restitution"], self.lo, self.L, self.periodic,
                                    self.nc_coll, start, order, self.dv)
        if self.th_type == "berendsen" and (self.th_until < 0 or self.t < self.th_until):
            self._berendsen()
        if self.th_type == "langevin":
            self._update_ou()

    def _berendsen(self):
        th = self.cfg["thermostat"]
        kx, ky, kz, na = kernels.kinetic(self.vel, self.m, self.state)
        if na == 0:
            return
        T = 2 * (kx + ky + kz) / (3 * na)
        if T <= 0:
            return
        lam = math.sqrt(max(1 + self.dt / th["tau"] * (th["temperature"] / T - 1), 0.0))
        act = self.state == kernels.ACTIVE
        self.vel[act] *= lam

    def diagnostics(self) -> dict:
        kx, ky, kz, na = kernels.kinetic(self.vel, self.m, self.state)
        kin = kx + ky + kz
        pe = float(self.pe.sum())
        nst = int((self.state == kernels.STUCK).sum())
        nrm = int((self.state == kernels.REMOVED).sum())
        nn = max(na, 1)
        return {"step": self.step, "time": self.t, "kinetic": kin, "kx": kx, "ky": ky, "kz": kz,
                "potential": pe, "total": kin + pe, "temperature": 2 * kin / (3 * nn),
                "tx": 2 * kx / nn, "ty": 2 * ky / nn, "tz": 2 * kz / nn,
                "n_active": na, "n_stuck": nst, "n_removed": nrm,
                "collected_charge": float(self.absorb_q.sum()), "collisions": self.ncoll}

    # ------------------------------------------------------------------ run
    def summary(self) -> str:
        sc = self.scales
        lines = [f"PPDyn: {self.n} particles in {len(self.cfg['species'])} species, {self.nsteps} steps of dt = {self.dt:.4g}",
                 f"  units: {sc.describe()}",
                 f"  interaction: {self.cfg['interaction']['model']}, kappa = {self.kappa:.4g}, "
                 f"cutoff = {self.rc:.4g}, method = {self.method}"
                 + (f" ({'x'.join(map(str, self.nc))} cells)" if self.method == "cells" else ""),
                 f"  domain: {self.cfg['domain']['shape']} {np.round(self.L, 4).tolist()}, "
                 f"periodic {['xyz'[d] for d in range(3) if self.periodic[d]] or 'none'}"]
        if self.use_ou or self.th_type != "none":
            lines.append(f"  drag rate (mean) = {self.gamma.mean():.4g}, thermostat = {self.th_type}")
        if len(self.prims):
            lines.append(f"  objects: {len(self.prims)} primitives")
        if self.dt * math.sqrt(3) > 0.2 and self.pair:
            lines.append("  warning: dt * omega_pd > 0.2; consider a smaller time step")
        import numba
        lines.append(f"  threads: {numba.get_num_threads()}")
        return "\n".join(lines)

    def run(self, writer: io.Writer | None = None, progress=True, callback=None):
        out = self.cfg["output"]
        self.log(self.summary())
        self.compute_forces()
        if writer and self.step == 0:
            writer.history(self.diagnostics())
            writer.snapshot(self)
        ckpt = os.path.join(out["dir"], "checkpoint.h5") if writer else None
        t0 = time.time()
        last = 0.0
        first = self.step
        while self.step < self.nsteps:
            self.advance()
            if writer:
                if self.step % out["history_every"] == 0:
                    writer.history(self.diagnostics())
                if self.step % out["every"] == 0 or self.step == self.nsteps:
                    writer.snapshot(self)
                if out["checkpoint_every"] and self.step % out["checkpoint_every"] == 0:
                    io.save_checkpoint(ckpt, self)
            if callback:
                callback(self)
            now = time.time()
            if progress and (now - last > 0.5 or self.step == self.nsteps):
                last = now
                done = self.step - first
                rate = done / max(now - t0, 1e-9)
                eta = (self.nsteps - self.step) / max(rate, 1e-9)
                frac = self.step / self.nsteps
                bar = "#" * int(30 * frac)
                sys.stderr.write(f"\r  [{bar:<30}] {100 * frac:5.1f}%  step {self.step}/{self.nsteps}  "
                                 f"{rate:8.1f} steps/s  ETA {eta:6.0f} s ")
                sys.stderr.flush()
        if progress:
            sys.stderr.write("\n")
        self.wall_time = time.time() - t0
        if writer:
            writer.finalize(self)
            io.save_checkpoint(ckpt, self)
        steps_done = self.step - first
        if steps_done:
            self.log(f"done: {steps_done} steps in {self.wall_time:.2f} s "
                     f"({steps_done * self.n / max(self.wall_time, 1e-9) / 1e6:.2f} M particle-steps/s)")
        return self


def run(cfg: dict, log=print, progress=None):
    """Convenience: build, run and write a simulation from a completed config."""
    out = cfg["output"]["dir"]
    if not os.path.isabs(out) and cfg.get("_source", "<string>") != "<string>":
        cfg = copy.deepcopy(cfg)
    if cfg["output"]["overwrite"] and not cfg["run"]["restart"] and os.path.isdir(out):
        shutil.rmtree(out)
    sim = Simulation(cfg, log)
    writer = io.Writer(sim.cfg, sim)
    try:
        sim.run(writer, progress=cfg["run"]["progress"] if progress is None else progress)
    finally:
        writer.close()
    return sim, writer.path
