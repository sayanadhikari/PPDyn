"""Numba kernels. All arrays are in reduced units.

Particle state: 0 = active, 1 = stuck on a wall/object (frozen), 2 = removed.
Pair forces and collisions skip removed particles; stuck particles keep acting
on others through their (possibly zeroed) charge but never move.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

ACTIVE, STUCK, REMOVED = 0, 1, 2
FASTMATH = {"nnan", "ninf", "nsz", "arcp", "contract", "afn", "reassoc"}

# boundary codes per face
BC_PERIODIC, BC_REFLECT, BC_ABSORB, BC_OPEN = 0, 1, 2, 3


# ======================================================================
# cell list
# ======================================================================
@njit(cache=True)
def build_cells(pos, state, lo, L, nc):
    """Counting sort of particles into cells. Returns (start, order).

    start has ncell+1 entries; particles of cell c are order[start[c]:start[c+1]].
    Particles outside non-periodic boxes are clamped into the edge cells.
    Removed particles are left out.
    """
    n = pos.shape[0]
    ncell = nc[0] * nc[1] * nc[2]
    cid = np.empty(n, np.int64)
    counts = np.zeros(ncell + 1, np.int64)
    for i in range(n):
        if state[i] == REMOVED:
            cid[i] = -1
            continue
        c3 = np.empty(3, np.int64)
        for d in range(3):
            s = (pos[i, d] - lo[d]) / L[d]
            s = s - math.floor(s) if s >= 1.0 or s < 0.0 else s
            k = int(s * nc[d])
            if k >= nc[d]:
                k = nc[d] - 1
            if k < 0:
                k = 0
            c3[d] = k
        c = (c3[0] * nc[1] + c3[1]) * nc[2] + c3[2]
        cid[i] = c
        counts[c + 1] += 1
    for c in range(ncell):
        counts[c + 1] += counts[c]
    start = counts.copy()
    fill = counts[:-1].copy()
    order = np.empty(counts[ncell], np.int64)
    for i in range(n):
        c = cid[i]
        if c >= 0:
            order[fill[c]] = i
            fill[c] += 1
    return start, order


@njit(cache=True, inline="always")
def _neigh_range(c, n, periodic, reach, out):
    """Unique neighbour cell indices along one axis (cells within +-reach)."""
    k = 0
    if periodic and 2 * reach + 1 >= n:
        for x in range(n):
            out[k] = x
            k += 1
        return k
    for d in range(-reach, reach + 1):
        x = c + d
        if periodic:
            x = x % n
        elif x < 0 or x >= n:
            continue
        out[k] = x
        k += 1
    return k


# ======================================================================
# pair forces
# ======================================================================
@njit(cache=True, inline="always")
def _min_image(d, L, per):
    # positions on periodic axes are always wrapped into the box, so |d| < L
    if per:
        if d > 0.5 * L:
            return d - L
        if d < -0.5 * L:
            return d + L
    return d


@njit(parallel=True, cache=True)
def _gather(pos, q, order):
    m = order.shape[0]
    ps = np.empty((m, 3))
    qs = np.empty(m)
    for p in prange(m):
        j = order[p]
        ps[p, 0] = pos[j, 0]
        ps[p, 1] = pos[j, 1]
        ps[p, 2] = pos[j, 2]
        qs[p] = q[j]
    return ps, qs


@njit(parallel=True, fastmath=FASTMATH, cache=True)
def _pair_sorted(ps, qs, order, force, pe, kappa, rc, lo, L, periodic, nc, reach, start):
    m = ps.shape[0]
    rc2 = rc * rc
    ushift = math.exp(-kappa * rc) / rc if rc < 1e300 else 0.0   # energy shift: U(rc) = 0
    h0 = 0.5 * L[0]
    h1 = 0.5 * L[1]
    h2 = 0.5 * L[2]
    for pi in prange(m):
        i = order[pi]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        u = 0.0
        qi = qs[pi]
        if qi != 0.0:
            xi = ps[pi, 0]
            yi = ps[pi, 1]
            zi = ps[pi, 2]
            c = np.empty(3, np.int64)
            for d in range(3):
                s = (ps[pi, d] - lo[d]) / L[d]
                if s >= 1.0 or s < 0.0:
                    s = s - math.floor(s)
                k = int(s * nc[d])
                c[d] = min(max(k, 0), nc[d] - 1)
            w = 2 * reach + 1
            ax = np.empty(max(w, nc[0]), np.int64)
            ay = np.empty(max(w, nc[1]), np.int64)
            az = np.empty(max(w, nc[2]), np.int64)
            nx = _neigh_range(c[0], nc[0], periodic[0], reach, ax)
            ny = _neigh_range(c[1], nc[1], periodic[1], reach, ay)
            nz = _neigh_range(c[2], nc[2], periodic[2], reach, az)
            px = periodic[0]
            py = periodic[1]
            pz = periodic[2]
            for a in range(nx):
                for b in range(ny):
                    base = (ax[a] * nc[1] + ay[b]) * nc[2]
                    for e in range(nz):
                        cc = base + az[e]
                        for pj in range(start[cc], start[cc + 1]):
                            if pj == pi:
                                continue
                            qj = qs[pj]
                            dx = xi - ps[pj, 0]
                            dy = yi - ps[pj, 1]
                            dz = zi - ps[pj, 2]
                            if px:
                                if dx > h0:
                                    dx -= L[0]
                                elif dx < -h0:
                                    dx += L[0]
                            if py:
                                if dy > h1:
                                    dy -= L[1]
                                elif dy < -h1:
                                    dy += L[1]
                            if pz:
                                if dz > h2:
                                    dz -= L[2]
                                elif dz < -h2:
                                    dz += L[2]
                            r2 = dx * dx + dy * dy + dz * dz
                            if r2 >= rc2 or r2 == 0.0 or qj == 0.0:
                                continue
                            r = math.sqrt(r2)
                            ex = math.exp(-kappa * r)
                            qq = qi * qj
                            uu = qq * ex / r
                            f = uu * (1.0 + kappa * r) / r2
                            fx += f * dx
                            fy += f * dy
                            fz += f * dz
                            u += uu - qq * ushift
        force[i, 0] = fx
        force[i, 1] = fy
        force[i, 2] = fz
        pe[i] = 0.5 * u


def pair_cells(pos, q, state, force, pe, kappa, rc, lo, L, periodic, nc, reach, start, order):
    """Screened-Coulomb forces with a linked-cell list. Particle data are first
    gathered in cell order (cache locality); each particle accumulates its own
    force, so the kernel is race-free. Removed particles are not in ``order``."""
    force[:] = 0.0
    pe[:] = 0.0
    ps, qs = _gather(pos, q, order)
    _pair_sorted(ps, qs, order, force, pe, kappa, rc, lo, L, periodic, nc, reach, start)


@njit(parallel=True, fastmath=FASTMATH, cache=True)
def pair_direct(pos, q, state, force, pe, kappa, rc, L, periodic):
    """All-pairs screened Coulomb with minimum image on periodic axes."""
    n = pos.shape[0]
    rc2 = rc * rc
    ushift = math.exp(-kappa * rc) / rc if rc < 1e300 else 0.0
    for i in prange(n):
        fx = 0.0
        fy = 0.0
        fz = 0.0
        u = 0.0
        if state[i] != REMOVED and q[i] != 0.0:
            for j in range(n):
                if j == i or state[j] == REMOVED or q[j] == 0.0:
                    continue
                dx = _min_image(pos[i, 0] - pos[j, 0], L[0], periodic[0])
                dy = _min_image(pos[i, 1] - pos[j, 1], L[1], periodic[1])
                dz = _min_image(pos[i, 2] - pos[j, 2], L[2], periodic[2])
                r2 = dx * dx + dy * dy + dz * dz
                if r2 >= rc2 or r2 == 0.0:
                    continue
                r = math.sqrt(r2)
                ex = math.exp(-kappa * r)
                qq = q[i] * q[j]
                f = qq * (1.0 + kappa * r) * ex / (r2 * r)
                fx += f * dx
                fy += f * dy
                fz += f * dz
                u += qq * (ex / r - ushift)
        force[i, 0] = fx
        force[i, 1] = fy
        force[i, 2] = fz
        pe[i] = 0.5 * u


# ======================================================================
# external forces
# ======================================================================
# layout of the float parameter vector ``ep`` (see forces.py)
EP_G, EP_E, EP_LE0, EP_LEG, EP_LEO, EP_HW, EP_HO, EP_CF = 0, 3, 6, 9, 12, 15, 18, 21
EP_VA, EP_VL, EP_VC, EP_VM, EP_VW = 24, 25, 26, 29, 30
EP_SIZE = 32
# int flags ``ef``
EF_G, EF_E, EF_LE, EF_H, EF_C, EF_CS, EF_V, EF_VS, EF_VCYL, EF_GRID, EF_GOUT = range(11)
EF_SIZE = 11


@njit(cache=True, inline="always")
def _scale(code, q, m, rr):
    if code == 0:
        return 1.0
    if code == 1:
        return q
    if code == 2:
        return m
    if code == 3:
        return rr
    return rr * rr


@njit(cache=True)
def _trilinear(px, py, pz, g0, gd, grid, outside):
    nx = grid.shape[0]
    ny = grid.shape[1]
    nz = grid.shape[2]
    s = np.empty(3)
    s[0] = (px - g0[0]) / gd[0]
    s[1] = (py - g0[1]) / gd[1]
    s[2] = (pz - g0[2]) / gd[2]
    nn = (nx, ny, nz)
    res = np.zeros(3)
    idx = np.empty(3, np.int64)
    w = np.empty(3)
    for d in range(3):
        n = nn[d]
        x = s[d]
        if outside == 2:  # periodic: grid spans n-1 intervals
            span = n - 1
            x = x - span * math.floor(x / span)
        elif x < 0.0 or x > n - 1:
            if outside == 0:
                return res
            x = min(max(x, 0.0), n - 1.0)
        k = int(math.floor(x))
        if k > n - 2:
            k = n - 2
        if k < 0:
            k = 0
        idx[d] = k
        w[d] = x - k
    for a in range(2):
        wa = w[0] if a else 1.0 - w[0]
        for b in range(2):
            wb = w[1] if b else 1.0 - w[1]
            for c in range(2):
                wc = w[2] if c else 1.0 - w[2]
                ww = wa * wb * wc
                for d in range(3):
                    res[d] += ww * grid[idx[0] + a, idx[1] + b, idx[2] + c, d]
    return res


@njit(parallel=True, fastmath=FASTMATH, cache=True)
def external_forces(pos, q, m, rr, state, force, t, ep, ef,
                    tab_r, tab_f, tab_off, tab_c, tab_s, grid_o, grid_d, grid):
    """Add all external forces (not accelerations) to ``force``."""
    n = pos.shape[0]
    ntab = tab_off.shape[0] - 1
    for i in prange(n):
        if state[i] != ACTIVE:
            continue
        x = pos[i, 0]
        y = pos[i, 1]
        z = pos[i, 2]
        qi = q[i]
        mi = m[i]
        ri = rr[i]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        if ef[EF_G]:
            fx += mi * ep[EP_G]
            fy += mi * ep[EP_G + 1]
            fz += mi * ep[EP_G + 2]
        if ef[EF_E]:
            fx += qi * ep[EP_E]
            fy += qi * ep[EP_E + 1]
            fz += qi * ep[EP_E + 2]
        if ef[EF_LE]:
            fx += qi * (ep[EP_LE0] + ep[EP_LEG] * (x - ep[EP_LEO]))
            fy += qi * (ep[EP_LE0 + 1] + ep[EP_LEG + 1] * (y - ep[EP_LEO + 1]))
            fz += qi * (ep[EP_LE0 + 2] + ep[EP_LEG + 2] * (z - ep[EP_LEO + 2]))
        if ef[EF_H]:
            fx -= mi * ep[EP_HW] * (x - ep[EP_HO])
            fy -= mi * ep[EP_HW + 1] * (y - ep[EP_HO + 1])
            fz -= mi * ep[EP_HW + 2] * (z - ep[EP_HO + 2])
        if ef[EF_C]:
            s = _scale(ef[EF_CS], qi, mi, ri)
            fx += s * ep[EP_CF]
            fy += s * ep[EP_CF + 1]
            fz += s * ep[EP_CF + 2]
        if ef[EF_V]:
            dx = x - ep[EP_VC]
            dy = y - ep[EP_VC + 1]
            dz = 0.0 if ef[EF_VCYL] else z - ep[EP_VC + 2]
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            if r > 0.0:
                amp = ep[EP_VA] * math.exp(-r / ep[EP_VL]) * (1.0 + ep[EP_VM] * math.sin(ep[EP_VW] * t))
                amp *= _scale(ef[EF_VS], qi, mi, ri) / r
                fx += amp * dx
                fy += amp * dy
                fz += amp * dz
        for k in range(ntab):
            dx = x - tab_c[k, 0]
            dy = y - tab_c[k, 1]
            dz = z - tab_c[k, 2]
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            if r > 0.0:
                a0 = tab_off[k]
                a1 = tab_off[k + 1]
                fr = np.interp(r, tab_r[a0:a1], tab_f[a0:a1])
                amp = fr * _scale(tab_s[k], qi, mi, ri) / r
                fx += amp * dx
                fy += amp * dy
                fz += amp * dz
        if ef[EF_GRID]:
            e = _trilinear(x, y, z, grid_o, grid_d, grid, ef[EF_GOUT])
            fx += qi * e[0]
            fy += qi * e[1]
            fz += qi * e[2]
        force[i, 0] += fx
        force[i, 1] += fy
        force[i, 2] += fz


# ======================================================================
# integrator pieces (BAOAB)
# ======================================================================
@njit(parallel=True, fastmath=FASTMATH, cache=True)
def kick(vel, force, m, state, h):
    for i in prange(vel.shape[0]):
        if state[i] == ACTIVE:
            s = h / m[i]
            vel[i, 0] += s * force[i, 0]
            vel[i, 1] += s * force[i, 1]
            vel[i, 2] += s * force[i, 2]


@njit(parallel=True, fastmath=FASTMATH, cache=True)
def drift(pos, vel, state, h):
    for i in prange(pos.shape[0]):
        if state[i] == ACTIVE:
            pos[i, 0] += h * vel[i, 0]
            pos[i, 1] += h * vel[i, 1]
            pos[i, 2] += h * vel[i, 2]


@njit(parallel=True, fastmath=FASTMATH, cache=True)
def ou_step(vel, state, c1, c2, xi):
    """Exact Ornstein-Uhlenbeck update v <- c1 v + c2 xi (per particle)."""
    for i in prange(vel.shape[0]):
        if state[i] == ACTIVE:
            for d in range(3):
                vel[i, d] = c1[i] * vel[i, d] + c2[i] * xi[i, d]


@njit(parallel=True, cache=True)
def kinetic(vel, m, state):
    """Return per-axis kinetic energy sums and number of active particles."""
    n = vel.shape[0]
    kx = 0.0
    ky = 0.0
    kz = 0.0
    na = 0
    for i in prange(n):
        if state[i] == ACTIVE:
            kx += 0.5 * m[i] * vel[i, 0] * vel[i, 0]
            ky += 0.5 * m[i] * vel[i, 1] * vel[i, 1]
            kz += 0.5 * m[i] * vel[i, 2] * vel[i, 2]
            na += 1
    return kx, ky, kz, na


# ======================================================================
# boundaries
# ======================================================================
@njit(cache=True, inline="always")
def _absorb(i, pos, vel, q, state, absorb_t, absorb_q, t, remove, neutralize):
    vel[i, 0] = 0.0
    vel[i, 1] = 0.0
    vel[i, 2] = 0.0
    absorb_t[i] = t
    absorb_q[i] = q[i]
    if neutralize:
        q[i] = 0.0
    state[i] = REMOVED if remove else STUCK


@njit(parallel=True, cache=True)
def apply_box(pos, vel, q, state, absorb_t, absorb_q, t, lo, hi, bc, rest, remove, neutralize):
    """bc is a (3,2) int array of face codes."""
    n = pos.shape[0]
    for i in prange(n):
        if state[i] != ACTIVE:
            continue
        for d in range(3):
            L = hi[d] - lo[d]
            x = pos[i, d]
            if bc[d, 0] == BC_PERIODIC:
                if x < lo[d] or x >= hi[d]:
                    pos[i, d] = x - L * math.floor((x - lo[d]) / L)
                continue
            if x < lo[d]:
                code = bc[d, 0]
                if code == BC_REFLECT:
                    pos[i, d] = lo[d] + (lo[d] - x) * rest
                    vel[i, d] = -rest * vel[i, d]
                    if pos[i, d] > hi[d]:
                        pos[i, d] = hi[d]
                elif code == BC_ABSORB:
                    pos[i, d] = lo[d]
                    _absorb(i, pos, vel, q, state, absorb_t, absorb_q, t, remove, neutralize)
                    break
                else:
                    state[i] = REMOVED
                    absorb_t[i] = t
                    break
            elif x > hi[d]:
                code = bc[d, 1]
                if code == BC_REFLECT:
                    pos[i, d] = hi[d] - (x - hi[d]) * rest
                    vel[i, d] = -rest * vel[i, d]
                    if pos[i, d] < lo[d]:
                        pos[i, d] = lo[d]
                elif code == BC_ABSORB:
                    pos[i, d] = hi[d]
                    _absorb(i, pos, vel, q, state, absorb_t, absorb_q, t, remove, neutralize)
                    break
                else:
                    state[i] = REMOVED
                    absorb_t[i] = t
                    break


@njit(parallel=True, cache=True)
def apply_sphere(pos, vel, q, state, absorb_t, absorb_q, t, c, R, code, rest, remove, neutralize):
    n = pos.shape[0]
    for i in prange(n):
        if state[i] != ACTIVE:
            continue
        dx = pos[i, 0] - c[0]
        dy = pos[i, 1] - c[1]
        dz = pos[i, 2] - c[2]
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        if r <= R:
            continue
        nx = dx / r
        ny = dy / r
        nz = dz / r
        if code == BC_REFLECT:
            rn = R - (r - R) * rest
            if rn < 0.0:
                rn = 0.0
            pos[i, 0] = c[0] + nx * rn
            pos[i, 1] = c[1] + ny * rn
            pos[i, 2] = c[2] + nz * rn
            vn = vel[i, 0] * nx + vel[i, 1] * ny + vel[i, 2] * nz
            if vn > 0.0:
                vel[i, 0] -= (1.0 + rest) * vn * nx
                vel[i, 1] -= (1.0 + rest) * vn * ny
                vel[i, 2] -= (1.0 + rest) * vn * nz
        elif code == BC_ABSORB:
            pos[i, 0] = c[0] + nx * R
            pos[i, 1] = c[1] + ny * R
            pos[i, 2] = c[2] + nz * R
            _absorb(i, pos, vel, q, state, absorb_t, absorb_q, t, remove, neutralize)
        else:
            state[i] = REMOVED
            absorb_t[i] = t


# ======================================================================
# solid objects (signed distance functions)
# ======================================================================
OBJ_BOX, OBJ_SPHERE, OBJ_CYL = 0, 1, 2


@njit(cache=True)
def sdf_one(kind, par, x, y, z):
    """par: cx, cy, cz, p0, p1, p2, axis. Box: half sizes. Sphere: radius.
    Cylinder: radius, half length, axis (0,1,2)."""
    px = x - par[0]
    py = y - par[1]
    pz = z - par[2]
    if kind == OBJ_SPHERE:
        return math.sqrt(px * px + py * py + pz * pz) - par[3]
    if kind == OBJ_BOX:
        qx = abs(px) - par[3]
        qy = abs(py) - par[4]
        qz = abs(pz) - par[5]
        ox = max(qx, 0.0)
        oy = max(qy, 0.0)
        oz = max(qz, 0.0)
        return math.sqrt(ox * ox + oy * oy + oz * oz) + min(max(qx, max(qy, qz)), 0.0)
    ax = int(par[6])
    if ax == 0:
        rad, h = math.sqrt(py * py + pz * pz), abs(px)
    elif ax == 1:
        rad, h = math.sqrt(px * px + pz * pz), abs(py)
    else:
        rad, h = math.sqrt(px * px + py * py), abs(pz)
    dr = rad - par[3]
    dh = h - par[4]
    o1 = max(dr, 0.0)
    o2 = max(dh, 0.0)
    return min(max(dr, dh), 0.0) + math.sqrt(o1 * o1 + o2 * o2)


@njit(cache=True)
def sdf_union(kinds, pars, x, y, z):
    best = 1e300
    who = -1
    for k in range(kinds.shape[0]):
        d = sdf_one(kinds[k], pars[k], x, y, z)
        if d < best:
            best = d
            who = k
    return best, who


@njit(cache=True)
def sdf_normal(kinds, pars, x, y, z, h):
    gx = sdf_union(kinds, pars, x + h, y, z)[0] - sdf_union(kinds, pars, x - h, y, z)[0]
    gy = sdf_union(kinds, pars, x, y + h, z)[0] - sdf_union(kinds, pars, x, y - h, z)[0]
    gz = sdf_union(kinds, pars, x, y, z + h)[0] - sdf_union(kinds, pars, x, y, z - h)[0]
    g = math.sqrt(gx * gx + gy * gy + gz * gz)
    if g == 0.0:
        return 0.0, 0.0, 1.0
    return gx / g, gy / g, gz / g


@njit(parallel=True, cache=True)
def apply_objects(pos, vel, q, state, absorb_t, absorb_q, hit_obj, t, kinds, pars, action, rest,
                  h, remove, neutralize):
    n = pos.shape[0]
    for i in prange(n):
        if state[i] != ACTIVE:
            continue
        x = pos[i, 0]
        y = pos[i, 1]
        z = pos[i, 2]
        d, k = sdf_union(kinds, pars, x, y, z)
        if d >= 0.0:
            continue
        nx, ny, nz = sdf_normal(kinds, pars, x, y, z, h)
        if action[k] == 0:  # reflect: mirror across the surface
            e = rest[k]
            pos[i, 0] = x - (1.0 + e) * d * nx
            pos[i, 1] = y - (1.0 + e) * d * ny
            pos[i, 2] = z - (1.0 + e) * d * nz
            vn = vel[i, 0] * nx + vel[i, 1] * ny + vel[i, 2] * nz
            if vn < 0.0:
                vel[i, 0] -= (1.0 + e) * vn * nx
                vel[i, 1] -= (1.0 + e) * vn * ny
                vel[i, 2] -= (1.0 + e) * vn * nz
        else:
            pos[i, 0] = x - d * nx
            pos[i, 1] = y - d * ny
            pos[i, 2] = z - d * nz
            hit_obj[i] = k
            _absorb(i, pos, vel, q, state, absorb_t, absorb_q, t, remove, neutralize)


@njit(parallel=True, cache=True)
def inside_objects(pts, kinds, pars):
    out = np.zeros(pts.shape[0], np.bool_)
    for i in prange(pts.shape[0]):
        out[i] = sdf_union(kinds, pars, pts[i, 0], pts[i, 1], pts[i, 2])[0] < 0.0
    return out


# ======================================================================
# contact collisions (Jacobi update: race-free and order independent)
# ======================================================================
@njit(parallel=True, cache=True)
def collide(pos, vel, m, rad, state, rest, lo, L, periodic, nc, start, order, dv):
    n = pos.shape[0]
    ncoll = 0
    for i in prange(n):
        dv[i, 0] = 0.0
        dv[i, 1] = 0.0
        dv[i, 2] = 0.0
        if state[i] == REMOVED or rad[i] <= 0.0:
            continue
        c = np.empty(3, np.int64)
        for d in range(3):
            s = (pos[i, d] - lo[d]) / L[d]
            if s >= 1.0 or s < 0.0:
                s = s - math.floor(s)
            c[d] = min(max(int(s * nc[d]), 0), nc[d] - 1)
        ax = np.empty(max(3, nc[0]), np.int64)
        ay = np.empty(max(3, nc[1]), np.int64)
        az = np.empty(max(3, nc[2]), np.int64)
        nx = _neigh_range(c[0], nc[0], periodic[0], 1, ax)
        ny = _neigh_range(c[1], nc[1], periodic[1], 1, ay)
        nz = _neigh_range(c[2], nc[2], periodic[2], 1, az)
        for a in range(nx):
            for b in range(ny):
                for e in range(nz):
                    cc = (ax[a] * nc[1] + ay[b]) * nc[2] + az[e]
                    for p in range(start[cc], start[cc + 1]):
                        j = order[p]
                        if j == i or rad[j] <= 0.0:
                            continue
                        dx = _min_image(pos[i, 0] - pos[j, 0], L[0], periodic[0])
                        dy = _min_image(pos[i, 1] - pos[j, 1], L[1], periodic[1])
                        dz = _min_image(pos[i, 2] - pos[j, 2], L[2], periodic[2])
                        r2 = dx * dx + dy * dy + dz * dz
                        s = rad[i] + rad[j]
                        if r2 >= s * s or r2 == 0.0:
                            continue
                        r = math.sqrt(r2)
                        ux = dx / r
                        uy = dy / r
                        uz = dz / r
                        vn = (vel[i, 0] - vel[j, 0]) * ux + (vel[i, 1] - vel[j, 1]) * uy + (vel[i, 2] - vel[j, 2]) * uz
                        if vn >= 0.0:
                            continue  # separating
                        # a stuck partner behaves as infinitely heavy
                        w = 1.0 if state[j] != ACTIVE else m[j] / (m[i] + m[j])
                        imp = -(1.0 + rest) * w * vn
                        dv[i, 0] += imp * ux
                        dv[i, 1] += imp * uy
                        dv[i, 2] += imp * uz
                        ncoll += 1
    for i in prange(n):
        if state[i] == ACTIVE:
            vel[i, 0] += dv[i, 0]
            vel[i, 1] += dv[i, 1]
            vel[i, 2] += dv[i, 2]
    return ncoll // 2
