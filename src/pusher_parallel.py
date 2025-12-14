from numba import jit, prange
import numpy as np
import h5py
import config
import sys

# ------------------------------------------------------------------
# Trilinear interpolator for E-field (JIT-accelerated)
# ------------------------------------------------------------------
@jit(nopython=True)
def interp_E(px, py, pz, xg, yg, zg, Ex, Ey, Ez):
    # enforce periodicity: wrap px, py, pz into [xg[0], xg[-1]), etc.
    width_x = xg[-1] - xg[0]
    width_y = yg[-1] - yg[0]
    width_z = zg[-1] - zg[0]
    px = (px - xg[0]) % width_x + xg[0]
    py = (py - yg[0]) % width_y + yg[0]
    pz = (pz - zg[0]) % width_z + zg[0]
    nx, ny, nz = len(xg), len(yg), len(zg)
    # locate base cell indices
    i = np.searchsorted(xg, px) - 1
    j = np.searchsorted(yg, py) - 1
    k = np.searchsorted(zg, pz) - 1
    # clamp to valid range
    if i < 0:
        i = 0
    elif i > nx - 2:
        i = nx - 2
    if j < 0:
        j = 0
    elif j > ny - 2:
        j = ny - 2
    if k < 0:
        k = 0
    elif k > nz - 2:
        k = nz - 2
    # fractions
    xd = (px - xg[i]) / (xg[i + 1] - xg[i])
    yd = (py - yg[j]) / (yg[j + 1] - yg[j])
    zd = (pz - zg[k]) / (zg[k + 1] - zg[k])

    # trilinear interpolation
    Epx = 0.0
    Epy = 0.0
    Epz = 0.0
    for di in (0, 1):
        wx = (1 - xd) if di == 0 else xd
        for dj in (0, 1):
            wy = (1 - yd) if dj == 0 else yd
            for dk in (0, 1):
                wz = (1 - zd) if dk == 0 else zd
                w = wx * wy * wz
                idx = i + di
                jdx = j + dj
                kdx = k + dk
                Epx += w * Ex[idx, jdx, kdx]
                Epy += w * Ey[idx, jdx, kdx]
                Epz += w * Ez[idx, jdx, kdx]
    return Epx, Epy, Epz


# ------------------------------------------------------------------
# Modified Verlet with PIC E-field coupling (periodic)
# ------------------------------------------------------------------
@jit(nopython=True, parallel=True)
def verlet_periodic(t, pos, vvel, uvel, acc, Q, M,
                    xg, yg, zg, Ex, Ey, Ez,
                    xg_obj, yg_obj, zg_obj, obj_mask,
                    particle_status, impact_heatmap):
    # half-kick
    for i in prange(config.N):
        uvel[i, :] = vvel[i, :] + acc[i, :] * config.dt / 2.0

    # drift + periodic wrap
    for i in prange(config.N):
        pos[i, :] = pos[i, :] + uvel[i, :] * config.dt
        pos[i, 0] -= int(pos[i, 0] / config.Lx) * 2.0 * config.Lx
        pos[i, 1] -= int(pos[i, 1] / config.Ly) * 2.0 * config.Ly
        pos[i, 2] -= int(pos[i, 2] / config.Lz) * 2.0 * config.Lz

        # object boundary interaction (absorb/reflect/attach)
        # Only process if particle is still active (status == 0)
        if particle_status[i] == 0:
            # map pos back to object-grid indices
            dx = xg_obj[1] - xg_obj[0]
            dy = yg_obj[1] - yg_obj[0]
            dz = zg_obj[1] - zg_obj[0]
            ix = int(round((pos[i,0] - xg_obj[0]) / dx))
            iy = int(round((pos[i,1] - yg_obj[0]) / dy))
            iz = int(round((pos[i,2] - zg_obj[0]) / dz))
            # clamp indices
            ix = max(0, min(ix, obj_mask.shape[0]-1))
            iy = max(0, min(iy, obj_mask.shape[1]-1))
            iz = max(0, min(iz, obj_mask.shape[2]-1))
            if obj_mask[ix, iy, iz] == 1:
                # Record impact in heatmap (atomic-like increment using simple addition)
                # Note: In parallel mode, this is approximate but acceptable for heatmap
                impact_heatmap[ix, iy, iz] += 1.0

                # Determine interaction type using deterministic hash-based selection
                # This avoids random number generation issues in parallel Numba
                hash_val = ((i * 73856093) ^ (t * 19349663) ^ (ix * 83492791)) % 1000
                rand_val = (hash_val / 1000.0)

                # Apply interaction probabilities
                if rand_val < config.absorb_prob:
                    # Absorb: mark particle as inactive
                    particle_status[i] = 1
                    # Move particle far away (effectively remove from simulation)
                    # Move absorbed particle to a position outside box that will wrap with periodic BC
                    # Use a position that's clearly outside but will be wrapped by periodic boundary
                    pos[i, 0] = config.Lx * 2.0 + 1.0
                    pos[i, 1] = config.Ly * 2.0 + 1.0
                    pos[i, 2] = config.Lz * 2.0 + 1.0
                    vvel[i, :] = 0.0
                    uvel[i, :] = 0.0
                elif rand_val < (config.absorb_prob + config.reflect_prob):
                    # Reflect: bounce particle back
                    pos[i, :] -= uvel[i, :] * config.dt
                    uvel[i, :] = -uvel[i, :]
                else:
                    # Attach: stop particle at surface
                    particle_status[i] = 2
                    pos[i, :] -= uvel[i, :] * config.dt
                    vvel[i, :] = 0.0
                    uvel[i, :] = 0.0

    # reset forces
    for i in prange(config.N):
        acc[i, 0] = 0.0
        acc[i, 1] = 0.0
        acc[i, 2] = 0.0

    # Coulomb / Yukawa forces
    for i in prange(config.N):
        for j in range(config.N):
            if i == j:
                continue
            # periodic pairwise displacement
            xdiff = (pos[i, 0] - pos[j, 0]) - round(
                (pos[i, 0] - pos[j, 0]) / (2.0 * config.Lx)
            ) * 2.0 * config.Lx
            ydiff = (pos[i, 1] - pos[j, 1]) - round(
                (pos[i, 1] - pos[j, 1]) / (2.0 * config.Ly)
            ) * 2.0 * config.Ly
            zdiff = (pos[i, 2] - pos[j, 2]) - round(
                (pos[i, 2] - pos[j, 2]) / (2.0 * config.Lz)
            ) * 2.0 * config.Lz
            r = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
            coef = (
                (1 + config.k * r) * np.exp(-config.k * r) * (Q[i] * Q[j]) / (r * r * r)
            )
            acc[i, 0] += xdiff * coef / M[i]
            acc[i, 1] += ydiff * coef / M[i]
            acc[i, 2] += zdiff * coef / M[i]

    # PIC-field coupling
    for i in prange(config.N):
        Ex_p, Ey_p, Ez_p = interp_E(
            pos[i, 0], pos[i, 1], pos[i, 2], xg, yg, zg, Ex, Ey, Ez
        )
        # Safety check: prevent NaN/Inf and clip extremely large values
        if np.isnan(Ex_p) or np.isinf(Ex_p):
            Ex_p = 0.0
        if np.isnan(Ey_p) or np.isinf(Ey_p):
            Ey_p = 0.0
        if np.isnan(Ez_p) or np.isinf(Ez_p):
            Ez_p = 0.0
        # Clip extremely large field values to prevent numerical overflow
        E_mag_p = np.sqrt(Ex_p*Ex_p + Ey_p*Ey_p + Ez_p*Ez_p)
        max_E_safe = 1e6  # 1 MV/m safety limit
        if E_mag_p > max_E_safe:
            scale = max_E_safe / E_mag_p
            Ex_p *= scale
            Ey_p *= scale
            Ez_p *= scale

        acc[i, 0] += Q[i] * Ex_p / M[i]
        acc[i, 1] += Q[i] * Ey_p / M[i]
        acc[i, 2] += Q[i] * Ez_p / M[i]

    # final half-kick
    for i in prange(config.N):
        vvel[i, :] = uvel[i, :] + acc[i, :] * config.dt / 2.0

    # hard-sphere collisions (periodic)
    for i in prange(config.N):
        for j in range(config.N):
            if i == j:
                continue
            # periodic pairwise displacement for collision
            xdiff = (pos[i,0] - pos[j,0]) - round((pos[i,0] - pos[j,0])/(2.0*config.Lx)) * 2.0*config.Lx
            ydiff = (pos[i,1] - pos[j,1]) - round((pos[i,1] - pos[j,1])/(2.0*config.Ly)) * 2.0*config.Ly
            zdiff = (pos[i,2] - pos[j,2]) - round((pos[i,2] - pos[j,2])/(2.0*config.Lz)) * 2.0*config.Lz
            # velocity differences
            vxdiff = vvel[i,0] - vvel[j,0]
            vydiff = vvel[i,1] - vvel[j,1]
            vzdiff = vvel[i,2] - vvel[j,2]
            r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
            if r < 2*config.rc:
                # elastic collision update
                mu_i = 2 * M[j] / (M[i] + M[j])
                factor = (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) / (r*r)
                vvel[i,0] -= mu_i * factor * xdiff
                vvel[i,1] -= mu_i * factor * ydiff
                vvel[i,2] -= mu_i * factor * zdiff
                mu_j = 2 * M[i] / (M[i] + M[j])
                vvel[j,0] += mu_j * factor * xdiff
                vvel[j,1] += mu_j * factor * ydiff
                vvel[j,2] += mu_j * factor * zdiff

    return pos, vvel, uvel, acc, Q, particle_status, impact_heatmap


# ------------------------------------------------------------------
# Modified Verlet-reflecting with PIC E-field coupling
# ------------------------------------------------------------------
@jit(nopython=True, parallel=True)
def verlet_reflecting(t, pos, vvel, uvel, acc, Q, M,
                      xg, yg, zg, Ex, Ey, Ez,
                      xg_obj, yg_obj, zg_obj, obj_mask,
                    particle_status, impact_heatmap):
    # half-kick
    for i in prange(config.N):
        uvel[i, :] = vvel[i, :] + acc[i, :] * config.dt / 2.0
        pos[i, :] = pos[i, :] + uvel[i, :] * config.dt
        # reflecting BC
        if pos[i, 0] > config.Lx or pos[i, 0] < -config.Lx:
            pos[i, 0] -= uvel[i, 0] * config.dt
            uvel[i, 0] = -uvel[i, 0]
        if pos[i, 1] > config.Ly or pos[i, 1] < -config.Ly:
            pos[i, 1] -= uvel[i, 1] * config.dt
            uvel[i, 1] = -uvel[i, 1]
        if pos[i, 2] > config.Lz:
            pos[i, 2] -= uvel[i, 2] * config.dt
            uvel[i, 2] = -uvel[i, 2]

        # object boundary interaction (absorb/reflect/attach)
        # Only process if particle is still active (status == 0)
        if particle_status[i] == 0:
            # map pos back to object-grid indices
            dx = xg_obj[1] - xg_obj[0]
            dy = yg_obj[1] - yg_obj[0]
            dz = zg_obj[1] - zg_obj[0]
            ix = int(round((pos[i,0] - xg_obj[0]) / dx))
            iy = int(round((pos[i,1] - yg_obj[0]) / dy))
            iz = int(round((pos[i,2] - zg_obj[0]) / dz))
            # clamp indices
            ix = max(0, min(ix, obj_mask.shape[0]-1))
            iy = max(0, min(iy, obj_mask.shape[1]-1))
            iz = max(0, min(iz, obj_mask.shape[2]-1))
            if obj_mask[ix, iy, iz] == 1:
                # Record impact in heatmap (atomic-like increment using simple addition)
                # Note: In parallel mode, this is approximate but acceptable for heatmap
                impact_heatmap[ix, iy, iz] += 1.0

                # Determine interaction type using deterministic hash-based selection
                # This avoids random number generation issues in parallel Numba
                hash_val = ((i * 73856093) ^ (t * 19349663) ^ (ix * 83492791)) % 1000
                rand_val = (hash_val / 1000.0)

                # Apply interaction probabilities
                if rand_val < config.absorb_prob:
                    # Absorb: mark particle as inactive
                    particle_status[i] = 1
                    # Move particle far away (effectively remove from simulation)
                    # Move absorbed particle to a position outside box that will wrap with periodic BC
                    # Use a position that's clearly outside but will be wrapped by periodic boundary
                    pos[i, 0] = config.Lx * 2.0 + 1.0
                    pos[i, 1] = config.Ly * 2.0 + 1.0
                    pos[i, 2] = config.Lz * 2.0 + 1.0
                    vvel[i, :] = 0.0
                    uvel[i, :] = 0.0
                elif rand_val < (config.absorb_prob + config.reflect_prob):
                    # Reflect: bounce particle back
                    pos[i, :] -= uvel[i, :] * config.dt
                    uvel[i, :] = -uvel[i, :]
                else:
                    # Attach: stop particle at surface
                    particle_status[i] = 2
                    pos[i, :] -= uvel[i, :] * config.dt
                    vvel[i, :] = 0.0
                    uvel[i, :] = 0.0

    # reset forces + gravity
    for i in prange(config.N):
        acc[i, 0] = 0.0
        acc[i, 1] = 0.0
        acc[i, 2] = -(pos[i, 2] + config.Lz) * config.g

    # Coulomb / Yukawa forces
    for i in prange(config.N):
        for j in range(config.N):
            if i == j:
                continue
            xdiff = pos[i, 0] - pos[j, 0]
            ydiff = pos[i, 1] - pos[j, 1]
            zdiff = pos[i, 2] - pos[j, 2]
            r = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
            coef = (
                (1 + config.k * r) * np.exp(-config.k * r) * (Q[i] * Q[j]) / (r * r * r)
            )
            acc[i, 0] += xdiff * coef / M[i]
            acc[i, 1] += ydiff * coef / M[i]
            acc[i, 2] += zdiff * coef / M[i]

    # PIC-field coupling
    for i in prange(config.N):
        Ex_p, Ey_p, Ez_p = interp_E(
            pos[i, 0], pos[i, 1], pos[i, 2], xg, yg, zg, Ex, Ey, Ez
        )
        acc[i, 0] += Q[i] * Ex_p / M[i]
        acc[i, 1] += Q[i] * Ey_p / M[i]
        acc[i, 2] += Q[i] * Ez_p / M[i]

    # final half-kick + reflecting collisions & KE
    for i in prange(config.N):
        vvel[i, :] = uvel[i, :] + acc[i, :] * config.dt / 2.0

    # optional hard-sphere collisions (unchanged)
    for i in prange(config.N):
        for j in range(config.N):
            if i == j:
                continue
            xdiff = pos[i, 0] - pos[j, 0]
            ydiff = pos[i, 1] - pos[j, 1]
            zdiff = pos[i, 2] - pos[j, 2]
            vxdiff = vvel[i, 0] - vvel[j, 0]
            vydiff = vvel[i, 1] - vvel[j, 1]
            vzdiff = vvel[i, 2] - vvel[j, 2]
            r = np.sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff)
            if r < 2 * config.rc:
                # elastic collision update
                mu = 2 * M[j] / (M[i] + M[j])
                factor = (vxdiff * xdiff + vydiff * ydiff + vzdiff * zdiff) / (r * r)
                vvel[i, 0] -= mu * factor * xdiff
                vvel[i, 1] -= mu * factor * ydiff
                vvel[i, 2] -= mu * factor * zdiff
                # symmetric for j
                mu_j = 2 * M[i] / (M[i] + M[j])
                vvel[j, 0] += mu_j * factor * xdiff
                vvel[j, 1] += mu_j * factor * ydiff
                vvel[j, 2] += mu_j * factor * zdiff

    # return without KE update
    return pos, vvel, uvel, acc, Q, particle_status, impact_heatmap
