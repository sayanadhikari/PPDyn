from numba import jit
import numpy as np

@jit(nopython=True)
def half_t1(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N):
    for i in range(N):
        ux[i] = vx[i] + ax[i] * dt/2.0
        uy[i] = vy[i] + ay[i] * dt/2.0
        uz[i] = vz[i] + az[i] * dt/2.0
        x[i] = x[i] + ux[i] * dt
        y[i] = y[i] + uy[i] * dt
        z[i] = z[i] + uz[i] * dt
        x[i] = x[i] - (int(x[i]/Lx)) * 2.0 * Lx      # Periodic Boundary Condition
        y[i] = y[i] - (int(y[i]/Ly)) * 2.0 * Ly      # Periodic Boundary Condition
        z[i] = z[i] - (int(z[i]/Lz)) * 2.0 * Lz      # Periodic Boundary Condition
    return x,y,z,ux,uy,uz

@jit(nopython=True)
def full(x,y,z,ax,ay,az,Lx,Ly,Lz,N):
    for i in range(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = 0.0
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] ) - np.round((x[i]-x[j])/(2.0*Lx)) * 2.0*Lx
                ydiff = ( y[i]-y[j] ) - np.round((y[i]-y[j])/(2.0*Ly)) * 2.0*Ly
                zdiff = ( z[i]-z[j] ) - np.round((z[i]-z[j])/(2.0*Lz)) * 2.0*Lz
                r = np.sqrt(np.square(xdiff*xdiff) + np.square(ydiff*ydiff) + np.square(zdiff*zdiff))
                fx = xdiff/(np.power(r,3))
                fy = ydiff/(np.power(r,3))
                fz = zdiff/(np.power(r,3))
                ax[i] += fx
                ay[i] += fy
                az[i] += fz
    return ax,ay,az

@jit(nopython=True)
def half_t2(ux,uy,uz,vx,vy,vz,ax,ay,az,dt,N,KE):
    for i in range(N):
        vx[i] = ux[i] + ax[i] * dt / 2.0
        vy[i] = uy[i] + ay[i] * dt / 2.0
        vz[i] = uz[i] + az[i] * dt / 2.0
        KE += ( np.square(vx[i]) + np.square(vy[i]) + np.square(vz[i]) ) / 2.0
    return vx,vy,vz,KE
