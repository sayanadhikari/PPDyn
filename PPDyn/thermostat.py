from numba import jit
import numpy as np
@jit(nopython=True)
def berendsen(vx,vy,vz,dt,Temp,KE,N,t,tmax):
    tau = 10.0*dt
    scl = np.sqrt(1.0 + (dt/tau) * ((Temp/(2.0*KE/(3.0*float(N)) )) -1.0))

    if (t <= tmax/2.0):
        for i in range(N):
            vx[i] = scl * vx[i]
            vy[i] = scl * vy[i]
            vz[i] = scl * vz[i]
    else:
        vx[i] = vx[i]
        vy[i] = vy[i]
        vz[i] = vz[i]
    return vx,vy,vz
