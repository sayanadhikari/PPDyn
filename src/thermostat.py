from numba import jit
import numpy as np
import config

@jit(nopython=True)
def berendsen(t,vvel,KE):
    tau = 10.0*config.dt
    scl = np.sqrt(1.0 + (config.dt/tau) * ((config.Temp/(2.0*KE/(3.0*float(config.N)) )) -1.0))

    if (t <= config.tmax/2.0):
        vvel[:,0] = scl * vvel[:,0]
        vvel[:,1] = scl * vvel[:,1]
        vvel[:,2] = scl * vvel[:,2]
    else:
        None
    return vvel
