import numpy as np
from os.path import join as pjoin
import config

def configSpace(f,dsetE,t,pos,vvel,KE):

    f["/%d"%int(t)+"/position/x"] = pos[:,0]
    f["/%d"%int(t)+"/position/y"] = pos[:,1]
    f["/%d"%int(t)+"/position/z"] = pos[:,2]
    f["/%d"%int(t)+"/velocity/vx"] = vvel[:,0]
    f["/%d"%int(t)+"/velocity/vy"] = vvel[:,1]
    f["/%d"%int(t)+"/velocity/vz"] = vvel[:,2]
    # f["/%d"%int(t)+"/energy"] = KE

    dsetE.resize(dsetE.shape[0]+1, axis=0)
    dsetE[-1:] = KE

    with open(pjoin(config.dataDir,'energy.txt'),"ab") as f:
        np.savetxt(f, np.column_stack([t, KE]))

    return 0

def attributes(f):
    f.attrs["dp"] = config.dumpPeriod
    f.attrs["dt"] = config.dt
    f.attrs["Nt"] = int(config.tmax/(config.dt*config.dumpPeriod))
    f.attrs["N"]  = config.N
    f.attrs["tmax"] = config.tmax
    f.attrs["Lx"] = config.Lx
    f.attrs["Ly"] = config.Ly
    f.attrs["Lz"] = config.Lz
    return 0
