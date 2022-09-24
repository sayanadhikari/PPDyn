import numpy as np
from os.path import join as pjoin
import config

def configSpace(dsetE,dsetPart,dsetVel,t,pos,vvel,KE):

    dsetE.resize(dsetE.shape[0]+1, axis=0)
    dsetE[-1:] = KE
    dsetPart[int(t/config.dumpPeriod),:,:] = pos
    dsetVel[int(t/config.dumpPeriod),:,:]  = vvel

    with open(pjoin(config.dataDir,'energy.txt'),"ab") as fenergy:
        np.savetxt(fenergy, np.column_stack([t, KE]))

    return 0

def attributes(f):
    f.attrs["dp"] = config.dumpPeriod
    f.attrs["dt"] = config.dt
    f.attrs["Nt"] = int(config.tmax/(config.dt))
    f.attrs["N"]  = config.N
    f.attrs["tmax"] = config.tmax
    f.attrs["Lx"] = config.Lx
    f.attrs["Ly"] = config.Ly
    f.attrs["Lz"] = config.Lz
    return 0
