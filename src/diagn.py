import numpy as np
from os.path import join as pjoin
def configSpace(f,dsetE,dsetQ,t,x,y,z,vx,vy,vz,KE,Qcollect,path):

    f["/%d"%int(t)+"/position/x"] = x
    f["/%d"%int(t)+"/position/y"] = y
    f["/%d"%int(t)+"/position/z"] = z
    f["/%d"%int(t)+"/velocity/vx"] = vx
    f["/%d"%int(t)+"/velocity/vy"] = vy
    f["/%d"%int(t)+"/velocity/vz"] = vz
    # f["/%d"%int(t)+"/energy"] = KE

    dsetE.resize(dsetE.shape[0]+1, axis=0)
    dsetE[-1:] = KE

    dsetQ.resize(dsetQ.shape[0]+1, axis=0)
    dsetQ[-1:] = Qcollect

    with open(pjoin(path,'energy.txt'),"ab") as f:
        np.savetxt(f, np.column_stack([t, KE]))

    return 0

def attributes(f,tmax,Lx,Ly,Lz,N,dt,dumpPeriod):
    f.attrs["dp"] = dumpPeriod
    f.attrs["dt"] = dt
    f.attrs["Nt"] = int(tmax/(dt*dumpPeriod))
    f.attrs["N"]  = N
    f.attrs["tmax"] = tmax
    f.attrs["Lx"] = Lx
    f.attrs["Ly"] = Ly
    f.attrs["Lz"] = Lz
    return 0

def dustDiagn(f,fduration):

    f["/fall_duration"] = fduration
    # f["/%d"%int(t)+"/energy"] = KE
    return 0
