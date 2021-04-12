import numpy as np
def configSpace(f,dset,t,x,y,z,KE):

    f["/%d"%int(t)+"/position/x"] = x
    f["/%d"%int(t)+"/position/y"] = y
    f["/%d"%int(t)+"/position/z"] = z
    # f["/%d"%int(t)+"/energy"] = KE
    dset.resize(dset.shape[0]+1, axis=0)
    dset[-1:] = KE

    with open("data/energy.txt","ab") as f:
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
