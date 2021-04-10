# import h5py
def configSpace(f,dset,t,x,y,z,KE):

    f["/%d"%int(t)+"/position/x"] = x
    f["/%d"%int(t)+"/position/y"] = y
    f["/%d"%int(t)+"/position/z"] = z
    # f["/%d"%int(t)+"/energy"] = KE
    dset.resize(dset.shape[0]+1, axis=0)
    dset[-1:] = KE
    return 0

def attributes(f,tmax,Lx,Ly,Lz,dt,dumpPeriod):
    f.attrs["dp"] = dumpPeriod
    f.attrs["dt"] = dt
    f.attrs["Nt"] = int(tmax/(dt*dumpPeriod))
    f.attrs["tmax"] = tmax
    f.attrs["Lx"] = Lx
    f.attrs["Ly"] = Ly
    f.attrs["Lz"] = Lz
    return 0
