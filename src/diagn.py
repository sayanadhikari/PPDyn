# import h5py
def configSpace(t,N,tmax,x,y,z,Lx,Ly,Lz,f,dt,dumpPeriod):
    # fp = open(path+'data%d'%int(t)+'.dat', 'w')
    # for i in range(N):
    #     fp.write("%f"%x[i]+" %f"%y[i]+" %f"%z[i]+"\n")
    # fp.close()
    f["/%d"%int(t)+"/position/x"] = x
    f["/%d"%int(t)+"/position/y"] = y
    f["/%d"%int(t)+"/position/z"] = z
    f.attrs["dp"] = dumpPeriod
    f.attrs["dt"] = dt
    f.attrs["Nt"] = int(tmax/(dt*dumpPeriod))
    f.attrs["tmax"] = tmax
    f.attrs["Lx"] = Lx
    f.attrs["Ly"] = Ly
    f.attrs["Lz"] = Lz
    return 0
