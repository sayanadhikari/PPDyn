# import h5py
def configSpace(t,N,Nt,x,y,z,path):
    print('TimeSteps = %d'%int(t)+' of %d'%Nt)
    f = open(path+'data%d'%int(t)+'.dat', 'w')
    for i in range(N):
        f.write("%f"%x[i]+" %f"%y[i]+" %f"%z[i]+"\n")
    f.close()
    return 0
