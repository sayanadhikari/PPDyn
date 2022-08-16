from numba import jit
import numpy as np
import random

import config

@jit(nopython=True)
def initial_periodic(Q,M):
    random.seed(99999999)
    pos   = np.empty((config.N,3), dtype=np.float64)
    uvel  = np.empty((config.N,3), dtype=np.float64)
    vvel  = np.empty((config.N,3), dtype=np.float64)
    acc   = np.empty((config.N,3), dtype=np.float64)
    fduration = np.zeros(config.N, dtype=np.float32)


    # svx  = 0.0  # velocity sum correction term in X
    # svy  = 0.0  # velocity sum correction term in Y
    # svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,config.tmax,config.Nt)
    data_num = np.arange(start=0, stop=config.Nt, step=config.dumpPeriod, dtype=np.int64)

    pos[:,0] = np.random.random(config.N)*2.0*config.Lx - config.Lx
    pos[:,1] = np.random.random(config.N)*2.0*config.Ly - config.Ly
    pos[:,2] = np.random.random(config.N)*2.0*config.Lz - config.Lz

    # Maxwellian
    vvel[:,0] = np.random.normal(0, config.Temp, config.N)
    vvel[:,1] = np.random.normal(0, config.Temp, config.N)
    vvel[:,2] = np.random.normal(0, config.Temp, config.N)

    # for i in range(N):
    #     x[i] = (random.random())*2.0*Lx - Lx
    #     y[i] = (random.random())*2.0*Ly - Ly
    #     z[i] = (random.random())*2.0*Lz - Lz
    #     vx[i] = (random.random())*Vxmax - Vxmax/2.0
    #     vy[i] = (random.random())*Vymax - Vymax/2.0
    #     vz[i] = (random.random())*Vzmax - Vzmax/2.0
    #     svx = svx + vx[i]
    #     svy = svy + vy[i]
    #     svz = svz + vz[i]

    # for i in range(N):
    #     vx[i] = vx[i] - svx/N
    #     vy[i] = vy[i] - svy/N
    #     vz[i] = vz[i] - svz/N

    # acc = 0.0*acc

    for i in range(config.N):
        acc[i,:] = 0.0
        # acc[i,1] = 0.0
        # acc[i,2] = 0.0
        f0 = 10
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] ) - round((pos[i,0]-pos[j,0])/(2.0*config.Lx)) * 2.0*config.Lx
                ydiff = ( pos[i,1]-pos[j,1] ) - round((pos[i,1]-pos[j,1])/(2.0*config.Ly)) * 2.0*config.Ly
                zdiff = ( pos[i,2]-pos[j,2] ) - round((pos[i,2]-pos[j,2])/(2.0*config.Lz)) * 2.0*config.Lz
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r) #+ zdiff*g + Lz*g  # zdiff/(r*r*r)
                acc[i,0] += fx/M[i]
                acc[i,1] += fy/M[i]
                acc[i,2] += fz/M[i]
        r1 = np.sqrt((pos[i,0]*pos[i,0] )+ ( pos[i,1]*pos[i,1]) + (pos[i,2]*pos[i,2]))
        acc[i,0] += (f0*np.exp(- r1 /30))/M[i]            #lambda_c/lambda_d =30
        acc[i,1] += (f0*np.exp(- r1 /30))/M[i]
        acc[i,2] += (f0*np.exp(- r1 /30))/M[i]
    return pos,vvel,uvel,acc,time,data_num,fduration


@jit(nopython=True)
def initial_reflecting(Q,M):
    random.seed(99999999)
    pos   = np.empty((config.N,3), dtype=np.float64)
    uvel  = np.empty((config.N,3), dtype=np.float64)
    vvel  = np.empty((config.N,3), dtype=np.float64)
    acc   = np.empty((config.N,3), dtype=np.float64)
    fduration = np.zeros(config.N, dtype=np.float32)

    # svx  = 0.0  # velocity sum correction term in X
    # svy  = 0.0  # velocity sum correction term in Y
    # svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,config.tmax,config.Nt)
    data_num = np.arange(start=0, stop=config.Nt, step=config.dumpPeriod, dtype=np.int32)

    # for i in range(N):
    #     x[i] = (random.random())*2.0*Lx - Lx
    #     y[i] = (random.random())*2.0*Ly - Ly
    #     z[i] =  Lz #(random.random())*2.0*Lz - Lz
    #     vx[i] = (random.random())*Vxmax - Vxmax/2.0
    #     vy[i] = (random.random())*Vymax - Vymax/2.0
    #     vz[i] = (random.random())*Vzmax - Vzmax/2.0
    #     svx = svx + vx[i]
    #     svy = svy + vy[i]
    #     svz = svz + vz[i]
    pos[:,0] = np.random.random(config.N)*2.0*config.Lx - config.Lx
    pos[:,1] = np.random.random(config.N)*2.0*config.Ly - config.Ly
    pos[:,2] = np.ones(config.N)*config.Lz #np.random.random(N)*2.0*Lz - Lz

    # Random velocity
    # vx = np.random.random(N)*Vxmax - Vxmax/2.0
    # vy = np.random.random(N)*Vymax - Vymax/2.0
    # vz = np.random.random(N)*Vzmax - Vzmax/2.0

    # Maxwellian
    vvel[:,0] = np.random.normal(0, config.Temp, config.N)
    vvel[:,1] = np.random.normal(0, config.Temp, config.N)
    vvel[:,2] = np.random.normal(0, config.Temp, config.N)

    # svx = svx + np.sum(vx)
    # svy = svy + np.sum(vy)
    # svz = svz + np.sum(vz)
    #
    # vx = vx - svx/N
    # vy = vy - svy/N
    # vz = vz - svz/N

    for i in range(config.N):
        acc[i,0] = 0.0
        acc[i,1] = 0.0
        acc[i,2] = -(pos[i,2]+config.Lz)*config.g
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] )
                ydiff = ( pos[i,1]-pos[j,1] )
                zdiff = ( pos[i,2]-pos[j,2] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+config.k*r)*np.exp(-config.k*r)*(Q[i]*Q[j])/(r*r*r) # + zdiff*g + Lz*g # zdiff/(r*r*r)
                acc[i,0] += fx/M[i]
                acc[i,1] += fy/M[i]
                acc[i,2] += fz/M[i]

    return pos,vvel,uvel,acc,time,data_num,fduration
