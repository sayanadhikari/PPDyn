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
    sv    = np.zeros((config.N,3), dtype=np.float64)


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
    if config.maxwell_load:
        vvel[:,0] = np.random.normal(0, config.Temp, config.N)
        vvel[:,1] = np.random.normal(0, config.Temp, config.N)
        vvel[:,2] = np.random.normal(0, config.Temp, config.N)
    # Random
    else:
        vvel[:,0] = np.random.random(config.N)*config.Vxmax - config.Vxmax/2.0
        vvel[:,1] = np.random.random(config.N)*config.Vymax - config.Vymax/2.0
        vvel[:,2] = np.random.random(config.N)*config.Vzmax - config.Vzmax/2.0

        sv[:,0] = sv[:,0] + vvel[:,0]
        sv[:,1] = sv[:,1] + vvel[:,1]
        sv[:,2] = sv[:,2] + vvel[:,2]

        vvel[:,0] = vvel[:,0] - sv[:,0]/config.N
        vvel[:,1] = vvel[:,1] - sv[:,1]/config.N
        vvel[:,2] = vvel[:,2] - sv[:,2]/config.N

    # acc = 0.0*acc

    for i in range(config.N):
        acc[i,:] = 0.0
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
    return pos,vvel,uvel,acc,time,data_num


@jit(nopython=True)
def initial_reflecting(Q,M):
    random.seed(99999999)
    pos   = np.empty((config.N,3), dtype=np.float64)
    uvel  = np.empty((config.N,3), dtype=np.float64)
    vvel  = np.empty((config.N,3), dtype=np.float64)
    acc   = np.empty((config.N,3), dtype=np.float64)
    sv    = np.zeros((config.N,3), dtype=np.float64)


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
    if config.maxwell_load:
        vvel[:,0] = np.random.normal(0, config.Temp, config.N)
        vvel[:,1] = np.random.normal(0, config.Temp, config.N)
        vvel[:,2] = np.random.normal(0, config.Temp, config.N)
    # Random
    else:
        vvel[:,0] = np.random.random(config.N)*config.Vxmax - config.Vxmax/2.0
        vvel[:,1] = np.random.random(config.N)*config.Vymax - config.Vymax/2.0
        vvel[:,2] = np.random.random(config.N)*config.Vzmax - config.Vzmax/2.0

        sv[:,0] = sv[:,0] + vvel[:,0]
        sv[:,1] = sv[:,1] + vvel[:,1]
        sv[:,2] = sv[:,2] + vvel[:,2]

        vvel[:,0] = vvel[:,0] - sv[:,0]/config.N
        vvel[:,1] = vvel[:,1] - sv[:,1]/config.N
        vvel[:,2] = vvel[:,2] - sv[:,2]/config.N

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
    return pos,vvel,uvel,acc,time,data_num
