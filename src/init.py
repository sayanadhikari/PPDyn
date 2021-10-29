from numba import jit
import numpy as np
import random

@jit(nopython=True)
def initial_periodic_radial(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,k,dumpPeriod,g,Q,M,Temp):
    random.seed(99999999)
    x  = np.empty(N, dtype=np.float64)
    y  = np.empty(N, dtype=np.float64)
    z  = np.empty(N, dtype=np.float64)
    vx = np.empty(N, dtype=np.float64)
    vy = np.empty(N, dtype=np.float64)
    vz = np.empty(N, dtype=np.float64)
    ux = np.empty(N, dtype=np.float64)
    uy = np.empty(N, dtype=np.float64)
    uz = np.empty(N, dtype=np.float64)
    ax = np.empty(N, dtype=np.float64)
    ay = np.empty(N, dtype=np.float64)
    az = np.empty(N, dtype=np.float64)
    fduration = np.zeros(N, dtype=np.float64)

    svx  = 0.0  # velocity sum correction term in X
    svy  = 0.0  # velocity sum correction term in Y
    svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,tmax,Nt)
    data_num = np.arange(start=0, stop=Nt, step=dumpPeriod, dtype=np.int64)

    for i in range(N):
        x[i] = (random.random())*2.0*Lx - Lx
        y[i] = (random.random())*2.0*Ly - Ly
        z[i] = (random.random())*2.0*Lz - Lz
        r = np.sqrt(x[i]*x[i] + y[i]*y[i]+z[i]*z[i])
        while r>Lx:
            x[i] = (random.random())*2.0*Lx - Lx
            y[i] = (random.random())*2.0*Ly - Ly
            z[i] = (random.random())*2.0*Lz - Lz
            r = np.sqrt(x[i]*x[i] + y[i]*y[i]+z[i]*z[i])

        vx[i] = (random.random())*Vxmax - Vxmax/2.0
        vy[i] = (random.random())*Vymax - Vymax/2.0
        vz[i] = (random.random())*Vzmax - Vzmax/2.0
        svx = svx + vx[i]
        svy = svy + vy[i]
        svz = svz + vz[i]

    for i in range(N):
        vx[i] = vx[i] - svx/N
        vy[i] = vy[i] - svy/N
        vz[i] = vz[i] - svz/N


    for i in range(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = 0.0
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] ) - round((x[i]-x[j])/(2.0*Lx)) * 2.0*Lx
                ydiff = ( y[i]-y[j] ) - round((y[i]-y[j])/(2.0*Ly)) * 2.0*Ly
                zdiff = ( z[i]-z[j] ) - round((z[i]-z[j])/(2.0*Lz)) * 2.0*Lz
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+k*r)*np.exp(-k*r)/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+k*r)*np.exp(-k*r)/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+k*r)*np.exp(-k*r)/(r*r*r) #+ zdiff*g + Lz*g  # zdiff/(r*r*r)
                ax[i] = ax[i] + fx
                ay[i] = ay[i] + fy
                az[i] = az[i] + fz
    return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num,fduration


@jit(nopython=True)
def initial_periodic(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,k,dumpPeriod,g,Q):
    random.seed(99999999)
    x  = np.empty(N, dtype=np.float64)
    y  = np.empty(N, dtype=np.float64)
    z  = np.empty(N, dtype=np.float64)
    vx = np.empty(N, dtype=np.float64)
    vy = np.empty(N, dtype=np.float64)
    vz = np.empty(N, dtype=np.float64)
    ux = np.empty(N, dtype=np.float64)
    uy = np.empty(N, dtype=np.float64)
    uz = np.empty(N, dtype=np.float64)
    ax = np.empty(N, dtype=np.float64)
    ay = np.empty(N, dtype=np.float64)
    az = np.empty(N, dtype=np.float64)

    svx  = 0.0  # velocity sum correction term in X
    svy  = 0.0  # velocity sum correction term in Y
    svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,tmax,Nt)
    data_num = np.arange(start=0, stop=Nt, step=dumpPeriod, dtype=np.int64)

    for i in range(N):
        x[i] = (random.random())*2.0*Lx - Lx
        y[i] = (random.random())*2.0*Ly - Ly
        z[i] = (random.random())*2.0*Lz - Lz
        vx[i] = (random.random())*Vxmax - Vxmax/2.0
        vy[i] = (random.random())*Vymax - Vymax/2.0
        vz[i] = (random.random())*Vzmax - Vzmax/2.0
        svx = svx + vx[i]
        svy = svy + vy[i]
        svz = svz + vz[i]

    for i in range(N):
        vx[i] = vx[i] - svx/N
        vy[i] = vy[i] - svy/N
        vz[i] = vz[i] - svz/N


    for i in range(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = 0.0
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] ) - round((x[i]-x[j])/(2.0*Lx)) * 2.0*Lx
                ydiff = ( y[i]-y[j] ) - round((y[i]-y[j])/(2.0*Ly)) * 2.0*Ly
                zdiff = ( z[i]-z[j] ) - round((z[i]-z[j])/(2.0*Lz)) * 2.0*Lz
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+k*r)*np.exp(-k*r)/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+k*r)*np.exp(-k*r)/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+k*r)*np.exp(-k*r)/(r*r*r) #+ zdiff*g + Lz*g  # zdiff/(r*r*r)
                ax[i] = ax[i] + fx
                ay[i] = ay[i] + fy
                az[i] = az[i] + fz
    return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num


@jit(nopython=True)
def initial_reflecting(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,k,dumpPeriod,g,Q,M,Temp):
    random.seed(99999999)
    x  = np.empty(N, dtype=np.float64)
    y  = np.empty(N, dtype=np.float64)
    z  = np.empty(N, dtype=np.float64)
    vx = np.empty(N, dtype=np.float64)
    vy = np.empty(N, dtype=np.float64)
    vz = np.empty(N, dtype=np.float64)
    ux = np.empty(N, dtype=np.float64)
    uy = np.empty(N, dtype=np.float64)
    uz = np.empty(N, dtype=np.float64)
    ax = np.empty(N, dtype=np.float64)
    ay = np.empty(N, dtype=np.float64)
    az = np.empty(N, dtype=np.float64)
    fduration = np.zeros(N, dtype=np.float64)

    svx  = 0.0  # velocity sum correction term in X
    svy  = 0.0  # velocity sum correction term in Y
    svz  = 0.0  # velocity sum correction term in Z

    ###### Initialize time array and data dump array ######
    time  = np.linspace(0,tmax,Nt)
    data_num = np.arange(start=0, stop=Nt, step=dumpPeriod, dtype=np.int64)

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
    x = np.random.random(N)*2.0*Lx - Lx
    y = np.random.random(N)*2.0*Ly - Ly
    z = np.zeros(N) #np.random.random(N)*2.0*Lz - Lz

    # Random velocity
    # vx = np.random.random(N)*Vxmax - Vxmax/2.0
    # vy = np.random.random(N)*Vymax - Vymax/2.0
    # vz = np.random.random(N)*Vzmax - Vzmax/2.0

    # Maxwellian
    vx = np.random.normal(0, Temp, N)
    vy = np.random.normal(0, Temp, N)
    vz = np.random.normal(0, Temp, N)

    # svx = svx + np.sum(vx)
    # svy = svy + np.sum(vy)
    # svz = svz + np.sum(vz)
    #
    # vx = vx - svx/N
    # vy = vy - svy/N
    # vz = vz - svz/N






    for i in range(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = -(z[i]+Lz)*g
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] )
                ydiff = ( y[i]-y[j] )
                zdiff = ( z[i]-z[j] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx = xdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy = ydiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz = zdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r) # + zdiff*g + Lz*g # zdiff/(r*r*r)
                ax[i] += fx/M[i]
                ay[i] += fy/M[i]
                az[i] += fz/M[i]
    return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num,fduration
