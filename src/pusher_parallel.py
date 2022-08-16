from numba import jit, prange
import numpy as np
import config

@jit(nopython=True, parallel=True)
def verlet_periodic(t,pos,vvel,uvel,acc,Q,M,KE,fduration,Qcollect):
    for i in prange(config.N):
        uvel[i,:] = vvel[i,:] + acc[i,:] * config.dt/2.0
        # uvel[i,0] = vvel[i,0] + acc[i,0] * config.dt/2.0
        # uvel[i,1] = vvel[i,1] + acc[i,1] * config.dt/2.0
        # uvel[i,2] = vvel[i,2] + acc[i,2] * config.dt/2.0
        pos[i,0] = pos[i,0] + uvel[i,0] * config.dt
        pos[i,1] = pos[i,1] + uvel[i,1] * config.dt
        pos[i,2] = pos[i,2] + uvel[i,2] * config.dt
        # r1 = np.sqrt((pos[i,0]*pos[i,0] )+ ( pos[i,1]*pos[i,1]) + (pos[i,2]*pos[i,2]))
        # pos[i,0] -= (int(r1/config.Lx))*pos[i,0]*2.0    # Periodic Boundary Condition
        # pos[i,1] -= (int(r1/config.Ly))*pos[i,1]*2.0    # Periodic Boundary Condition
        # pos[i,2] -= (int(r1/config.Lz))*pos[i,2]*2.0    # Periodic Boundary Condition
        pos[i,0] = pos[i,0] - (int(pos[i,0]/config.Lx)) * 2.0 * config.Lx      # Periodic Boundary Condition
        pos[i,1] = pos[i,1] - (int(pos[i,1]/config.Ly)) * 2.0 * config.Ly      # Periodic Boundary Condition
        pos[i,2] = pos[i,2] - (int(pos[i,2]/config.Lz)) * 2.0 * config.Lz      # Periodic Boundary Condition

    for i in prange(config.N):
        acc[i,:] = 0.0
        # acc[i,0] = 0.0
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
        acc[i,0] += ((f0*np.exp(- r1 /30))*(pos[i,0]/r1))/M[i]            #lambda_c/lambda_d =30
        acc[i,1] += ((f0*np.exp(- r1 /30))*(pos[i,1]/r1))/M[i]
        acc[i,2] += ((f0*np.exp(- r1 /30))*(pos[i,2]/r1))/M[i] 

    for i in prange(config.N):
        vvel[i,:] = uvel[i,:] + acc[i,:] * config.dt / 2.0
        # vvel[i,0] = uvel[i,0] + acc[i,0] * config.dt / 2.0
        # vvel[i,1] = uvel[i,1] + acc[i,1] * config.dt / 2.0
        # vvel[i,2] = uvel[i,2] + acc[i,2] * config.dt / 2.0
        KE += ((vvel[i,0]*vvel[i,0]) + (vvel[i,1]*vvel[i,1]) + (vvel[i,2]*vvel[i,2]) ) / 2.0
    return pos,vvel,uvel,acc,Q,KE,fduration,Qcollect

@jit(nopython=True, parallel=True)
def verlet_reflecting(t,pos,vvel,uvel,acc,Q,M,KE,fduration,Qcollect):
    for i in prange(config.N):
        uvel[i,:] = vvel[i,:] + acc[i,:] * config.dt/2.0
        # uvel[i,0] = vvel[i,0] + acc[i,0] * config.dt/2.0
        # uvel[i,1] = vvel[i,1] + acc[i,1] * config.dt/2.0
        # uvel[i,2] = vvel[i,2] + acc[i,2] * config.dt/2.0
        if (pos[i,2] <= -config.Lz):
            pos[i,2] = -config.Lz
            uvel[i,:] = 0.0
            # uvel[i,2] = 0.0
            # uvel[i,1] = 0.0
            # uvel[i,0] = 0.0
            if Q[i] > 0.0:
                fduration[i] = t
                Qcollect += Q[i]
            # print(i,fduration[i],t)
            Q[i]  = 0.0      #Make charges zero as they hit the ground
        pos[i,:] = pos[i,:] + uvel[i,:] * config.dt
        # pos[i,0] = pos[i,0] + uvel[i,0] * config.dt
        # pos[i,1] = pos[i,1] + uvel[i,1] * config.dt
        # pos[i,2] = pos[i,2] + uvel[i,2] * config.dt
        # Reflecting boundary
        if (pos[i,0] > config.Lx or pos[i,0] < -config.Lx):
            pos[i,0] -= uvel[i,0] * config.dt
            uvel[i,0] = -uvel[i,0]
        if (pos[i,1] > config.Ly or pos[i,1] < -config.Ly):
            pos[i,1] -= uvel[i,1] * config.dt
            uvel[i,1] = -uvel[i,1]
        if (pos[i,2] > config.Lz):
            pos[i,2] -= uvel[i,2] * config.dt
            uvel[i,2] = -uvel[i,2]

    for i in prange(config.N):
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

    for i in prange(config.N):
        vvel[i,:] = uvel[i,:] + acc[i,:] * config.dt / 2.0
        # vvel[i,0] = uvel[i,0] + acc[i,0] * config.dt / 2.0
        # vvel[i,1] = uvel[i,1] + acc[i,1] * config.dt / 2.0
        # vvel[i,2] = uvel[i,2] + acc[i,2] * config.dt / 2.0


    for i in prange(config.N):
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] )
                ydiff = ( pos[i,1]-pos[j,1] )
                zdiff = ( pos[i,2]-pos[j,2] )
                vxdiff = ( vvel[i,0]-vvel[j,0] )
                vydiff = ( vvel[i,1]-vvel[j,1] )
                vzdiff = ( vvel[i,2]-vvel[j,2] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                if (r < 2*config.rc):
                    vvel[i,0] = vvel[i,0] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vvel[i,1] = vvel[i,1] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vvel[i,2] = vvel[i,2] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)
                    vvel[j,0] = vvel[j,0] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vvel[j,1] = vvel[j,1] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vvel[j,2] = vvel[j,2] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)

    for i in prange(config.N):
        KE += ((vvel[i,0]*vvel[i,0]) + (vvel[i,1]*vvel[i,1]) + (vvel[i,2]*vvel[i,2]) ) / 2.0

    return pos,vvel,uvel,acc,Q,KE,fduration,Qcollect


@jit(nopython=True, parallel=True)
def verlet_mixed(t,pos,vvel,uvel,acc,Q,M,KE,fduration,Qcollect):
    for i in prange(config.N):
        uvel[i,:] = vvel[i,:] + acc[i,:] * config.dt/2.0
        if (pos[i,2] <= -config.Lz):
            pos[i,2] = -config.Lz
            uvel[i,:] = 0.0
            # uvel[i,2] = 0.0
            # uvel[i,1] = 0.0
            # uvel[i,0] = 0.0
            if Q[i] > 0.0:
                fduration[i] = t
                Qcollect += Q[i]
            # print(i,fduration[i],t)
            Q[i]  = 0.0      #Make charges zero as they hit the ground
        pos[i,:] = pos[i,:] + uvel[i,:] * config.dt
        # pos[i,0] = pos[i,0] + uvel[i,0] * config.dt
        # pos[i,1] = pos[i,1] + uvel[i,1] * config.dt
        # pos[i,2] = pos[i,2] + uvel[i,2] * config.dt

        # Periodic Boundary Condition on sides
        pos[i,0] = pos[i,0] - (int(pos[i,0]/config.Lx)) * 2.0 * config.Lx      # Periodic Boundary Condition
        pos[i,1] = pos[i,1] - (int(pos[i,1]/config.Ly)) * 2.0 * config.Ly      # Periodic Boundary Condition
        # Reflecting boundary (Will not work due to charge neutralization)
        if (pos[i,2] > config.Lz):
            pos[i,2] -= uvel[i,2] * config.dt
            uvel[i,2] = -uvel[i,2]

    for i in prange(config.N):
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
                acc[i,0] += fx /M[i]
                acc[i,1] += fy /M[i]
                acc[i,2] += fz /M[i]

    for i in prange(config.N):
        vvel[i,:] = uvel[i,:] + acc[i,:] * config.dt / 2.0
        # vvel[i,0] = uvel[i,0] + acc[i,0] * config.dt / 2.0
        # vvel[i,1] = uvel[i,1] + acc[i,1] * config.dt / 2.0
        # vvel[i,2] = uvel[i,2] + acc[i,2] * config.dt / 2.0


    for i in prange(config.N):
        for j in range(config.N):
            if (i != j):
                xdiff = ( pos[i,0]-pos[j,0] )
                ydiff = ( pos[i,1]-pos[j,1] )
                zdiff = ( pos[i,2]-pos[j,2] )
                vxdiff = ( vvel[i,0]-vvel[j,0] )
                vydiff = ( vvel[i,1]-vvel[j,1] )
                vzdiff = ( vvel[i,2]-vvel[j,2] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                if (r < 2*config.rc):
                    vvel[i,0] = vvel[i,0] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vvel[i,1] = vvel[i,1] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vvel[i,2] = vvel[i,2] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)
                    vvel[j,0] = vvel[j,0] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vvel[j,1] = vvel[j,1] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vvel[j,2] = vvel[j,2] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)

    for i in prange(config.N):
        KE += ((vvel[i,0]*vvel[i,0]) + (vvel[i,1]*vvel[i,1]) + (vvel[i,2]*vvel[i,2]) ) / 2.0

    return pos,vvel,uvel,acc,Q,KE,fduration,Qcollect
