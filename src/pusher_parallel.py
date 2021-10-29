from numba import jit, prange
import numpy as np

@jit(nopython=True, parallel=True)
def verlet_periodic_radial(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE,k,g,Q,M,a,Fr,Fel,Fdrag,fduration,t,Qcollect):
    for i in prange(N):
        ux[i] = vx[i] + ax[i] * dt/2.0
        uy[i] = vy[i] + ay[i] * dt/2.0
        uz[i] = vz[i] + az[i] * dt/2.0
        x[i] = x[i] + ux[i] * dt
        y[i] = y[i] + uy[i] * dt
        z[i] = z[i] + uz[i] * dt
        r = np.sqrt(x[i]*x[i] + y[i]*y[i]+z[i]*z[i])
        #r_boundary = (int(r/Lx)) * 2.0 * Lx
        x[i] -= (int(r/Lx))*x[i]*2.0    # Periodic Boundary Condition
        y[i] -= (int(r/Lx))*y[i]*2.0    # Periodic Boundary Condition
        z[i] -= (int(r/Lx))*z[i]*2.0    # Periodic Boundary Condition
        
        #x[i] -= r_boundary * np.sign(x[i])   # Periodic Boundary Condition
        #y[i] -= r_boundary * np.sign(y[i])     # Periodic Boundary Condition
        #z[i] -= r_boundary * np.sign(z[i])     # Periodic Boundary Condition


    for i in prange(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = 0.0
        fx=0
        fy=0
        fz=0
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] ) - round((x[i]-x[j])/(2.0*Lx)) * 2.0*Lx
                ydiff = ( y[i]-y[j] ) - round((y[i]-y[j])/(2.0*Ly)) * 2.0*Ly
                zdiff = ( z[i]-z[j] ) - round((z[i]-z[j])/(2.0*Lz)) * 2.0*Lz
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                fx += xdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fy += ydiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                fz += zdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r) #+ zdiff*g + Lz*g # zdiff/(r*r*r)
                #az[i] += fcy/M[i]
        #Electric field
        xdiff = x[i] # diff to origin
        ydiff = y[i] #- Ly/2.0
        zdiff = z[i]
        r = np.sqrt(xdiff*xdiff + ydiff*ydiff +zdiff*zdiff)

        fel = np.interp(r,Fr,Fel)
        fx += a[i]*fel * xdiff/r
        fy += a[i]*fel * ydiff/r
        fz += a[i]*fel * zdiff/r
        fdrag = np.interp(r,Fr,Fdrag)
        fx += a[i]*a[i]*fdrag * xdiff/r
        fy += a[i]*a[i]*fdrag * ydiff/r
        fz += a[i]*a[i]*fdrag * zdiff/r
        
        ax[i] += fx/M[i]
        ay[i] += fy/M[i]
        az[i] += fz/M[i]

    for i in prange(N):
        vx[i] = ux[i] + ax[i] * dt / 2.0
        vy[i] = uy[i] + ay[i] * dt / 2.0
        vz[i] = uz[i] + az[i] * dt / 2.0
        KE += ((vx[i]*vx[i]) + (vy[i]*vy[i]) + (vz[i]*vz[i]) ) / 2.0
    
    #z = np.zeros(len(x))   #For 2D
    return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE,Q,fduration,Qcollect
    #return x,y,z,vx,vy,z,ux,uy,z,ax,ay,z,KE,Q,fduration,Qcollect



@jit(nopython=True, parallel=True)
def verlet_periodic(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE,k,g,Q,M,E):
    for i in prange(N):
        ux[i] = vx[i] + ax[i] * dt/2.0
        uy[i] = vy[i] + ay[i] * dt/2.0
        uz[i] = vz[i] + az[i] * dt/2.0
        x[i] = x[i] + ux[i] * dt
        y[i] = y[i] + uy[i] * dt
        z[i] = z[i] + uz[i] * dt
        x[i] = x[i] - (int(x[i]/Lx)) * 2.0 * Lx      # Periodic Boundary Condition
        y[i] = y[i] - (int(y[i]/Ly)) * 2.0 * Ly      # Periodic Boundary Condition
        z[i] = z[i] - (int(z[i]/Lz)) * 2.0 * Lz      # Periodic Boundary Condition

    for i in prange(N):
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
                fz = zdiff*(1+k*r)*np.exp(-k*r)/(r*r*r) #+zdiff*g  +Lz*g  # zdiff/(r*r*r)
                ax[i] += fx
                ay[i] += fy
                az[i] += fz

    for i in prange(N):
        vx[i] = ux[i] + ax[i] * dt / 2.0
        vy[i] = uy[i] + ay[i] * dt / 2.0
        vz[i] = uz[i] + az[i] * dt / 2.0
        KE += ((vx[i]*vx[i]) + (vy[i]*vy[i]) + (vz[i]*vz[i]) ) / 2.0
    return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE

@jit(nopython=True, parallel=True)
def verlet_reflecting(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE,k,g,Q,M,E,drag,fduration,t,Qcollect):
    for i in prange(N):
        ux[i] = vx[i] + ax[i] * dt/2.0
        uy[i] = vy[i] + ay[i] * dt/2.0
        uz[i] = vz[i] + az[i] * dt/2.0
        if (z[i] <= -Lz):
            z[i] = -Lz
            uz[i] = 0.0
            ux[i] = 0.0
            uy[i] = 0.0
            if Q[i] > 0.0:
                fduration[i] = t
                Qcollect += Q[i]
            # print(i,fduration[i],t)
            Q[i]  = 0.0      #Make charges zero as they hit the ground
        x[i] = x[i] + ux[i] * dt
        y[i] = y[i] + uy[i] * dt
        z[i] = z[i] + uz[i] * dt
        # Reflecting boundary
        if (x[i] > Lx or x[i] < -Lx):
            x[i] -= ux[i] * dt
            ux[i] = -ux[i]
        if (y[i] > Ly or y[i] < -Ly):
            y[i] -= uy[i] * dt
            uy[i] = -uy[i]
        if (z[i] > Lz):
            z[i] -= uz[i] * dt
            uz[i] = -uz[i]

    for i in prange(N):
        ax[i] = 0.0
        ay[i] = 0.0
        az[i] = -(z[i]+Lz)*g
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] )
                ydiff = ( y[i]-y[j] )
                zdiff = ( z[i]-z[j] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                #Coloumb
                fcx = xdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # xdiff/(r*r*r)
                fcy = ydiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r)    # ydiff/(r*r*r)
                #fz = zdiff*(1+k*r)*np.exp(-k*r)*(Q[i]*Q[j])/(r*r*r) #+ zdiff*g + Lz*g # zdiff/(r*r*r)

                ax[i] += fcx/M[i]
                ay[i] += fcy/M[i]

        #r = np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i])
        #r2 = r*r
        a = 0.01
        b = 0.1
        #Electric field
        xdiff = x[i] # diff to origin
        ydiff = y[i] #- Ly/2.0
        r = np.sqrt(xdiff*xdiff + ydiff*ydiff)
        
        fx = Q[i] * E * xdiff / (0.2+np.exp(3.0-a*r*r)+a*0.2*r*r)   #1.5/((a*0.3)+np.exp(a*x[i]*x[i])+(0.2*a*x[i]*x[i]))
        fy = Q[i] * E * ydiff / (0.2+np.exp(3.0-a*r*r)+a*0.2*r*r)
        
        fx += drag[i] * xdiff / (0.2+np.exp(3.0-b*r*r)+b*0.2*r*r)   #1.5/((a*0.3)+np.exp(a*x[i]*x[i])+(0.2*a*x[i]*x[i]))
        fy += drag[i] * ydiff / (0.2+np.exp(3.0-b*r*r)+b*0.2*r*r)   #1.5/((a*0.3)+np.exp(a*y[i]*y[i])+(0.2*a*y[i]*y[i]))
        #fz = -(z[i]/np.abs(z[i]))*(Q[i]*E)/z[i]    #1.5/((a*0.3)+np.exp(a*r2)+(0.2*a*r2))

        #Ion drag force
        """
        fx += x[i] * drag[i] * 1/(b*0.3+np.exp(b*r2)+0.2*b*r2)
        fy += y[i] * drag[i] * 1/(b*0.3+np.exp(b*r2)+0.2*b*r2)
        fz += z[i] * drag[i] * 1/(b*0.3+np.exp(b*r2)+0.2*b*r2)
        """
        ax[i] += fx/M[i]
        ay[i] += fy/M[i]

    for i in prange(N):
        vx[i] = ux[i] + ax[i] * dt / 2.0
        vy[i] = uy[i] + ay[i] * dt / 2.0
        vz[i] = uz[i] + az[i] * dt / 2.0


    for i in prange(N):
        for j in range(N):
            if (i != j):
                xdiff = ( x[i]-x[j] )
                ydiff = ( y[i]-y[j] )
                zdiff = ( z[i]-z[j] )
                vxdiff = ( vx[i]-vx[j] )
                vydiff = ( vy[i]-vy[j] )
                vzdiff = ( vz[i]-vz[j] )
                r = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                rc= 1e-5
                if (r < 2*rc):
                    vx[i] = vx[i] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vy[i] = vy[i] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vz[i] = vz[i] - ( 2*M[j] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)
                    vx[j] = vx[j] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * xdiff / (r*r)
                    vy[j] = vy[j] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * ydiff / (r*r)
                    vz[j] = vz[j] + ( 2*M[i] / (M[i]+M[j]) ) * (vxdiff*xdiff + vydiff*ydiff + vzdiff*zdiff) * zdiff / (r*r)

    for i in prange(N):
        KE += ((vx[i]*vx[i]) + (vy[i]*vy[i]) + (vz[i]*vz[i]) ) / 2.0

    z = np.zeros(len(x))   
    #return x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE,Q,fduration,Qcollect
    return x,y,z,vx,vy,z,ux,uy,z,ax,ay,z,KE,Q,fduration,Qcollect
