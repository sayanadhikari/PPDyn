#!/usr/bin/ python

"""
** Plasma Particle Dynamics (PPDyn) **

@file                main.py
@authors             Sayan Adhikari <sayan.adhikari@fys.uio.no>
                     Rupak Mukherjee <rupakm@princeton.edu>
@date                17.03.2021

"""


import numpy as np
from numba import jit
import time
import os
import ini

## User defined functions
from initial import initial
from thermostat import berendsen
from pusher import verlet
import diagn



def main():
    """ PPDyn main() function """
    params = ini.parse(open('input.ini').read())
    #========== Input Parameters ===========

    Lx      = float(params['simbox']['Lx'])  # System length in X
    Ly      = float(params['simbox']['Ly'])   # System length in Y
    Lz      = float(params['simbox']['Lz'])   # System length in Z

    N       = int(params['particles']['N'])    # Number of particles

    Vxmax   = float(params['particles']['Vxmax']) # Maximum velocity in X
    Vymax   = float(params['particles']['Vymax']) # Maximum velocity in Y
    Vzmax   = float(params['particles']['Vzmax']) # Maximum velocity in Z

    Temp    = float(params['particles']['Temp'])

    tmax    = float(params['time']['tmax'])  # Final time
    dt      = float(params['time']['dt']) # time step size
    Nt      = round(tmax/dt) #number of time steps

    #========= Boundary ==========
    btype   = str(params['boundary']['btype']) # Type of boundary

    #========= Diagnostics =======
    dumpPeriod = int(params['diagnostics']['dumpPeriod'])
    path ="data/"  # DO NOT CHANGE THE PATH

    #========= Initialize ========
    x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num = initial(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,dumpPeriod)

    #========= Time Loop =========

    for t in range(len(time)):
        KE = 0.0   # Reset KE
        x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE = verlet(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE)

        #============ Diagnostics Write ===================
        if t%dumpPeriod==0:
#             diagn.configSpace(t,N,Nt,x,y,z,path)
        #============  Thermostat =========================
        vx,vy,vz = berendsen(vx,vy,vz,dt,Temp,KE,N,t,tmax)
    return 0
    #========== End of Time Loop ======

if __name__== "__main__":
	start = time.time()
	main()
	end = time.time()
	print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
