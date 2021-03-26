#!/usr/bin/env python

"""
** Plasma Particle Dynamics (PPDyn) **

@file                main.py
@authors             Sayan Adhikari <sayan.adhikari@fys.uio.no>
                     Rupak Mukherjee <rupakm@princeton.edu>
@date                17.03.2021

"""


import numpy as np
from numba import jit
import h5py
import time
import os
import ini
import sys
import argparse

## User defined functions
from PPDyn.initial import initial
from PPDyn.thermostat import berendsen
from PPDyn.diagn import configSpace
# import PPDyn.diagn




def ppdyn(argv):
    print("===================================================================")
    print("Running PPDyn (Plasma Particle Dynamics)")
    print("Author: Dr. Sayan Adhikari, PostDoc @ UiO, Norway")
    print("::::::: Dr. Rupak Mukherjee, Associate Research Physicist @ PPPL, NJ")
    print("Input: Edit input.ini file to change the parameters for simulation")
    print("===================================================================")
    """ PPDyn main() function """
    parser = argparse.ArgumentParser(description='Plasma Particle Dynamics (PPDyn)')
    parser.add_argument('--i', default='input.ini', type=str, help='Input file name')
    args        = parser.parse_args()
    inputFile   = args.i

    params = ini.parse(open(inputFile).read())
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
    dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
    path        = "data/"  # DO NOT CHANGE THE PATH
    #========== Data Directory Setup =============
    if os.path.exists(path):
        print("Data directory exists. Current data will be replaced after the run.")
    else:
        os.mkdir(path)
    dumpData    = bool(params['diagnostics']['dumpData'])
    f           = h5py.File(path+"particle.hdf5","w")

    #========== Options ============
    parallelMode    = bool(params['options']['parallelMode'])
    if parallelMode:
        from PPDyn.pusher_parallel import verlet
        print("Running in Parallel Mode")
    else:
        from PPDyn.pusher_serial import verlet
        print("Running in Serial Mode")
    #========= Initialize ========
    x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num = initial(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,dumpPeriod)

    #========= Time Loop =========

    for t in range(len(time)):
        KE = 0.0   # Reset KE
        x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE = verlet(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE)
        #============ Diagnostics Write ===================
        if dumpData:
            if t%dumpPeriod==0:
                configSpace(t,N,tmax,x,y,z,Lx,Ly,Lz,f,dt,dumpPeriod)
                print('TimeSteps = %d'%int(t)+' of %d'%Nt+' Energy: %e'%KE)
        #============  Thermostat =========================
        vx,vy,vz = berendsen(vx,vy,vz,dt,Temp,KE,N,t,tmax)
    return "All done"
    #========== End of Time Loop ======

if __name__== "__ppdyn__":
	# start = time.time()
	ppdyn(sys.argv[1:])
	# end = time.time()
	# print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
