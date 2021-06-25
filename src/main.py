#!/usr/env/ python

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
from os.path import join as pjoin
import ini
import sys
import argparse

## User defined functions

from thermostat import berendsen
import diagn


path        = "data/"

def main(argv):
    """ PPDyn main() function """
    parser = argparse.ArgumentParser(description='Plasma Particle Dynamics (PPDyn)')
    parser.add_argument('-i','--input', default='input.ini', type=str, help='Input file name')
    args        = parser.parse_args()
    inputFile   = args.input

    params = ini.parse(open(inputFile).read())
    #========== Input Parameters ===========

    Lx      = float(params['simbox']['Lx'])  # System length in X
    Ly      = float(params['simbox']['Ly'])   # System length in Y
    Lz      = float(params['simbox']['Lz'])   # System length in Z

    N       = int(params['particles']['N'])    # Number of particles

    Vxmax   = float(params['particles']['Vxmax']) # Maximum velocity in X
    Vymax   = float(params['particles']['Vymax']) # Maximum velocity in Y
    Vzmax   = float(params['particles']['Vzmax']) # Maximum velocity in Z

    k       = float(params['screening']['k'])

    Temp    = float(params['particles']['Temp'])

    tmax    = float(params['time']['tmax'])  # Final time
    dt      = float(params['time']['dt']) # time step size
    Nt      = round(tmax/dt) #number of time steps

    #========= Boundary ==========
    btype   = str(params['boundary']['btype']) # Type of boundary

    #========= Diagnostics =======
    dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
    path        = "data/"  # DO NOT CHANGE THE PATH
    os.makedirs(path, exist_ok=True)

    dumpData    = bool(params['diagnostics']['dumpData'])
    f           = h5py.File(path+"particle.hdf5","w")
    if dumpData:
        diagn.attributes(f,tmax,Lx,Ly,Lz,N,dt,dumpPeriod)
        dset = f.create_dataset('energy', (1,), maxshape=(None,), dtype='float64', chunks=(1,))

    vtkData     = bool(params['diagnostics']['vtkData'])
    realTime    = bool(params['diagnostics']['realTime'])
    #========== Options ============
    parallelMode    = bool(params['options']['parallelMode'])
    if parallelMode:
        if btype == 'periodic':
            from pusher_parallel import verlet_periodic as verlet
            from init import initial_periodic as initial
            print("Running in Parallel Mode (Periodic boundary)")
        elif btype == 'reflecting':
            from pusher_parallel import verlet_reflecting as verlet
            from init import initial_reflecting as initial
            print("Running in Parallel Mode (Reflecting boundary)")
    else:
        if btype == 'periodic':
            from pusher_serial import verlet_periodic as verlet
            from init import initial_periodic as initial
            print("Running in Serial Mode (Periodic boundary)")
        elif btype == 'reflecting':
            from pusher_serial import verlet_reflecting as verlet
            from init import initial_reflecting as initial
            print("Running in Serial Mode (Reflecting boundary)")
    #========= Initialize ========
    x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num = initial(Lx,Ly,Lz,Vxmax,Vymax,Vzmax,N,tmax,Nt,k,dumpPeriod)

    #========= Time Loop =========

    for t in range(len(time)):
        KE = 0.0   # Reset KE
        x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE = verlet(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE,k)
        #============  Thermostat =========================
        vx,vy,vz = berendsen(vx,vy,vz,dt,Temp,KE,N,t,tmax)

        #============ Diagnostics Write ===================
        if dumpData:
            if t%dumpPeriod==0:
                diagn.configSpace(f,dset,t,x,y,z,KE)
                print('TimeSteps = %d'%int(t)+' of %d'%Nt+' Energy: %e'%KE)

    if vtkData:
        from vtk_data import vtkwrite
        print('Writing VTK files for Paraview visualization ...')
        vtkwrite()
    return 0
    #========== End of Time Loop ======

if __name__== "__main__":
	start = time.time()
	main(sys.argv[1:])
	end = time.time()
	print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
	os.remove(pjoin(path,'energy.txt'))
