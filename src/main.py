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
from tasktimer import TaskTimer

## User defined functions
from thermostat import berendsen
import diagn
import config


path        = "data/"

def main():
    """ PPDyn main() function """
    timer = TaskTimer()
    timer.task('Step: Reading Variables')

    # parser = argparse.ArgumentParser(description='Plasma Particle Dynamics (PPDyn)')
    # parser.add_argument('-i','--input', default='input.ini', type=str, help='Input file name')
    # args        = parser.parse_args()
    # inputFile   = args.input
    #
    # params = ini.parse(open(inputFile).read())
    # #========== Input Parameters ===========
    #
    # Lx      = float(params['simbox']['Lx'])  # System length in X
    # Ly      = float(params['simbox']['Ly'])   # System length in Y
    # Lz      = float(params['simbox']['Lz'])   # System length in Z
    #
    # N       = int(params['particles']['N'])    # Number of particles
    #
    # Vxmax   = float(params['particles']['Vxmax']) # Maximum velocity in X
    # Vymax   = float(params['particles']['Vymax']) # Maximum velocity in Y
    # Vzmax   = float(params['particles']['Vzmax']) # Maximum velocity in Z
    #
    # k       = float(params['screening']['k'])
    #
    # g       = float(params['gravity']['g_0'])
    # #rc      = float(params['cutoff radius']['r_c'])
    # Temp    = float(params['particles']['Temp'])
    #
    # tmax    = float(params['time']['tmax'])  # Final time
    # dt      = float(params['time']['dt']) # time step size
    # Nt      = round(tmax/dt) #number of time steps
    #
    # dist    = bool(params['particles']['dist'])

    #========= Charge and Mass distribution ========
    if config.dist:
        mean = config.mean
        stdDev = config.stdDev
        M  = np.random.normal(loc=mean,scale=stdDev,size=config.N) #mass of particles (Gaussian)
        Q = M**(2/3) # charge of particles
    else:
        M = np.ones(config.N)
        Q = M

    # #========= Boundary ==========
    # btype   = str(params['boundary']['btype']) # Type of boundary
    #
    # #========= Diagnostics =======
    # dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
    # pathName    = str(params['directory']['path'])
    # path = pjoin(config.dataDir)
    # path        = "data/"  # DO NOT CHANGE THE PATH
    if  os.path.exists(config.dataDir)== False:
        os.mkdir(config.dataDir)
    # dumpData    = bool(params['diagnostics']['dumpData'])
    f  = h5py.File(pjoin(config.dataDir,"particle.hdf5"),"w")
    if config.dumpData:
        diagn.attributes(f)
        dsetE = f.create_dataset('energy', (1,), maxshape=(None,), dtype='float64', chunks=(1,))
        # dsetQ = f.create_dataset('Qcollect', (1,), maxshape=(None,), dtype='float64', chunks=(1,))

    # vtkData     = bool(params['diagnostics']['vtkData'])
    # realTime    = bool(params['diagnostics']['realTime'])
    # #========== Options ============
    # parallelMode    = bool(params['options']['parallelMode'])
    if config.parallelMode:
        if config.btype == 'periodic':
            from pusher_parallel import verlet_periodic as verlet
            from init import initial_periodic as initial
            print("Running in Parallel Mode (Periodic boundary)\nWarning!! Make sure gravity is zero. Otherwise may lead to unexpected result.")
        elif config.btype == 'reflecting':
            from pusher_parallel import verlet_reflecting as verlet
            from init import initial_reflecting as initial
            print("Running in Parallel Mode (Reflecting boundary)")
        elif config.btype == 'mixed':
            from pusher_parallel import verlet_mixed as verlet
            from init import initial_reflecting as initial
            print("Running in Parallel Mode (Mixed boundary)")
    else:
        print("Serial version not supported anymore!!\nChange parallelMode  = True")
        exit()
        # if config.btype == 'periodic':
        #     from pusher_serial import verlet_periodic as verlet
        #     from init import initial_periodic as initial
        #     print("Running in Serial Mode (Periodic boundary)")
        # elif config.btype == 'reflecting':
        #     from pusher_serial import verlet_reflecting as verlet
        #     from init import initial_reflecting as initial
        #     print("Running in Serial Mode (Reflecting boundary)")
    #========= Initialize ========
    timer.task('Step: Initialization')
    pos,vvel,uvel,acc,time,data_num,fdist = initial(Q,M)
    #========= Time Loop =========
    timer.task('Step: Time Solution')
    for t in range(len(time)):
        KE = 0.0   # Reset KE
        pos,vvel,uvel,acc,Q,KE,fdist = verlet(t,pos,vvel,uvel,acc,Q,M,KE,fdist)
        #============  Thermostat =========================
        vvel = berendsen(t,vvel,KE)
        #============ Diagnostics Write ===================
        if config.dumpData:
            if t%config.dumpPeriod==0:
                diagn.configSpace(f,dsetE,t,pos,vvel,KE,fdist)
                print('TimeSteps = %d'%int(t)+' of %d'%config.Nt+' Energy: %e'%KE)

    timer.task('Step: Diagnostics')

    if config.vtkData:
        from vtk_data import vtkwrite
        print('Writing VTK files for Paraview visualization ...')
        vtkwrite(config.dataDir)
    os.remove(pjoin(config.dataDir,'energy.txt'))

    timer.task(None)
    print(timer)
    #========== End of Time Loop ======

if __name__== "__main__":
	start = time.time()
	main()
	end = time.time()
	print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
