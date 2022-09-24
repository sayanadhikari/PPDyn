#!/usr/env/ python

"""
** Plasma Particle Dynamics (PPDyn) **

@file                main.py
@authors             Sayan Adhikari <sayan.adhikari@fys.uio.no>

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
import shutil

## User defined functions
from thermostat import berendsen
import diagn
import config


path        = "data/"

def main():
    """ PPDyn main() function """
    timer = TaskTimer()
    timer.task('Step: Reading Variables')

    #========= Charge and Mass distribution ========
    if config.dist:
        mean = config.mean
        stdDev = config.stdDev
        M  = np.random.normal(loc=mean,scale=stdDev,size=config.N) #mass of particles (Gaussian)
        Q = M**(2/3) # charge of particles
    else:
        M = np.ones(N)
        Q = M
    #======== Diagnostics and data management =======
    # if  os.path.exists(config.dataDir)== False:
    #     os.rmdir(config.dataDir)
    shutil.rmtree(config.dataDir, ignore_errors=True)
    os.mkdir(config.dataDir)
    # dumpData    = bool(params['diagnostics']['dumpData'])
    f  = h5py.File(pjoin(config.dataDir,"particle.hdf5"),"w")
    if config.dumpData:
        diagn.attributes(f)
        dsetE = f.create_dataset('energy', (1,), maxshape=(None,), dtype='float64', chunks=(1,))
        dsetPart = f.create_dataset("position", (config.dumpNt, config.N, 3), dtype='float64', compression="gzip", compression_opts=9)
        dsetVel = f.create_dataset("velocity", (config.dumpNt, config.N, 3), dtype='float64', compression="gzip", compression_opts=9)

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
    pos,vvel,uvel,acc,time,data_num = initial(Q,M)
    #========= Time Loop =========
    timer.task('Step: Time Solution')
    for t in range(len(time)):
        KE = 0.0   # Reset KE

        pos,vvel,uvel,acc,Q,KE = verlet(t,pos,vvel,uvel,acc,Q,M,KE)
        #============  Thermostat =========================
        # vvel = berendsen(t,vvel,KE)
        #============ Diagnostics Write ===================
        if config.dumpData:
            if t%config.dumpPeriod==0:
                diagn.configSpace(dsetE,dsetPart,dsetVel,t,pos,vvel,KE)
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
