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
from init import load_object_grid, copy_pic_data


path        = "data/"

def main():
    """ PPDyn main() function """
    timer = TaskTimer()
    timer.task('Step: Reading Variables')

    #========= Charge and Mass distribution ========
    # if config.dist == 'gaussian':
    #     mean = config.mean
    #     stdDev = config.stdDev
    #     radii = np.random.normal(loc=mean, scale=stdDev, size=config.N)
    #     M  = np.random.normal(loc=mean,scale=stdDev,size=config.N) #mass of particles (Gaussian)
    #     Q = M**(2/3) # charge of particles
    # else:
    #     M = np.ones(N)
    #     Q = M

    #========= Charge and Mass distribution ========
    if config.dist == 'gaussian':
        # Gaussian distribution of particle radii
        mean = config.mean
        stdDev = config.stdDev
        # radii drawn from a normal distribution
        radii = np.random.normal(loc=mean, scale=stdDev, size=config.N)
        # normalize radii to unit mean to maintain earlier normalization
        radii_mean = np.mean(radii)
        radii /= radii_mean
        # compute mass from radii (assumes density provided in config)
        density = config.density
        M = density * (4.0/3.0) * np.pi * radii**3
        # Q = M ** (2/3)  # charge proportional to volume^(2/3), i.e., surface area (commented out)
        # normalize masses to unit mean to maintain earlier normalization
        M_mean = np.mean(M)
        M /= M_mean
        # adjust charge accordingly if Q was computed before normalization
        Q = M ** (2/3)
    elif config.dist == 'uniform_radius':
        # Uniform distribution of radius between r_min and r_max
        radii = np.random.uniform(low=config.r_min, high=config.r_max, size=config.N)
        # normalize radii to unit mean to maintain earlier normalization
        radii_mean = np.mean(radii)
        radii /= radii_mean
        density = config.density
        M = density * (4.0/3.0) * np.pi * radii ** 3
        # Q = M ** (2/3)   # Or Q = 4 * np.pi * radii**2 if you want Q ∝ surface area (commented out)
        # normalize masses to unit mean to maintain earlier normalization
        M_mean = np.mean(M)
        M /= M_mean
        # adjust charge accordingly if Q was computed before normalization
        Q = M ** (2/3)
    else:
        M = np.ones(config.N)
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
        # store particle radii for VTK export
        f.create_dataset('radius', data=radii)

    # vtkData     = bool(params['diagnostics']['vtkData'])
    # realTime    = bool(params['diagnostics']['realTime'])
    # #========== Options ============
    # parallelMode    = bool(params['options']['parallelMode'])
    if config.parallelMode:
        if config.btype == 'periodic':
            from pusher_parallel import verlet_periodic as verlet
            from init import initial_periodic as initial
            from init import load_pic_field
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
    # ========= Initialize ========
    timer.task('Step: Initialization')
    copy_pic_data()
    if config.object_data:
        xg_obj, yg_obj, zg_obj, obj_mask = load_object_grid()
        pos, vvel, uvel, acc, time, data_num = initial(Q, M,
                                                       xg_obj, yg_obj, zg_obj, obj_mask)
    else:
        pos, vvel, uvel, acc, time, data_num = initial(Q, M)
    if config.PIC_data:
        xg, yg, zg, Ex, Ey, Ez = load_pic_field()
        # if config.object_data:
            # xg_obj etc already set
    # ========= Time Loop =========
    timer.task('Step: Time Solution')
    # Ensure PIC_data is enabled before entering time loop
    if not config.PIC_data:
        raise ValueError("PIC_data must be set to True for this simulation version.")
    for t in range(len(time)):
        if config.object_data:
            pos, vvel, uvel, acc, Q = verlet(
                t, pos, vvel, uvel, acc, Q, M,
                xg, yg, zg, Ex, Ey, Ez,
                xg_obj, yg_obj, zg_obj, obj_mask
            )
        else:
            pos, vvel, uvel, acc, Q = verlet(
                t, pos, vvel, uvel, acc, Q, M,
                xg, yg, zg, Ex, Ey, Ez
            )

        # Compute total kinetic energy including particle mass in a serial loop to avoid race conditions
        KE = 0.0
        # compute kinetic energy including particle mass
        for i in range(config.N):
            KE += 0.5 * M[i] * (vvel[i,0]**2 + vvel[i,1]**2 + vvel[i,2]**2)
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
