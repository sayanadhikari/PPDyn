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
from pathlib import Path
import ini
import sys
import argparse
try:
    from tasktimer import TaskTimer
except ImportError:
    # Simple dummy TaskTimer if not available
    class TaskTimer:
        def __init__(self): pass
        def task(self, name): pass
        def __str__(self): return "TaskTimer not available"
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
    # Remove existing data directory/file and recreate
    data_dir = config.dataDir.strip()
    if os.path.exists(data_dir):
        if os.path.isdir(data_dir):
            shutil.rmtree(data_dir, ignore_errors=True)
        else:
            os.remove(data_dir)
    os.makedirs(data_dir, exist_ok=True)
    # dumpData    = bool(params['diagnostics']['dumpData'])
    f  = h5py.File(pjoin(data_dir,"particle.hdf5"),"w")
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
        # Initialize particle status: 0=active, 1=absorbed, 2=attached
        particle_status = np.zeros(config.N, dtype=np.int32)
        # Initialize impact heatmap (same shape as object mask)
        impact_heatmap = np.zeros_like(obj_mask, dtype=np.float64)
    else:
        pos, vvel, uvel, acc, time, data_num = initial(Q, M)
        particle_status = None
        impact_heatmap = None
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
            pos, vvel, uvel, acc, Q, particle_status, impact_heatmap = verlet(
                t, pos, vvel, uvel, acc, Q, M,
                xg, yg, zg, Ex, Ey, Ez,
                xg_obj, yg_obj, zg_obj, obj_mask,
                particle_status, impact_heatmap
            )
        else:
            # For non-object case, verlet doesn't need status/heatmap
            # This case shouldn't happen with current code, but keep for compatibility
            result = verlet(
                t, pos, vvel, uvel, acc, Q, M,
                xg, yg, zg, Ex, Ey, Ez
            )
            if len(result) > 5:
                pos, vvel, uvel, acc, Q, particle_status, impact_heatmap = result
            else:
                pos, vvel, uvel, acc, Q = result

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
                diagn.configSpace(dsetE,dsetPart,dsetVel,t,pos,vvel,KE,particle_status)
                print('TimeSteps = %d'%int(t)+' of %d'%config.Nt+' Energy: %e'%KE)

    timer.task('Step: Diagnostics')
    # Save impact heatmap and particle status if object_data is enabled
    if config.object_data and impact_heatmap is not None:
        f.create_dataset('impact_heatmap', data=impact_heatmap)
        f.create_dataset('particle_status', data=particle_status)
        print(f'Impact heatmap saved. Total impacts: {np.sum(impact_heatmap):.0f}')
        print(f'  Absorbed particles: {np.sum(particle_status == 1)}')
        print(f'  Attached particles: {np.sum(particle_status == 2)}')
        print(f'  Active particles: {np.sum(particle_status == 0)}')
    if config.vtkData:
        from vtk_data import vtkwrite
        print('Writing VTK files for Paraview visualization ...')
        vtkwrite(data_dir)
    if os.path.exists(pjoin(data_dir,'energy.txt')):
        os.remove(pjoin(data_dir,'energy.txt'))

    timer.task(None)
    print(timer)
    #========== End of Time Loop ======

if __name__== "__main__":
	start = time.time()
	main()
	end = time.time()
	print("Elapsed (after compilation) = %s"%(end - start)+" seconds")

