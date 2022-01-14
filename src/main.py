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

## User defined constants
from constants import K_zd, q_e, u


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
    r_min = float(params['simbox']['r_min']) 
    r_max = float(params['simbox']['r_max']) 

    N       = int(params['particles']['N'])    # Number of particles

    Vxmax   = float(params['particles']['Vxmax']) # Maximum velocity in X
    Vymax   = float(params['particles']['Vymax']) # Maximum velocity in Y
    Vzmax   = float(params['particles']['Vzmax']) # Maximum velocity in Z

    k       = float(params['screening']['k'])

    g       = float(params['gravity']['g_0'])
    #rc      = float(params['cutoff radius']['r_c'])
    Temp    = float(params['particles']['Temp'])

    tmax    = float(params['time']['tmax'])  # Final time
    dt      = float(params['time']['dt']) # time step size
    Nt      = round(tmax/dt) #number of time steps

    dist    = bool(params['particles']['dist'])

    #========= Charge and Mass distribution ========
    if dist:
        mean = float(params['particles']['mean'])
        stdDev = float(params['particles']['stdDev'])
        M  = np.random.normal(loc=mean,scale=stdDev,size=(N)) #mass of particles (Gaussian)
        Q = M**(2/3) # charge of particles
    else:
        M = np.ones(N)
        Q = M


    density = float(params['particles']['density'])*np.ones(N)

    #========= a scaling =========
    a_set = str(params['particles']['a']).split(',')
    a_set = [float(a) for a in a_set]
    a = np.random.choice(a_set,(N,))
    a_scale = a[:]/np.mean(a_set)


    M = a[:]*a[:]*a[:]*density
    E = float(params['plasma']['E'])

    #===== Neutral drag  =========
    R = 8.31446261815324 #Gas const
    k_b = 1.38064852e-23 #BOltzman const
    delta = float(params['neutrals']['delta'])
    mn = float(params['neutrals']['mn']) * u
    Tn = float(params['neutrals']['Tn'])
    P = float(params['neutrals']['pressure'])
    V = Lx*Ly*Lz
    v_th_neutrals = np.sqrt((8*k_b*Tn)/(np.pi*mn))
    N_neutrals = (P*V*6.02214076e23)/(R*Tn)
    print('N neutrals: ', N_neutrals)
    Kn_drag = delta*1.33*np.pi*N_neutrals*v_th_neutrals*mn
    print('Kn_drag: ', Kn_drag)


    drag = float(params['plasma']['drag']) #Maybe delete

    #===== Particle charge =======
    Te = float(params['plasma']['Te'])
    phi_hat = float(params['plasma']['phi_hat'])
    Zd = a[:]*Te*K_zd*phi_hat
    Q = Zd[:]*q_e

    #========= Force Field =======
    Fpath = str(params['plasma']['Ffield'])
    Ffield = np.load(Fpath)
    if Ffield['arr_0'][-1] - Lx > 0.00001: 
        Ff_rmax = Ffield["arr_0"][-1]; raise ValueError(f'Force field npz must have same range in r {Ff_rmax} as Lx {Lx}')
    Fr = Ffield['arr_0']
    Fel = Ffield['arr_1']
    Fion = Ffield['arr_2']
    Ftot = Ffield['arr_3']
    #Fel[:] = Fel[:]*E
    #Fdrag[:] = Fdrag[:]*drag

    #========= Boundary ==========
    btype   = str(params['boundary']['btype']) # Type of boundary
    geometry = str(params['boundary']['geometry'])

    #========= Diagnostics =======
    dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
    pathName    = str(params['directory']['path'])
    path        = pjoin(pathName)
    # path        = "data/"  # DO NOT CHANGE THE PATH
    if  os.path.exists(path)== False:
        os.mkdir(path)
    dumpData    = bool(params['diagnostics']['dumpData'])
    f           = h5py.File(pjoin(path,"particle.hdf5"),"w")
    if dumpData:
        diagn.attributes(f,tmax,Lx,Ly,Lz,a,Q,M,N,dt,dumpPeriod)
        dsetE = f.create_dataset('energy', (1,), maxshape=(None,), dtype='float64', chunks=(1,))
        dsetQ = f.create_dataset('Qcollect', (1,), maxshape=(None,), dtype='float64', chunks=(1,))

    vtkData     = bool(params['diagnostics']['vtkData'])
    realTime    = bool(params['diagnostics']['realTime'])
    #========== Options ============
    parallelMode    = bool(params['options']['parallelMode'])
    if parallelMode:
        if btype == 'periodic' and geometry == 'radial':
            from pusher_parallel import verlet_periodic_radial as verlet
            from init import initial_periodic_radial as initial
            print("Running in Parallel Mode (Periodic boundary, radial geometry)")
        elif btype == 'periodic' and geometry == 'box':
            from pusher_parallel import verlet_periodic as verlet
            from init import initial_periodic as initial
            print("Running in Parallel Mode (Periodic boundary)")
        elif btype == 'reflecting' and geometry == 'box':
            from pusher_parallel import verlet_reflecting as verlet
            from init import initial_reflecting as initial
            print("Running in Parallel Mode (Reflecting boundary)")
    else:
        if btype == 'periodic' and geometry == 'box':
            from pusher_serial import verlet_periodic as verlet
            from init import initial_periodic as initial
            print("Running in Serial Mode (Periodic boundary)")
        elif btype == 'reflecting' and geometry == 'box':
            from pusher_serial import verlet_reflecting as verlet
            from init import initial_reflecting as initial
            print("Running in Serial Mode (Reflecting boundary)")
    #========= Initialize ========
    x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,time,data_num,fduration = initial(Lx,Ly,Lz,r_min,r_max,Vxmax,Vymax,Vzmax,N,tmax,Nt,k,dumpPeriod,g,Q,M,Temp)

    #========= Time Loop =========

    for t in range(len(time)):
        KE = 0.0   # Reset KE
        Qcollect = 0.0 # Initialize Q_collect
        x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,KE,Q,fduration,Qcollect = verlet(x,y,z,vx,vy,vz,ux,uy,uz,ax,ay,az,dt,Lx,Ly,Lz,N,KE,k,g,Q,M,a,a_scale,Fr,Fel,Fion,Kn_drag,fduration,t,Qcollect)
        #============  Thermostat =========================
        # vx,vy,vz = berendsen(vx,vy,vz,dt,Temp,KE,N,t,tmax)

        #============ Diagnostics Write ===================
        if dumpData:
            if t%dumpPeriod==0:
                diagn.configSpace(f,dsetE,dsetQ,t,x,y,z,vx,vy,vz,KE,Qcollect,path)
                print('TimeSteps = %d'%int(t)+' of %d'%Nt+' Energy: %e'%KE)

    diagn.dustDiagn(f,fduration)
    if vtkData:
        from vtk_data import vtkwrite
        print('Writing VTK files for Paraview visualization ...')
        vtkwrite(path)
    os.remove(pjoin(path,'energy.txt'))
    return 0
    #========== End of Time Loop ======

if __name__== "__main__":
	start = time.time()
	main(sys.argv[1:])
	end = time.time()
	print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
