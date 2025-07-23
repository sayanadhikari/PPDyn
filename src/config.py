#! /usr/bin/env python
import numpy as np
import ini
from os.path import join as pjoin
import sys
import os
import argparse
import h5py


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
maxwell_load = bool(params['particles']['maxwell_load']) # Maximum velocity in Z

k       = float(params['screening']['k'])

g       = float(params['gravity']['g_0'])
rc      = float(params['cutoff radius']['rc'])
Temp    = float(params['particles']['Temp'])

tmax    = float(params['time']['tmax'])  # Final time
dt      = float(params['time']['dt']) # time step size
Nt      = round(tmax/dt) #number of time steps

dist = str(params['particles']['dist'])   # Now can be 'gaussian' or 'uniform_radius'

mean = float(params['particles']['mean'])
stdDev = float(params['particles']['stdDev'])
density = float(params['particles'].get('density', 1.0))  # default density 1
r_min = float(params['particles'].get('r_min', 0.0))      # default to 0
r_max = float(params['particles'].get('r_max', 0.0))      # default to 0

#========= Boundary ==========
btype   = str(params['boundary']['btype']) # Type of boundary

#========= Diagnostics =======
dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
dataDir    = str(params['directory']['dataDir'])
picDir    = str(params['directory']['picDir'])

dumpData    = bool(params['diagnostics']['dumpData'])

vtkData     = bool(params['diagnostics']['vtkData'])
realTime    = bool(params['diagnostics']['realTime'])
#========== Options ============
parallelMode    = bool(params['options']['parallelMode'])
PIC_data = bool(params['options']['PIC_data'])
object_data = bool(params['options']['object_data'])

#====== Additional =======
dumpNt = round(Nt/dumpPeriod)

# Electric field normalization method: 'mean' or 'max'
E_norm_method = str(params['options'].get('E_norm_method', 'mean'))


# PINC coupling
with h5py.File("PIC_data/E.grid.h5", "r") as f:
    denorm = f.attrs["Axis denormalization factor"][0]
    # Get the last n entry (sorted by key)
    keys = sorted(f.keys())
    last_key = keys[-1]
    grid_shape = f[last_key].shape
    Lx = 0.5 * grid_shape[0] * denorm
    Ly = 0.5 * grid_shape[1] * denorm
    Lz = 0.5 * grid_shape[2] * denorm

