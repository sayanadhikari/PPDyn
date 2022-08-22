#! /usr/bin/env python
import numpy as np
import ini
from os.path import join as pjoin
import sys
import os
import argparse


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

f0      = float(params['particles']['f0'])

k       = float(params['screening']['k'])
lambda_c = float(params['screening']['lambda_c'])
nu = float(params['neutral']['nu'])

g       = float(params['gravity']['g_0'])
rc      = float(params['cutoff radius']['rc'])
Temp    = float(params['particles']['Temp'])

tmax    = float(params['time']['tmax'])  # Final time
dt      = float(params['time']['dt']) # time step size
Nt      = round(tmax/dt) #number of time steps

dist    = bool(params['particles']['dist'])

mean = float(params['particles']['mean'])
stdDev = float(params['particles']['stdDev'])

#========= Boundary ==========
btype   = str(params['boundary']['btype']) # Type of boundary

#========= Diagnostics =======
dumpPeriod  = int(params['diagnostics']['dumpPeriod'])
dataDir    = str(params['directory']['dataDir'])

dumpData    = bool(params['diagnostics']['dumpData'])

vtkData     = bool(params['diagnostics']['vtkData'])
realTime    = bool(params['diagnostics']['realTime'])
#========== Options ============
parallelMode    = bool(params['options']['parallelMode'])
