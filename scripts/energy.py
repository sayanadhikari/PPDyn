#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
import matplotlib as mp
# import plotly.graph_objects as go
#========= Configuration ===========

try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')

file_name = "particle"#"rhoNeutral" #"P"

h5 = h5py.File(pjoin(DIR,file_name+'.hdf5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
Lz = h5.attrs["Lz"]

N = h5.attrs["N"]

dp   = h5.attrs["dp"]
Nt   = h5.attrs["Nt"]


data_num = np.arange(start=0, stop=Nt, step=1, dtype=int)

time = data_num*dp
energy = h5["/energy"]
energy = 3*(np.array(energy[:-1]))/N

#### FIG SIZE CALC ############
figsize = np.array([77,77/1.618]) #Figure size in mm
dpi = 300                         #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution
#########################
mp.rc('text', usetex=True)
mp.rc('font', family='sans-serif', size=14, serif='Computer Modern Roman')
mp.rc('axes', titlesize=14)
mp.rc('axes', labelsize=14)
mp.rc('xtick', labelsize=14)
mp.rc('ytick', labelsize=14)
mp.rc('legend', fontsize=14)

fig,ax1 = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax1.plot(time[10:],energy[10:])
ax1.set_xlabel("$timestep$")
ax1.set_ylabel("$Energy$")

# ax2.plot(N,energy[10:])
# ax2.set_xlabel("$timestep$")
# ax2.set_ylabel("$Energy$")




plt.show()
