#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mp
from os.path import join as pjoin
import sys
#import pandas as pd

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
fduration = h5["/fall_duration"]
Qcollect = h5["/Qcollect"]

Qcurrent = np.zeros(Qcollect.shape)
Qcurrent[:] = Qcollect[:]
Qcurrent[Qcurrent == 0.0] = np.nan


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

fig,ax1 = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

sns.kdeplot(fduration,fill=True,linewidth=1, ax=ax1)
# sns.kdeplot(fduration,fill=True,common_norm=False, palette="crest",
#    alpha=.5, linewidth=0, ax=ax1)

ax1.set_xlabel("$Fall~time \mathrm[steps]$")
ax1.set_ylabel("A.U.")

fig,ax2 = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

ax2.plot(time,Qcurrent[:-1],'.')

ax2.set_xlabel("$\mathrm[timesteps]$")
ax2.set_ylabel("Dust Current")




plt.show()
