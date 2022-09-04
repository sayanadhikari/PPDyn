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
Nt   = dp*h5.attrs["Nt"]


data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

time = data_num*dp
i = len(data_num)-1

x = h5["/%d"%data_num[i]+"/position/x"]
y = h5["/%d"%data_num[i]+"/position/y"]
z = h5["/%d"%data_num[i]+"/position/z"]
x = np.array(x)
y = np.array(y)
z = np.array(z)

vx = h5["/%d"%data_num[i]+"/velocity/vx"]
vy = h5["/%d"%data_num[i]+"/velocity/vy"]
vz = h5["/%d"%data_num[i]+"/velocity/vz"]
vx = np.array(vx)
vy = np.array(vy)
vz = np.array(vz)

#### FIG SIZE CALC ############
# figsize = np.array([77,77/1.618]) #Figure size in mm
figsize = np.array([100,100/1.618]) #Figure size in mm
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
#
fig,ax = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)

## Uncomment for 3D
# ax = plt.figure().add_subplot(projection='3d')
# ax.quiver(x, y, vx, vy, normalize=False)
ax.quiver(x, y, vx, vy) #2D

ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

# ax2.plot(N,energy[10:])
# ax2.set_xlabel("$timestep$")
# ax2.set_ylabel("$Energy$")

plt.show()
