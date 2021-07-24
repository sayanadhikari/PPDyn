#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt

# import plotly.graph_objects as go
#========= Configuration ===========

DIR ="../data"

file_name = "particle"#"rhoNeutral" #"P"

h5 = h5py.File('../data/'+file_name+'.hdf5','r')

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

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8, 6))
ax1.plot(time[10:],energy[10:])
ax1.set_xlabel("$timestep$")
ax1.set_ylabel("$Energy$")

# ax2.plot(N,energy[10:])
# ax2.set_xlabel("$timestep$")
# ax2.set_ylabel("$Energy$")




plt.show()
