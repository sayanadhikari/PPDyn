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

dp   = h5.attrs["dp"]
Nt   = h5.attrs["Nt"]


data_num = np.arange(start=0, stop=Nt, step=1, dtype=int)

time = data_num
energy = h5["/energy"]
energy = np.array(energy[:-1])

fig,ax = plt.subplots(figsize=(6, 6))
plt.plot(time[10:],energy[10:])
ax.set_xlabel("$time [s]$")
ax.set_ylabel("$Energy$")




plt.show()
