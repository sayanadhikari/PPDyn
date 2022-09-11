#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib as mp
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
#matplotlib.use('MacOSX')
# import plotly.graph_objects as go
#========= Configuration ===========
show_anim = True
save_anim = False
interval = 10 #in milliseconds
try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')
# folder = sys.argv[1]

file_name = "particle"#"rhoNeutral" #"P"

h5 = h5py.File(pjoin(DIR,file_name+'.hdf5'),'r')

dataforce = h5["/force"]
datapos = h5["/position"]

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
Lz = h5.attrs["Lz"]

dp   = h5.attrs["dp"]
Nt   = h5.attrs["Nt"]
N   = h5.attrs["N"]

data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

F_col = np.zeros((len(data_num),N), dtype=np.float64)
F_rep = np.zeros((len(data_num),N), dtype=np.float64)
F_flow = np.zeros((len(data_num),N), dtype=np.float64)
F_rand = np.zeros((len(data_num),N), dtype=np.float64)
F_ndrag = np.zeros((len(data_num),N), dtype=np.float64)
posx = np.zeros((len(data_num),N), dtype=np.float64)
posy = np.zeros((len(data_num),N), dtype=np.float64)

for i in range(len(data_num)):
    F_col[i,:] = dataforce[i,:,0] #h5["/%d"%data_num[i]+"/F_col"]
    F_rep[i,:] = dataforce[i,:,1] #h5["/%d"%data_num[i]+"/F_rep"]
    F_flow[i,:] = dataforce[i,:,2] #h5["/%d"%data_num[i]+"/F_flow"]
    F_rand[i,:] = dataforce[i,:,3] #h5["/%d"%data_num[i]+"/F_rand"]
    F_ndrag[i,:] = dataforce[i,:,4] #h5["/%d"%data_num[i]+"/F_ndrag"]
    posx[i,:] = datapos[i,:,0] #h5["/%d"%data_num[i]+"/position/x"]
    posy[i,:] = datapos[i,:,1] #h5["/%d"%data_num[i]+"/position/y"]

time = data_num*dp

F_col_sum = np.sum(F_col,axis=1)
F_rep_sum = np.sum(F_rep,axis=1)
F_flow_sum = np.sum(F_flow,axis=1)
F_rand_sum = np.sum(F_rand,axis=1)
F_ndrag_sum = np.sum(F_ndrag,axis=1)

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

part_num = 300

fig,ax = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax.plot(time,F_col[:,part_num],label='$F_{col}$')
ax.plot(time,F_rep[:,part_num],label='$F_{rep}$')
ax.plot(time,F_flow[:,part_num],label='$F_{flow}$')
ax.legend()
ax.set_xlabel("$timestep$")
ax.set_ylabel("$Force$")

fig,ax1 = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# ax1.plot(time,posx[:,part_num],label='$x$')
# ax1.plot(time,posy[:,part_num],label='$y$')
ax1.plot(posx[:,part_num],posy[:,part_num],label='$pos$')
# ax1.plot(time,F_rep[:,100],label='$F_{rep}$')

ax1.legend()
ax1.set_xlabel("$x$")
ax1.set_ylabel("$y$")

fig,ax2 = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax2.plot(time,F_col_sum/N,label='$F_{col}$')
ax2.plot(time,F_rep_sum/N,label='$F_{rep}$')
ax2.plot(time,F_flow_sum/N,label='$F_{flow}$')
ax2.legend()
ax2.set_xlabel("$timestep$")
ax2.set_ylabel("$Force$")




plt.show()
