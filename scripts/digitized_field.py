import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace, logspace


def zero_section(r,f,radius):
    N = len(f)
    for i in range(N):
        if np.abs(r[i]) < radius:
            f[i] = 0

    return f

def mirror(arr):
    N=len(arr)
    middle = int(len(arr)/2)
    new_arr = np.zeros((N,))
    new_arr[0:middle] = arr[0:middle]
    new_arr[middle:-1] = -arr[middle:0:-1]
    return new_arr


savefile = True
doplot = True

DIR = '../staticfields/'
filename = 'F_digitized'

df_fEl = pd.read_csv('../staticfields/digitized/Fel.csv', header=None)
df_Fion = pd.read_csv('../staticfields/digitized/Fion.csv', header=None)
radius = 70 #mm
r = linspace(-radius,radius,6*radius+1)
scale = 1e-12

fEl = np.interp(r,df_fEl[0].to_list(),df_fEl[1].to_list())
fIon = np.interp(r,df_Fion[0].to_list(),df_Fion[1].to_list())
fEl = mirror(fEl)
fIon = mirror(fIon)
#Scaling, note that fIon is multiplied with 2
#fEl = zero_section(r,fEl,15)
fEl[:] = fEl[:]*scale
fIon[:] = fIon[:]*scale*2
r[:] = 1e-3*r[:]
ftot = fEl[:] + fIon[:]

if doplot:
    fig,ax = plt.subplots()
    ax.plot(r,fEl,label='$F_{EL}$',color='blue')
    ax.plot(r,fIon,label='$F_{ion}$',color='r')
    ax.plot(r,ftot,label='$F_{total}$',color='orange')
    ax.set_xlabel('r')
    ax.set_ylabel('F')
    fig.legend()
    plt.grid()
    plt.show()

if savefile:
    np.savez(DIR+filename,r,fEl,fIon,ftot)