import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace, logspace

savefile = True
doplot = True

DIR = '../staticfields/'
filename = 'Ftot'

fEl_max = 1
fdrag_max = 1

r_max = 150
x_scale = 50e-3
r = np.linspace(-r_max,r_max,r_max*2+1)
a = 0.0007
b = 0.01

fEl = -np.sign(r[:]) / (0.0+np.exp(3.0-a*r[:]*r[:])+a*0.1*r[:]*r[:])
fdrag = np.sign(r[:]) / (0.2+0.1*np.exp(5.0-b*r[:]*r[:])+b*0.1*r[:]*r[:]) 
fdrag[:] = fdrag[:]+  0.5*np.sign(r[:])/(r[:]*10+(np.sign(r[:])*0.01+0.001))
fEl = fEl_max*fEl[:]/np.max(fEl)
fdrag[r_max:-1] = -fdrag[r_max:0:-1]
fdrag = fdrag_max*fdrag[:]/np.max(fdrag)

#zero at origin
zero_radius = 25
#fdrag[r_max-zero_radius:r_max+zero_radius]=linspace(fdrag[r_max-zero_radius],fdrag[r_max+zero_radius],2*zero_radius)
fEl[r_max-zero_radius:r_max+zero_radius]=linspace(fEl[r_max-zero_radius],fEl[r_max+zero_radius],2*zero_radius)

ftot = fdrag[:] + fEl[:]
r[:] *= x_scale/r_max

if doplot:
    fig,ax = plt.subplots()
    ax.plot(r,fEl,label='$F_{EL}$',color='blue')
    ax.plot(r,fdrag,label='$F_{drag}$',color='r')
    ax.plot(r,ftot,label='$F_{total}$',color='orange')
    ax.set_xlabel('r')
    ax.set_ylabel('F')
    fig.legend()
    plt.grid()
    plt.show()

if savefile:
    np.savez(DIR+filename,r,fEl,fdrag)