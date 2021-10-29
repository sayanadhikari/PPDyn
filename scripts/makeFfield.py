import numpy as np
import matplotlib.pyplot as plt

savefile = False
doplot = True

DIR = '../staticfields/'
filename = 'Ftot'

fEl_max = 1
fdrag_max = 1

r_max = 150
r = np.linspace(-r_max,r_max,r_max*2+1)
a = 0.001
b = 0.01

fEl = -np.sign(r[:]) / (0.0+np.exp(3.0-a*r[:]*r[:])+a*0.1*r[:]*r[:])
fdrag = np.sign(r[:]) / (0.2+np.exp(3.0-b*r[:]*r[:])+b*0.2*r[:]*r[:]) 
fdrag[:] = fdrag[:]+  0.5*np.sign(r[:])/(r[:]*10+(np.sign(r[:])*0.01+0.001))
fEl = fEl_max*fEl[:]/np.max(fEl)
fdrag = fdrag_max*fdrag[:]/np.max(fdrag)
ftot = fdrag[:] + fEl[:]

if doplot:
    fig,ax = plt.subplots()
    ax.plot(r,fEl,label='$F_{EL}$',color='blue')
    ax.plot(r,fdrag,label='$F_{drag}$',color='r')
    ax.plot(r,ftot,label='$F_{total}$',color='orange')
    ax.set_xlabel('r')
    ax.set_ylabel('F')
    fig.legend()
    plt.show()

if savefile:
    np.savez(DIR+filename,r,fEl,fdrag)