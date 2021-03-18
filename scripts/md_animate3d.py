#!/usr/bin/ python

import numpy as np
import matplotlib.pyplot as plt
# import random
# import math
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
# import plotly.graph_objects as go
#========= Configuration ===========
show_anim = True
save_anim = True
interval = 0.01#in seconds
datadt = 10
DIR ="../data"
Lx = 10.0
Ly = 10.0
Lz = 10.0

tmax = 50.0
dt   = 0.010
Nt   = round(tmax/dt)


data_num = np.arange(start=0, stop=Nt, step=datadt, dtype=int)

if (show_anim == True):
    def animate(i):
        file=DIR+'/data%d'%data_num[i]+'.dat'
        datax,datay,dataz = np.loadtxt(file, unpack=True)
        ax1.cla()
        img1 = ax1.scatter(datax,datay,dataz,marker='o',color='b',alpha=1.0,s=10)
        ax1.set_title('TimeSteps = %d'%i+'\n Phase Space')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$y$")
        ax1.set_xlim([-Lx, Lx])
        ax1.set_ylim([-Ly, Ly])
        ax1.set_ylim([-Lz, Lz])



if (show_anim == True):
    # fig,ax1 = plt.subplots(1,1,1, projection='3d')
    fig = plt.figure(figsize=(6, 6))
    ax1 = plt.axes(projection ="3d")
    ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval*1e+3,blit=False)
    # ani.save('phase_space.gif',writer='imagemagick')
    plt.show()
    if(save_anim == True):
        try:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=(1/interval), metadata=dict(artist='Me'), bitrate=1800)
        except RuntimeError:
            print("ffmpeg not available trying ImageMagickWriter")
            writer = animation.ImageMagickWriter(fps=(1/interval))
        ani.save('mdanimation3d.mp4')
