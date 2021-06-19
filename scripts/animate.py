#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
#matplotlib.use('MacOSX')
# import plotly.graph_objects as go
#========= Configuration ===========
show_anim = True
save_anim = False
interval = 10 #in milliseconds

DIR ="../data"

file_name = "particle"#"rhoNeutral" #"P"

h5 = h5py.File('../data/'+file_name+'.hdf5','r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
Lz = h5.attrs["Lz"]

dp   = h5.attrs["dp"]
Nt   = h5.attrs["Nt"]


data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

if (show_anim == True):
    def animate(i):
        # file=DIR+'/data%d'%data_num[i]+'.dat'
        datax = h5["/%d"%data_num[i]+"/position/x"]
        datay = h5["/%d"%data_num[i]+"/position/y"]
        dataz = h5["/%d"%data_num[i]+"/position/z"]
        ax1.cla()
        img1 = ax1.scatter(datax,datay,dataz,marker='o',color='b',alpha=1.0,s=10)
        ax1.set_title('TimeSteps = %d'%(i*dp)+'\n Phase Space')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$y$")
        ax1.set_zlabel("$z$")
        ax1.set_xlim([-Lx, Lx])
        ax1.set_ylim([-Ly, Ly])
        ax1.set_zlim([-Lz, Lz])



if (show_anim == True):
    # fig,ax1 = plt.subplots(1,1,1, projection='3d')
    fig = plt.figure(figsize=(6, 6))
    ax1 = plt.axes(projection ="3d")
    ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval,blit=False)
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
