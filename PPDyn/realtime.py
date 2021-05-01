#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import ini
import argparse
import os.path
from os.path import join as pjoin

def eview(argv):
#     parser = argparse.ArgumentParser(description='Plasma Particle Dynamics (PPDyn)')
#     parser.add_argument('-i','--input', default='input.ini', type=str, help='Input file name')
#     parser.add_argument('-s','--show', action='store_true', help='Show Animation')
#     args        = parser.parse_args()
#     inputFile   = args.input
#     show_anim = args.show

    params  = ini.parse(open(argv).read())
    N       = int(params['particles']['N'])    # Number of particles
    tmax    = float(params['time']['tmax'])
    realTime    = bool(params['diagnostics']['realTime'])
    #========= Configuration ===========
    def energy(N):
        path ='data'
        # N = 100
        def animate(i):
            if os.path.exists(pjoin(path,'energy.txt')):
                time,energy = np.loadtxt(pjoin(path,'energy.txt'),unpack=True)
                ax.clear()
                ax.plot(time, energy/N)
                ax.set_xlabel("$timestep$")
                ax.set_ylabel("$Energy$")
                ax.set_title("Timestep: %d"%time[-1])

        fig,ax = plt.subplots(figsize=(6, 6))

        ani = animation.FuncAnimation(fig, animate, interval=1000)

        plt.show()

    if realTime:
        energy(N)
        plt.close()

if __name__== "__pyview__":
	eview(sys.argv[1:])
