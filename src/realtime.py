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


parser = argparse.ArgumentParser(description='Plasma Particle Dynamics (PPDyn)')
parser.add_argument('-i','--input', default='input.ini', type=str, help='Input file name')
args        = parser.parse_args()
inputFile   = args.input

params = ini.parse(open(inputFile).read())
N       = int(params['particles']['N'])    # Number of particles
tmax    = float(params['time']['tmax'])

energy(N)
