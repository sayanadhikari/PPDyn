#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin
import sys
import matplotlib as mp
# import plotly.graph_objects as go
#========= Configuration ===========

try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')

file_name = "particle"#"rhoNeutral" #"P"

h5 = h5py.File(pjoin(DIR,file_name+'.hdf5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
Lz = h5.attrs["Lz"]

N = h5.attrs["N"]

dp   = h5.attrs["dp"]
Nt   = h5.attrs["Nt"]
nx,ny,nz = 256,256,50
dx,dy,dz = (2*Lx/nx),(2*Ly/ny), (2*Lz/nz)

data_num = np.arange(start=0, stop=Nt, step=1, dtype=int)

time = data_num*dp
i = 1000
x = h5["/%d"%data_num[i]+"/position/x"]
y = h5["/%d"%data_num[i]+"/position/y"]
z = h5["/%d"%data_num[i]+"/position/z"]
x = np.array(x)
y = np.array(y)
z = np.array(z)

vx = h5["/%d"%data_num[i]+"/velocity/vx"]
vy = h5["/%d"%data_num[i]+"/velocity/vy"]
vz = h5["/%d"%data_num[i]+"/velocity/vz"]
vx = np.array(vx)
vy = np.array(vy)
vz = np.array(vz)
print(len(x))
with open(pjoin(DIR,'velocity.vti'), 'a') as file:
    file.write('{0:s}\n'.format("<VTKFile type=\"ImageData\">"))
    file.write('{0:s}{1:f} {2:f} {3:f}{4:s}\n'.format("<ImageData Origin= \"",-Lx, -Ly, -Lz,"\""))
    file.write('{0:s}{1:f} {2:f} {3:f}{4:s}\n'.format("Spacing= ",dx, dy, dz,"\""))
    file.write('{0:s}{1:f} {2:s} {3:f} {4:s} {5:f}{6:s}\n'.format("WholeExtent= \"0 ",(nx-1), "0", (ny-1), "0", (nz-1),"\">"))
    file.write('{0:s}\n'.format("<PointData>"))
    file.write('{0:s}\n'.format("<DataArray Name=\"object\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">"))
    for k in range(len(x)):
        for j in range(len(y)):
            for i in range(len(z)):
                print(i,j,k)
                file.write('{0:f} {1:f} {2:f}\n'.format(x[i],y[j],z[k]))

    file.write('{0:s}\n'.format("</DataArray>"))
    file.write('{0:s}\n'.format("</PointData>"))
    file.write('{0:s}\n'.format("</ImageData>"))
    file.write('{0:s}\n'.format("</VTKFile>"))
