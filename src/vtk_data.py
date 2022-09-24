import numpy as np
from pyevtk.hl import pointsToVTK
from os.path import join as pjoin
import h5py
import os

def vtkwrite(path):
    file_name = "particle"#"rhoNeutral" #"P"
    if os.path.exists(pjoin(path,'vtkdata')) == False:
        os.mkdir(pjoin(path,'vtkdata'))
    h5 = h5py.File(pjoin(path,file_name+'.hdf5'),'r')

    Lx = h5.attrs["Lx"]
    Ly = h5.attrs["Ly"]
    Lz = h5.attrs["Lz"]

    dp   = h5.attrs["dp"]
    Nt   = h5.attrs["Nt"]

    data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

    datapos = h5["/position"]

    for i in range(len(data_num)):
        datax = np.array(datapos[i,:,0])
        datay = np.array(datapos[i,:,1])
        dataz = np.array(datapos[i,:,2])
        pointsToVTK(pjoin(path,'vtkdata','points_%d'%i), datax, datay, dataz)
