import numpy as np
from pyevtk.hl import pointsToVTK
from os.path import join as pjoin
import h5py
import os

def vtkwrite():
    file_name = "particle"#"rhoNeutral" #"P"
    if os.path.exists(pjoin('data','vtkdata')) == False:
        os.mkdir(pjoin('data','vtkdata'))
    h5 = h5py.File(pjoin('data',file_name+'.hdf5'),'r')

    Lx = h5.attrs["Lx"]
    Ly = h5.attrs["Ly"]
    Lz = h5.attrs["Lz"]

    dp   = h5.attrs["dp"]
    Nt   = h5.attrs["Nt"]

    data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

    for i in range(len(data_num)):
        datax = h5["/%d"%data_num[i]+"/position/x"]
        datay = h5["/%d"%data_num[i]+"/position/y"]
        dataz = h5["/%d"%data_num[i]+"/position/z"]
        datax = np.array(datax)
        datay = np.array(datay)
        dataz = np.array(dataz)
        pointsToVTK(pjoin('data','vtkdata','points_%d'%i), datax, datay, dataz)
