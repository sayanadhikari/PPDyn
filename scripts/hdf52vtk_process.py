#!/usr/bin/env python3

import numpy as np
from pyevtk.hl import pointsToVTK
from os.path import join as pjoin
import h5py
import os
import sys
sys.path.insert(1, '../src/')
from vtk_data import vtkwrite

try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')

print('Writing VTK files for Paraview visualization ...')
vtkwrite(DIR)
