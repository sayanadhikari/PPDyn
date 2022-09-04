#!/usr/bin/env python3

import numpy as np
from os.path import join as pjoin
import h5py
import os
import sys


try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')

file_name = "particle"

file = h5py.File(pjoin(DIR,file_name+'.hdf5'),'w')
# print(file)
# exit()
file.attrs['Nt']=100000
# file.attrs.modify('Nt',"100000");
#del file.attrs['attr5'];
file.close()
