# import sys
# sys.path.insert(1, '../PPDyn')
# from __init__ import *
from PPDyn import ppdyn
from PPDyn.ppdplot import pyview
import time
from multiprocessing import Process
# import multiprocessing, logging
# from PPDyn.realtime import eview
from PPDyn.dash_ppdyn import dash_energy
import subprocess


input = 'input.ini'

# p = Process(target=dash_energy(input))
# p.start()
def dash_run():
    subprocess.call('xterm -e python3 ediag.py', cwd='./', shell =True)

p = Process(target=dash_run())
p.start()

start = time.time()
ppdyn(input)
end = time.time()
# p.join()

print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
# p.join()
# p.close()
pyview(input)
