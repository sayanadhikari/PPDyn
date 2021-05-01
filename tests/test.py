from PPDyn import ppdyn
from PPDyn.ppdplot import pyview
import time
from multiprocessing import Process
from PPDyn.realtime import eview


input = 'input.ini'
realTime = True


if realTime:
    p = Process(target=eview(input))
    p.start()
    # p.join()

start = time.time()
ppdyn(input)
end = time.time()

print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
pyview(input)
