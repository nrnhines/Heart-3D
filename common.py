from neuron import h
pc = h. ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())
pr_ = (rank == 0)

#rank = 0
#nhost = 8192

from time import time
firsttime = time()
lasttime = time()
def timeit(mes=None):
  global lasttime
  t = time()
  if mes and pr_:
    print ("%s time %g of %g" %(mes, t - lasttime, t - firsttime))
  lasttime = t

def pr(str):
  if (pr_):
    print (str)
