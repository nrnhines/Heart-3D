from neuron import h, gui
import sys
this_module = sys.modules[__name__]
from common import timeit
import cellorg
from cellorg import xyz, gid2org, npts
timeit('import cellorg')

def raster(fname):
  f = open(fname)
  ras = []
  for line in f:
    [t,gid] = [float(x) for x in line.split()]
    gid = int(gid)
    ras.append((t, gid))
  f.close()
  return ras

dtt = 2.0
grp = 0
sgrp = 0

pras = None
def pras_(ras):
  global tt
  p = []
  tt = 0.0
  #segregate into groups of duration dtt
  for t, gid in ras:
    ilayer, icircle, ipt = gid2org(gid)
    if ilayer == 0:
      while t >= tt:
        tt += dtt
        pp = []
        p.append(pp)
      pp.append((t, (ipt, icircle)))
  return p

def pl(pp):
    global tt, grp
    grp = round(grp, 0)
    if len(pp) > 0:
      tt = int(pp[0][0])
    else:
      tt = grp*dtt
    ttstr[0]='tt=%g'%tt
    for p in pp:
      ipt = p[1][0]
      icircle = p[1][1]
      g.mark(ipt, icircle, "S", 10, 1, 1)

def ttcallback():
  global sgrp
  g.erase()
  n = len(pras)
  sgrp = grp/n*100
  ig = int(grp)
  pl(pras[ig if ig < n else n-1])

def ttscallback():
  global grp
  grp = round(sgrp*len(pras)/100, 0)
  ttcallback()

def movie(pras):
  global grp
  for grp in range(len(pras)):
    ttcallback()

def getdat(fname):
  import pickle
  pfile = fname + ".pkl"
  try:
    1/0 # need to check dates between fname and fname.pkl
    p = pickle.load(open(pfile, "rb"))
    timeit("pickle load")
  except:
    ras = raster(fname)
    timeit("input raster")
    p = pras_(ras)
    timeit ("construct pras")
    pickle.dump(p, open(pfile, "wb"))
    timeit("pickle dump")
  return p

def bld(x, y):
  global g, tt, ttstr, grp, sgrp, pras
  tt = 0.0
  ttstr = h.ref('tt=0.000000000')
  v = h.VBox()
  v.intercept(1)
  g = h.Graph()
  g.size(0, x, 0, y)
  h.xpanel("", 1)
  h.xvarlabel(ttstr)
  h.xvalue("group", (this_module, "grp"), 1, ttcallback)
  h.xslider((this_module, "sgrp"), 0, 100, ttscallback)
  h.xpanel()
  v.intercept(0)
  v.map("wave front", 100, 100, 700, 500)
  return v

if __name__ == "__main__":
  vb = bld(npts[0][-1], len(npts[0]))
  fname="spkpurkmean1500.dat"
  pras = getdat(fname)
  movie(pras)
