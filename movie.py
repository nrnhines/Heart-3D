from neuron import h, gui
import time
import cellorg
from cellorg import xyz, gid2org

def raster(fname):
  f = open(fname)
  ras = []
  for line in f:
    [t,gid] = [float(x) for x in line.split()]
    gid = int(gid)
    ras.append((t, gid))
  f.close()
  return ras

def movie(r, g1, g2):
  tt = 0.0
  c = 1.
  for t, gid in r:
    x, y, z = xyz(*gid2org(gid))
    if t > tt:
      time.sleep(0.1)
      tt += 1.0
      c = 1 + tt/10.
      #if int(tt)%50 == 0: g2.erase(); g1.erase()
    g1.mark(x, z, "S", 10, c, 1)
    g2.mark(y, z, "S", 10, c, 1)

if __name__ == "__main__":
  ras = raster("spk000.dat")
  g1 = h.Graph()
  g1.size(-14000, 15000, 48000, 50000)
  g2 = h.Graph()
  g2.size(-2000, 2000, 48000, 50000)
  movie(ras, g1, g2)
