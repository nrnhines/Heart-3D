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

gid2orgmap = {}
def write_morphfile(outlayer, morphfile):
  mf = open(morphfile, "w")
  for icircle in range(cellorg.ncircle[outlayer]):
    for ipt in range(cellorg.npts[outlayer][icircle]):
      x,y,z = xyz(outlayer, icircle, ipt)
      gid = cellorg.org2gid(outlayer, icircle, ipt)
      assert (gid not in gid2orgmap)
      gid2orgmap[gid] = (outlayer, icircle, ipt)
      mf.write("%d %g %g %g\n"%(gid, x, y, z))
  mf.close()

def write_spkfile(fname, spkfile): # depends on gid2orgmap from write_morphfile
  f = open(fname, "r")
  sf = open(spkfile, "w")

  for line in f:
    [t,gid] = [float(x) for x in line.split()]
    gid = int(gid)
    if gid in gid2orgmap:
      sf.write(line)

  f.close()
  sf.close()

def write_files(fname, outlayer):
  spkfile="layer%d.spk"%outlayer
  morphfile="morphology_layer%d.txt"%outlayer
  write_morphfile(outlayer, morphfile)
  timeit ("wrote %s" % (morphfile))
  print("wait...writing %s will take a while"%spkfile)
  write_spkfile(fname, spkfile)
  timeit ("wrote %s" % (spkfile))

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("spike_filename")
  parser.add_argument("layer")
  args = parser.parse_args()
  fname=args.spike_filename
  outlayer=int(args.layer)
  write_files(fname, outlayer)
