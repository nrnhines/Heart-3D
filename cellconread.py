from common import h, pc, rank, nhost
from math import pi
import mkgap

gidinfo = {}
connections = {}
ncell = 0
ncon = 0

def gid2rank(gid): # RR distribution
  return gid%nhost

def cellconread():
  # new Heart-3D paraboloid organization
  global ncon, ncell, connections
  import cellorg, mkgap
  from cellorg import sim_layers, sim_circles
  ncell = cellorg.ngid
  #old way iterating over all possible cells takes 5.4 seconds
  for gid in range(rank, ncell, nhost):
    ilayer, icircle, ipt = cellorg.gid2org(gid)
    if icircle < cellorg.ncircle[ilayer] - 1:
      if cellorg.is_simulated(ilayer, icircle, ipt):
        xyz = cellorg.xyz(ilayer, icircle, ipt)
        gidinfo[gid] = xyz
  '''
  #new way iterating only over cells that exist takes
  import param as p
  for ilayer in sim_layers:
    for icircle in sim_circles[ilayer]:
      i0 = cellorg.angle2ipt(p.simulation_angledeg[0]*2*pi/360, ilayer, icircle)
      i1 = cellorg.angle2ipt(p.simulation_angledeg[1]*2*pi/360, ilayer, icircle)
      for ipt in range(i0, i1+1):
        if cellorg.is_simulated(ilayer, icircle, ipt):
          gid = cellorg.org2gid(ilayer, icircle, ipt)
          if gid%nhost == rank:
            gidinfo[gid] = cellorg.xyz(ilayer, icircle, ipt)
  '''

  for gid in gidinfo:
    # because of floating round-off error which may or may not create
    # a gap with area close to 0, guarantee gap pairs by only creating
    # gaps where gid1 < gid2
    mkgap.gaps_for_gid(gid)
  # for parallel, copy gid2 gaps to ranks that need them
  mkgap.gaps_gid2_copy()
  connections = mkgap.gaps

if __name__ == '__main__':
  cellconread()
  print ("%d %d %d" % (rank, len(gidinfo), len(mkgap.gaps)))















