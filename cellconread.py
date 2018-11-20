from common import h, pc, rank, nhost

gidinfo = {}
connections = {}
ncell = 0
ncon = 0

def gid2rank(gid): # RR distribution
  return gid%nhost

def cellconread():
  # new Heart-3D paraboloid organization
  global ncon, ncell
  import cellorg, mkgap

  ncell = cellorg.ngid
  for gid in range(rank, ncell, nhost):
    ilayer, icircle, ipt = cellorg.gid2org(gid)
    if icircle < cellorg.ncircle[ilayer] - 1:
      xyz = cellorg.xyz(ilayer, icircle, ipt)
      if cellorg.is_simulated(xyz):
        gidinfo[gid] = xyz

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
  print ("%d %d" % (len(gidinfo), len(mkgap.gaps)))















