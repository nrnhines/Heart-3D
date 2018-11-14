from common import h, pc, rank, nhost

gidinfo = {}
connections = {}
ncell = 0
ncon = 0

def gid2rank(gid): # RR distribution
  return gid%nhost

def readcells(fname): #read and save positions (according to gid2rank distribution)
  global ncell
  f = open(fname)
  ncell = 0
  for line in f:
    a = line.split()
    gid = int(a[0])
    ncell += 1
    if gid2rank(gid) == rank:
      gidinfo[gid] = (float(a[1]), float(a[2]), float(a[3]))
  f.close()
  
def readconnections(fname):
  global ncon
  f = open(fname)
  ncon = 0
  for line in f:
    a= line.split()
    gid1, gid2 = int(a[0]), int(a[1])
    ncon += 1
    if gid1 in gidinfo or gid2 in gidinfo:
      connections[ncon] = (gid1, gid2)
  f.close()

def cellconread():
  #readcells('../heart3d/morphology.txt')
  #readconnections('../heart3d/connections.txt')

  # new Heart-3D paraboloid organization
  global ncon, ncell
  import cellorg, mkgap

  ncell = cellorg.ngid
  for gid in range(rank, ncell, nhost):
    gidinfo[gid] = cellorg.xyz(*cellorg.gid2org(gid))

  for gid in gidinfo:
    mkgap.gaps_for_gid(gid)
  ncon = 0 # unfortunately, no longer global number. Just for this rank
  for g in mkgap.gaps:
    if g not in connections:
      connections[ncon] = g
      ncon += 1

if __name__ == '__main__':
  cellconread()
  print ("%d %d" % (len(gidinfo), len(connections)))















