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
  readcells('../heart3d/morphology.txt')
  readconnections('../heart3d/connections.txt')

if __name__ == '__main__':
  cellconread()
  print len(gidinfo), len(connections)















