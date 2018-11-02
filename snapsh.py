# snapshots of all the soma voltages beginning at snapshot_begin and at
# snapshot_interval thereafter. That data is stored in file snapsh.txt

snapshot_begin = 100
snapshot_interval = 20

from common import h, pc, rank, nhost
from net3d import gidinfo

somav_ptrvec = None
somav_vec = None

vcnts = None
def send_to_master(srcvec):
  global vcnts
  if vcnts is None:
    vcnts = h.Vector(nhost)
    vcnts.x[0] = len(srcvec)
  destvec = h.Vector()
  pc.alltoall(srcvec, vcnts, destvec)
  pc.barrier()
  return destvec

def snapsh_init():
  global somav_ptrvec, somav_vec, allgid
  n = len(gidinfo) #or else count the segments
  somav_ptrvec = h.PtrVector(n)
  somav_vec = h.Vector(n)
  gidvec = h.Vector(n)
  for i, (gid, cellinfo) in enumerate(gidinfo.items()):
    seg = cellinfo.cell.soma(0.5)
    somav_ptrvec.pset(i, seg._ref_v)
    gidvec.x[i] = gid
  allgid = send_to_master(gidvec)
  h.cvode.event(h.t + snapshot_begin, snapshout)

def snapsh():
  # print 'snapsh t=%g' % pc.t(0)
  somav_ptrvec.gather(somav_vec)

def snapsh_setup():
  global fih
  fih = h.FInitializeHandler(1, snapsh_init)

def snapshout():
  somav_ptrvec.gather(somav_vec)
  allsomav = send_to_master(somav_vec)
  if rank is 0:
    f = open('snapsh.txt', 'a')
    f.write('%d %g\n'%(len(allsomav), h.t))
    for i in range(len(allsomav)):
      line = '%d'%i
      line += ' %g' % allsomav.x[i]
      f.write('%d %g\n'%(allgid.x[i], allsomav.x[i]))
    f.close()
  h.cvode.event(h.t + snapshot_interval, snapshout)

if __name__ == '__main__':
  # test() invalid test
  for i in range(n):
    print i,
    print somav_vec.x[i],
    print ""

