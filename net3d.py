from common import h, pc, rank, nhost, timeit, pr
from cellconread import cellconread, gidinfo, connections, ncell, ncon
from purkspec import is_purkinje_gap

import snapsh
import ecg
from gjrecord import gj_record, gj_out
from cellorg import gid2org, xyz, nlayer, ncircle, npts
import mkgap
from math import pi
from util import isclose
import param

h.load_file("verifygap.hoc")

class CellInfo:
  def __init__(self, cell):
    self.cell = cell
    self.gaps = {} # srcgid:HalfGap
    self.is_purk = False

def mkcells(gidinfo):
  timeit()
  for gid in gidinfo:
    x,y,z = gidinfo[gid]
    cell = h.Cell()
    gidinfo[gid] = CellInfo(cell)

    # cell shape is actually an arc and area is fastidious with respect
    # to all 6 sides. But length
    # treated as line distance between org points (interior corners
    # in circumferential direction. Set diam so area is correct including
    # end areas.
    cell.soma.pt3dclear()
    cell.soma.pt3dadd(x,y,z, 1.)
    ilayer, icircle, ipt = gid2org(gid)
    x1,y1,z1 = xyz(ilayer, icircle, ipt + 1)
    cell.soma.pt3dadd(x1, y1, z1, 1.)
    length = cell.soma.L
    area = sum(mkgap.cell_side_areas(gid))
    diam = area/pi/length
    cell.soma.diam = diam
    assert(isclose(cell.soma(.5).area(), area, abs_tol=area*1e-5))

    cell.position(x,y,z)
    pc.set_gid2node(gid, rank)
    nc = cell.connect2target(None)
    pc.cell(gid, nc)
  x = pc.allreduce(len(gidinfo), 1)
  pr("Global number of real cells is %d"%x)
  timeit("mkcells")

def mkgaps(gidinfo, gaps):
  timeit()
  mark = set()
  for gapinfo in gaps.values():
    gg = (gapinfo.gid1, gapinfo.gid2)
    id = gapinfo.id
    mkhalfgap(gg[0], gg[1], id, gidinfo, mark)
    mkhalfgap(gg[1], gg[0], -id, gidinfo, mark)
  pc.setup_transfer()

  x = 0
  for cell in gidinfo.values():
    x += len(cell.gaps)
  x = pc.allreduce(x, 1)
  pr("Global number of halfgap is %d"%x)

  timeit("mkgaps")

def mkhalfgap(gid1, gid2, id, gidinfo, mark):
  # sgid is the gid for the voltage since single compartment
  global ncon

  if gid1 in gidinfo:
    cell1 = gidinfo[gid1]
    cell = cell1.cell
    gap = h.HalfGap(cell.soma(.5))
    gap.id = id
    gap.id1 = gid1
    gap.id2 = gid2
    gap.g = 0.0
    pc.target_var(gap, gap._ref_vgap, gid2)
    assert(gid2 not in cell1.gaps)
    cell1.gaps[gid2] = gap
    if gid1 not in mark:
      pc.source_var(cell.soma(.5)._ref_v, gid1, sec=cell.soma)
      mark.add(gid1)

def special_gap_params():
  try:
    f = open('Connections_Other_Information.txt')
  except:
    pr('No special gap parameters')
    return
  gidpair2ncon = {}
  for ncon in connections:
    gidpair2ncon[connections[ncon]] = ncon
  for line in f:
    #gid1, gid2, gmin, gmax, gvar, tc, tcvar, drift
    info = [float(x) for x in line.split()]
    for i in range(2):
      info[i] = int(info[i])
    for i in range(2):
      special_gap_params2(i, info, gidpair2ncon)
  f.close()

def special_gap_params2(i, info, gidpair2ncon):
  gid = info[i]
  if gid in gidinfo:
    gaps = gidinfo[gid].gaps
    # note assumption that pair order is same in info and connections
    pair = (info[0], info[1])
    if pair in gidpair2ncon:
      ncon = gidpair2ncon[pair]
    elif (info[1], info[0]) in gidpair2ncon:
      ncon = gidpair2ncon[(info[1], info[0])]
    else:
      print ("%d pair %s not found"%(rank, pair.__str__()))
      return
    #print ("%d %s"%(rank, connections[ncon].__str__()))
    gap = gaps[ncon]
    gmin, gmax, gvar, tc, tcvar, drift = info[2:]
    gap.meang = gmax
    gap.gmax = gmax
    gap.gmin = gmin
    gap.g = gmax
    gap.rg = gvar
    gap.meant = tc
    gap.rt = tcvar
    gap.drift=drift

def setallgaps(meang, interval, drift):
  npurkgap = 0
  for gid1, cellinfo in gidinfo.items():
    for gid2, gap in cellinfo.gaps.items():
      gapinfo = mkgap.gaps[(gid1, gid2) if gid1 < gid2 else (gid2, gid1)]
      area = gapinfo.area
      g = mkgap.abscond(area, meang)
      if is_purkinje_gap(gapinfo.gid1, gapinfo.gid2):
        g *= param.purkinje_gap_factor
        cellinfo.is_purk = True
        npurkgap += 1
      gap.meang = g
      gap.gmax = g
      gap.gmin = g
      gap.g = g
      gap.rg = interval
      gap.drift=drift
  npurkgap = pc.allreduce(npurkgap, 1)
  pr("number of purkinje gaps is %d" % (npurkgap/2))

def mknet():
  h.load_file("cell.hoc")

  timeit()
  cellconread()
  timeit("cellconread makes gapinfo")

  mkcells(gidinfo)
  mkgaps(gidinfo, mkgap.gaps)
  setallgaps(param.meang, 1000.0, 0.0)
  #special_gap_params()
  h.verifyHalfGap()


def mkmodel():
  mknet()
  #snapsh.snapsh_setup()
  ecg.ecg_setup()

# all cells in last circle of layer 0
def circlestim():
  #return vertstim()
  from cellorg import nlayer, ncircle, npts, org2gid
  r = []
  ilayer = 0
  for icircle in range(ncircle[ilayer])[-2:]:
    npt = npts[ilayer][icircle]
    for ipt in range(npt):
      r.append(org2gid(ilayer, icircle, ipt))
    if rank == 0:
      print ("circlestim (%d %d [0:%d])" % (ilayer, icircle, npt))
  r = h.Vector(r)
  return r

# all ipt=0 cells in all circles of layer 0
def vertstim():
  from cellorg import nlayer, ncircle, npts, org2gid
  r = []
  ilayer = 0
  for icircle in range(ncircle[ilayer]):
    for ipt in range(2):
      r.append(org2gid(ilayer, icircle, ipt))
  if rank == 0:
    print ("vertstim (%d [0:%d] [0,1])" % (ilayer, ncircle[ilayer]))
  r = h.Vector(r)
  return r

# first gid in gidinfo of rank 1 and all its adjacent gids.
def purkstim():
  r = []
  if rank == 0:
    for gid in gidinfo: # break, so just the first and all its connections
      r.append(gid)
      for gid2 in gidinfo[gid].gaps:
        r.append(gid2)
      break
  pr("purkstim" + str(r))
  r = h.Vector(r)
  pc.broadcast(r, 0)
  return r
  

def test1():
  pc.barrier()
  for r in range(nhost):
    if r == rank:
      for gid1 in gidinfo:
        cellinfo = gidinfo[gid1]
        for gid2, gap in cellinfo.gaps.items():
          print ("%d %d %g" %(gid1, gid2, gap.g))
    pc.barrier()

def test2():
  from matplotlib import pyplot
  pyplot.hist([c.cell.soma(.5).area() for c in gidinfo.values()])
  pyplot.show()
  pyplot.hist([c.cell.soma.L for c in gidinfo.values()])
  pyplot.show()
  pyplot.hist([gap.g for cellinfo in gidinfo.values() for gap in cellinfo.gaps.values()])
  pyplot.show()

def showpurk(): # only for nhost=1
  global g
  g = h.Graph()
  g.size(0, npts[0][-1], 0, ncircle[0])
  for gid, ci in gidinfo.items():
    if ci.is_purk:
      ilayer, icircle, ipt = gid2org(gid)
      g.mark(ipt, icircle, "S", 10)

if __name__ == '__main__':
  mknet()
  #test1()
  #test2()
  #purkstim()
  circlestim()
  if pc.nhost() > 1:
    pc.barrier()
    h.quit()
