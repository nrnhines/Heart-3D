from common import h, pc, rank, nhost, timeit, pr
from cellconread import cellconread, gidinfo, connections, ncell, ncon
import snapsh
import ecg
from gjrecord import gj_record, gj_out
import mkgap

class CellInfo:
  def __init__(self, cell):
    self.cell = cell
    self.gaps = {} # srcgid:HalfGap

def mkcells(gidinfo):
  timeit()
  for gid in gidinfo:
    x,y,z = gidinfo[gid]
    cell = h.Cell()
    gidinfo[gid] = CellInfo(cell)
    cell.position(x,y,z)
    pc.set_gid2node(gid, rank)
    nc = cell.connect2target(None)
    pc.cell(gid, nc)
  x = pc.allreduce(len(gidinfo), 1)
  pr("Global number of real cells is %d"%x)
  timeit("mkcells")

def mkgaps(gidinfo, connections):
  timeit()
  mark = set()
  for cid in connections:
    gid1, gid2 = connections[cid]
    gapinfo = mkgap.get_gap(gid1, gid2)
    mkhalfgap(gid1, gid2, 1, gapinfo, gidinfo, mark)
    mkhalfgap(gid2, gid1, -1, gapinfo, gidinfo, mark)
  pc.setup_transfer()

  x = 0
  for cell in gidinfo.values():
    x += len(cell.gaps)
  x = pc.allreduce(x, 1)
  pr("Global number of halfgap is %d"%x)

  timeit("mkgaps")

def mkhalfgap(gid1, gid2, gid2_polarity, gapinfo, gidinfo, mark):
  # sgid is the gid for the voltage since single compartment
  global ncon

  if gid1 in gidinfo:
    cell1 = gidinfo[gid1]
    cell = cell1.cell
    gap = h.HalfGap(cell.soma(.5))
    gap.id = gid2_polarity
    gap.g = gapinfo.g
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
  for cellinfo in gidinfo.values():
    for gap in cellinfo.gaps.values():
      gap.meang = meang
      gap.gmax = meang
      gap.gmin = meang
      gap.g = meang
      gap.rg = interval
      gap.drift=drift


def mknet():
  h.load_file("cell.hoc")

  timeit()
  cellconread()
  timeit("cellconread makes gapinfo")

  mkcells(gidinfo)
  mkgaps(gidinfo, connections)
  setallgaps(30.0, 1000.0, 0.0)
  special_gap_params()


def mkmodel():
  mknet()
  #snapsh.snapsh_setup()
  ecg.ecg_setup()

def test1():
  pc.barrier()
  for gid1 in gidinfo:
    cellinfo = gidinfo[gid1]
    for gid2, gap in cellinfo.gaps.items():
      print ("%d %d %g" %(gid1, gid2, gap.g))
  pc.barrier()

if __name__ == '__main__':
  mknet()
  #test1()
  if pc.nhost() > 1:
    pc.barrier()
    h.quit()
