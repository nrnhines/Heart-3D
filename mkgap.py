import param
from morphdef import distance
from cellorg import gid2org, org2gid, overlap, nlayer, ncircle, npts, xyz


gaps = {}

def isclose(a, b, abs_tol=1e-9):
  return abs(a - b) < abs_tol

class GapInfo:
  def __init__(self, gid1, gid2, f1, f2, n):
    self.gid1, self.gid2, self.f1, self.f2 = (gid1, gid2, f1, f2) if gid1 < gid2 else (gid2, gid1)
    self.ncellsep = n

  def same(self, gid1, gid2, f1, f2, n):
    g1, g2, f1, f2 = (gid1, gid2, f1, f2) if gid1 < gid2 else (gid2, gid1, f2, f1)
    if  gid1 != self.gid1 or gid2 != self.gid2 or n != self.ncellsep:
      return False
    if isclose(f1, self.f1) and isclose(f2, self.f2):
      return True
    return False

  def __repr__(self):
    return "(%d %d %g %g %d)\n" %(self.gid1, self.gid2, self.f1, self.f2, self.ncellsep)

def set_gap(g1, g2, f1, f2, n):
    g1, g2, f1, f2 = (g1, g2, f1, f2) if g1 < g2 else (g2, g1, f2, f1)
    key = (g1, g2)
    if key in gaps:
      g = gaps[key]
      assert(g.same(g1, g2, f1, f2, n))
      return g
    g = GapInfo(g1, g2, f1, f2, n)
    gaps[key] = g
    return g;

def get_gap(gid1, gid2):
    g = (gid1, gid2) if gid1 < gid2 else (gid2, gid1)
    if g in gaps:
      return gaps[g]
    return None

def gaps_for_gid(gid):
  o = gid2org(gid)
  ilayer, icircle, ipt = o
  gs = []

  #circum coordinate, end to end
  npt = npts[ilayer][icircle]
  for jpt in [(ipt - 1)%npt, (ipt + 1)%npt]:
    g2 = org2gid(ilayer, icircle, jpt)
    gs.append(set_gap(gid, g2, 1.0, 1.0, 0))

  # between layers (assume icircle the same)
  n = int(param.thickness/(nlayer - 1)/param.cell_diameter)
  for jlayer in [ilayer - 1, ilayer + 1]:
    if jlayer >= 0 and jlayer < nlayer:
      a = overlap(o, jlayer, icircle)
      for b in a:
        g2 = org2gid(jlayer, icircle, b[0])
        gs.append(set_gap(gid, g2, b[1], b[2], n))

  # between circles
  for jcircle in [icircle - 1, icircle + 1]:
    if jcircle >= 0 and jcircle < ncircle[ilayer]:
      # how many cells between icircle and jcircle
      d = distance(xyz(ilayer, icircle, 0), xyz(ilayer, jcircle, 0))
      n = int(d/param.cell_diameter)
      a = overlap(o, ilayer, jcircle)
      for b in a:
        g2 = org2gid(ilayer, jcircle, b[0])
        gs.append(set_gap(gid, g2, b[1], b[2], n))

  return gs

def test1(ilayer, icircle, ipt):
  print (gaps_for_gid(org2gid(ilayer, icircle, ipt)))

if __name__ == "__main__":
  test1(0,0,0)
