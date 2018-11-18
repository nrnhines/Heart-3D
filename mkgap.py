import param
from morphdef import distance
from cellorg import gid2org, org2gid, angle_overlap, nlayer, ncircle, npts, xyz
from cellorg import paraboloid, ipt2angle, gid_is_simulated
from util import side_area, area3pt

def conductance_density_circum(ilayer, icircle, ipt):
  return 1.0

def conductance_density_layer(ilayer, icircle, ipt):
  return 1.0

def conductance_density_parabola(ilayer, icircle, ipt):
  return 1.0

def abscond(a, g):
  return a*g

def area_circum(ilayer, icircle):
  assert(icircle < ncircle[ilayer] - 1)
  pinfo = paraboloid[ilayer][icircle]
  pinfo1 = paraboloid[ilayer][icircle+1]
  b = pinfo[2].p1b[1] if pinfo[2].p1b else []
  top = [pinfo[1]] + b + [pinfo1[1]]
  b = pinfo[2].p0b[1] if pinfo[2].p0b else []
  bottom = [pinfo[0]] + b + [pinfo1[0]]
  area = 0.0
  for i in range(len(top) - 1):
    area += area3pt([bottom[0], top[i], top[i+1]])
  for i in range(len(bottom) - 1):
    area += area3pt([top[-1], bottom[i], bottom[i+1]])
  return area

gaps = {}

def isclose(a, b, abs_tol=1e-7):
  return abs(a - b) < abs_tol

class GapInfo:
  def __init__(self, gid1, gid2, g):
    self.gid1, self.gid2 = (gid1, gid2) if gid1 < gid2 else (gid2, gid1)
    self.g = g

  def same(self, gid1, gid2, g):
    g1, g2 = (gid1, gid2) if gid1 < gid2 else (gid2, gid1)
    if  g1 != self.gid1 or g2 != self.gid2:
      return False
    if isclose(g, self.g):
      return True
    else:
      print ("same (%d,%d) %.20g   (%d, %d) %.20g" %(gid1, gid2, g, self.gid1, self.gid2, self.g))
    return False

  def __repr__(self):
    return "(%d %d %g)\n" %(self.gid1, self.gid2, self.g)

def set_gap(g1, g2, g):
    g1, g2 = (g1, g2) if g1 < g2 else (g2, g1)
    key = (g1, g2)
    if key in gaps:
      gap = gaps[key]
      assert(gap.same(g1, g2, g))
      return gap
    gap = GapInfo(g1, g2, g)
    gaps[key] = gap
    return gap;

def get_gap(gid1, gid2):
    key = (gid1, gid2) if gid1 < gid2 else (gid2, gid1)
    if key in gaps:
      return gaps[key]
    return None

def gaps_for_gid(gid):
  if not gid_is_simulated(gid):
      return None
  o = gid2org(gid)
  ilayer, icircle, ipt = o
  pinfo = paraboloid[ilayer][icircle]
  rf = pinfo[2] # RegionFace
  a = ipt2angle(1, ilayer, icircle)
  afirst= a*ipt
  alast = a*(ipt+1)
  gs = []

  #circum coordinate, end to end
  if icircle < ncircle[ilayer] - 1:
    npt = npts[ilayer][icircle]
    area = area_circum(ilayer, icircle)
    for jpt in [ipt - 1, ipt + 1]:
      dens_ipt = ipt if jpt < ipt else jpt%npts[ilayer][icircle]
      g = conductance_density_circum(ilayer, icircle, dens_ipt)
      g2 = org2gid(ilayer, icircle, jpt)
      if gid_is_simulated(g2): gs.append(set_gap(gid, g2, abscond(area, g)))

  # between layers
  if icircle < ncircle[ilayer] - 1:
    pinfo1 = paraboloid[ilayer][icircle+1]
    for jlayer in [ilayer - 1, ilayer + 1]:
      if jlayer >= 0 and jlayer < nlayer:
        # usually 2, sometimes 1, and rarely 3, circles in jlayer overlap icircle
        # n = int(param.layer_thickness/param.cell_diameter)
        n = 1

        #jcircle, b = rf.p0b if jlayer < ilayer else rf.p1b
        p0, plast, (jcircle, b) = (pinfo[0], pinfo1[0], rf.p0b) if jlayer < ilayer else (pinfo[1], pinfo1[1], rf.p1b)

        for p1 in b + [plast]:
          jpt, angles = angle_overlap(o, jlayer, jcircle)
          a0 = afirst
          for a1 in angles + [alast]:
            area = side_area(p0, p1, a1 - a0)
            g2 = org2gid(jlayer, jcircle, jpt)
            dens_layer, dens_circle, dens_ipt = (ilayer, icircle, ipt) if jlayer < ilayer else (jlayer, jcircle, jpt)
            g = conductance_density_layer(dens_layer, dens_circle, dens_ipt)
            if gid_is_simulated(g2) and jcircle < ncircle[jlayer] - 1: gs.append(set_gap(gid, g2, abscond(area, g)))
            jpt += 1
            a0 = a1
          jcircle += 1
          p0 = p1

  # between circles in same layer
  jlayer = ilayer
  for jcircle in [icircle - 1, icircle + 1]:
    if jcircle >= 0 and jcircle < ncircle[jlayer]:
      # how many cells between icircle and jcircle
      d = distance(xyz(ilayer, icircle, 0), xyz(jlayer, jcircle, 0))
      n = int(d/param.cell_diameter)

      jpt, angles = angle_overlap(o, jlayer, jcircle)
      a0 = afirst
      pinfo1 = paraboloid[ilayer][jcircle]
      p0, p1, dens_circle, dens_ipt = (pinfo[0], pinfo[1], icircle, ipt) if jcircle < icircle else (pinfo1[0], pinfo1[1], jcircle, jpt)
      for a1 in angles + [alast]:
        area = side_area(p0, p1, a1 - a0)
        dens_circle, dens_ipt = (icircle, ipt) if jcircle < icircle else (jcircle, jpt)
        g = conductance_density_parabola(ilayer, dens_circle, dens_ipt)
        g2 = org2gid(jlayer, jcircle, jpt)
        if gid_is_simulated(g2) and jcircle < ncircle[jlayer] - 1: gs.append(set_gap(gid, g2, abscond(area, g)))
        jpt += 1
        a0 = a1

  return gs

def test1(ilayer, icircle, ipt):
  print (gaps_for_gid(org2gid(ilayer, icircle, ipt)))

if __name__ == "__main__":
  test1(0,0,0)
