# Cell position p000 = xyz(ilayer, icircle, ipt) refers to the cell corner
# of least r (ilayer), least z (icircle), least angle (ipt).
# The cell volume is defined by 4 arcs having the same subtend angle
# The primary arc is
# p000 , p001 = xyz(ilayer, icircle, ipt+1)
# from which one calculates the subtend angle (a1, a2).
# The other three arcs are the same subtend angle for the circles defined by
# points
# p100 = p000 + dlayer*gradnorm(p000)
# p010 = (x cos(a1), x sin(a1), z) where x,z from xyz(ilayer, icircle+1, 0)
# p110 = p010 + dlayer*gradnorm(p010)
from math import cos, sin
from cellorg import xyz, ipt2angle, paraboloid
from morphdef import normgrad, addmul
from neuron import h
import param

from p100from import paranormal0, maxiter

class Corners:
  def __repr__(s):
    a = ""
    for i in [0, 1]:
      for j in [0, 1]:
        for k in [0, 1]:
          p = eval("s.p%d%d%d"%(i,j,k))
          a += "p%d%d%d = (%g, %g, %g)\n" % (i, j, k, p[0], p[1], p[2])
    return a          

def cellcorners(ilayer, icircle, ipt):
  c = Corners()
  c.p000 = xyz(ilayer, icircle, ipt)
  c.p001 = xyz(ilayer, icircle, ipt + 1)

  a1 = ipt2angle(ipt, ilayer, icircle)
  a2 = ipt2angle(ipt+1, ilayer, icircle)
  pt = xyz(ilayer, icircle, 0)
  pn = paranormal0(pt)
  dlayer = param.layer_thickness
  pt = addmul(pt, dlayer, pn)
  c.p100 = (cos(a1)*pt[0], sin(a1)*pt[0], pt[2])
  c.p101 = (cos(a2)*pt[0], sin(a2)*pt[0], pt[2])
  
  pt = xyz(ilayer, icircle+1, 0)
  c.p010 = (cos(a1)*pt[0], sin(a1)*pt[0], pt[2])
  c.p011 = (cos(a2)*pt[0], sin(a2)*pt[0], pt[2])

  pn = paranormal0(pt)
  pt = addmul(pt, dlayer, pn)
  c.p110 = (cos(a1)*pt[0], sin(a1)*pt[0], pt[2])
  c.p111 = (cos(a2)*pt[0], sin(a2)*pt[0], pt[2])

  return c

def test1(ilayer, icircle, ipt):
  print ("cellcorners(%d, %d, %d)" % (ilayer, icircle, ipt))
  c = cellcorners(ilayer, icircle, ipt)
  print (c)

if __name__ == "__main__":
  test1(0, 0, 0)
  test1(1, 0, 0)
  from cellorg import nlayer, ncircle
  ilayer = nlayer-2
  for ilayer in range(nlayer-1):
    for icircle in range(ncircle[ilayer]-1):
      cellcorners(ilayer, icircle, 10)
  print ("maxiter = ", maxiter)
