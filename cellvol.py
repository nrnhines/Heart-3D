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
from cellorg import ipt2angle, paraboloid
from morphdef import normgrad, addmul
import param


class Corners:
  def __repr__(s):
    a = ""
    for i in [0, 1]:
      for j in [0, 1]:
        for k in [0, 1]:
          p = eval("s.p%d%d%d"%(i,j,k))
          a += "p%d%d%d = (%g, %g, %g)\n" % (i, j, k, p[0], p[1], p[2])
    return a          

def xyz(pt, a):
  return cos(a)*pt[0], sin(a)*pt[0], pt[2]

def cellcorners(ilayer, icircle, ipt):
  c = Corners()
  o0 = paraboloid[ilayer][icircle]
  o1 = paraboloid[ilayer][icircle+1]
  a0 = ipt2angle(ipt, ilayer, icircle)
  a1 = ipt2angle(ipt+1, ilayer, icircle)

  c.p000 = xyz(o0[0], a0)
  c.p001 = xyz(o0[0], a1)
  c.p100 = xyz(o0[1], a0)
  c.p101 = xyz(o0[1], a1)

  c.p010 = xyz(o1[0], a0)
  c.p011 = xyz(o1[0], a1)
  c.p110 = xyz(o1[1], a0)
  c.p111 = xyz(o1[1], a1)

  return c

def test1(ilayer, icircle, ipt):
  print ("cellcorners(%d, %d, %d)" % (ilayer, icircle, ipt))
  c = cellcorners(ilayer, icircle, ipt)
  print (c)

def test2(): # x, z slice through ipt = 0
  from neuron import h, gui
  g = h.Graph(0)
  g.view(2)
  for ilayer in range(nlayer-1):
    for icircle in range(ncircle[ilayer]-1):
      c = cellcorners(ilayer, icircle, 0)
      g.beginline(ilayer+1, 1)
      g.line(c.p000[0], c.p000[2])
      g.line(c.p010[0], c.p010[2])
      g.line(c.p110[0], c.p110[2])
      g.line(c.p100[0], c.p100[2])
      g.line(c.p000[0], c.p000[2])
      g.flush()
  return g
if __name__ == "__main__":
  from p100from import maxiter
  test1(0, 0, 0)
  test1(0, 1, 0)
  test1(1, 0, 0)
  from cellorg import nlayer, ncircle
  ilayer = nlayer-2
  for ilayer in range(nlayer-1):
    for icircle in range(ncircle[ilayer]-1):
      cellcorners(ilayer, icircle, 10)
  print ("maxiter = ", maxiter)
  test2()
