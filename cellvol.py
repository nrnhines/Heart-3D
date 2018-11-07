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
from cellorg import xyz, ipt2angle
from morphdef import normgrad, addmul
from neuron import h
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

maxiter = 0

def paranormal0(pt):
  #print ("paranormal0 ", pt)
  # return the normal direction from the inner paraboloid that contains pt
  assert(pt[1] == 0.0)
  # two dimensional problem due to assertion.
  x = pt[0]
  z = pt[2]
  lam = 0.0
  # three equations in three unknowns (x, z) on parboloid (x', z') is pt
  # unknowns x, z, lam
  # a*x^2 - c*z = 0
  # x + lam*2*a*x - x' = 0
  # z - lam*c  - z' = 0
  # can reduce to single cubic equation in x but might as well solve
  # by newton method. Jacobian is
  # [2*a*x      , -c, 0    ]
  # [1 + 2*a*lam, 0 , 2*a*x]
  # [0          , 1 , -c   ]
  b = h.Vector(3)
  m = h.Matrix(3,3)
  dx = h.Vector(3)
  iter = 0
  for i in range(100):
    mset(m, x, z, lam)
    bset(b, x, z, lam, pt[0], pt[2])
    m.solv(b, dx)
    #dx.printf()
    x += dx.x[0]
    z += dx.x[1]
    lam += dx.x[2]
    iter += 1
    if dx.sumsq() < 1e-14:
      break
  assert (iter < 100)
  global maxiter
  maxiter = iter if iter > maxiter else maxiter
  return normgrad((x, 0, z))

def mset(m, x, z, lam):
  a, _, c = param.abc
  m.zero()
  m.x[0][0] = 2.0*a*x
  m.x[0][1] = -c
  m.x[1][0] = 1.0 + 2.0*a*lam
  m.x[1][2] = 2.0*a*x
  m.x[2][1] = 1.0
  m.x[2][2] = -c

def bset(b, x, z, lam, x0, z0):
  a, _, c = param.abc
  b.x[0] = -(a*x*x - c*z)
  b.x[1] = -(x + lam*2.0*a*x - x0)
  b.x[2] = -(z - lam*c - z0)


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
