#see cellvol.py for the meaning of p100from

from neuron import h
from math import cos, sin
from morphdef import normgrad, addmul
import param

def p100from(p000, dlayer):
  pn = paranormal0(p000)
  p100 = addmul(p000, dlayer, pn)
  return p100

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
