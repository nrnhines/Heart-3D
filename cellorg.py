import param as p
from morphdef import circle_origins, outer_layer_origins
from math import pi, cos, sin

# structure that allows fast calculation of cell position and neighbors
# paraboloid contains layers contains circle origins

# return number of points on circle
def circle0_n(pt):
  c = 2.*pi*pt[0]
  n = int(c/p.nominal_cell_length)
  return n

# return n points around circle
def circle_discrete(n, pt):
  z = pt[2]
  r = pt[0]
  pts = []
  for i in range(n):
    a = 2.*pi*i/n
    pts.append((r*sin(a), r*cos(a), z))
  return pts

def morphorg():
  sep = p.internal_surface_circle_distance
  rstart = p.hole_radius
  zend = p.nominal_height
  paraboloid = [circle_origins(sep, rstart, zend)]
  for i in range(1, p.n_layer):
    paraboloid.append(outer_layer_origins(i, paraboloid[0]))
  npts = [circle0_n(pt) for pt in paraboloid[0]]
  return paraboloid, npts

def test1():
  #iterate over all points
  paraboloid, npts = morphorg()
  for ilayer, layer in enumerate(paraboloid):
    for icircle, circle in enumerate(layer):
      n = npts[icircle]
      for ipt, pt in enumerate(circle_discrete(n, circle)):
        # do something with ilayer, icircle, ipt
        pass

if __name__ == "__main__":
  test1()
