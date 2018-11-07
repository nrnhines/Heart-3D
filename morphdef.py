import param as p
from math import sqrt, pi, fabs

#distance between two points
def distance(p1, p2):
  return sqrt(sum([(x-y)**2 for x, y in zip(p1, p2)]))

# circle origins [(x, y=0, z)] along p.abc paraboloid that are sep apart
# as measured on chord surface of paraboloid, starting with radius rstart
# and ending just after zend. These circles lie on layer 0.
def circle_origins(sep, rstart, zend):
  origins = []
  pt = (rstart, 0.0, p.abc[0]*rstart**2/p.abc[2])
  while True:
    origins.append(pt)
    if pt[2] > zend:
      break
    pt = next_circle_origin(sep, pt)
  return origins

def addmul(p1, d, p2):
  return tuple([a + d*b for a, b in zip(p1, p2)])

def normgrad(pt): # pt on layer 0 with y=0
  g = (2.*p.abc[0]*pt[0], 0.0, -p.abc[2])
  norm = distance(g, (0., 0., 0.))
  return (g[0]/norm, 0., g[2]/norm)

def add_help_(pt):
  # unit vector in direction of increasing z on parabola from point p
  # g is normal to parabola
  g = (2.*p.abc[0]*pt[0], 0.0, -p.abc[2])
  norm = distance(g, (0., 0., 0.))
  g = (-g[2]/norm, 0.0, g[0]/norm) # up the parabola
  return g

def add_(pt, d):
  # add distance d to pt in direction of increasing z on parabola
  # direction up the parabola (increasing x and z)
  v = add_help_(pt)
  return tuple([a + b*d for a, b in zip(pt, v)])

def const_sep_layer_origins(i_layer, sep, p0, zend):
  d = p.layer_thickness*i_layer
  o = []
  pt = addmul(p0, d, normgrad(p0))
  pt0 = p0
  while True:
    o.append(pt)
    if pt[2] > zend:
      break
    pt, pt0 = next_const_sep_layer_origin(pt, pt0, d, sep)
  return o

def next_circle_origin(sep, p1):
  # next point on the p.abc (y=0) parabola such the chord distance is sep
  # want p2[2] > p1[2], distance(p1, p2) = sep, and
  # p.abc[0]*p2[0]**2 - p.abc[2]*p2[2] = 0.0
  # First approx along gradient normal
  p2 = add_(p1, sep)
  iter=1
  while True:
    p2 = (p2[0], 0.0, p.abc[0]*p2[0]**2/p.abc[2]) # onto the parabola again
    d = distance(p1, p2)
    if fabs(d - sep) < 1e-9:
      return p2 # on the parabola and distance is sep
    p2 = add_(p2, sep - d)
    iter += 1
    assert(iter < 20)

def next_const_sep_layer_origin(pt, pt0, d, sep):
  # pt current point on the layer
  # pt0 point on inner layer paraboloid corresponding to pt (distance d)
  # note that pt is normal to the inner paraboloid at pt0
  # d separation in layer dimension (normal to layer0)
  # sep desired separation between pt and next pt along the outer layer surface

  # first approx along the gradient normal on the inner surface, overestimates.
  p0 = add_(pt0, sep)
  iter = 1
  while True:
    p0 = (p0[0], 0.0, p.abc[0]*p0[0]**2/p.abc[2]) # onto the parabola again
    d0 = distance(p0, pt0)
    pt1 = addmul(p0, d, normgrad(p0))
    d1 = distance(pt, pt1)
    x = (d1 - sep)/d1
    if fabs(d1 - sep) < 1e-9:
      return pt1, p0
    p0 = add_(p0, -x*d0)
    iter += 1
    assert(iter < 20)
    
def outer_layer_origins(i_layer, origins):
  o = []
  d = i_layer * p.layer_thickness
  for pt in origins:
    g = (2.*p.abc[0]*pt[0], 0.0, -p.abc[2])
    norm = distance(g, (0., 0., 0.))
    p2 = tuple([a + b*d/norm for a, b in zip(pt, g)])
    o.append(p2)
  return o


def test1():
  # Each layer has same number of circles
  origins = circle_origins(p.layer_surface_circle_distance, p.hole_radius, p.nominal_height)
  nsec = 0
  o4 = outer_layer_origins(4, origins)
  for p0, p4 in zip(origins, o4):
     c = 2.*pi*p0[0]
     c4 = 2.*pi*p4[0]
     n = int(c/p.nominal_region_length)
     print ("%d sections of length %g (layer 4 length %g)" %(n, c/n, c4/n))
     nsec += n
  print ("%d circles in one layer with total %d sections"%(len(origins), nsec))
  return origins

def test2():
  # Circles in each layer have constant separation
  origins = circle_origins(p.layer_surface_circle_distance, p.hole_radius, p.nominal_height)
  o4 = const_sep_layer_origins(4, p.layer_surface_circle_distance, origins[0], p.nominal_height)
  for o in [origins, o4]:
    nsec = 0
    for p0 in o:
      c = 2.*pi*p0[0]
      n = int(c/p.nominal_region_length)
      nsec += n
    print ("%d circles in layer with total %d sections"%(len(o), nsec))
  
  return origins, o4


if __name__ == "__main__":
  origins = test1()
  a = test2()
