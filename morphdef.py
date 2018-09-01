import param as p
from math import sqrt, pi, fabs

#distance between two points
def distance(p1, p2):
  return sqrt(sum([(x-y)**2 for x, y in zip(p1, p2)]))

# circle origins [(x, y=0, z)] along p.abc paraboloid that are sep apart
# as measured on chord surface of paraboloid, starting with radious rstart
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

def outer_layer_origins(i_layer, origins):
  o = []
  d = i_layer * p.layer_width
  for pt in origins:
    g = (2.*p.abc[0]*pt[0], 0.0, -p.abc[2])
    norm = distance(g, (0., 0., 0.))
    p2 = tuple([a + b*d/norm for a, b in zip(pt, g)])
    o.append(p2)
  return o


def test1():
  origins = circle_origins(p.internal_surface_circle_distance, p.hole_radius, p.nominal_height)
  nsec = 0
  o4 = outer_layer_origins(4, origins)
  for p0, p4 in zip(origins, o4):
     c = 2.*pi*p0[0]
     c4 = 2.*pi*p4[0]
     n = int(c/p.nominal_cell_length)
     print ("%d sections of length %g (layer 4 length %g)" %(n, c/n, c4/n))
     nsec += n
  print ("%d circles in one layer with total %d sections"%(len(origins), nsec))
  return origins

if __name__ == "__main__":
  origins = test1()
