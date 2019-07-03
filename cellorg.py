from common import timeit, pr, pc, rank
timeit()
import param as p
from morphdef import circle_origins
from morphdef import distance, const_sep_layer_origins, addmul, normgrad
from p100from import p100from
from math import pi, cos, sin, tan


# structure that allows fast calculation of cell position and neighbors
# paraboloid contains layers contains circle origins
# The origin is a pair of (r, 0.0, z) tuples specifying the p000 and p100
# corners of the cell/region. These corners are shared by the cell/region
# (ilayer, icircle-1, 0) corners p010 and p110
# We really need to save space by using (r, z) instead of (r, 0.0, z)
# Since each of the four parabola edges is a piecewise linear function
# between the two corner points and corner points in the adjacent layer,
# it is useful to also supply those extra interior (r, 0.0, z) points.

# Instead of origin being (p000, p100). It is now
# ((r, 0.0, z), (r, 0.0, z), RegionFace)
class RegionFace:
  def __init__(self):
    #self.p0 = None # p000
    #self.p1 = None # p100
    self.p0b = None  # (jcircle, [p110] of layer-1 where p110 points are between p000 and p010
    self.p1b = None  # (jcircle, [p010] of layer+1 where p010 points are between p100 and p110
    # the b stands for breakpoint. This misses the point that there may
    # be no breakpoints and jcircle is the first circle in the adjacent layer
    # that is relevant to gap connectivity.

  def __str__(self):
    return "<RegionFace %s %s>" % (str(self.p0b), str(self.p1b))

# return number of points on circle
def circle0_n(pt):
  c = 2.*pi*pt[0]
  n = int(c/p.nominal_region_length)
  return n

# return n points around circle
def circle_discrete(n, pt):
  z = pt[2]
  r = pt[0]
  pts = []
  for i in range(n):
    a = (2.*pi/n)*i
    pts.append((r*cos(a), r*sin(a), z))
  return pts


def morphorg():
  global paraboloid, ncircle
  csep = p.layer_surface_circle_distance
  dsep = p.layer_thickness
  rstart = p.hole_radius
  zend = p.nominal_height
  #try first from file
  paraboloid = paraboloid_from_file()
  if paraboloid is None:
    # start with p000 corners
    paraboloid = [circle_origins(csep, rstart, zend)]
    for i in range(1, p.n_layer):
      paraboloid.append(const_sep_layer_origins(i, csep, paraboloid[0][0], zend))

    # replace p000 corners with (p000, p100) corner pairs.
    for i in range(p.n_layer):
      for j, p000 in enumerate(paraboloid[i]):
        paraboloid[i][j] = (p000, p100from(p000, dsep), RegionFace())

    # Fill in the RegionFace info of paraboloid
    ncircle = [len(paraboloid[i]) for i in range(p.n_layer)]
    for ilayer in range(p.n_layer):
      for icircle in range(ncircle[ilayer]):
        fill_interlayer_overlap(paraboloid, ncircle, ilayer, icircle)

    paraboloid_to_file(paraboloid)

  else:
    ncircle = [len(paraboloid[i]) for i in range(p.n_layer)]

  npts = [[circle0_n(o[0]) for o in paraboloid[i]] for i in range(p.n_layer)]
  return paraboloid, npts, ncircle

def paraboloid_filename():
  parm=(p.n_layer, p.abc, p.hole_radius, p.nominal_thickness, p.layer_surface_circle_distance, p.layer_thickness, p.cell_length, p.cell_diameter)
  return "paraboloid."+str(parm.__hash__()%0xffffffff)

def paraboloid_from_file():
  par = None
  if rank == 0:
    try:
      import pickle
      fname = paraboloid_filename()
      f = open(fname, "rb")
      par = pickle.load(f)
      print ("paraboloid read from " + fname)
    except:
      print (fname + " does not exist. Will be created")
  par = pc.py_broadcast(par, 0)
  pc.barrier()
  return par

def paraboloid_to_file(par):
  if rank == 0:
    try:
      import pickle
      fname = paraboloid_filename()
      f = open(fname, "wb")
      pickle.dump(par, f)
      print ("dumped paraboloid to " + fname)
    except:
      print ("could not dump " + fname)
      pass
  pc.barrier()

pt2circle_maxiter = 0

def pt2circle(ilayer, pt, use=0):
  global pt2circle_maxiter
  # given a point in layer, what is the icircle that contains the point
  # Note that an icircle domain is from 0 to sep in the parabola dimension
  assert(pt[1] == 0.0) # assume x, z on the layer surface.
  circles = paraboloid[ilayer]
  sep = p.layer_surface_circle_distance
  # kind of a discrete newton method
  # the dz between circles is an increasing function of z so start at end
  i = len(circles) - 1
  z = circles[i][use][2]
  if pt[2] >= z:
    #assert(pt[2] <= z + sep) # bug since sep is larger for greater layers
    return i
  maxiter = 0
  while i > 0 and pt[2] < z:
    dz = z - circles[i-1][use][2]
    j = int((pt[2] - z)/dz)
    i += j if j else -1
    z = circles[i][use][2]
    assert(pt[2] < circles[i+1][use][2])
    maxiter += 1
    assert (maxiter < 20)
  if maxiter > pt2circle_maxiter:
    pt2circle_maxiter = maxiter
  if False and pt[2] < z:
    print ("ilayer=", ilayer, " pt= ", pt)
    print ("i=", i)
    print("pt[2]= ", pt[2], " z=", z)
  assert(pt[2] >= z - 1e-8)
  return i

# Interlayer overlap between ilayer, icircle, anything)
# and the relevant jcircle in jlayer (must be ilayer + 1 or ilayer - 1)
# Fill in the RegionFace info for both parabola edges.
# with jcircle, [p, ..] where jcircle identifies
# the first circle that overlaps icircle and the p are the (r, 0, z)p000
# for jcircle+1, etc, that are interior to parabola edge.
# Note that the list is empty if there are no interior points.
def fill_interlayer_overlap(paraboloid, ncircle, ilayer, icircle):
  nlayer = p.n_layer
  if icircle >= ncircle[ilayer] - 1:
    return # the last circle has no RegionFace
  o0i = paraboloid[ilayer][icircle]
  o1i = paraboloid[ilayer][icircle + 1]
  rf = o0i[2]

  # (i=0,j=1 for jlayer = ilayer - 1) and (i=1,j=0 for jlayer = ilayer + 1)
  for jlayer in range(ilayer-1, ilayer+2, 2):
    if jlayer >= 0 and jlayer < nlayer:
      i,j = (0,1) if jlayer < ilayer else (1,0)
      nj = ncircle[jlayer]
      a = o0i[i] # p100i for ilayer
      b = o1i[i] # p110i
      # what is jcircle
      jcircle = pt2circle(jlayer, a, use=j)
      c = paraboloid[jlayer][jcircle][j] # p000j
      if c[2] > a[2]:
        jcircle -= 1
      result = (jcircle, [])
      while True:
        jcircle += 1
        if (jcircle >= nj):
          break
        c = paraboloid[jlayer][jcircle][j]
        if c[2] >= b[2]:
          break
        result[1].append(c)
      if i == 1:
        rf.p1b = result
      else:
        rf.p0b = result

paraboloid, npts, ncircle = morphorg()
nlayer = len(paraboloid)
nlayerpts = [sum(npts[i]) for i in range(p.n_layer)]
ngid = sum(nlayerpts)
circle_offset = [[0] for i in range(p.n_layer)]
layer_offset = [0]
for ilayer in range(p.n_layer):
  for npt in npts[ilayer]:
    circle_offset[ilayer].append(circle_offset[ilayer][-1] + npt)
  layer_offset.append(layer_offset[-1] + circle_offset[ilayer][-1])

timeit("abstract model definition")
pr("ngid = %d"%ngid)

def org2gid(ilayer, icircle, ipt):
  gid = layer_offset[ilayer] + circle_offset[ilayer][icircle] + ipt%npts[ilayer][icircle]
  return gid
  
def gid2org(gid):
  ilayer = gid2layer(gid)
  r = gid - layer_offset[ilayer]
  icircle = gid2org_help_circle(ilayer, r)
  ipt = r - circle_offset[ilayer][icircle]
  return ilayer, icircle, ipt

itermax=0

def gid2layer(gid): # return last i where layer_offset[i] <= gid
  # nlayerpts is concave so use discrete newton method iwth decreasing indices
  i = nlayer - 1
  iter = 0
  while True:
    if layer_offset[i] <= gid:
      break
    i -= int((layer_offset[i] - gid)/nlayerpts[i]) + 1
    iter += 1
  global itermax
  if iter > itermax: itermax = iter
  return i

def gid2org_help_circle(ilayer, r): # return last i where circle_offset[i] <= r
  # note npts is concave
  # something like a discrete newton method can work with decreasing indices
  # No more than 7 iterations.
  i = len(npts[ilayer]) - 1
  iter = 0
  while True:
    if circle_offset[ilayer][i] <= r:
      break
    j = int((circle_offset[ilayer][i] - r)/npts[ilayer][i])
    i -= j if j else 1
    iter += 1
  global itermax
  if iter > itermax: itermax = iter
  return i

def gid_is_simulated(gid):
  return is_simulated(xyz(*gid2org(gid)))

def is_simulated(xyz):
    return distance(xyz, p.simulation_center) < p.simulation_region

def xyz(ilayer, icircle, ipt): # note that ipt refers to the proximal point on the section
  n = npts[ilayer][icircle]
  ipt = ipt%n # wrap around
  pt = paraboloid[ilayer][icircle][0]
  a = (2*pi/n)*ipt
  r = pt[0]
  return r*cos(a), r*sin(a), pt[2]

def xyz_center(ilayer, icircle, ipt):
 x1, y1, z1 = xyz(ilayer, icircle, ipt)
 x2, y2, z2 = xyz(ilayer, icircle, ipt+1)
 return (x1 + x2)/2., (y1 + y2)/2., (z1 + z2)/2.

# Angle overlap in the circumferential direction
# between pt=(ilayer, icircle, ipt) and  the relevant points
# in layer, circle.
# Return jpt, [angle1, angle2, ...] where jpt is the
# identifies the first point on circle that overlaps ipt, and the
# angles are the distal angles that ovelap ipt. Note that
# the list is empty if the the jpt distal angle is larger than ipt distal
# angle. The most common overlap is jpt, [angle1], which means that
# the distal angle of jpt is interior to ipt and the distal angle of
# jpt+1%npt[layer][circle] is exterior to ipt
# The overlapping circle sections
def angle_overlap(pt, jlayer, jcircle):
  ilayer = pt[0]
  icircle = pt[1]
  ipt = pt[2]
  ai0 = ipt2angle(ipt, ilayer, icircle) # proximal angle of pt
  ai1 = ipt2angle(ipt+1, ilayer, icircle) # distal angle of pt
  if ai1 == 0:
    ai1 = 2.*pi

  jpt = angle2ipt(ai0, jlayer, jcircle) # jpt contains ai0
  result = (jpt, [])
  while True:
    jpt += 1
    aj0 = ipt2angle(jpt, jlayer, jcircle)
    if aj0 >= ai1 or aj0 == 0.:
      break
    result[1].append(aj0)
    
  return result

def ipt2angle(ipt, ilayer, icircle):
  n = npts[ilayer][icircle]
  ipt = ipt%n # occasionally convenient for ipt = -1
  return (2*pi/n)*ipt

def angle2ipt(angle, ilayer, icircle): #ipt that contains angle (note that ipt refers to the proximal point on the section)
  n = npts[ilayer][icircle]
  return int(angle*(n/(2*pi)))%n

def test1():
  from time import time
  #iterate over all points
  t = time()
  paraboloid, npts = morphorg()
  print ("morphorg time %g" % (time() - t,))
  t = time()
  cnt=0
  for ilayer, layer in enumerate(paraboloid):
    for icircle, circle in enumerate(layer):
      n = npts[ilayer][icircle]
      for ipt, pt in enumerate(circle_discrete(n, circle)):
        # do something with ilayer, icircle, ipt
        gid = org2gid(ilayer, icircle, ipt)
        org = gid2org(gid)
        assert(org == (ilayer, icircle, ipt))
        cnt += 1
        pass
  print ("visit all %d points time %g" % (cnt, (time() - t)))

  t = time()
  cnt = 0
  for ilayer in range(nlayer):
    for icircle in range(ncircle - 1):
      for ipt in range(npts[ilayer][icircle]):
        o = overlap((ilayer, icircle, ipt), icircle + 1)
        cnt += len(o)
  print ("all %d circle to circle connections time %g" % (cnt, (time() - t)))


def test2(layer, circle, ipt):
  print(layer, circle, ipt)
  print("npts ", npts[ilayer][circle], npts[ilayer][circle+1])
  a = overlap((layer, circle, ipt), layer, circle + 1)
  print(a)
  for x in a:
    b = overlap((layer, circle + 1, x[0]), layer, circle)
    print(b)

  print("length d_angle of circle ", circle, distance(xyz(layer, circle, 0), xyz(layer, circle, 1)), ipt2angle(1, layer, circle))
  print("length d_angle of circle ", circle+1, distance(xyz(layer, circle+1, 0), xyz(layer, circle+1, 1)), ipt2angle(1, layer, circle+1))
  
  for i, y in enumerate([a, b]):
    c = circle+(1-i)
    for x in y:
      print([layer, c, x[0]], xyz(layer, c, x[0]))
  print ()

def test3(l1, c1, c2):
  # sum of all fractions around each circle should total npts for the circle
  n1 = npts[l1][c1]
  n2 = npts[l1][c2]
  s1 = 0.
  s2 = 0.
  for i in range(n1):
    a = overlap((l1, c1, i), l1, c2)
    for o in a:
      s1 += o[1]
      s2 += o[2]
  print("ll=%d c1=%d c2=%d n1=%d n2=%d s1=%g s2=%g" % (l1, c1, c2, n1, n2, s1, s2))

def test4():
  from neuron import h, gui
  # distance between circles as function of surface distance for
  # inner and outemost layers. Also length of cells.
  g1 = h.Graph()
  #for j, layer in enumerate([0, nlayer-1]):
  for j, layer in enumerate(range(nlayer)):
    p0 = xyz(layer, 0, 0)
    d_surf = [distance(xyz(layer, i, 0), p0) for i in range(1, ncircle[layer])]
    d_circle = [distance(xyz(layer, i, 0), xyz(layer, i-1, 0)) for i in range(1, ncircle[layer])]
    print(d_circle, d_surf)
    h.Vector(d_circle).line(g1, h.Vector(d_surf), j+1, 2)
  return g1

def test5(): #pt2circle
  from neuron import h, gui
  o0 = paraboloid[0]
  n = 4
  o1 = paraboloid[n]
  d = distance(o1[0][0], o0[0][0])
  g = h.Graph()
  for p in o1:
    g.mark(p[0][0], p[0][2], "|", 10, 1, 1)
  for o in o0:
    p = o[0]
    p1 = addmul(p, d, normgrad(p))
    i = pt2circle(n, p1)
    g.mark(p1[0], p1[2], "|", 10, 2, 1)
    g.mark(o1[i][0][0], o1[i][0][2], "|", 10, 3, 1)
    d1 = distance(p1, o1[i][0])
    if i < len(o1)-1:
      zi = o1[i][0][2]
      z = p1[2]
      zip = o1[i+1][0][2]
      print("zi=%g z=%g zip=%g"%(zi, z, zip))
      assert (zi <= z+1e-9 and z <= zip+1e-9)
  return g

if __name__ == "__main__":
 
  #test1()
  if False:
    c = 1
    n = npts[0][c]
    test2(1, c, 0)
    test2(1, c, n-1)
    print ()
    test2(1, c, 1)
    test2(1, c, n-2)

  #test3(1, 0, 1)
  #test3(4, 101, 100)

  #a = test4()
  #a = test5()
