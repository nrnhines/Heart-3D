import param as p
from morphdef import circle_origins, outer_layer_origins
from morphdef import distance, const_sep_layer_origins, addmul, normgrad
from math import pi, cos, sin, tan

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
    #paraboloid.append(outer_layer_origins(i, paraboloid[0]))
    paraboloid.append(const_sep_layer_origins(i, sep, paraboloid[0][0], zend))
  npts = [[circle0_n(pt) for pt in paraboloid[i]] for i in range(p.n_layer)]
  return paraboloid, npts

pt2circle_maxiter = 0

def pt2circle(ilayer, pt):
  global pt2circle_maxiter
  # given a point in a layer, what is the icircle that contains the point
  # Note that an icircle domain is from -sep/2 to +sep/2.
  assert(pt[1] == 0.0) # assume x, z on the layer surface.
  circles = paraboloid[ilayer]
  sep2 = p.internal_surface_circle_distance * .5
  # kind of a discrete newton method
  # the dz between circles is an increasing function of z so start at end
  i = len(circles) - 1
  z = circles[i][2]
  if pt[2] >= z:
    assert(pt[2] <= z + sep2)
    return i
  maxiter = 0
  while i > 0 and pt[2] < z:
    dz = z - circles[i-1][2]
    j = int((pt[2] - z)/dz)
    i += j if j else -1
    z = circles[i][2]
    assert(pt[2] < circles[i+1][2])
    maxiter += 1
    assert (maxiter < 20)
  if maxiter > pt2circle_maxiter:
    pt2circle_maxiter = maxiter
  if i == 0 and pt[2] < z :
    assert(pt[2] >= z - sep2)
    return i
  d1 = distance(pt, circles[i])
  d2 = distance(pt, circles[i+1])
  if (d1 > d2):
    i += 1
  return i

paraboloid, npts = morphorg()
nlayer = len(paraboloid)
nlayerpts = [sum(npts[i]) for i in range(p.n_layer)]
ncircle = [len(paraboloid[i]) for i in range(p.n_layer)]
ngid = sum(nlayerpts)
circle_offset = [[0] for i in range(p.n_layer)]
layer_offset = [0]
for ilayer in range(p.n_layer):
  for npt in npts[ilayer]:
    circle_offset[ilayer].append(circle_offset[ilayer][-1] + npt)
  layer_offset.append(layer_offset[-1] + circle_offset[ilayer][-1])

def org2gid(ilayer, icircle, ipt):
  gid = layer_offset[ilayer] + circle_offset[ilayer][icircle] + ipt
  return gid
  
def gid2org(gid):
  ilayer = gid2layer(gid)
  r = gid - layer_offset[ilayer]
  icircle = gid2org_help_circle(ilayer, r)
  ipt = r - circle_offset[ilayer][icircle]
  return ilayer, icircle, ipt

itermax=0

def gid2layer(gid): # return last i where layer_offset[i] <= gid
  # ugh linear search
  for i in range(1, p.n_layer + 1):
    if gid < layer_offset[i]:
      return i - 1
  assert(False)

def gid2org_help_circle(ilayer, r): # return last i where circle_offset[i] <= r
  # note npts is concave
  # something like a discrete newton method can work with decreasing indices
  # No more than 7 iterations.
  i = len(npts[ilayer]) - 1
  iter = 0
  while True:
    if circle_offset[ilayer][i] <= r:
      break
    i -= int((circle_offset[ilayer][i] - r)/npts[ilayer][i]) + 1
    iter += 1
  global itermax
  if iter > itermax: itermax = iter
  return i

def xyz(ilayer, icircle, ipt): # note that ipt refers to the proximal point on the section
  n = npts[ilayer][icircle]
  ipt = ipt%n # wrap around
  pt = paraboloid[ilayer][icircle]
  a = 2*pi*ipt/n
  r = pt[0]
  return r*sin(a), r*cos(a), pt[2]

def xyz_center(ilayer, icircle, ipt):
 x1, y1, z1 = xyz(ilayer, icircle, ipt)
 x2, y2, z2 = xyz(ilayer, icircle, ipt+1)
 return (x1 + x2)/2., (y1 + y2)/2., (z1 + z2)/2.

#fractional edge overlap between (ilayer, icircle, ipt) pt1 and
# the relevant points in layer, circle
# the overlapping circle sections
# return (overlap1, overlap2) where overlaps are the
# fractional length of sections pt1 and pt2. Note, because circles have
# different circumference a section in icircle can subtend a one or more
# sections in circle. Return a list of 3tuples where each tuple is
# i on circle, fractional overlap of ipt on i, fractional overlap of i on ipt)
def overlap(pt, layer, circle):
  deb=False
  if deb: print("\nenter overlap ", pt, circle)
  result = []
  ilayer = pt[0]
  icircle = pt[1]
  ipt = pt[2]
  aesc = ipt2angle(1, layer, circle) # angle of each section on circle
  aesp = ipt2angle(1, ilayer, icircle) # angle of each section on pcircle

  amin = ipt2angle(ipt, ilayer, icircle) # angle of pt on it's own circle is the angle on circle
  imin = angle2ipt(amin, layer, circle) # imin contains amin
  aimin = ipt2angle(imin, layer, circle)# angle of imin

  amax = ipt2angle(ipt + 1, ilayer, icircle)
  imax = angle2ipt(amax, layer, circle) # imax contain amax
  aimax = ipt2angle(imax, layer, circle)# angle of imax
  if imax < imin:
   imax = npts[layer][circle]
   amax = 2*pi
   aimax = 2*pi

  if deb: print("ipt=%d amin=%g amax=%g" %(ipt, amin, amax))
  if deb: print("imin=%d aimin=%g imax=%d aimax=%g" %(imin, aimin, imax, aimax))

  if aimax < amax: #ipt subtends some of the current imax section
    imax += 1 # could be npt and therefore refer to 0
    aimax = ipt2angle(imax, layer, circle) # could be 2pi less than true angle
    if imax == npts[layer][circle]:
      aimax = 2*pi
  if deb: print("imin=%d aimin=%g imax=%d aimax=%g" %(imin, aimin, imax, aimax))

  for i in range(imin, imax): # the sections on circle that are involved
    if deb: print("begin i=%d"%i)
    # ipt section subtend fraction on circle sections
    aicprox = ipt2angle(i, layer, circle)
    fcprox = 0.0 # beginning overlap distance of ipt on i section of circle
    if aicprox < amin: # fcprox > 0
      fcprox = fracangle(amin - aicprox, aesc)
    aicdist = ipt2angle(i + 1, layer, circle)
    if aicdist < aicprox:
      aicdist += 2*pi
    fcdist = 1.0 # ending overlap distance of ipt on i section of circle
    if aicdist > amax:
     fcdist = fracangle(amax - aicprox, aesc)
    if deb: print("i=%d aicprox=%g aicdist=%g" %(i, aicprox, aicdist))

    # circle section subtend fraction on ipt section
    # Note, the relevant angles are still the same, ie.
    # amin, amax, aicprox, aicdist, but use aescp instead of aesc
    # and the angle relationship is opposite
    fpcprox = 0.0
    if aicprox > amin:
      fpcprox = fracangle(aicprox - amin, aesp)
    fpcdist = 1.0
    if aicdist < amax:
      fpcdist = fracangle(aicdist - amin, aesp)

    result.append((i, fcdist - fcprox, fpcdist - fpcprox))
    if deb: print("end i=%d"%i)

  if deb: print("leave overlap", pt, circle, "\n")
  return result

def fracangle(a, a0):
  # isoceles triangle with vertex angle a0. Return fractional distance
  # along base from angle a. Note a=0 return 0, a=a0 return 1
  theta = 0.5*a0
  phi = a - theta
  return 0.5*(1.0 + tan(phi)/tan(theta))

def ipt2angle(ipt, ilayer, icircle):
  n = npts[ilayer][icircle]
  ipt = ipt%n # occasionally convenient for ipt = -1
  return 2*pi*ipt/n

def angle2ipt(angle, ilayer, icircle): #ipt that contains angle (note that ipt refers to the proximal point on the section)
  n = npts[ilayer][icircle]
  return int(angle*n/(2*pi))%n

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
  g2 = h.Graph()
  #for j, layer in enumerate([0, nlayer-1]):
  for j, layer in enumerate(range(nlayer)):
    p0 = xyz(layer, 0, 0)
    d_surf = [distance(xyz(layer, i, 0), p0) for i in range(1, ncircle[layer])]
    d_circle = [distance(xyz(layer, i, 0), xyz(layer, i-1, 0)) for i in range(1, ncircle[layer])]
    print(d_circle, d_surf)
    h.Vector(d_circle).line(g1, h.Vector(d_surf), j+1, 2)
  return g1, g2

def test5(): #pt2circle
  from neuron import h, gui
  o0 = paraboloid[0]
  n = 4
  o1 = paraboloid[n]
  d = distance(o1[0], o0[0])
  g = h.Graph()
  for p in o1:
    g.mark(p[0], p[2], "|", 10, 1, 1)
  for p in o0:
    p1 = addmul(p, d, normgrad(p))
    i = pt2circle(n, p1)
    g.mark(p1[0], p1[2], "|", 10, 2, 1)
    g.mark(o1[i][0], o1[i][2], "|", 10, 3, 1)
    d1 = distance(p1, o1[i])
    if i > 0 and p1[2] < o1[i][2]:
      d2 = distance(p1, o1[i-1])
      print(i, " < ", distance(p1, o1[i]), distance(p1, o1[i-1]))
      assert (d1 <= d2)
    if i < len(o1)-1 and p1[2] > o1[i][2]:
      d2 = distance(p1, o1[i+1])
      print(i, " > ", distance(p1, o1[i]), distance(p1, o1[i+1]))
      assert (d1 <= d2)

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
  a = test5()
