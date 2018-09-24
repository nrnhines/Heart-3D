import param as p
from morphdef import circle_origins, outer_layer_origins, distance
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
    paraboloid.append(outer_layer_origins(i, paraboloid[0]))
  npts = [circle0_n(pt) for pt in paraboloid[0]]
  return paraboloid, npts

paraboloid, npts = morphorg()
nlayer = len(paraboloid)
nlayerpts = sum(npts)
ncircle = len(paraboloid[0])
ngid = nlayerpts*nlayer
circle_offset = [0]
for npt in npts:
  circle_offset.append(circle_offset[-1] + npt)

def org2gid(ilayer, icircle, ipt):
  gid = ilayer*nlayerpts + circle_offset[icircle] + ipt
  return gid
  
def gid2org(gid):
  ilayer = int(gid/nlayerpts)
  r = gid - ilayer*nlayerpts
  icircle = gid2org_help_(r)
  ipt = r - circle_offset[icircle]
  return ilayer, icircle, ipt

itermax=0

def gid2org_help_(r): # return last i where circle_offset[i] <= r
  # note npts is concave
  # something like a discrete newton method can work with decreasing indices
  # No more than 7 iterations.
  i = len(npts) - 1
  iter = 0
  while True:
    if circle_offset[i] <= r:
      break
    i -= int((circle_offset[i] - r)/npts[i]) + 1
    iter += 1
  global itermax
  if iter > itermax: itermax = iter
  return i

def xyz(ilayer, icircle, ipt): # note that ipt refers to the proximal point on the section
  n = npts[icircle]
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
# the overlapping circle sections
# assume same layer, return (overlap1, overlap2) where overlaps are the
# fractional length of sections pt1 and pt2. Note, because circles have
# different circumference a section in pcircle can subtend a one or more
# sections in circle. Return a list of 3tuples where each tuple is
# i on circle, fractional overlap of ipt on i, fractional overlap of i on ipt)
def overlap(pt, circle):
  deb=False
  if deb: print("\nenter overlap ", pt, circle)
  result = []
  pcircle = pt[1]
  ipt = pt[2]
  aesc = ipt2angle(1, circle) # angle of each section on circle
  aesp = ipt2angle(1, pcircle) # angle of each section on pcircle

  amin = ipt2angle(ipt, pcircle) # angle of pt on it's own circle is the angle on circle
  imin = angle2ipt(amin, circle) # imin contains amin
  aimin = ipt2angle(imin, circle)# angle of imin

  amax = ipt2angle(ipt + 1, pcircle)
  imax = angle2ipt(amax, circle) # imax contain amax
  aimax = ipt2angle(imax, circle)# angle of imax
  if imax < imin:
   imax = npts[circle]
   amax = 2*pi
   aimax = 2*pi

  if deb: print("ipt=%d amin=%g amax=%g" %(ipt, amin, amax))
  if deb: print("imin=%d aimin=%g imax=%d aimax=%g" %(imin, aimin, imax, aimax))

  if aimax < amax: #ipt subtends some of the current imax section
    imax += 1 # could be npt and therefore refer to 0
    aimax = ipt2angle(imax, circle) # could be 2pi less than true angle
    if imax == npts[circle]:
      aimax = 2*pi
  if deb: print("imin=%d aimin=%g imax=%d aimax=%g" %(imin, aimin, imax, aimax))

  for i in range(imin, imax): # the sections on circle that are involved
    if deb: print("begin i=%d"%i)
    # ipt section subtend fraction on circle sections
    aicprox = ipt2angle(i, circle)
    fcprox = 0.0 # beginning overlap distance of ipt on i section of circle
    if aicprox < amin: # fcprox > 0
      fcprox = fracangle(amin - aicprox, aesc)
    aicdist = ipt2angle(i + 1, circle)
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

def ipt2angle(ipt, icircle):
  n = npts[icircle]
  ipt = ipt%n # occasionally convenient for ipt = -1
  return 2*pi*ipt/n

def angle2ipt(angle, icircle): #ipt that contains angle (note that ipt refers to the proximal point on the section)
  n = npts[icircle]
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
      n = npts[icircle]
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
      for ipt in range(npts[icircle]):
        o = overlap((ilayer, icircle, ipt), icircle + 1)
        cnt += len(o)
  print ("all %d circle to circle connections time %g" % (cnt, (time() - t)))


def test2(layer, circle, ipt):
  print(layer, circle, ipt)
  print("npts ", npts[circle], npts[circle+1])
  a = overlap((layer, circle, ipt), circle + 1)
  print(a)
  for x in a:
    b = overlap((layer, circle + 1, x[0]), circle)
    print(b)

  print("length d_angle of circle ", circle, distance(xyz(layer, circle, 0), xyz(layer, circle, 1)), ipt2angle(1, circle))
  print("length d_angle of circle ", circle+1, distance(xyz(layer, circle+1, 0), xyz(layer, circle+1, 1)), ipt2angle(1, circle+1))
  
  for i, y in enumerate([a, b]):
    c = circle+(1-i)
    for x in y:
      print([layer, c, x[0]], xyz(layer, c, x[0]))
  print ()

def test3(l1, c1, c2):
  # sum of all fractions around each circle should total npts for the circle
  n1 = npts[c1]
  n2 = npts[c2]
  s1 = 0.
  s2 = 0.
  for i in range(n1):
    a = overlap((l1, c1, i), c2)
    for o in a:
      s1 += o[1]
      s2 += o[2]
  print("ll=%d c1=%d c2=%d n1=%d n2=%d s1=%g s2=%g" % (l1, c1, c2, n1, n2, s1, s2))

if __name__ == "__main__":
  test1()
  c = 1
  n = npts[c]
  test2(1, c, 0)
  test2(1, c, n-1)
  print ()
  test2(1, c, 1)
  test2(1, c, n-2)

  test3(1, 0, 1)
  test3(4, 101, 100)
