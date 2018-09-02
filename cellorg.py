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

def xyz(ilayer, icircle, ipt):
  n = npts[icircle]
  ipt = ipt%n # wrap around
  pt = paraboloid[ilayer][icircle]
  a = 2*pi*ipt/n
  r = pt[0]
  return r*sin(a), r*cos(a), pt[2]

def test1():
  #iterate over all points
  paraboloid, npts = morphorg()
  for ilayer, layer in enumerate(paraboloid):
    for icircle, circle in enumerate(layer):
      n = npts[icircle]
      for ipt, pt in enumerate(circle_discrete(n, circle)):
        # do something with ilayer, icircle, ipt
        gid = org2gid(ilayer, icircle, ipt)
        org = gid2org(gid)
        assert(org == (ilayer, icircle, ipt))
        pass

if __name__ == "__main__":
  test1()
