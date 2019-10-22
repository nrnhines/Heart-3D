from neuron import h
from cellorg import gid2org, nlayer, ncircle, npts
from math import pi

'''
The purkinje spec is entirely in layer 0
The spec is a list of "rectangles" defined as
((angle, n), (bottom, top))
where angle specifies an ipt depending on icircle and n
is the number of cells, range(ipt, ipt+n). Note that if angle is
close to 360, the n may wrap over crossing point 0
Bottom, top are values
of icircle (positive means relative to the 0 bottom,  and negative means
relative to ncircle. First and last numbers are inclusive.
Thus, the list of (ilayer, icircle, ipt) cells in
((0, 2), (20, -1)) are
[(0, icircle, ipt) for icircle in range(20,ncircle[0]+1-10) \
  for ipt in range(3, npt[0][icircle]+1-1)]
'''
from param import purkinje_spec

# setup the pspec based on purkinje_spec to allow fast checking of inside
# by converting to actual (min, max) for icircle and ipt. Since
# (min, max) for ipt can differ depending on icircle, we need a list
# of them of length circlemax-circlemin. (So during is_purkinje_gap we
# do not have to deal with negative numbers.)
# We have a problem with regard to the points on circles where our
# desired region crosses 0. So if min > max, we reverse the sense of
# the inside and return True if ipt >= min or ipt <= max

pspec = []
for p in purkinje_spec:
  # first the circles which come from p[1]
  c = p[1]
  circles = c if c[1] >= 0 else (c[0], ncircle[0] + c[1])

  # now the pts which come from p[0]. Need circles[1] - circles[0] + 1
  # of them in case one or both of the p[0][0] or p[0][1] are negative
  angle = p[0][0]
  n = p[0][1]
  pts = []

  for icircle in range(circles[0], circles[1]+1):
    npoint = npts[0][icircle]
    first = int(angle/360*npoint)
    last = (first + n - 1)%npoint
    pts.append((first, last)) # last may be less than first, that is ok.

  pspec.append([circles, pts])

def is_purkinje_gap(gid1, gid2):
  return in_purkspec(gid2org(gid1)) and in_purkspec(gid2org(gid2))

def in_purkspec(ilayer, icircle, ipt):
  if ilayer != 0:
    return False
  for rec in pspec:
    if inside(icircle, rec[0]):
      if inside(ipt, rec[1][icircle-rec[0][0]]):
        return True
  return False

def inside(i, pair):
  if (pair[0] < pair[1]):
    return i >= pair[0] and i <= pair[1]
  else:
    return i >= pair[0] or i <= pair[1]

  
def show():
  g = h.Graph()
  g.size(0, npts[0][-1], 0, ncircle[0])
  for icircle in range(ncircle[0]):
    for ipt in range(npts[0][icircle]):
      if in_purkspec(0, icircle, ipt):
        draw_cell(g, ipt, icircle)
  return g
def draw_cell(g, x, y):
  g.beginline()
  g.line(x, y)
  g.line(x+1, y)
  g.line(x+1, y+1)
  g.line(x, y+1)
  g.line(x, y)
  g.flush()

if __name__ == "__main__":
  print("purkinje_spec ", purkinje_spec)
  g = show()
  #print("pspec", pspec)

