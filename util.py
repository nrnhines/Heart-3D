from math import sqrt, pi
from morphdef import distance

def isclose(a, b, abs_tol=1e-7):
  r = abs(a - b) < abs_tol
  if r == False:
    print("isclose %20g %20g %g" % (a, b, abs_tol))
  return r

# area of 4 xyz pts (with pt[1] = 0.0)
def end_area(pts):
  d = [distance(pts[i], pts[(i+1)%4]) for i in range(4)]
  x = distance(pts[0], pts[2])
  return accurate_triangle_area(d[0], d[1], x) + accurate_triangle_area(x, d[2], d[3])

# 2 xyz pts with pt[1] = 0.0 so trhat x == r and a subtend angle in radians
def side_area(p1, p2, angle):
  x = frustum_area(p2[2] - p1[2], p1[0], p2[0])
  return angle*x/(2.*pi)

def area3pt(pts):
  d = [distance(pts[i], pts[(i+1)%3]) for i in range(3)]
  return accurate_triangle_area(*d)

n_triang_zero_ = 0
def accurate_triangle_area(x, y, z):
  global n_triang_zero_
  # x,y,z sides of triangle
  # from http://http.cs.berkeley.edu/~wkahan/Triangle.pdf
  # W. Kahan
  a = [x, y, z]
  a.sort()
  if (a[0] - (a[2] - a[1])) > 0:
    x = .25*sqrt((a[2]+(a[1]+a[0])) * (a[0]-(a[2]-a[1])) \
      * (a[0]+(a[2]-a[1])) * (a[2]+(a[1]-a[0])))
    return x
  #print("accurate_triangle_area: not a triangle: " + str(a))
  #assert((a[0] - (a[2] - a[1])) > 0)
  n_triang_zero_ += 1
  return 0.0

def n_triang_zero():
  return n_triang_zero_

def frustum_area(h, r1, r2):
  return pi*(r1 + r2)*sqrt((r1 - r2)**2 + h**2)

if __name__ == "__main__":
  print(accurate_triangle_area(3., 4., 5.))
  print(accurate_triangle_area(4., 5., 3.))
  print(accurate_triangle_area(5., 3., 4.))

  print(frustum_area(1., 2., 3.))
  print(frustum_area(1., 3., 2.))
  print(frustum_area(0., sqrt(.5), 1.))
  print(frustum_area(0., 1., sqrt(.5)))
  print(side_area((0., 0., 0.), (1., 0., 0.), pi))

  print(end_area([(0., 0., 0.), (1., 0., 0.), (1., 0., 1.), (0., 0., 1.)]))

