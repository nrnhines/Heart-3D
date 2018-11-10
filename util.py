from math import sqrt, pi

def accurate_triangle_area(x, y, z):
  # x,y,z sides of triangle
  # from http://http.cs.berkeley.edu/~wkahan/Triangle.pdf
  # W. Kahan
  a = [x, y, z]
  a.sort()
  assert ((a[0] - (a[2] - a[1])) > 0)
  x = .25*sqrt((a[2]+(a[1]+a[0])) * (a[0]-(a[2]-a[1])) \
    * (a[0]+(a[2]-a[1])) * (a[2]+(a[1]-a[0])))
  return x

def frustum_area(h, r1, r2):
  return pi*(r1 + r2)*sqrt((r1 - r2)**2 + h**2)

if __name__ == "__main__":
  print(accurate_triangle_area(3., 4., 5.))
  print(accurate_triangle_area(4., 5., 3.))
  print(accurate_triangle_area(5., 3., 4.))

  print(frustum_area(1., 2., 3.))
  print(frustum_area(1., 3., 2.))

