from neuron import h, gui
from cellorg import org2gid, xyz, p, nlayer, ncircle, npts

class View():
  def __init__(self, pts):
    diam = p.cell_diameter
    self.pts = pts
    self.sections = []
    self.sl =h. SectionList()
    for pt in pts:
      gid = org2gid(*pt)
      sec = h.Section(name=str(gid))
      sec.pt3dclear()
      x,y,z = xyz(*pt)
      sec.pt3dadd(x, y, z, diam)
      x,y,z = xyz(pt[0], pt[1], pt[2] + 1)
      sec.pt3dadd(x, y, z, diam)
      self.sections.append(sec)
      self.sl.append(sec=sec)
    self.sh = h.Shape(self.sl)


def test1():
  pts = []
  for ilayer in range(nlayer):
    for icircle in list(range(10)) + list(range(ncircle-5, ncircle)):
      for ipt in range(npts[icircle]):
        pts.append((ilayer, icircle, ipt))
  return View(pts)

if __name__ == "__main__":
  v = test1()
