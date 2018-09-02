from neuron import h, gui
from cellorg import org2gid, xyz, p, nlayer, ncircle, npts, gid2org

class View():
  def __init__(self, pts):
    diam = p.cell_diameter
    self.pts = pts
    self.sections = {}
    self.sl =h.SectionList()
    for pt in pts:
      gid = org2gid(*pt)
      sec = h.Section(name=str(gid))
      sec.pt3dclear()
      x,y,z = xyz(*pt)
      sec.pt3dadd(x, y, z, diam)
      x,y,z = xyz(pt[0], pt[1], pt[2] + 1)
      sec.pt3dadd(x, y, z, diam)
      self.sections[sec] = gid
      self.sl.append(sec=sec)
    self.sh = h.Shape(self.sl)
    self.sh.menu_tool("print info", self.callback)

  def callback(self, type, x, y, keystate):
    #info about nearest section to mouse
    if type == 2:
      s = self.sh
      d = s.nearest(x, y)
      arc = s.push_selected()
      if arc >= 0:
        s.select()
        sec = h.cas()
        gid = self.sections[sec]
        ilayer, icircle, ipt = gid2org(gid)
        x,y,z = xyz(ilayer, icircle, ipt)
        print ("gid %d   id (%d, %d, %d) prox pt at (%g, %g, %g) length %g"%(gid, ilayer, icircle, ipt, x, y, z, sec.L))
        h.pop_section()

def test1():
  pts = []
  for ilayer in range(nlayer):
    for icircle in list(range(10)) + list(range(ncircle-5, ncircle)):
      for ipt in range(npts[icircle]):
        pts.append((ilayer, icircle, ipt))
  return View(pts)

if __name__ == "__main__":
  v = test1()
