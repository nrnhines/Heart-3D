from neuron import h, gui
from cellorg import org2gid, xyz, p, nlayer, ncircle, npts, gid2org
from cellorg import angle2ipt, ipt2angle
from cellorg import gid_is_simulated
from morphdef import distance

class View():
  def __init__(self, pts, master=True):
    self.master = master
    diam = p.cell_diameter
    self.pts = pts
    self.sections = {}
    self.sl =h.SectionList()
    for pt in pts:
      gid = org2gid(*pt)
      sec = h.Section(name=str(gid))
      sec.pt3dclear()
      #draw the line a little shorter so we can see the junctions
      p1 = xyz(*pt)
      p2 = xyz(pt[0], pt[1], pt[2] + 1)
      dp = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
      x,y,z = p1[0] + .05*dp[0], p1[1] + .05*dp[1], p1[2] + .05*dp[2]
      sec.pt3dadd(x, y, z, diam-1)
      x,y,z = p2[0] - .05*dp[0], p2[1] - .05*dp[1], p2[2] - .05*dp[2]
      sec.pt3dadd(x, y, z, diam-1)
      self.sections[sec] = gid
      self.sl.append(sec=sec)
    self.sh = h.Shape(self.sl)
    self.sh.menu_tool("print info", self.callback)
    if not master:
      for sec in self.sections:
        self.sh.color(gid2org(self.sections[sec])[0] + 1, sec=sec)

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
        if self.master:
          neighborhood(ilayer, icircle, ipt)

def neighborhood(ilayer, icircle, ipt):
  global focusedview
  x,y,z = xyz(ilayer, icircle, ipt)
  gid = org2gid(ilayer, icircle, ipt)  
  angle = ipt2angle(ipt, ilayer, icircle)
  pts = []
  for jlayer in range(nlayer):
    for jcircle in range(max([icircle - 5, 0]), min([icircle + 1 + 5, ncircle[jlayer]])):
      npt = npts[jlayer][jcircle]
      n = min([int(npts[jlayer][jcircle]/2), 5])
      kpt = angle2ipt(angle, jlayer, jcircle)
      for jpt in [i%npt for i in range(kpt-n, kpt+1 + n)]:
        pts.append((jlayer, jcircle, jpt))
  focusedview = View(pts, master=False)
	
def test1(): # tip and base
  pts = []
  for ilayer in range(nlayer):
    for icircle in list(range(10)) + list(range(ncircle[ilayer]-5, ncircle[ilayer])):
      for ipt in range(npts[ilayer][icircle]):
        pts.append((ilayer, icircle, ipt))
  return View(pts)

def test2():
  pts = []
  for ilayer in range(1):
    for icircle in range(ncircle[ilayer]):
      for ipt in range(0, npts[ilayer][icircle], 5):
        if gid_is_simulated(org2gid(ilayer, icircle, ipt)):
          pts.append((ilayer, icircle, ipt))
  return View(pts)

if __name__ == "__main__":
  v = test2()
