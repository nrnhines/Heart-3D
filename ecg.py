from common import h, pc, rank, nhost
from net3d import gidinfo
from math import sqrt

resistivity_of_cytoplasm = 1000 #ohm-cm

nelectrode = 6

ecg_t = h.Vector()
ecg_v = [h.Vector() for _ in range(nelectrode)]
rx = [h.Vector() for _ in range(nelectrode)]
vx = [h.Vector() for _ in range(nelectrode)]
e_coord = [[0, 4000, 6000], [0, 4000, 5500], [0, 4000, 6500], [0, -4000, 6000], [0, -4000, 5500], [0, -4000, 6500]]

imem_ptrvec = None
imem_vec = None

def transfer_resistance(cell, exyz):
  xyz = [cell.x, cell.y, cell.z]
  for i, x in enumerate(exyz):
    xyz[i] -= x
    xyz[i] *= xyz[i]
  d = sqrt(sum(xyz))
  rx = resistivity_of_cytoplasm/d * (0.01)
  return rx

def ecg_init():
  global imem_ptrvec, imem_vec, rx, vx
  n = len(gidinfo) #or else count the segments
  imem_ptrvec = h.PtrVector(n)
  imem_vec = h.Vector(n)
  for i, cellinfo in enumerate(gidinfo.values()):
    seg = cellinfo.cell.soma(0.5)
    imem_ptrvec.pset(i, seg._ref_i_membrane_)
  rx = h.Matrix(nelectrode, n)
  vx = h.Vector(nelectrode)
  for i in range(nelectrode):
    for j, cellinfo in enumerate(gidinfo.values()):
      rx.setval(i, j, transfer_resistance(cellinfo.cell, e_coord[i]))
    #rx.setval(i,1,1.0)

def ecg():
  # print 'ecg t=%g' % pc.t(0)
  imem_ptrvec.gather(imem_vec)
  #s = pc.allreduce(imem_vec.sum(), 1) #verify sum i_membrane_ == stimulus
  #if rank == 0: print pc.t(0), s

  #sum up the weighted i_membrane_. Result in vx
  rx.mulv(imem_vec, vx)
    
  # append to Vector
  ecg_t.append(pc.t(0))
  for i in range(nelectrode): 
    ecg_v[i].append(vx.x[i])

def ecg_setup():
  global bscallback, fih
  h.cvode.use_fast_imem(1)
  bscallback = h.beforestep_callback(h.cas()(.5))
  bscallback.set_callback(ecg)
  fih = h.FInitializeHandler(1, ecg_init)

def ecg_final():
  #return
  for i in range(nelectrode):
    pc.allreduce(ecg_v[i], 1)
    #return

def ecgout():
  ecg_final()
  if rank is not 0: return
  f = open('ecg.txt', 'w')
  for i in range(len(ecg_t)):
    line = '%g'%ecg_t.x[i]
    for j in range(nelectrode):
      line += ' %g' % ecg_v[j].x[i]
    f.write(line + '\n')
  f.close()

def test():
  h.load_file("stdgui.hoc")
  #h.cvode_active(1)
  h('create soma[5]')

  h.tstop=.2
  ecg_setup()
  ecg_init()
  h.run()
  ecg_final()
  

if __name__ == '__main__':
  # test() invalid test
  for i in range(len(ecg_t)):
    print ecg_t.x[i],
    for j in range(nelectrode):
      print ecg_v[j].x[i],
    print ""

