from common import h, rank, nhost
from cellconread import gidinfo, connections
  
try:
  import cPickle as pickle
except:
  import pickle

#record GJ conductance specified by cid argument.                 
#Note cellconread.connection[cid] = (gid1, gid2) and cid is the line
# number for connections.txt. On the rank where gid1 exists, we find 
# gidinfo[gid1].gaps[cid] = HalfGap instance

tvec = None
gj_trajec_dict = {} # {cid : gvec}

def gj_record(cid):
  global tvec
  if cid not in connections: return

  gid = connections[cid][0]
  # Note. if all HalfGap associated with gid are wanted and you want the
  # arg to mean gid, instead of cid, then consider
  # for cid, halfgap in gidinfo[gid].gaps.items():
  if gid in gidinfo:
    halfgap = gidinfo[gid].gaps[cid]
    gvec = h.Vector()
    gvec.record(halfgap, halfgap._ref_gmax)
    gj_trajec_dict[cid] = gvec
    if tvec is None:
      tvec = h.Vector()
      tvec.record(halfgap, h._ref_t)

def gj_out(prefix):
  if tvec is None: return

  fname = "%s_%04d.dat" % (prefix, rank)
  f = open(fname, "w")
  pickle.dump((tvec, gj_trajec_dict), f)
  f.close()

def gj_in(prefix, i, j):
  fname = "%s_%04d.dat" % (prefix, rank)
  f.open(fname, "r")
  result = pickle.load(f)
  f.close()
  return result
