from neuron import h
pc = h. ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())

rank = 0
nhost = 1000
