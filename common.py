from neuron import h
pc = h. ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())

