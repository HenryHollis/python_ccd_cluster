import leidenalg
import numpy as np
import igraph as ig
import os
pid = os.getpid()
print(pid)
emat = np.array([[1, 2, 3], [1, 2, 34]], dtype = np.float64)
# G = ig.Graph.Erdos_Renyi(100, 0.1)
# part = leidenalg.find_partition(G, leidenalg.ccdModularityVertexPartition, emat)
print(emat)
leidenalg._c_leiden.incmat(emat, emat)
print(emat)

G = ig.Graph.Erdos_Renyi(100, 0.1)
part = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)