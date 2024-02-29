import leidenalg
import numpy as np
import igraph as ig
import random
random.seed(42) #in this case only needed so that Erdos_renyi graph stays const
# import os
# pid = os.getpid()
# print(pid)
#emat = np.array([[1, 2, 3], [1, 2, 34], [5, 6, 7], [6, 7, 8], [6, 7, 8]], dtype = np.float64)
emat  = np.random.rand(12, 200)

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(emat.T, rowvar=False)
print(correlation_matrix.shape)
print(emat.shape)

G = ig.Graph.Erdos_Renyi(200, 0.01)
part = leidenalg.find_partition(G, leidenalg.ccdModularityVertexPartition, emat = emat, refmat=correlation_matrix)
print(part._membership[1:10])

# part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)
# print(part2._membership[1:10])