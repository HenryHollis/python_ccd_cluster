import leidenalg
import numpy as np
import igraph as ig
# import os
# pid = os.getpid()
# print(pid)
#emat = np.array([[1, 2, 3], [1, 2, 34], [5, 6, 7], [6, 7, 8], [6, 7, 8]], dtype = np.float64)
emat  = np.random.rand(12, 5)

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(emat.T, rowvar=False)

print(correlation_matrix)
#print(emat)

G = ig.Graph.Erdos_Renyi(5, 0.1)
part = leidenalg.find_partition(G, leidenalg.ccdModularityVertexPartition, emat = emat, refmat=correlation_matrix)
print(part._membership)

# part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)
# print(part2._membership)