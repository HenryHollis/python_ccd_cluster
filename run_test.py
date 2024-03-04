#import leidenalg
import louvain
import leidenalg
import numpy as np
import igraph as ig
import random
from sklearn.metrics.cluster import adjusted_rand_score
random.seed(42) #in this case only needed so that Erdos_renyi graph stays const
# import os
# pid = os.getpid()
# print(pid)
#emat = np.array([[1, 2, 3], [1, 2, 34], [5, 6, 7], [6, 7, 8], [6, 7, 8]], dtype = np.float64)
emat  = np.random.rand(12,1000)

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(emat.T, rowvar=False)
print(correlation_matrix.shape)
print(emat.shape)

G = ig.Graph.Erdos_Renyi(1000, 0.2)
part = leidenalg.find_partition(G, leidenalg.ccdModularityVertexPartition, emat, correlation_matrix, seed = 42)
part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition, seed = 42)
# part3 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition, seed = 42)
#part2 = louvain.find_partition(G, louvain.ccdModularityVertexPartition)





# part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)
# print(part2._membership[1:10])
# Calculate Adjusted Rand Index
ari = adjusted_rand_score(part2._membership, part._membership)

# Print the ARI
print("Adjusted Rand Index:", ari)