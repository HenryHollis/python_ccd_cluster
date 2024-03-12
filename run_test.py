#import leidenalg
import louvain
import leidenalg
import numpy as np
import igraph as ig
import random
import time
from sklearn.metrics.cluster import adjusted_rand_score
random.seed(42) #in this case only needed so that Erdos_renyi graph stays const
np.random.seed(42)
# emat  = np.random.rand(12,3) #creates bug
emat  = np.random.rand(12,50) 
# Calculate the correlation matrix
refcor = np.random.rand(12, 100)
refcor = np.corrcoef(refcor.T, rowvar=False)
print(refcor)
# ccd = leidenalg.calcCCD(correlation_matrix, emat)
# print("CCD: {:.6f}".format(ccd))

#G = ig.Graph(n=50, edges=[[0, 1], [1,2], [2,3],[3,4], [4,0], [4,1], [4,2]])
G = ig.Graph.Erdos_Renyi(n=50, p=.2) #creates bug



#ig.plot(G)
t0 = time.time()
part = louvain.find_partition(G, louvain.ccdModularityVertexPartition, emat, refmat=refcor, seed = 42)
t1 = time.time()
print("time: {}".format(t1-t0))
print(part._membership)
#part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition, seed = 42)
# part3 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition, seed = 42)
#part2 = louvain.find_partition(G, louvain.ccdModularityVertexPartition)





# part2 = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition)
# print(part2._membership[1:10])
# Calculate Adjusted Rand Index
# ari = adjusted_rand_score(part2._membership, part._membership)

# Print the ARI
# print("Adjusted Rand Index:", ari)