#import leidenalg
import louvain
import leidenalg
import numpy as np
import igraph as ig
import random
import time
from sklearn.metrics.cluster import adjusted_rand_score
random.seed(142) #in this case only needed so that Erdos_renyi graph stays const
np.random.seed(142)


# emat  = np.random.rand(12,3) #creates bug
emat  = np.random.rand(12,2000) 
# Calculate the correlation matrix
refcor = np.random.rand(12, 100)
refcor = np.corrcoef(refcor.T, rowvar=False)
print(refcor)
# ccd = leidenalg.calcCCD(correlation_matrix, emat)
# print("CCD: {:.6f}".format(ccd))

# G = ig.Graph(n=6, edges=[[0, 1], [1,2], [2,0], [2,3],[3,4], [4,5], [5,3]])
G = ig.Graph.Erdos_Renyi(n=2000, p=.05) 


#ig.plot(G)
t0 = time.time()
part = louvain.find_partition(G, louvain.ccdModularityVertexPartition, emat, refmat=refcor, seed = 42)
t1 = time.time()
print("time: {}".format(t1-t0))

# _plot(G, "/Users/henryhollis/Desktop/ccd_clustering.png", part._membership)
print("louvain ccd: # unique clusters:")
print((part._membership))

t0 = time.time()
part2 = louvain.find_partition(G, louvain.ModularityVertexPartition, emat, refmat=refcor, seed = 42)
t1 = time.time()
print("time2: {}".format(t1-t0))

# with open('/Users/henryhollis/Desktop/run_test_membership.txt', 'a') as file:
#     for item in part._membership:
#         file.write(f"{item} ")
#     file.write("-------------------")
# part2 = louvain.find_partition(G, leidenalg.ModularityVertexPartition, seed = 42)
# print("louvain stock: # unique clusters:")
# print(np.unique(part2._membership))
# _plot(G, "/Users/henryhollis/Desktop/stock_louvain.png", part2._membership)

# # Calculate Adjusted Rand Index
# ari = adjusted_rand_score(part2._membership, part._membership)

# # Print the ARI
# print("Adjusted Rand Index:", ari)