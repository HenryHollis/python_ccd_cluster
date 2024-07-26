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

n = 2000
emat  = np.random.rand(12,n) 
# Calculate the correlation matrix
refcor = np.random.rand(12, 100)
refcor = np.corrcoef(refcor.T, rowvar=False)

k = 10  # Define the upper bound for the random numbers.
vec = np.random.randint(0, k, size=n)  # Generate 100 random integers in the range [0, k).
vec = np.array(vec).astype(np.int32).reshape(1,-1)

G = ig.Graph.Erdos_Renyi(n=n, p=.05) 


#ig.plot(G)
t0 = time.time()
part = louvain.find_partition(G, louvain.ccdModularityVertexPartition, emat, refmat=refcor, subject_info=vec, seed = 42)
t1 = time.time()
print("time: {}".format(t1-t0))

# _plot(G, "/Users/henryhollis/Desktop/ccd_clustering.png", part._membership)
print("louvain ccd: # unique clusters:")
print((part._membership))

# t0 = time.time()
# part2 = louvain.find_partition(G, louvain.ModularityVertexPartition, emat, refmat=refcor, seed = 42)
# t1 = time.time()
# print("time2: {}".format(t1-t0))
