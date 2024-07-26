#import leidenalg
import louvain
import leidenalg
import numpy as np
import random
import scipy
import math
random.seed(42) #in this case only needed so that Erdos_renyi graph stays const
np.random.seed(42)

num_cells = 100
num_subjects = 5
# emat  = np.random.rand(12,3) #creates bug
emat  = np.random.rand(12,num_cells) 
subjects = np.repeat(range(num_subjects), math.ceil(num_cells/num_subjects)).reshape(1,-1).astype(np.int32)
subjects = subjects[:, 0:num_cells]
print(subjects.shape)
# Calculate the correlation matrix
refcor = np.random.rand(12, 100)
refcor = np.corrcoef(refcor.T, rowvar=False)
print(refcor)

#True calcCCD formula from Hughey
def calcDist(r1, r2):
    tmp = r1-r2
    tmp = tmp ** 2
    return(math.sqrt(np.sum(tmp)))
def calcDistNull(r1):
    np.fill_diagonal(r1, 0.)
    tmp = r1 ** 2
    return(math.sqrt(np.sum(tmp)))

def calcCCDsimple(ref, emat):
    upper_ref = np.triu(ref)
    corr_mat = np.array(scipy.stats.spearmanr(emat.T))[0,:,:]
    upper_corrmat = np.triu(corr_mat)
    return(calcDist(upper_corrmat, upper_ref))

def sumByGroup(matrix, groups):
    # Convert inputs to numpy arrays for easier manipulation
    matrix = np.array(matrix)
    groups = np.array(groups).flatten()
    # Get the unique groups
    unique_groups = np.unique(groups)
    # Initialize the result matrix with zeros
    result = np.zeros((matrix.shape[0], len(unique_groups)))
    
    # Sum the columns of the matrix according to the groups
    for i, group in enumerate(unique_groups):
        result[:, i] = matrix[:, groups == group].sum(axis=1)
    return(result)

def calcCCS(ref, emat):
    upper_ref = np.triu(ref)
    corr_mat = np.array(scipy.stats.spearmanr(emat.T))[0,:,:]
    upper_corrmat = np.triu(corr_mat)
    ccd = calcDist(upper_corrmat, upper_ref)
    ncd = calcDistNull(upper_corrmat)
    return(ncd-ccd)

print("CCD expected: {:.4f}".format(calcCCDsimple(refcor, emat)))
print("CCD leiden: {:.4f}".format(leidenalg.calcCCD(refcor, emat)))
print("CCD louvain: {:.4f}".format(louvain.calcCCD(refcor, emat)))
grouped_emat = sumByGroup(emat, subjects)
print("CCS expected: {:.4f}".format(calcCCS(refcor, grouped_emat)))
print("CCS louvain: {:.4f}".format(louvain.calcCCS(refcor, emat, subjects)))
print("CCS leiden: {:.4f}".format(leidenalg.calcCCS(refcor, emat, subjects)))
