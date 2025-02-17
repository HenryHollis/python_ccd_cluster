# %%
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import scipy

import random
import leidenalg
import louvain
import time
import math

# %%
adata = sc.read_h5ad("Greenham/Data/Greenham.h5ad")

# %%
# Copy the 'counts' layer to 'X'
adata.X = adata.layers['counts'].copy()

# %%
np.any(np.isnan(adata.X.data))

# %%
adata.obs['phase'] = adata.obs['orig.ident'].str.extract(r'(\d+)')[0]

# %% [markdown]
# The Seurat Umap from Greenham with her clustering

# %%
sc.pl.umap(adata, color='seurat_clusters')  # replace 'gene1', 'gene2' with the names of genes or metadata you want to color by

# %% [markdown]
# Calculate my own UMAP, color by Greenham's clusters

# %%
print(type(adata.X))
print(adata.X.shape)

# %%
adata.var['highly_variable'] = adata.var['highly_variable'] == 1

# %%
# If UMAP has not been computed, compute it (optional)
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# # Identify highly variable genes (optional)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# # Subset to highly variable genes (optional)
# adata = adata[:, adata.var.highly_variable]

# Scale the data to unit variance and zero mean
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=30)

sc.pp.neighbors(adata, use_rep = 'X_pca', n_neighbors= 20)
sc.tl.umap(adata)

# Plot the UMAP
sc.pl.umap(adata, color='seurat_clusters')  


# %% [markdown]
# Create a Pseudobulk

# %%
num_subjects = len(np.unique(adata.obs['orig.ident']))

pseudobulk_samples = np.zeros((num_subjects, len(adata.var.index)))

# Group by subjects and sum the expression counts
grouped = adata.obs.groupby('orig.ident').indices

for sub, cell_indices in grouped.items():
    sub_idx = list(grouped.keys()).index(sub)
    pseudobulk_samples[sub_idx, :] = adata.layers['counts'][cell_indices, :].sum(axis=0)

# Create a new AnnData object for pseudobulk samples
adata_pseudobulk = AnnData(X=pseudobulk_samples)
adata_pseudobulk.obs['subject'] = np.unique(adata.obs['orig.ident'])
adata_pseudobulk.var = adata.var
# Print the shape of the pseudobulk samples
print("Pseudobulk samples shape:", adata_pseudobulk.X.shape)

# Display the pseudobulk data
print("Pseudobulk data:")
print(adata_pseudobulk.X)
adata

# %%
subjects = [list(grouped.keys()).index(i) for i in adata.obs['orig.ident']]
subjects_vect = np.array(subjects).astype(np.int32).reshape(1, -1)

# %%
# List of core clock genes
genes_of_interest = ['AT3G09600','AT5G17300','AT1G01060','AT2G46830','AT3G54500','AT5G06980',
                     'AT3G12320','AT5G02810','AT5G37260','AT2G25930','AT5G61380','AT4G39260', 'AT2G21660'] 
#took out AT1G01520 with low expression

# Ensure that the gene names exist in the AnnData object
genes_of_interest = [gene for gene in genes_of_interest if gene in adata_pseudobulk.var_names]

# Subset the AnnData object
adata_subset = adata_pseudobulk[:, genes_of_interest]

emat = adata_subset.X
refmat = np.array(scipy.stats.spearmanr(emat))[0,:,:]

# Extract the gene symbols from 'adata.vars'
gene_symbols = adata.var.index  # Assuming the gene symbols are the index of 'adata.var'

# Create a dictionary mapping gene symbols to their corresponding indices
gene_indices = {gene: idx for idx, gene in enumerate(gene_symbols)}

# Get the indices of the genes of interest in the same order as 'gene_list'
indices_list = [gene_indices[gene] for gene in genes_of_interest if gene in gene_indices]

print(indices_list)


# %%
import matplotlib.pyplot as plt

plt.imshow( refmat, cmap = "RdBu" ,vmin=-1, vmax=1)
plt.colorbar()
plt.title( "Reference Matrix" )
plt.show()

# %%
from scanpy import _utils

connectivities = adata.obsp['connectivities']
g = _utils.get_igraph_from_adjacency(connectivities, directed=False)
print(type(g))


# %%
emat = adata.X.T
emat = emat[indices_list,:]
print(emat.shape)
print(subjects_vect.shape)


# %%
part = louvain.find_partition(
                    g,
                    louvain.ModularityVertexPartition, emat, refmat, seed=4
                )
membership_louvainStock= part._membership

# %%
part_ccd = louvain.find_partition(
                    g,
                    louvain.ccdModularityVertexPartition, emat, refmat,
                    subject_info=subjects_vect, seed=4, ccs_weight=10.
                )
membership = part_ccd._membership



# %%
from sklearn.metrics.cluster import adjusted_rand_score


import warnings
# Suppress specific warnings
# Suppress all warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    adata.obs['louvainStock'] = [str(i) for i in membership_louvainStock]
    adata.obs['louvainccd'] = [str(i) for i in membership]
    # Calculate UMAP
    sc.tl.umap(adata)
    pcs_to_plot = ['1,2', '2,3', '1,3', '1, 4']
    # Plot UMAP with Louvain clusters
    sc.pl.umap(adata, color='louvainccd', legend_loc='on data')
    sc.pl.umap(adata, color='louvainStock', legend_loc='on data')
    sc.pl.pca(adata, color= 'louvainccd' , components = pcs_to_plot, show=True)

# ari = adjusted_rand_score(membership,adata.obs['seurat_clusters'].tolist())
ari = adjusted_rand_score(membership,membership_louvainStock)
# Print the ARI
print("Adjusted Rand Index:", ari)

# %%
#Read in JTK results
aryth = pd.read_csv("Greenham/Data/GreenhamJTKNonCyclersBHQmoreThan2.csv")

#Select the first column with the gene names
no_cycle_genes = aryth.CycID

# Ensure the genes in the list are present in the anndata object
gene_list = [gene for gene in no_cycle_genes if gene in adata.var_names]

# Subset the anndata object
adata_noCycle = adata[:, gene_list].copy()

# %%
# Copy the 'counts' layer to 'X'
adata_noCycle.X = adata_noCycle.layers['counts'].copy()

# %%
# # If UMAP has not been computed, compute it (optional)
# sc.pp.normalize_total(adata_noCycle, target_sum=1e4)

# # Log-transform the data
# sc.pp.log1p(adata_noCycle)

# # Scale the data to unit variance and zero mean
# sc.pp.scale(adata_noCycle, max_value=10)
# sc.pp.pca(adata_noCycle, n_comps=30)

# sc.pp.neighbors(adata_noCycle, use_rep = 'X_pca', n_neighbors= 20)
# sc.tl.umap(adata_noCycle)

# Plot the UMAP
sc.pl.umap(adata_noCycle, color='seurat_clusters')  # replace 'gene1', 'gene2' with the names of genes or metadata you want to color by
sc.pl.umap(adata_noCycle, color='louvainccd')  # replace 'gene1', 'gene2' with the names of genes or metadata you want to color by


# %%
import matplotlib.pyplot as plt

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
    
    return result

def get_corgram(membership, clusterID):
    cluster_A_cells = [i for i in range(len(membership)) if membership[i] == clusterID]
    cluster_A_emat = adata.X.T[np.ix_(indices_list, cluster_A_cells)]
    clusterA_groups = [subjects[i] for i in cluster_A_cells]
    pseudoBulk = sumByGroup(cluster_A_emat, clusterA_groups)
    print("Cluster {} contains {} samples".format(clusterID, len(np.unique(clusterA_groups))))
    print(np.array(clusterA_groups).astype(np.int32).reshape(1, -1))
    ccs = louvain.calcCCS(refmat, cluster_A_emat, np.array(clusterA_groups).astype(np.int32).reshape(1, -1))
    # Get correlation matrix:
    corr_mat = np.array(scipy.stats.spearmanr(pseudoBulk.T))[0,:,:]
    return(corr_mat, ccs)

def plot_corgrams(membership):
    num_clusters = len(np.unique(membership))
    cols = 2
    rows = math.ceil((num_clusters+1) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6, rows*2))
    # Flatten the axes array
    axes = axes.flatten()

    # Iterate over the list and plot each item in a subplot
    for i in range(num_clusters):
        corgram, ccs = get_corgram(membership, i)
        heatmap = axes[i].imshow(corgram, cmap = "RdBu")
        axes[i].set_title('Cluster: {}, CCS: {:.4f}'.format(i, ccs))
        fig.colorbar(heatmap, ax=axes[i])

    heatmap_ref = axes[-1].imshow(refmat, cmap = "RdBu")
    axes[-1].set_title('Reference')
    fig.colorbar(heatmap_ref, ax=axes[-1])

    # # Hide any unused subplots
    # for j in range(i + 1, len(axes)):
    #     fig.delaxes(axes[j])

    # Display the plot
    plt.tight_layout()
    plt.show()

plot_corgrams(membership)
plot_corgrams(membership_louvainStock)


