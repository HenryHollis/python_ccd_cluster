# %%
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData

import random
import leidenalg
import louvain
import time
import math

# %%
import scipy.sparse as sp

def simulate_sparse_counts(num_genes, num_cells, sparsity_level=0.9, mean_expression=10, dispersion=1):
    """
    Simulate sparse single-cell RNA-seq counts.

    Parameters:
    - num_genes: Number of genes.
    - num_cells: Number of cells.
    - sparsity_level: Proportion of zero counts (default is 0.9).
    - mean_expression: Mean expression level (lambda parameter of Poisson).
    - dispersion: Dispersion parameter for negative binomial distribution.

    Returns:
    - sparse_matrix: Sparse matrix of gene expression counts.
    """
    
    # Calculate size parameter for negative binomial (equivalent to number of trials)
    size = (mean_expression ** 2) / (dispersion - mean_expression) if dispersion > mean_expression else 1
    
    # Generate negative binomial counts
    counts = np.random.negative_binomial(size, size / (size + mean_expression), ( num_cells, num_genes))
    counts = counts.astype(np.float64)
    # Convert to sparse matrix
    # sparse_matrix = sp.csr_matrix(counts)
    
    return counts



# %% [markdown]
# **Create Non_rhythmic Genes**

# %%
import matplotlib.pyplot as plt

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Create an AnnData object
adata = sc.AnnData()

# Define the number of cells and genes
n_cells = 5000
n_genes = 100
offset_multiplier = 5 #How different are the cell types
num_cell_types = 3
num_subjects = 12
sparsity_level = 0.  # percent of the counts that are zeros
mean_expression = 10  # Mean expression level
dispersion = 2  # Dispersion parameter for the negative binomial distribution
subject_variability = .1  # Variability in expression levels between subjects
amplitude_scaling_factor = np.random.uniform(0, 3) #amp will be this number * "mean_expression" from above


# Simulate the sparse counts
sparse_counts = simulate_sparse_counts(n_genes, n_cells, sparsity_level, mean_expression, dispersion)

# Repeat each string n times
offsets = random.choices( range(num_cell_types), k=n_cells)

# Simulate gene expression data
gene_expression_data = sparse_counts  # Zipf distribution for sparsity
for i in range(n_genes):
    gene_expression_data[:,i] += (np.array(offsets)*offset_multiplier)
# Create an AnnData object with the simulated data
adata = sc.AnnData(X=gene_expression_data)

# Plot the expression levels of the first gene
plt.figure(figsize=(10, 6))
plt.plot(adata.X[:,0], marker='o', linestyle='-', color='b')
plt.xlabel('Cell Index')
plt.ylabel('Expression Level')
plt.title('Expression Levels of the First Gene')
plt.grid(True)
plt.show()

# %% [markdown]
# **Create Rhythmic Genes**

# %%
np.random.seed(42)

# Generate gene expression data for 10 genes from a cosine wave
n_cells = adata.shape[0]
n_new_genes = 12
new_gene_names = [f'NewGene_{i}' for i in range(n_genes+1, n_genes+n_new_genes + 1)]


# Create a cosine wave for each gene
gene_expression_data = np.zeros((n_cells, n_new_genes))
for i in range(n_new_genes):
    # phase = np.random.uniform(0, 2 * np.pi)
    if i < 3:
        phase = np.random.uniform(-.5, .5)
    elif (i in [3, 4, 5]):
        phase =  np.random.uniform(3*np.pi/2-.5, 3*np.pi/2+ .5)

    else:
        phase = np.random.uniform(np.pi-.5, np.pi+.5)

    noise = np.random.normal(0, 2, n_cells)
    freq = np.pi*2/24
    subjects = np.repeat(range(num_subjects), math.ceil(n_cells/num_subjects))
    subjects = subjects[0:n_cells]
    timepoints = subjects/num_subjects*24
    gene_expression_data[:, i] = amplitude_scaling_factor *mean_expression  *  np.cos(freq * timepoints - phase ) + (np.array(offsets)*offset_multiplier) + noise
    


# Create a new AnnData object for the new genes
new_genes_adata = ad.AnnData(X=gene_expression_data, var=pd.DataFrame(index=new_gene_names))

# Concatenate the new genes with the existing adata
adata2 = ad.concat([adata, new_genes_adata], axis = 1)

# %% [markdown]
# Introduce Sparsity with Zero Masking

# %%
# Introduce sparsity
zero_mask = np.random.rand(n_cells, n_genes+n_new_genes) < sparsity_level
adata2.X[zero_mask] = 0
assert(zero_mask.shape == adata2.X.shape)

# %% [markdown]
# **Plot Synthetic Genes**

# %%
import matplotlib.pyplot as plt
# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot rhythmic genes
axes[0].set_title('Rhythmic Genes')
for i in range(n_genes, n_genes+1):
    axes[0].scatter(timepoints, adata2.X[:, i], label=f'Gene_{i + 1}')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Expression Level')

# Plot non-rhythmic genes
axes[1].set_title('First 10 Non-Rhythmic Genes')
for i in range(0, n_genes):
    axes[1].scatter(timepoints, adata2.X[:, i], label=f'Gene_{i + 1}')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Expression Level')
#axes[1].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()


# %% [markdown]
# **Add some housekeeping variables to the adata2 objects**

# %%
cell_identities = [f'batch_{i:.2f}' for i in offsets]
subjects_obs = [f's_{i}' for i in subjects]

# Add the 'Cell_Identity' column to adata.obs
adata2.obs['Cell_Identity'] = cell_identities
adata2.obs['subject'] = subjects_obs
adata2.obs['subject_num'] = subjects
adata2.obs['timepoints'] = timepoints

assert(len(cell_identities) == n_cells)
assert(len(subjects_obs) == n_cells)

# %%
# Perform PCA
sc.tl.pca(adata2)

pcs_to_plot = ['1,2', '2,3', '1,3', '1, 4']

# Plot PCA for the selected PCs
sc.pl.pca(adata2, color='Cell_Identity', components=pcs_to_plot, show=True)
sc.pl.pca(adata2, color='timepoints', components=pcs_to_plot, show=True)


# %% [markdown]
# **Add between subject varition**

# %%
# Add subject-specific variability
subject_effects = np.random.normal(1.0, subject_variability, (adata2.shape[1], num_subjects)) #num genes by num_subjects
for cell in range(n_cells):
    subject = subjects[cell]
    adata2.X[cell,:] = adata2.X[cell,:] * subject_effects[:, subject]
    # base_counts[:, cell] = base_counts[:, cell] * subject_effects[:, subject]


# %% [markdown]
# **Plot synthetic data again, with subject variation**

# %%
from matplotlib.colors import ListedColormap, BoundaryNorm
# Get the number of unique membership values
unique_memberships = np.unique(subjects)
num_colors = len(unique_memberships)

# Generate a random permutation of color indices
color_indices = np.random.permutation(num_colors)

# Create a discrete colormap with randomly ordered colors
cmap = plt.cm.get_cmap('hsv', num_colors)
colors = cmap(color_indices)
random_cmap = ListedColormap(colors)

# Define the boundaries and norm for the colormap
bounds = np.arange(-0.5, num_colors + 0.5, 1)
norm = BoundaryNorm(bounds, random_cmap.N)


# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot rhythmic genes
axes[0].set_title('Example Rhythmic Gene, showing samples by color')
for i in range(n_genes, n_genes+1):
   scatter =  axes[0].scatter(timepoints, adata2.X[:, i], c = subjects,cmap=random_cmap, label=f'Gene_{i + 1}')
# Adding a colorbar
cbar = plt.colorbar(scatter, ticks=np.arange(num_colors))
cbar.set_label('Membership')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Expression Level')

# Plot non-rhythmic genes
axes[1].set_title('First 10 Non-Rhythmic Genes')
for i in range(0, n_genes):
    axes[1].scatter(timepoints, adata2.X[:, i],label=f'Gene_{i + 1}')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Expression Level')
#axes[1].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()

# %%
# Perform PCA
sc.tl.pca(adata2)

# Specify the PCs you want to plot (e.g., PC2, PC3)
pcs_to_plot = ['1,2', '2,3', '1,3', '1, 4']

# Plot PCA for the selected PCs
sc.pl.pca(adata2, color='subject', components=pcs_to_plot, show=True)
sc.pl.pca(adata2, color='Cell_Identity', components=pcs_to_plot, show=True)

# Preprocess the data (e.g., log-transform and scale)
sc.pp.scale(adata2)

# Calculate the neighborhood graph
sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=10)  # Adjust parameters as needed


# %% [markdown]
# **Create All Cell Pseudobulk**

# %%
num_subjects = len(np.unique(adata2.obs['subject']))
pseudobulk_samples = np.zeros((num_subjects, n_genes+n_new_genes))

# Group by subjects and sum the expression counts
grouped = adata2.obs.groupby('subject_num').indices

for sub, cell_indices in grouped.items():
    pseudobulk_samples[sub, :] = adata2.X[cell_indices, :].sum(axis=0)

# Create a new AnnData object for pseudobulk samples
adata_pseudobulk = AnnData(X=pseudobulk_samples)
adata_pseudobulk.obs['subject'] = np.unique(subjects)
adata_pseudobulk.obs['timepoints'] = np.unique(timepoints)
# Print the shape of the pseudobulk samples
print("Pseudobulk samples shape:", adata_pseudobulk.X.shape)

# Display the pseudobulk data
print("Pseudobulk data:")
print(adata_pseudobulk.X)

# %%

# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Get the number of unique membership values
unique_memberships = np.unique(subjects)
num_colors = len(unique_memberships)

# Generate a random permutation of color indices
color_indices = np.random.permutation(num_colors)

# Create a discrete colormap with randomly ordered colors
cmap = plt.cm.get_cmap('hsv', num_colors)
colors = cmap(color_indices)
random_cmap = ListedColormap(colors)

# Define the boundaries and norm for the colormap
bounds = np.arange(-0.5, num_colors + 0.5, 1)
norm = BoundaryNorm(bounds, random_cmap.N)

# Plot rhythmic genes
axes[0].set_title('Example Rhtymic Gene, all-cell pseudobulk, colored by subject (1 color per dot)')
for i in range(n_genes, n_genes+1):
    scatter = axes[0].scatter( adata_pseudobulk.obs['timepoints'], adata_pseudobulk.X[:, i],c = adata_pseudobulk.obs['subject'], cmap= random_cmap, label=f'Gene_{i + 1}')
# Adding a colorbar
cbar = plt.colorbar(scatter, ticks=np.arange(num_colors))
cbar.set_label('Membership')

axes[0].set_xlabel('Time')
axes[0].set_ylabel('Expression Level')

# Plot non-rhythmic genes
axes[1].set_title('First 10 Non-Rhythmic Genes')
for i in range(0, n_genes):
    axes[1].scatter(adata_pseudobulk.obs['timepoints'], adata_pseudobulk.X[:, i], label=f'Gene_{i + 1}')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Expression Level')
#axes[1].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()
plt.hist(adata_pseudobulk.X[:, 1])

# %%
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import igraph
def _plot(g, membership=None, draw = "kk"):
    layouts = ["circle", "drl", "fr", "kk", "fr3d", "kk3d", "large", "lgl", "random", "rt", "tree", "rt_circular", "sphere"]
    if draw not in layouts:
        print("sorry, {} not viable layout option.".format(draw))
        return
    if membership is not None:
        gcopy = g.copy()
        edges = []
        edges_colors = []
        for edge in g.es():
            if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
                edges.append(edge)
                edges_colors.append("gray")
            else:
                edges_colors.append("black")
        gcopy.delete_edges(edges)
        layout = gcopy.layout(draw)
        g.es["color"] = edges_colors
    else:
        layout = g.layout(draw)
        g.es["color"] = "gray"
    visual_style = {}
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["edge_color"] = g.es["color"]
    # visual_style["bbox"] = (4000, 2500)
    visual_style["vertex_size"] = 15
    visual_style["layout"] = gcopy.layout(draw)
    visual_style["bbox"] = (1024, 768)
    visual_style["margin"] = 40
    #visual_style["edge_label"] = g.es["weight"]
    for vertex in g.vs():
        vertex["label"] = vertex.index
    if membership is not None:
        colors = []
        for i in range(0, max(membership)+1):
            colors.append('%06X' % randint(0, 0xFFFFFF))
        for vertex in g.vs():
            vertex["color"] = str('#') + colors[membership[vertex.index]]
        visual_style["vertex_color"] = g.vs["color"]
    fig, ax = plt.subplots()
    igraph.plot(g, target = ax, **visual_style)

# %%
from scanpy import _utils
import importlib
import warnings
from scanpy import logging as logg
from scipy import sparse
# from scanpy.tools._compat import old_positionals
from scanpy._utils import _choose_graph
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from natsort import natsorted
from collections.abc import Sequence
from types import MappingProxyType
from scipy.sparse import spmatrix
from collections.abc import Mapping, Sequence

try:
    from leidenalg.VertexPartition import MutableVertexPartition
except ImportError:

    class MutableVertexPartition:
        pass

    MutableVertexPartition.__module__ = "leidenalg.VertexPartition"
try:
    from louvain.VertexPartition import MutableVertexPartition
except ImportError:

    class MutableVertexPartition:
        pass

    MutableVertexPartition.__module__ = "louvain.VertexPartition"

from typing import TYPE_CHECKING, Literal, Any

#interface my custom clustering with scanpy:
def ccdCluster_leiden(
    adata: AnnData,
    emat,
    refmat,
    resolution: float = 1,
    *,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    random_state: _utils.AnyRandom = 0,
    key_added: str = "louvainccd",
    adjacency: sparse.spmatrix | None = None,
    directed: bool | None = None,
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: type[MutableVertexPartition] | None = None,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    copy: bool = False,
    flavor: Literal["leidenalg", "ipgraph"] = "leidenalg",
    **clustering_args,
) -> AnnData | None:
    """\
    Cluster cells into subgroups [Traag18]_.

    Cluster cells using the Leiden algorithm [Traag18]_,
    an improved version of the Louvain algorithm [Blondel08]_.
    It has been proposed for single-cell analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.
        Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain `(obs_key, list_of_categories)`.
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
        2 is faster and the default for underlying packages.
    partition_type
        Type of partition to use.
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, leiden looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, leiden looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Whether to copy `adata` or modify it inplace.
    flavor
        Which package's implementation to use.
    **clustering_args
        Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
        or :meth:`igraph.Graph.community_leiden` from `igraph`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['leiden']['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """
    if flavor not in {"igraph", "leidenalg"}:
        raise ValueError(
            f"flavor must be either 'igraph' or 'leidenalg', but '{flavor}' was passed"
        )
    igraph_spec = importlib.util.find_spec("igraph")
    if igraph_spec is None:
        raise ImportError(
            "Please install the igraph package: `conda install -c conda-forge igraph` or `pip3 install igraph`."
        )
    if flavor == "igraph":
        if directed:
            raise ValueError(
                "Cannot use igraph's leiden implemntation with a directed graph."
            )
        if partition_type is not None:
            raise ValueError(
                "Do not pass in partition_type argument when using igraph."
            )
    else:
        try:
            import leidenalg

            msg = 'Use of leidenalg is discouraged and will be deprecated in the future.  Please use `flavor="igraph"` `n_iterations=2` to achieve similar results.  `directed` must also be `False` to work with `igraph`\'s implementation.'
            warnings.warn(msg, FutureWarning)
        except ImportError:
            raise ImportError(
                "Please install the leiden algorithm: `conda install -c conda-forge leidenalg` or `pip3 install leidenalg`."
            )
    clustering_args = dict(clustering_args)

    start = logg.info("running Leiden ccd clustering")
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        adjacency = _utils._choose_graph(adata, obsp, neighbors_key)
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories=restrict_categories,
            adjacency=adjacency,
        )
    # Prepare find_partition arguments as a dictionary,
    # appending to whatever the user provided. It needs to be this way
    # as this allows for the accounting of a None resolution
    # (in the case of a partition variant that doesn't take it on input)
    clustering_args["n_iterations"] = n_iterations
    #if resolution is not None:
        #clustering_args["resolution_parameter"] = resolution
    if flavor == "leidenalg":
        directed = True if directed is None else directed
        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if partition_type is None:
            partition_type = leidenalg.ccdModularityVertexPartition
        if use_weights:
            clustering_args["weights"] = np.array(g.es["weight"]).astype(np.float64)
        clustering_args["seed"] = random_state
        part = leidenalg.find_partition(g, partition_type, emat, refmat, **clustering_args)
    else:
        g = _utils.get_igraph_from_adjacency(adjacency, directed=False)
        if use_weights:
            clustering_args["weights"] = "weight"
        clustering_args.setdefault("objective_function", "modularity")
        with _utils.set_igraph_random_state(random_state):
            part = g.community_leiden(**clustering_args)
    # store output into adata.obs
    groups = np.array(part.membership)
    if restrict_to is not None:
        if key_added == "leidenccd":
            key_added += "_R"
        groups = rename_groups(
            adata,
            key_added=key_added,
            restrict_key=restrict_key,
            restrict_categories=restrict_categories,
            restrict_indices=restrict_indices,
            groups=groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    # store information on the clustering parameters
    adata.uns["louvainccd"] = {}
    adata.uns["leidenccd"]["params"] = dict(
        resolution=resolution,
        random_state=random_state,
        n_iterations=n_iterations,
    )
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key_added!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return (adata, g )

def cluster_louvain(
    adata: AnnData,
    emat,
    refmat,
    sample_ids = None,
    resolution: float | None = None,
    *,
    random_state: _utils.AnyRandom = 0,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    key_added: str = "louvainccd",
    adjacency: spmatrix | None = None,
    flavor: Literal["vtraag", "igraph", "rapids"] = "vtraag",
    directed: bool = False,
    use_weights: bool = False,
    partition_type: type[MutableVertexPartition] | None = None,
    partition_kwargs: Mapping[str, Any] = MappingProxyType({}),
    neighbors_key: str | None = None,
    obsp: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """\
    Cluster cells into subgroups [Blondel08]_ [Levine15]_ [Traag17]_.

    Cluster cells using the Louvain algorithm [Blondel08]_ in the implementation
    of [Traag17]_. The Louvain algorithm has been proposed for single-cell
    analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first,
    or explicitly passing a ``adjacency`` matrix.

    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        For the default flavor (``'vtraag'``) or for ```RAPIDS```, you can provide a
        resolution (higher resolution means finding more and smaller clusters),
        which defaults to 1.0.
        See “Time as a resolution parameter” in [Lambiotte09]_.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain ``(obs_key, list_of_categories)``.
    key_added
        Key under which to add the cluster labels. (default: ``'louvain'``)
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    flavor
        Choose between to packages for computing the clustering.

        ``'vtraag'``
            Much more powerful than ``'igraph'``, and the default.
        ``'igraph'``
            Built in ``igraph`` method.
        ``'rapids'``
            GPU accelerated implementation.

            .. deprecated:: 1.10.0
                Use :func:`rapids_singlecell.tl.louvain` instead.
    directed
        Interpret the ``adjacency`` matrix as directed graph?
    use_weights
        Use weights from knn graph.
    partition_type
        Type of partition to use.
        Only a valid argument if ``flavor`` is ``'vtraag'``.
    partition_kwargs
        Key word arguments to pass to partitioning,
        if ``vtraag`` method is being used.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, louvain looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, louvain looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    copy
        Copy adata or modify it inplace.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obs['louvain' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['louvain']['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """
    partition_kwargs = dict(partition_kwargs)
    start = logg.info("running Louvain clustering")
    if (flavor != "vtraag") and (partition_type is not None):
        raise ValueError(
            "`partition_type` is only a valid argument " 'when `flavour` is "vtraag"'
        )
    adata = adata.copy() if copy else adata
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories=restrict_categories,
            adjacency=adjacency,
        )
    if flavor in {"vtraag", "igraph"}:
        if flavor == "igraph" and resolution is not None:
            logg.warning('`resolution` parameter has no effect for flavor "igraph"')
        if directed and flavor == "igraph":
            directed = False
        if not directed:
            logg.debug("    using the undirected graph")
        g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
        if use_weights:
            weights = np.array(g.es["weight"]).astype(np.float64)
        else:
            weights = None
        if flavor == "vtraag":
            import louvain

            if partition_type is None:
                partition_type = louvain.ccdModularityVertexPartition
            if resolution is not None:
                partition_kwargs["resolution_parameter"] = resolution
            if use_weights:
                partition_kwargs["weights"] = weights
                louvain.set_rng_seed(random_state)
            else:
                partition_kwargs["seed"] = random_state
            logg.info('    using the "louvain" package of Traag (2017)')
            part = louvain.find_partition(
                graph = g,
                partition_type = partition_type,
                emat = emat,
                refmat = refmat,
                subject_info= sample_ids,
                **partition_kwargs
            )
            # adata.uns['louvain_quality'] = part.quality()
        else:
            part = g.community_multilevel(weights=weights)
        groups = np.array(part.membership)
    elif flavor == "rapids":
        msg = (
            "`flavor='rapids'` is deprecated. "
            "Use `rapids_singlecell.tl.louvain` instead."
        )
        warnings.warn(msg, FutureWarning)
        # nvLouvain only works with undirected graphs,
        # and `adjacency` must have a directed edge in both directions
        import cudf
        import cugraph

        offsets = cudf.Series(adjacency.indptr)
        indices = cudf.Series(adjacency.indices)
        if use_weights:
            sources, targets = adjacency.nonzero()
            weights = adjacency[sources, targets]
            if isinstance(weights, np.matrix):
                weights = weights.A1
            weights = cudf.Series(weights)
        else:
            weights = None
        g = cugraph.Graph()

        if hasattr(g, "add_adj_list"):
            g.add_adj_list(offsets, indices, weights)
        else:
            g.from_cudf_adjlist(offsets, indices, weights)

        logg.info('    using the "louvain" package of rapids')
        if resolution is not None:
            louvain_parts, _ = cugraph.louvain(g, resolution=resolution)
        else:
            louvain_parts, _ = cugraph.louvain(g)
        groups = (
            louvain_parts.to_pandas()
            .sort_values("vertex")[["partition"]]
            .to_numpy()
            .ravel()
        )
    else:
        raise ValueError('`flavor` needs to be "vtraag" or "igraph" or "taynaud".')
    if restrict_to is not None:
        if key_added == "louvainccd":
            key_added += "_R"
        groups = rename_groups(
            adata,
            key_added=key_added,
            restrict_key=restrict_key,
            restrict_categories=restrict_categories,
            restrict_indices=restrict_indices,
            groups=groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )
    adata.uns["louvainccd"] = {}
    adata.uns["louvainccd"]["params"] = dict(
        resolution=resolution,
        random_state=random_state,
    )
    logg.info(
        "    finished",
        time=start,
        deep=(
            f"found {len(np.unique(groups))} clusters and added\n"
            f"    {key_added!r}, the cluster labels (adata.obs, categorical)"
        ),
    )
    return (adata, g )


# correlation_matrix = np.corrcoef(emat.T, rowvar=False)


# %% [markdown]
# **Show Expression Matrix**

# %%
import scipy
from matplotlib.colors import Normalize

emat = adata2.X.T[n_genes:n_genes + n_new_genes,:]

# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot rhythmic genes
axes[0].set_title('Rhythmic Genes')
for i in range(0,1):
    axes[0].scatter(adata2.obs['timepoints'], emat[i, :], label=f'Gene_{i + 1}')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Expression Level')

# Plot correlation matrix:
corr_mat = np.array(scipy.stats.spearmanr(emat.T))[0,:,:]

axes[1].set_title('correlation matrix')
im = axes[1].imshow(corr_mat, cmap = "RdBu",  norm=Normalize(vmin=-1, vmax=1))
# Add the color bar
cbar = axes[1].figure.colorbar(im, ax = axes[1])
#axes[1].legend()
cbar.ax.set_ylabel("Spearman", rotation = -90, va = "bottom")

# Add text annotations to each box
for i in range(corr_mat.shape[0]):
    for j in range(corr_mat.shape[1]):
        text = axes[1].text(j, i, f'{corr_mat[i, j]:.2f}',
                           ha="center", va="center", color="black")
        
plt.tight_layout()

# Show the plot
plt.show()


# %% [markdown]
# Pseudobulk Matrix

# %%
import scipy
emat = adata_pseudobulk.X.T[n_genes:n_genes + n_new_genes,:]

# Create a figure with two subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot rhythmic genes
axes[0].set_title('Rhythmic Genes')
for i in range(0,1):
    axes[0].scatter(adata_pseudobulk.obs['timepoints'], emat[i, :], label=f'Gene_{i + 1}')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Expression Level')

# Plot correlation matrix:
corr_mat = np.array(scipy.stats.spearmanr(emat.T))[0,:,:]

axes[1].set_title('correlation matrix')
im = axes[1].imshow(corr_mat, cmap = "RdBu",  norm=Normalize(vmin=-1, vmax=1))
# Add the color bar
cbar = axes[1].figure.colorbar(im, ax = axes[1])
#axes[1].legend()
cbar.ax.set_ylabel("Spearman", rotation = -90, va = "bottom")

# Add text annotations to each box
for i in range(corr_mat.shape[0]):
    for j in range(corr_mat.shape[1]):
        text = axes[1].text(j, i, f'{corr_mat[i, j]:.2f}',
                           ha="center", va="center", color="black")
        
plt.tight_layout()

# Show the plot
plt.show()

# %% [markdown]
# **Show Reference Matrix**

# %%
import math
refmat = np.array([[1.0000000,  0.77547090,  0.72492855,  0.27817942, -0.63637681, -0.60375141, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.8473980],
                        [0.7754709,  1.00000000,  0.63439613,  0.07402797, -0.62632300, -0.34987550, -0.7461844, -0.6450780, -0.70865725, -0.7845410, -0.7654845, -0.7983427],
                        [0.7249286,  0.63439613,  1.00000000,  0.06541974, -0.59727560, -0.30024636, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101],
                        [0.2781794,  0.07402797,  0.06541974,  1.00000000, -0.01245765, -0.72253596, -0.4099044, -0.1411756,  0.25538496, -0.0252816, -0.3401805, -0.0781101],
                        [-0.6363768, -0.62632300, -0.59727560, -0.01245765,  1.00000000,  0.28367324,  0.6234166,  0.6454257,  0.59510653,  0.6712806,  0.6618797,  0.7597038],
                        [-0.6037514, -0.34987550, -0.30024636, -0.72253596,  0.28367324,  1.00000000,  0.6772739,  0.4242223, -0.06776682,  0.3366267,  0.6955807,  0.3810191],
                        [-0.8614806, -0.74618443, -0.60317949, -0.40990436,  0.62341661,  0.67727389,  1.0000000,  0.7132144,  0.52923596,  0.7673822,  0.9111478,  0.7487607],
                        [-0.7471112, -0.64507795, -0.63649530, -0.14117556,  0.64542570,  0.42422234,  0.7132144,  1.0000000,  0.60794410,  0.7467579,  0.7732704,  0.7756198],
                        [-0.5945529, -0.70865725, -0.56958405,  0.25538496,  0.59510653, -0.06776682,  0.5292360,  0.6079441,  1.00000000,  0.7868302,  0.5543211,  0.7530874],
                        [-0.8234182, -0.78454102, -0.71446119, -0.02528160,  0.67128060,  0.33662668,  0.7673822,  0.7467579,  0.78683019,  1.0000000,  0.8117621,  0.8738338],
                        [-0.9146447, -0.76548454, -0.64551113, -0.34018047,  0.66187971,  0.69558073,  0.9111478,  0.7732704,  0.55432112,  0.8117621,  1.0000000,  0.8443479],
                        [-0.8473980, -0.79834269, -0.75951011, -0.07811010,  0.75970381,  0.38101906,  0.7487607,  0.7756198,  0.75308740,  0.8738338,  0.8443479,  1.0000000]])
refmat = np.array([[1,	0.7754709,	0.72492855,	0.27817942,	0.3,	-0.60375141,	-0.8614806,	-0.7471112,	-0.59455286,	-0.8234182,	-0.9146447,	-0.847398],
                    [0.7754709,	1,	0.63439613,	0.07402797,	0.07,	-0.3498755,	-0.7461844,	-0.645078,	-0.70865725,	-0.784541,	-0.7654845,	-0.7983427],
                    [0.7249286,	0.63439613,	1,	0.06541974,	0.05,	-0.30024636,	-0.6031795,	-0.6364953,	-0.56958405,	-0.7144612,	-0.6455111,	-0.7595101],
                    [0.2781794,	0.07402797,	0.06541974,	1,	0.8,	-0.72253596,	-0.4099044,	-0.1411756,	0.25538496,	-0.0252816,	-0.3401805,	-0.0781101],
                    [0.3,	0.07,	0.05,	0.8,	1,	-0.7,	-0.5,	-0.1,	0.05,	0.1,	-0.3,	0.01],
                    [-0.6037514,	-0.3498755,	-0.30024636,	-0.72253596,	-0.7,	1,	0.6772739,	0.4242223,	-0.06776682,	0.3366267,	0.6955807,	0.3810191],
                    [-0.8614806,	-0.74618443,	-0.60317949,	-0.40990436,	-0.5,	0.67727389,	1,	0.7132144,	0.52923596,	0.7673822,	0.9111478,	0.7487607],
                    [-0.7471112,	-0.64507795,	-0.6364953,	-0.14117556,	-0.1,	0.42422234,	0.7132144,	1,	0.6079441,	0.7467579,	0.7732704,	0.7756198],
                    [-0.5945529,	-0.70865725,	-0.56958405,	0.25538496,	0.05,	-0.06776682,	0.529236,	0.6079441,	1,	0.7868302,	0.5543211,	0.7530874],
                    [-0.8234182,	-0.78454102,	-0.71446119,	-0.0252816,	0.1,	0.33662668,	0.7673822,	0.7467579,	0.78683019,	1,	0.8117621,	0.8738338],
                    [-0.9146447,	-0.76548454,	-0.64551113,	-0.34018047,	-0.3,	0.69558073,	0.9111478,	0.7732704,	0.55432112,	0.8117621,	1,	0.8443479],
                    [-0.847398,	-0.79834269,	-0.75951011,	-0.0781101,	0.01,	0.38101906,	0.7487607,	0.7756198,	0.7530874,	0.8738338,	0.8443479,	1]])
refmat = np.array([[1, 0.7754709, 0.72492855, 0.27817942, 0.3, 0.38, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.847398],
                    [0.7754709, 1, 0.63439613, 0.07402797, 0.07, 0.1, -0.7461844, -0.645078, -0.70865725, -0.784541, -0.7654845, -0.7983427],
                    [0.7249286, 0.63439613, 1, 0.06541974, 0.05, 0.02, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101],
                    [0.2781794, 0.07402797, 0.06541974, 1, 0.8, 0.77, -0.4099044, -0.1411756, 0.25538496, -0.0252816, -0.3401805, -0.0781101],
                    [0.3, 0.07, 0.05, 0.8, 1, 0.66, -0.5, -0.1, 0.05, 0.1, -0.3, 0.01],
                    [0.38, 0.1, 0.02, 0.77, 0.66, 1, -0.3, 0.01, -0.04, 0.2, -0.2, -0.1],
                    [-0.8614806, -0.74618443, -0.60317949, -0.40990436, -0.5, -0.3, 1, 0.7132144, 0.52923596, 0.7673822, 0.9111478, 0.7487607],
                    [-0.7471112, -0.64507795, -0.6364953, -0.14117556, -0.1, 0.01, 0.7132144, 1, 0.6079441, 0.7467579, 0.7732704, 0.7756198],
                    [-0.5945529, -0.70865725, -0.56958405, 0.25538496, 0.05, -0.04, 0.529236, 0.6079441, 1, 0.7868302, 0.5543211, 0.7530874],
                    [-0.8234182, -0.78454102, -0.71446119, -0.0252816, 0.1, 0.2, 0.7673822, 0.7467579, 0.78683019, 1, 0.8117621, 0.8738338],
                    [-0.9146447, -0.76548454, -0.64551113, -0.34018047, -0.3, -0.2, 0.9111478, 0.7732704, 0.55432112, 0.8117621, 1, 0.8443479],
                    [-0.847398, -0.79834269, -0.75951011, -0.0781101, 0.01, -0.1, 0.7487607, 0.7756198, 0.7530874, 0.8738338, 0.8443479, 1]])
plt.imshow( refmat, cmap = "RdBu" )
plt.colorbar()
plt.title( "Reference Matrix" )
plt.show()
# print(emat)
print("CCD leiden: {}".format(leidenalg.calcCCD(refmat, emat)))
print("CCD louvain: {}".format(louvain.calcCCD(refmat, emat)))

def calcDist(r1, r2):
    tmp = r1-r2
    tmp = tmp ** 2
    return(math.sqrt(np.sum(tmp)))
def calcCCDsimple(ref, emat):
    upper_ref = np.triu(ref)
    corr_mat = np.array(scipy.stats.spearmanr(emat.T))[0,:,:]
    upper_corrmat = np.triu(corr_mat)
    return(calcDist(upper_corrmat, upper_ref))

print(calcCCDsimple(refmat, emat))

# %%
import igraph as ig
emat = adata2.X.T[n_genes:n_genes + n_new_genes,:]
subjects = subjects.astype(np.int32)
subjects = subjects.reshape(1, -1)
print(subjects.shape)

print("CCS louvain: {:.2f}".format(louvain.calcCCS(refmat, emat, subjects)))


# %%
_, G2 = cluster_louvain(adata2, emat, refmat, subjects,  partition_type= louvain.ModularityVertexPartition)  # You can adjust the 'resolution' parameter
membership_louvainStock = [int(i) for i in adata2.obs['louvainccd'].to_list()]
# _plot(G2, membership_louvainStock, draw = "kk")

# %%
# Perform Louvain clustering
t0 = time.time()
_, G = cluster_louvain(adata2, emat, refmat, sample_ids= subjects, partition_type= louvain.ccdModularityVertexPartition)  # You can adjust the 'resolution' parameter
membership = [int(i) for i in adata2.obs['louvainccd'].to_list()]
t1 = time.time()
print("time: {}".format(t1-t0))
# _plot(G, membership, draw = "kk")


# %%
# Assuming you have an 'adata' object with Louvain cluster assignments
import warnings
# Suppress specific warnings
# Suppress all warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # Calculate UMAP
    sc.tl.umap(adata2)
    pcs_to_plot = ['1,2', '2,3', '1,3', '1, 4']
    # Plot UMAP with Louvain clusters
    adata2.obs['louvainStock'] = [str(i) for i in membership_louvainStock]
    sc.pl.umap(adata2, color='louvainStock', legend_loc='on data')
    sc.pl.umap(adata2, color='louvainccd', legend_loc='on data')
    sc.pl.umap(adata2, color='Cell_Identity', legend_loc='on data')
    sc.pl.umap(adata2, color = 'subject', legend_loc = 'on data')
    # sc.pl.pca(adata2, color= 'louvainccd' , components = pcs_to_plot, show=True)
    # adata2.obs['louvainStock'] = [str(i) for i in membership_louvainStock]
    # sc.pl.pca(adata2, color = 'louvainStock', components = pcs_to_plot,  show = True)
    # sc.pl.pca(adata2, color='Cell_Identity', components=pcs_to_plot, show=True)
    # sc.pl.pca(adata2, color='timepoints', components=pcs_to_plot, show=True)


# %%
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


clusterA = 9
cluster_A_cells = [i for i in range(len(membership)) if membership[i] == clusterA]
cluster_A_emat = adata2.X.T[n_genes:n_genes + 12,cluster_A_cells]
clusterA_groups = subjects[0,cluster_A_cells] #subject groups
pseudoBulk = sumByGroup(cluster_A_emat, clusterA_groups)
print("Cluster contains {} groups".format(len(np.unique(clusterA_groups))))

# Plot correlation matrix:
corr_mat = np.array(scipy.stats.spearmanr(pseudoBulk.T))[0,:,:]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

heatmap1 = ax1.imshow(corr_mat, cmap = "RdBu",  norm=Normalize(vmin=-1, vmax=1))
ax1.set_title('Cell Population')
fig.colorbar(heatmap1, ax=ax1)

heatmap2 = ax2.imshow(refmat, cmap = "RdBu",  norm=Normalize(vmin=-1, vmax=1))
ax2.set_title('Reference')
fig.colorbar(heatmap2, ax=ax2)
# Display the plot
plt.show()


# %%
clusterA = 0
clusterB = 2
cluster_A_cells = [i for i in range(len(membership)) if membership[i] == clusterA]
print(cluster_A_cells)  

cluster_B_cells = [i for i in range(len(membership)) if membership[i] == clusterB]
print(cluster_B_cells)  
cluster_A_emat = adata2.X.T[n_genes:n_genes + n_new_genes,cluster_A_cells]
cluster_B_emat = adata2.X.T[n_genes:n_genes + n_new_genes,cluster_B_cells]

cluster_A_ccd = louvain.calcCCD(refmat, cluster_A_emat)
cluster_B_ccd = louvain.calcCCD(refmat, cluster_B_emat)

combined_emat = np.concatenate(( cluster_A_emat, cluster_B_emat), axis = 1)
combined_ccd = louvain.calcCCD(refmat, combined_emat)
print("cluster {} ccd {:.3f}".format(clusterA, cluster_A_ccd))
print("cluster {} ccd {:.3f}".format(clusterB, cluster_B_ccd))
print("{} & {} ccd {:.3f}".format(clusterA, clusterB, combined_ccd))


# %%
# from sklearn.metrics import adjusted_rand_score
# # Extract true labels
# true_labels = adata2.obs['Cell_Identity'].values

# # Extract Louvain clusters
# louvain_clusters = adata2.obs['louvain'].astype(str).values

# # Calculate Adjusted Rand Index
# ari = adjusted_rand_score(true_labels, louvain_clusters)

# # Print the ARI
# print("Adjusted Rand Index:", ari)


