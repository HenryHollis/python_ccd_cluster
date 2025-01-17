
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import argparse
import json
from scipy.sparse import csr_matrix
from sklearn.metrics.cluster import adjusted_rand_score


import random
import leidenalg
import louvain
import time

random.seed(42)
np.random.seed(42)

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
# def ccdCluster_leiden(
#     adata: AnnData,
#     emat,
#     refmat,
#     resolution: float = 1,
#     *,
#     restrict_to: tuple[str, Sequence[str]] | None = None,
#     random_state: _utils.AnyRandom = 0,
#     key_added: str = "louvainccd",
#     adjacency: sparse.spmatrix | None = None,
#     directed: bool | None = None,
#     use_weights: bool = True,
#     n_iterations: int = -1,
#     partition_type: type[MutableVertexPartition] | None = None,
#     neighbors_key: str | None = None,
#     obsp: str | None = None,
#     copy: bool = False,
#     flavor: Literal["leidenalg", "ipgraph"] = "leidenalg",
#     **clustering_args,
# ) -> AnnData | None:
#     """\
#     Cluster cells into subgroups [Traag18]_.

#     Cluster cells using the Leiden algorithm [Traag18]_,
#     an improved version of the Louvain algorithm [Blondel08]_.
#     It has been proposed for single-cell analysis by [Levine15]_.

#     This requires having ran :func:`~scanpy.pp.neighbors` or
#     :func:`~scanpy.external.pp.bbknn` first.

#     Parameters
#     ----------
#     adata
#         The annotated data matrix.
#     resolution
#         A parameter value controlling the coarseness of the clustering.
#         Higher values lead to more clusters.
#         Set to `None` if overriding `partition_type`
#         to one that doesn’t accept a `resolution_parameter`.
#     random_state
#         Change the initialization of the optimization.
#     restrict_to
#         Restrict the clustering to the categories within the key for sample
#         annotation, tuple needs to contain `(obs_key, list_of_categories)`.
#     key_added
#         `adata.obs` key under which to add the cluster labels.
#     adjacency
#         Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
#     directed
#         Whether to treat the graph as directed or undirected.
#     use_weights
#         If `True`, edge weights from the graph are used in the computation
#         (placing more emphasis on stronger edges).
#     n_iterations
#         How many iterations of the Leiden clustering algorithm to perform.
#         Positive values above 2 define the total number of iterations to perform,
#         -1 has the algorithm run until it reaches its optimal clustering.
#         2 is faster and the default for underlying packages.
#     partition_type
#         Type of partition to use.
#         Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
#         For the available options, consult the documentation for
#         :func:`~leidenalg.find_partition`.
#     neighbors_key
#         Use neighbors connectivities as adjacency.
#         If not specified, leiden looks .obsp['connectivities'] for connectivities
#         (default storage place for pp.neighbors).
#         If specified, leiden looks
#         .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
#     obsp
#         Use .obsp[obsp] as adjacency. You can't specify both
#         `obsp` and `neighbors_key` at the same time.
#     copy
#         Whether to copy `adata` or modify it inplace.
#     flavor
#         Which package's implementation to use.
#     **clustering_args
#         Any further arguments to pass to :func:`~leidenalg.find_partition` (which in turn passes arguments to the `partition_type`)
#         or :meth:`igraph.Graph.community_leiden` from `igraph`.

#     Returns
#     -------
#     Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

#     `adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
#         Array of dim (number of samples) that stores the subgroup id
#         (``'0'``, ``'1'``, ...) for each cell.

#     `adata.uns['leiden']['params']` : :class:`dict`
#         A dict with the values for the parameters `resolution`, `random_state`,
#         and `n_iterations`.
#     """
#     if flavor not in {"igraph", "leidenalg"}:
#         raise ValueError(
#             f"flavor must be either 'igraph' or 'leidenalg', but '{flavor}' was passed"
#         )
#     igraph_spec = importlib.util.find_spec("igraph")
#     if igraph_spec is None:
#         raise ImportError(
#             "Please install the igraph package: `conda install -c conda-forge igraph` or `pip3 install igraph`."
#         )
#     if flavor == "igraph":
#         if directed:
#             raise ValueError(
#                 "Cannot use igraph's leiden implemntation with a directed graph."
#             )
#         if partition_type is not None:
#             raise ValueError(
#                 "Do not pass in partition_type argument when using igraph."
#             )
#     else:
#         try:
#             import leidenalg

#             msg = 'Use of leidenalg is discouraged and will be deprecated in the future.  Please use `flavor="igraph"` `n_iterations=2` to achieve similar results.  `directed` must also be `False` to work with `igraph`\'s implementation.'
#             warnings.warn(msg, FutureWarning)
#         except ImportError:
#             raise ImportError(
#                 "Please install the leiden algorithm: `conda install -c conda-forge leidenalg` or `pip3 install leidenalg`."
#             )
#     clustering_args = dict(clustering_args)

#     start = logg.info("running Leiden ccd clustering")
#     adata = adata.copy() if copy else adata
#     # are we clustering a user-provided graph or the default AnnData one?
#     if adjacency is None:
#         adjacency = _utils._choose_graph(adata, obsp, neighbors_key)
#     if restrict_to is not None:
#         restrict_key, restrict_categories = restrict_to
#         adjacency, restrict_indices = restrict_adjacency(
#             adata,
#             restrict_key,
#             restrict_categories=restrict_categories,
#             adjacency=adjacency,
#         )
#     # Prepare find_partition arguments as a dictionary,
#     # appending to whatever the user provided. It needs to be this way
#     # as this allows for the accounting of a None resolution
#     # (in the case of a partition variant that doesn't take it on input)
#     clustering_args["n_iterations"] = n_iterations
#     #if resolution is not None:
#         #clustering_args["resolution_parameter"] = resolution
#     if flavor == "leidenalg":
#         directed = True if directed is None else directed
#         g = _utils.get_igraph_from_adjacency(adjacency, directed=directed)
#         if partition_type is None:
#             partition_type = leidenalg.ccdModularityVertexPartition
#         if use_weights:
#             clustering_args["weights"] = np.array(g.es["weight"]).astype(np.float64)
#         clustering_args["seed"] = random_state
#         part = leidenalg.find_partition(g, partition_type, emat, refmat, **clustering_args)
#     else:
#         g = _utils.get_igraph_from_adjacency(adjacency, directed=False)
#         if use_weights:
#             clustering_args["weights"] = "weight"
#         clustering_args.setdefault("objective_function", "modularity")
#         with _utils.set_igraph_random_state(random_state):
#             part = g.community_leiden(**clustering_args)
#     # store output into adata.obs
#     groups = np.array(part.membership)
#     if restrict_to is not None:
#         if key_added == "leidenccd":
#             key_added += "_R"
#         groups = rename_groups(
#             adata,
#             key_added=key_added,
#             restrict_key=restrict_key,
#             restrict_categories=restrict_categories,
#             restrict_indices=restrict_indices,
#             groups=groups,
#         )
#     adata.obs[key_added] = pd.Categorical(
#         values=groups.astype("U"),
#         categories=natsorted(map(str, np.unique(groups))),
#     )
#     # store information on the clustering parameters
#     adata.uns["louvainccd"] = {}
#     adata.uns["leidenccd"]["params"] = dict(
#         resolution=resolution,
#         random_state=random_state,
#         n_iterations=n_iterations,
#     )
#     logg.info(
#         "    finished",
#         time=start,
#         deep=(
#             f"found {len(np.unique(groups))} clusters and added\n"
#             f"    {key_added!r}, the cluster labels (adata.obs, categorical)"
#         ),
#     )
#     return (adata, g )

def cluster_louvain(
    adata: AnnData,
    emat,
    refmat,
    resolution: float | None = None,
    *,
    random_state: _utils.AnyRandom = 0,
    restrict_to: tuple[str, Sequence[str]] | None = None,
    key_added: str = "louvainccd",
    adjacency: spmatrix | None = None,
    flavor: Literal["vtraag", "igraph", "rapids"] = "vtraag",
    directed: bool = True,
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
                g,
                partition_type, emat, refmat,seed = 40
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


def process_data(emat_path, ref_path, graph_path):
    # Load matrices (assuming CSV format)
    reference = pd.read_csv(ref_path, index_col=0)  # Adjust as needed
    emat = pd.read_csv(emat_path, index_col=0)  # Adjust as needed

    # Example processing (you can replace this with your actual processing logic)
    print(f"Expression matrix shape: {emat.shape}")
    print(f"Reference matrix shape: {reference.shape}")

    # Load the neighbors graph
    with open(graph_path, "r") as f:
        graph = json.load(f)

    # Convert matrices back to sparse format
    connectivities = csr_matrix(graph["connectivities_key"])
    distances = csr_matrix(graph["distances_key"])
    params = graph["params"]

    # Add neighbors graph to an AnnData object
    adj = {
        "connectivities_key": connectivities,
        "distances_key": distances,
        "params": params,
    }
    
    # Return some result
    return emat, reference, adj

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process single-cell data with matrices.")
    parser.add_argument('emat_path', type=str, help="Path to the AnnData file")
    parser.add_argument('ref_path', type=str, help="Path to the expression matrix file")
    parser.add_argument('graph_path', type=str, help="Path to the JSON file containing connectivities")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with the file paths
    emat, refmat, adj = process_data(args.emat_path, args.ref_path, args.graph_path)

    print("Expression matrix shape:", emat.shape)
    print("Ref matrix shape:", refmat.shape)

    g = _utils.get_igraph_from_adjacency(adj["connectivities_key"], directed=False)
    partition = louvain.find_partition(
                    g,
                    louvain.ccdModularityVertexPartition, emat, refmat
                )

    print(partition._membership)
