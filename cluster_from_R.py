
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import argparse
import json
from scipy.sparse import csr_matrix
from sklearn.metrics.cluster import adjusted_rand_score
from scipy.io import mmread

import tempfile

import random
# import leidenalg
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

def process_data(emat_path, ref_path, graph_path, alg):
    # Load matrices (assuming CSV format)
    reference = pd.read_csv(ref_path, index_col=0)  # Adjust as needed
    emat = pd.read_csv(emat_path, index_col=0)  # Adjust as needed

    # Example processing (you can replace this with your actual processing logic)
    print(f"Expression matrix shape: {emat.shape}")
    print(f"Reference matrix shape: {reference.shape}")

    # Load the neighbors graph
    # Read sparse matrices
    print("reading graph")
    connectivities = mmread(graph_path).tocsr()
    # print(f"Number of nonzero elements: {connectivities.nnz}")
    if(alg == "louvain"):
        print("using original Louvain clustering")
        algorithm = louvain.ModularityVertexPartition
    elif (alg == "ccd"):
        print("using custom CCD louvain clustering")
        algorithm = louvain.ccdModularityVertexPartition
    else:
        print("No valid clustering alg has been selected, assuming Louvain")
        algorithm = louvain.ModularityVertexPartition
    # Convert matrices back to sparse format
    # connectivities = csr_matrix(graph["connectivities_key"])
    # distances = csr_matrix(graph["distances_key"])
    # params = graph["params"]

    # Add neighbors graph to an AnnData object
    # adj = {
    #     "connectivities_key": connectivities,
    #     "distances_key": distances,
    #     "params": params,
    # }
    
    # Return some result
    return emat, reference, connectivities, algorithm

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process single-cell data with matrices.")
    parser.add_argument('emat_path', type=str, help="Path to the AnnData file")
    parser.add_argument('ref_path', type=str, help="Path to the expression matrix file")
    parser.add_argument('graph_path', type=str, help="Path to the JSON file containing connectivities")
    parser.add_argument('alg', type = str, help = "clustering algorithm used, 'louvain' or 'ccd'")
    # Parse arguments
    args = parser.parse_args()

    # Call the function with the file paths
    emat, refmat, connectivities, alg = process_data(args.emat_path, args.ref_path, args.graph_path, args.alg)

    print("Expression matrix shape:", emat.shape)
    print("Ref matrix shape:", refmat.shape)

    g = _utils.get_igraph_from_adjacency(connectivities, directed=False)
    print(type(g))
    part = louvain.find_partition(
                    g,
                    alg, emat, refmat
                )

    # print(partition._membership)
    # Save result to a temporary file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w") as tmpfile:
        np.savetxt(tmpfile.name, part._membership, delimiter=",", header="ClusterID", comments="")
        print(tmpfile.name)  # Print the file path so R can capture it