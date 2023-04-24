from typing import Optional, Union, Dict, List
import anndata
import pandas as pd
import numpy as np
import scipy
import pickle
import networkx as nx

from .. import _utils as ucdutils
from .._data import metadata

def get_base_celltypes(
    root : Optional[str] = None,
    category : Optional[str] = None,
    as_digraph : bool = False,
    ) -> List[str]:
    """\
    Get UCDBase CellTypes
    
    Return a list of UCDbase celltypes.
    
    Params
    -------
    root
        Optional root cell type, if set, returns all decendants
        of that celltype along the UCD celltype hierarchy.
    category
        If category is set, overrides root behavior and returns
        a list of all celltypes in the subset category. Must be
        either 'primary', 'lines', or 'cancer'.
    as_digraph
        If true, return the underlying networkX graph that can
        be used to visualize result directly when calling a 
        root node. Call root node "cell" to return all nodes
        in graph.
    Returns
    -------
    celltypes : List[str]
        A list of celltypes, either all or subsets.
    """
    
    # If category set, check and return subset
    if category:
        assert category in ('primary','lines','cancer'), "Category must be one of 'primary', 'lines', or 'cancer'."
        return metadata[f"celltypes_{category}"]
    
    elif root:
        # If root is set, get hiearachy
        G = pickle.loads(ucdutils.decompress_b64_gzip(metadata['cellgraph']))
        if not as_digraph:
            return list(map(str.lower, nx.descendants(G, root)))
        else:
            return G.subgraph(set.union(nx.descendants(G, root),{root}))
    else:
        # Otherwise return all
        return metadata["celltypes_all"]
        
def assign_top_celltypes(
    adata : anndata.AnnData,
    key : str = "ucdbase",
    category : Optional[str] = None,
    groupby : Optional[str] = None,
    inplace : bool = True,
    key_added : str = 'pred_celltype',
    knnsmooth_neighbors: Optional[int] = None,
    knnsmooth_cycles : int = 1,
    ) -> Union[List[str], Dict[str, str]]:
    """\
    
    Gets top deconvolution predictions by cluster.
    
    Params
    -------
    adata
        Annotated dataset with deconvolution results stored
    groupby
        Optional variable, if not passed then the top celltype
        is predicted for each individual sample by row.
    category
        Which split category to use, defaults if none.
    key
        Key for deconvolution results, default key is 'ucdbase'
    inplace
        Whether or not to attach the result to the anndata object 'obs'
        as a column
    key_added
        The key to add as a column.
    knnsmooth_neighbors
        Optional smoothing for predictions, uses neighbors graph calculated
        using sc.tl.neighbors. If not none, passed integer referring to number
        of nearest neighbors to consider for smoothing, reccomended 3.
    knnsmooth_cycles
        Number of cycles to repeat smoothing, default 1.
        
    Returns
    -------
    celltypes : Union[List[str], Dict[str, str]]
        If no groupby is passed return a list of strings corresponding
        to the top celltype that can be used for annotation. If a group
        is passed, returns a dicitonary mapping group label to celltype.
    """
    
    # Get results
    preds = ucdutils.read_results(adata, key = key, category = category)
    
    # Optional smoothing
    if knnsmooth_neighbors:
        topk = _get_top_k(adata.obsp['distances'], knnsmooth_neighbors)
        X_sm = preds.values.copy()
        for cycle in range(knnsmooth_cycles):
            for i in range(len(X_sm)):
                X_sm[i] = X_sm[topk[i]].mean(0)
        preds = pd.DataFrame(X_sm, index = preds.index, columns = preds.columns)
        
        # Modify key name for smoothing
        key = f"{key}_sm{knnsmooth_neighbors}"
    
    if not groupby:
        celltypes = list(preds.idxmax(1))
    else:
        celltypes = preds.groupby(adata.obs[groupby]).mean().idxmax(1).to_dict()
    
    if inplace:
        if groupby:
            adata.obs[f"{key_added}_{key}"] = adata.obs[groupby].map(celltypes)
        else:
            adata.obs[f"{key_added}_{key}"] = celltypes
    else:
        return celltypes

    
def _get_top_k(m : scipy.sparse.csr_matrix, k : int = 3) -> np.ndarray:
    """\
    Get topk nearest neighbors from distance matrix
    """
    topk = np.zeros((m.shape[0], k)).astype(np.int32)
    for i in range(m.shape[0]):
        topk[i] = m.indices[m.indptr[i]:m.indptr[i+1]][np.argsort(m.data[m.indptr[i]:m.indptr[i+1]])[0:k]].astype(np.int32)
    return topk


