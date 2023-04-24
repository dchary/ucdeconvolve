####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

"""

======================
Deconvolve Base Plotting
======================

Description: Some basic plotting utilities for UC deconvolve to visualize base deconvolution results.
Note: Requires 'scanpy' library to be installed as these functions are mainly designed to sit on
top of 'scanpy.plotting' module.

"""

import anndata

from .. import _compat as ucdtools
from .. import _data as ucddata

from typing import Union, Optional, Tuple, List, Dict
import scanpy as sc
import seaborn as sns

def base_clustermap(
    adata : anndata.AnnData,
    groupby : str = 'leiden',
    category : Optional[str] = None,
    key : str = "ucdbase",
    n_top_celltypes : int = 30,
    max_filter : float = 0.1,
    **kwargs,
) -> Optional:
    """\
    
    Plot Clustered heatmap of top celltype predictions
    grouped by a column in 'adata.obs'
    
    Params
    -------
    adata
        The annotated dataset with deconvolution data
    groupby
        What column in 'adata.obs' to group celltype
        predictions by (i.e. 'leiden').
    category
        Which category of prediction data to use if
        split, or all of not split.
    key
        Key for deconvolution results, default is 'ucdbase'
    n_top_celltypes
        Number of top celltypes per category to take 
        and plot. Smaller means only the most common types.
    kwargs
        Keyword attributes for clustermap. See
        seaborn.clustermap for details.
    Returns
    -------
    A clustermap
    
    """
    
    # Get deconvolution predictions
    predictions = ucdtools.read_results(adata, category, key)
        
    # Get top predictions
    top_celltypes = predictions.mean(0).sort_values(
                        ascending = False).head(n_top_celltypes).index.tolist()
    
    # Subset predictions to only top celltypes
    predictions = predictions[top_celltypes]
    
    # Group predictions by column in 'adata.obs'
    predictions = predictions.groupby(adata.obs[groupby]).mean().T            
    predictions = predictions[predictions.max(1) > max_filter]

    # Append defaults to kwargs if not overwritten
    if 'cmap' not in kwargs.keys(): kwargs['cmap'] = 'viridis'
    if 'row_cluster' not in kwargs.keys(): kwargs['row_cluster'] = False

    # Generate the plot
    return sns.clustermap(predictions, **kwargs)
