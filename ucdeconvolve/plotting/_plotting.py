####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

"""

======================
Deconvolve Plotting
======================

Description: Some basic plotting utilities for UC deconvolve to visualize deconvolution results.
Note: Requires 'scanpy' library to be installed as these functions are mainly designed to sit on
top of 'scanpy.plotting' module.

"""

import anndata

from .. import tools as ucdtools
from .. import _data as ucddata
from typing import Union, Optional, Tuple, List, Dict
import scanpy as sc
import seaborn as sns

def deconvolve(
    adata : anndata.AnnData,
    basis : str,
    color : Optional[Union[str, List[str]]] = None,
    key : str = ucddata.metadata['default_key'],
    category : Optional[str] = None,
    **kwargs
) -> Optional:
    """\
    Plot Deconvolution
    
    Wrapper for scanpy function 'sc.pl.embedding'
    to help plot deconvolution results. Follows 
    the parameter conventions of its wrapped function
    with some exceptions noted below.
    
    Functions to read the results from the deconvolution
    run given by key, subset to category and then appends
    them to the 'adata.obs' dataframe of a copy of the
    passed adata object, allowing standard plotting
    module to the visualize the results.
    
    Params
    ------
    adata
        anndata object to plot
    basis
        The embedding to plot using, for example
        'X_pca' or 'X_umap' if calculated and present.
    color
        Refers to the cell type we want to plot contained
        within the category of split and result specificed by
        key. Can be one or more.
    key
        location of data in obsm and uns to plot containing
        numerical data and headers, respectively.
    category
        if the data results are split, indicate which split
        to use for plotting. defaults to 'all' assuming that
        we did not split the output. valid categories
        are 'all', 'primary', 'cell_lines', and 'cancer'.
    kwargs
        attributes to pass along to 'sc.pl.embedding', see
        documentation for details.
    
    Returns
    -------
    Plot(s)
    
    """
    
    # Read the prediciton results we want
    predictions = ucdtools.read_results(adata, category, key)
    
    # Create a copy of adata
    adata = adata.copy()
    
    # Attach predictions so we can use them
    adata.obs = adata.obs.merge(predictions, how = 'left', left_index = True, right_index = True)

    # If the request is spatial use a different function
    if basis == 'spatial':
        return sc.pl.spatial(adata, color = color, **kwargs)
    
    # Return the plot
    return sc.pl.embedding(adata, color = color, basis = basis, **kwargs)
