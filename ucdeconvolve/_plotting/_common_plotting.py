####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

# Common Plotting Utilities Shared Across Functions
from typing import Union, Optional, Tuple, List, Dict

from .. import _utils as ucdutils
from .. import _data as ucddata
from ._colormaps import CM
from ._plotutils import append_custom_cbar

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl

import scanpy as sc
import seaborn as sns
import anndata
import warnings
import numpy as np
import math

def spatial(
    adata : anndata.AnnData,
    color : Optional[Union[str, List[str]]] = None,
    key : str = "ucdbase",
    category : Optional[str] = None,
    labels : Optional[List[str]] = None,
    colormaps : Optional[List[ListedColormap]] = None,
    cbar_nrows : int = 4,
    title : str = "",
    **kwargs
    ) -> Optional[object]:
    """\
    Plot Spatial
    
    Wrapper for scanpy function 'sc.pl.spatial'
    to help plot deconvolution results on spatial data.
    Follows parameter conventions of wrapped function
    with some exceptions.
    
    Functions to read the results from the deconvolution
    run given by key, subset to category and then appends
    them to the 'adata.obs' dataframe of a copy of the
    passed adata object, allowing standard plotting
    module to the visualize the results.
    
    Params
    ------
    adata
        anndata object to plot
    color
        Refers to the cell type(s) we want to plot contained
        within the category of split and result specificed by
        key. Can be one or more. If more than one string is
        passed we try to plot an overlapped plot.
    key
        location of data in obsm and uns to plot containing
        numerical data and headers, respectively. Can be
        either 'ucdbase' or 'ucdselect'.
    category
        if the data results are split, indicate which split
        to use for plotting. defaults to 'all' assuming that
        we did not split the output. valid categories
        are 'all', 'primary', 'cell_lines', and 'cancer'.
    labels
        Labels for each color being plotted when using the
        overlapping colormap spatial function.
    colormaps
        Optional custom colormaps to use for each color.
    cbar_nrows
        Number of rows to spread cbars across, default is 3.
    kwargs
        attributes to pass along to 'sc.pl.spatial', see
        documentation for details.
     
    Returns
    -------
    Plot(s)
    
    """
    
    # Catch warnings
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=Warning)
    
        # Read the prediciton results we want
        predictions = ucdutils.read_results(adata, key, category)

        # Create a copy of adata
        adata = adata.copy()

        # Attach predictions so we can use them
        adata.obs[predictions.columns] = predictions.values
        
        # Now check if color is a list of more than one strings, if it is, 
        # we attempt to plot a single plot with overlayed colormaps
        if (isinstance(color, list) and len(color) > 1):
            
            # Check that we have enough colors to plot the colormap
            assert len(color) <= len(CM), "Not enough unique colors to plot"
            
            cbar_ncols = math.ceil(len(color) / cbar_nrows)
            
            # Create a new figure / ax instance using rc params
            fig = plt.figure(constrained_layout=False)
            
            gs = fig.add_gridspec(cbar_nrows, 1 + cbar_ncols, 
                height_ratios = [1] * cbar_nrows, 
                width_ratios = [1] + [0.025] * cbar_ncols)
            
            # Get ax for the spatial plot
            ax = fig.add_subplot(gs[:, 0])

            # Remove show if given and set to false
            show = kwargs.pop("show", None)
            
            # Remove other params that we may need if they are passed
            alpha_img = kwargs.pop("alpha_img", None)
            cmap = kwargs.pop("cmap", None)
            alpha_img = kwargs.pop("alpha_img", 1.0)
            
            
            # Get oclormaps
            colormaps = colormaps if colormaps else list(CM)
            
            # Iterate each celltype 
            for i, (cmap, celltype) in enumerate(zip(colormaps, color)):
                
                sc.pl.spatial(adata, color = celltype, ax = ax, show = False, 
                              alpha_img = alpha_img if i == 0 else 0.0, 
                              cmap = cmap.value, vmax = "p98", vmin = 0.0, title = title, **kwargs)
            
            coords = []
            # Create coords
            for j in range(1, cbar_ncols + 1):
                for i in range(0, cbar_nrows):
                    coords.append((i,j))
                    
            for i, (collection, cmap, name, coord) in enumerate(zip(ax.collections, colormaps[0:len(color)], color, coords)):
                collection.colorbar.ax.set_visible(False)
                
                cax = fig.add_subplot(gs[coord[0], coord[1]])
                
                ylims = collection.colorbar.ax.get_ylim()
    
                pc = mpl.colorbar.ColorbarBase(cax, 
                    orientation="vertical", 
                    cmap=cmap.value,
                    norm=mpl.colors.Normalize(ylims[0],ylims[1]),
                    extend="both")
                
                                
                cax.yaxis.set_ticks_position('left')
                cax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
                cax.yaxis.set_tick_params(which='both', labelsize='x-small')
                if labels: cax.set_ylabel(labels[i], fontsize = 'xx-small')
                
            plt.subplots_adjust(wspace = 0.35)
            plt.tight_layout()
            
            return ax

        return sc.pl.spatial(adata, color = color, **kwargs)
    
    
def embedding(
    adata : anndata.AnnData,
    basis : str = "X_umap",
    color : Optional[Union[str, List[str]]] = None,
    key : str = "ucdbase",
    category : Optional[str] = None,
    **kwargs
    ) -> Optional[object]:
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
        numerical data and headers, respectively. Can be
        either 'ucdbase' or 'ucdselect'.
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
    predictions = ucdutils.read_results(adata, key, category)
        
    # Create a copy of adata
    adata = adata.copy()
    
    # Attach predictions so we can use them
    adata.obs[predictions.columns] = predictions.values
    
    # Return the plot
    return sc.pl.embedding(adata, color = color, basis = basis, **kwargs)