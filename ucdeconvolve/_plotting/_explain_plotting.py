####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Optional, Union, Tuple, Dict, List
import anndata
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import numpy as np
import pandas as pd
import seaborn as sns
import textwrap3

from .._data import metadata as ucdmetadata
from .. import _utils as ucdutils

def explain_clustermap(
    adata : anndata.AnnData,
    key : Union[str, List[str]] = "ucdexplain",
    n_top_genes : int = 64,
    **kwargs
    ) -> Optional[plt.Axes]:
    """\
    
    Plot Explanation Results as Clustermap
    
    Plot Clustered heatmap of top feature attribution predictions
    grouped by the celltypes passed to the ``ucd.tl.explain`` function.
    
    Params
    -------
    adata
        The annotated dataset with deconvolution data
    key
        Key for deconvolution results, default is 'ucdexplain'.
    n_top_genes
        Number of top feature attributes (genes) per celltype
    kwargs
        Keyword attributes for clustermap. See
        seaborn.clustermap for details.

    Returns
    -------
    A clustermap
    
    """

    # First read results out
    df = ucdutils.read_results(
        adata = adata, 
        key = key,
        explain_n_top_genes = n_top_genes)
    
    # Groupby celltype
    df = df.groupby('celltype').mean()
    
    # Set defaults for kwargs if not set
    defaults = {"figsize" : (12,3), "cmap" : "RdBu_r", 
        "z_score" : 1, "vmin" : 0, "vmax" : 2}
    for k,v in defaults.items():
        if k not in kwargs: kwargs[k] = v
    return sns.clustermap(data = df, **kwargs)
    

def explain_boxplot(
    adata : anndata.AnnData,
    key : str = "ucdexplain",
    celltypes : Optional[Union[str, List[str]]] = None,
    n_top_genes : int = 16,
    ncols : int = 5,
    figsize : Tuple[int,int] = (3,3),
    dpi : int = 150,
    titlewidth : int = 24,
    barcolor : str = "lightblue",
    ax : Optional[plt.Axes] = None,
    return_fig : bool = False
    ) -> Optional[plt.Axes]:
    """\
    
    Plot Boxplots of Feature Attributions By Gene
    
    Params
    -------
    adata
        Annotated dataset with ucdexplain results.
    key
        UCDExplain results key, default is 'ucdexplain'
    celltypes
        The celltypes from the given run to plot. if none
        then plots all.
    n_top_genes
        Number of top attribution genes to plot.
    ncols
        Number of columns to plot for multiple celltypes before
        creating a new row
    figsize
        Size of individual subplot figure
    dpi
        Pixel density of plot
    titlewidth
        Width of subplot title before newline
    barcolor
        Color of bars
    ax
        Optional axes to plot on.
    return_fig
        Return figure or not
        
    Returns
    -------
    
    fig : plt.Figure
        Figure with underlying subplots
    
    """
    
    # Get celltype and mask keys from results key
    key_celltype = f"{key}_celltype"
    key_mask = f"{key}_mask"

    # Subset anndata based on mask
    adata_mask = adata[adata.obs[key_mask]]

    # Get unique celltypes if manually not passed
    if not celltypes:
        celltypes = tuple(adata_mask.obs[key_celltype].unique())
    else:
        # Process conditions for celltypes param
        if not celltypes:

            # Get all celltypes
            celltypes = tuple(adata_mask.obs[key_celltype].unique())

        elif isinstance(celltypes, str):

            # If just one celltype passed, then make tuple
            celltypes = (celltypes,)
        
    n_celltypes = len(celltypes)
    
    # Calculate nrows based on n_celltypes and ncols
    nrows = math.ceil(n_celltypes / ncols)
    
    # Calculate figsize
    figsize_full = (figsize[0] * min(ncols, n_celltypes),
                    figsize[1] * nrows)

    # Create subplots
    fig, axes = plt.subplots(nrows = nrows, 
                             ncols = min(ncols, n_celltypes),
                             figsize = figsize_full, dpi = dpi) if not ax else (ax.get_figure(), ax)

    # Flatten axes object
    axes = axes.flatten() if n_celltypes > 1 else [axes]

    # Iterate each celltype
    for ax, celltype in zip(axes, celltypes):

        # Get susbet from masked anndata
        adata_mask_ct = adata_mask[adata_mask.obs[key_celltype].eq(celltype)]

        # Get explanations data 
        expl = adata_mask_ct.obsm[key]

        # Get top gene indices in descending order
        top_gene_indices = np.argsort(np.array(expl.mean(0)).flatten())[::-1][0:n_top_genes]

        # Get top gene names
        top_gene_names = np.array(ucdmetadata['target_genes'])[top_gene_indices]

        # Make dataframe
        df = pd.DataFrame(data = expl[:, top_gene_indices].toarray(),
                index = adata_mask_ct.obs_names,
                columns = top_gene_names)

        # Format for plotting
        df = df.stack()\
                .rename_axis(['cell','gene'])\
                .to_frame("attribution")\
                .reset_index()\
                .drop(columns = ['cell'])

        # Actually plot it
        sns.boxplot(data = df, y = "gene", x = "attribution", color = barcolor,
                   linewidth = 0.5, ax = ax, fliersize = 0)
        sns.stripplot(data = df, y = "gene", x = "attribution", color = "k",
                     size = 2, ax = ax)
        ax.set_title(textwrap3.fill(celltype, width = titlewidth))
        sns.despine(ax = ax, trim = True, offset = 5)

        # Create tight_layout
        fig.tight_layout()
        
        # Return the figure
        if return_fig: return fig
    
    
    