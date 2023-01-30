####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List, Iterable, Dict
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import issparse, spmatrix
import scipy
import logging
from datetime import datetime

from ._data import metadata
from . import _utils as ucdutils
from ._propagator import UCDeconvolvePredPropagator
from ._metadata import __version__ as ucdversion

def postprocess_predictions(
    predictions : np.ndarray,
    split : bool = True,
    normalize_split_totals : bool = True,
    sort : bool = True,
    propagate : bool = True,
    ):
    """\
    Postprocess Deconvolution Results
    
    Apply a number of postprocessing routines
    to the underlying data.
    
    Params
    -------
    predictions
        The raw results from predictive model.
    split
        Whether or not to split underlying data into
        three categories, primary, cancer cell_line.
        Helps with interpretability downstream, default
        is True.
    normalize_split_totals
        If data is split, fractiosn no longer sum to 1. In order
        to get a better sense of cell states across data, it can help
        to normalize the row fractions to 1 for primary cells, while
        adding "primary" and "normal" cell columns to cell lines and 
        cancers, respectively, such that each split sums to 1.
    sort
        Sort columns of results by mean predictions. Default True.
    propagate
        whether or not to perform belief propagation and
        pass predictions up a cell-type heiarchy. helpful
        in interpreting overall deconvolution results.
        default is True.
        
    Returns
    -------
    predictions_pp : Union[np.ndarray, Dict[str, np.ndarray]
        Predictions either as a normalize numpy array or a set of them
        if split.
        
    """
    
    # Get logging object
    ucdlogger = logging.getLogger("UCD")
    
    # Convert predictions to float32
    predictions = predictions.astype(np.float32)
    
    # Invalid predictions are those where -1 is returned
    # Save a record here to set to NaN at the end
    invalid_responses = predictions.min(1) == -1
    
    # If more than one invalid, report it to the user
    if invalid_responses.sum() > 0:
        pct_invalid = (1.0 * invalid_responses.sum()) / invalid_responses.shape[0]
        ucdlogger.warning(f"{pct_invalid:.0%} of returned predictions invalid. Consider re-running.")

    # Report status
    ucdlogger.debug("Splitting predictions.") if split else ucdlogger.debug("Skipping split predictions.")
    
    # Either split predictions using celltype headers or keep all together
    predictions = split_predictions(predictions) if split else {'celltypes_all' : predictions}
    
    # Correct predictions if necessary
    predictions = correct_categories(predictions, is_split = split)
    
    # Normalize rows
    if normalize_split_totals and split:
        
        # Add a column to cancer split with normal probability
        predictions['celltypes_cancer']['normal'] = 1 - predictions['celltypes_cancer'].sum(1)
    
        # Add a column to lines split with primary probability
        predictions['celltypes_lines']['primary'] = 1 - predictions['celltypes_lines'].sum(1)
        
        # Normalize primary to sum to 1 across samples
        predictions['celltypes_primary'] = predictions['celltypes_primary'].divide(predictions['celltypes_primary'].sum(1), 0)
    
    # Preload belief propagator
    propagator = UCDeconvolvePredPropagator(sort_output = False) if propagate else None
    
    # Iterate each key in predictions and peform the following operations
    for key in predictions.keys():
    
        # Belief propagation
        if propagate:
            predictions[key] = propagator.fit_transform(predictions[key])

        # Sort columns by most common predicted cell type
        if sort:
            
            # Get sort order
            order = predictions[key].mean(0).sort_values(ascending = False).index.tolist()

            # Reindex based on sorted order
            predictions[key] = predictions[key][order]
            
        # If a prediction was invalid (i.e. server time-out, etc...) then we need to save that
        # prediction as NaN here
        predictions[key].loc[invalid_responses, :] = np.nan
                
    # Return the preprocessed prediction data
    return predictions


def correct_categories(
    predictions : Dict[str, np.ndarray],
    is_split : bool = False,
) -> Dict[str, pd.DataFrame]:
    """\
    Quick-fix correction for misannotated headers.
    Will remove in later build from header directly.
    
    Params
    -------
    predictions 
        Set of predictions split or not
    is_split
        Flag if split occured.
    Returns
    -------
    predictions
        Set of corrected split predictions as dictionaries.
        
    """
    
    # If we didn't split no need for correcton just return a dataframe
    if not is_split:
        
        # Create an all attribute
        all = pd.DataFrame(predictions['celltypes_all'], columns = metadata['celltypes_all'])
        
        # Drop the empty column
        all = all.drop(columns = [''])
        
        return { 'celltypes_all' : all }
    
    # Create dataframes from split data
    cancer = pd.DataFrame(predictions['celltypes_cancer'], columns = metadata['celltypes_cancer'])
    lines = pd.DataFrame(predictions['celltypes_lines'], columns = metadata['celltypes_lines'])
    primary = pd.DataFrame(predictions['celltypes_primary'], columns = metadata['celltypes_primary'])
    
    # Swap the misannotated cell line
    lines['hcc1428'] = primary['hcc1428']
    primary = primary.drop(columns = ['hcc1428'])
    
    return { 'celltypes_cancer' : cancer, 'celltypes_lines' : lines, 'celltypes_primary' : primary }


def split_predictions(
    predictions : np.ndarray
) -> Dict[str, np.ndarray]:
    """\
    
    Split raw predictions by predefined segments
    
    Params
    -------
    predictions
        The raw predictions array
    
    Returns
    -------
    splitpredictions : dict
        A dictionary holding all split predictions
    """
    
    # Buffer to hold all split predictions
    splitpredictions = dict()

    # Iterate each segment
    for targetkey in list(filter(lambda k : k.endswith('segments'), metadata.keys())):

        # Get actual target indices
        target = metadata[targetkey]

        # Get header information with each cell type
        header = metadata[targetkey.rsplit("_",1)[0]]

        # Create a buffer to hold the split predictions
        splitprediction = np.zeros((predictions.shape[0], len(header)))

        # Iterate each column index
        for i in range(predictions.shape[1]):

            # Ignore -1 as it doesn't apply in this split
            if target[i] == -1: continue

            # Sum up predictions for target columns in split
            splitprediction[:, target[i]] += predictions[:, i]

        # Assign to dict buffer holding all three splits
        splitpredictions[targetkey.rsplit("_",1)[0]] = splitprediction
    
    # Return buffer
    return splitpredictions


def append_predictions_to_anndata(
    predictions : Dict[str, pd.DataFrame],
    adata : anndata.AnnData,
    key_added : str = metadata['default_key'],
    additional_runinfo : Optional[dict] = dict()
):
    """\
    Append deconvolution results to anndata
    
    Params
    -------
    predictions
        Dictionary of prediction(s) as dataframes.
    adata
        Annotated dataset object to append predictions to
    key_added
        Key to use as the results key, default to 'ucd_results'
    additional_runinfo
        Optional dictionary one can pass to add specific runinfo params
    Returns
    -------
    adata
        Annotated dataset with appended results
        
    """
    
    # Create a dict to hold results metadata
    adata.uns[key_added] = dict()

    # Create a dict to hold headers and run information 
    headers = dict()
    runinfo = dict()

    # Add current version and timestamp of run
    runinfo['ucdversion'] = ucdversion
    runinfo['timestamp'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    
    # Append any additional run informations
    runinfo.update(additional_runinfo)
    
    # Append each prediction to anndata object
    for key, vals in predictions.items():
        
        # Get stem for each name (i.e. cancer, lines, primary)
        key_stem = key.split("celltypes_",1)[-1]
        
        # Put each prediction raw numpy array in the 'adata.obsm' dict
        adata.obsm[f"{key_added}_{key_stem}"] = predictions[key].values
        
        # Put header information into headers dict
        headers[key_stem] = predictions[key].columns.str.lower().tolist()
        
    # Create a dicts to hold headers and run information 
    adata.uns[key_added]['headers'] = headers
    adata.uns[key_added]['runinfo'] = runinfo
    
    # Return appended adata
    return adata