####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List, Dict, Iterable
import pandas as pd
import numpy as np
import anndata
import mudata
import logging
from more_itertools import collapse
from scipy.sparse import coo_matrix

from .. import _utils as ucdutils
from .._settings import settings
from .._api import _api as ucdapi
from .._api import _apiutils as ucdapiutils
from .._api import _statewatchers as ucdstatewatchers
from .._data import metadata as ucdmetadata

def deconvolve_explain(
    data : anndata.AnnData,
    celltypes : Union[str, List[str], Dict[Union[int, str], str]],
    groupby : Optional[str] = None,
    group_n : int = 16,
    group_frac : Optional[float] = None,
    token : Optional[str] = None,
    return_results : bool = False,
    key_added : str = 'ucdexplain',
    use_raw : Union[bool, Tuple[bool, bool]] = True,
    verbosity : Optional[int] = None,
    ) -> Optional[anndata.AnnData]:
    """\
    
    UniCell Deconvolve: Explain
    
    Explains cell type fraction prediction for provided transcriptomic data.
    
    Params
    -------
    data
        Transcriptomic data (obs x genes) to predict cell type fractions. Can
        be either a dataframe or annotated dataset. Note that in any case
        data will be converted to an annotated datset object before proceeding.
    celltypes
        Name of cell type(s) to get explanations for. If a single string is passed, this
        celltype is used for all samples. If a list of strings is passed, the list must
        be the same length as the dataset and each entry corresponds to which celltype
        to get explanatons for in the whole dataset. If a dictionary is passed, the key
        should corresponding to an 'adata.obs' column defined by 'groupby', alliowing for
        celltype expalantions to be generated specific to different clusters or conditions.
    groupby
        Groupby key in 'adata.obs' to arrange search for celltypes. If celltypes is given
        as a dict, this must be defined.
    group_n
        The number of samples to subsample from each group for explanations, as this is an
        expensive operation and most cells in a cluster will yield similar results.
    token
        UCDeconvolve API access token. If None, defaults to settings parameter.
    return_results
        Whether or not to return the predictions dict from the function,
        default to false as all data is written to anndata object either
        passed in, or created when passing in a dataframe, which will
        in that case be returned by default. Also returns the underlying
        anndata if it is a view as copying can destroy context internally.
    use_raw
        Use counts in 'adata.raw'. Default True, as by convention in 
        single cell analysis, log1p scaled counts before HVG filter are kept here.
    verbosity
        Level of verbosity for function information. Default is taken from package,
        set to 'logging.DEBUG' for more detailed information.
        
    Returns
    -------
    adata_mixture_orig : anndata.AnnData
        Results appended to anndata object if return_results or if original input was dataframe.
    
    """

    # Get verbosity
    verbosity = verbosity if verbosity else settings.verbosity
    
    # Get token
    token = token if token else settings.token
    
    # Convert dataframe to annotated dataset if one is passed
    if isinstance(data, pd.DataFrame):
        data = anndata.AnnData(data)
        orig_df = True
    else:
        orig_df = False

    # Wrap the entire lifecycle in a logging context
    with ucdutils.log(verbosity) as ucdlogger:
        
        # Wrap function in a timer
        with ucdutils.timer(ucdlogger, "Starting UCDeconvolveEXPLAIN Run.", "Run Complete", logging.INFO):
            
            # Wrap the script in a try/except to catch keyboard interrupts
            try:
                
                # Ensure inputs are correct type
                assert isinstance(data, (anndata.AnnData)), "data must either be a annotated dataset"

                # Assign vars
                adata_mixture = data.copy()

                # Keep references to originals
                adata_mixture_orig = data
                
                # Drop all obs data from mixture
                columns_to_drop = list(set(adata_mixture.obs.columns) - (set([groupby]) if groupby else set()))
                adata_mixture.obs = adata_mixture.obs.drop(columns = columns_to_drop)

                # Build explanation index column "explain_idx"
                adata_mixture = _build_celltypes_index(adata_mixture, 
                                    celltypes, groupby, group_n, group_frac)

                # Use raw from adata if specified, keep original 'adata' as reference to append to and return
                adata_mixture = ucdutils.try_get_raw(adata_mixture, use_raw = use_raw)

                # Ensure scaled data is not provided
                assert adata_mixture.X.min() >= 0, "Detected expression < 0, make sure counts are provided"
                
                # Standardize gene names
                adata_mixture.var_names = ucdutils.match_to_gene(adata_mixture.var_names.tolist())

                # Drop NA's and make copy for implicit view
                adata_mixture = adata_mixture[:,
                    ~adata_mixture.var_names.isna() & ~adata_mixture.var_names.duplicated()].copy()

                # Format inputs for streaming to UCD API
                adata_mixture = ucdutils.get_preprocessed_anndata(adata_mixture, desc = "Dataset", 
                                                                  progress = verbosity <= logging.INFO)
            
                # Package both objects into mudata object
                mdata = mudata.MuData({"mixture" : adata_mixture})

                # Format params 
                params = {"celltypes_key" : "explain_idx" }

                # Get signed URL
                metadata, url = ucdapi.start_run_request("explain", params = params)

                # Upload to URL
                with ucdutils.timer(ucdlogger, "Uploading Data", "Upload Complete", logging.INFO):
                    response = ucdapiutils.upload_mudata_to_cloud_storage_url(mdata, url, metadata)

                # Wait for run submission
                ucdstatewatchers.wait_for_submission(
                    metadata['run_id'], metadata['runtype'],
                    token = token,
                    showprogress = verbosity < logging.INFO)

                # Wait for run completion and results
                results = ucdstatewatchers.wait_for_completion(
                    metadata['run_id'], token = token,
                    showprogress = verbosity < logging.INFO,
                desc = "Waiting For Completion")

                # Download and attach results
                with ucdutils.timer(ucdlogger, "Download Results", "Download Complete", logging.INFO):
                    
                    # Include method
                    params['method'] = 'explain'
                    
                    # We attach results to the potentially subset mixture file first
                    # This behavior isn't consistent in this method between functions
                    # and needs to be updated
                    explanations = ucdapiutils.download_and_attach_results(
                            adata_mixture_orig, results['download_url'], metadata['runtype'],
                            key_added, params)

                    # Get explanations from mixture object
                    #explanations = adata_mixture_orig.obsm[key_added]

                    # Convert to coo
                    explanations = explanations.tocoo()

                    # Use indexer to get reindex positions for each value 
                    # in the subset matrix
                    reindexes = adata_mixture_orig.obs.index.get_indexer(
                        adata_mixture.obs.index)

                    # Expand explanations while maintaining sparsity constraint
                    # to match the original mixture object
                    explanations = coo_matrix(
                        (explanations.data,
                        (reindexes[explanations.row],explanations.col)),
                        shape = (adata_mixture_orig.shape[0], 28867)).tocsr()

                    # Attach to original anndata along with mask
                    adata_mixture_orig.obsm[key_added] = explanations
                    adata_mixture_orig.obs[f"{key_added}_mask"] = \
                        np.isin(np.arange(len(adata_mixture_orig)), reindexes)
        
                    # Create dictionary mapping numerical indexes to celltypes
                    num_to_celltype = dict(enumerate(ucdmetadata['celltypes_all']))

                    # Convert explain index to celltype names
                    celltypenames = adata_mixture.obs["explain_idx"].map(
                            num_to_celltype).values.astype(str)

                    # Include target celltype, ensure length of the string array
                    # is length of longest value in celltypenames
                    explain_celltype = np.array(
                        ["n/a"]*adata_mixture_orig.shape[0]).astype(celltypenames.dtype)

                    np.put_along_axis(explain_celltype, 
                        reindexes, adata_mixture.obs["explain_idx"].map(
                            num_to_celltype).values, 0)

                    # Attach to adata
                    adata_mixture_orig.obs[f"{key_added}_celltype"] = explain_celltype

                # Remove and clean up run
                ucdapi.remove_run(metadata['run_id'], token)

                # Return original data with attached results
                # if return is set or input was dataframe.
                if return_results or orig_df:
                    return adata_mixture_orig

            # If we have a premature keyboard interrupt,
            # we will want to send a kill signal and remove the run
            except KeyboardInterrupt:
                
                # Print newline before message
                ucdlogger.error('\n')
                
                # Notify user we are trying to have a graceful shutdown
                ucdlogger.error("KeyboardInterrupt. Attempting Graceful Shutdown")
                
                # Check if we have metadata object yet, if so we can try to kill the run
                # in case it is running and then remove it from records.
                if metadata:
                    ucdapi.kill_run(metadata['run_id'], token)
                    ucdapi.remove_run(metadata['run_id'], token)                 
            

def _build_celltypes_index(
    adata : anndata.AnnData,
    celltypes : Union[str, List[str], Dict[Union[int, str], str]],
    groupby : Optional[str] = None,
    group_n : int = 16,
    group_frac : Optional[float] = None,
    ) -> anndata.AnnData:
    """
    Build CellTypes Index
    
    Utility for building celltype index column, can accept several
    different formats. See details in main function above.

    
    """
    
    # If celltypes is a string, make list length of adata
    if isinstance(celltypes, str):
        
        # Make sure celltype is in list of available celltypes
        assert celltypes in ucdmetadata['celltypes_all'], "Celltype name mismatch"
        
        # Extend as list to match length of dataset and get index position
        celltypes = [ucdmetadata['celltypes_all'].index(celltypes)] * len(adata)
        
        # Set explain index
        adata.obs["explain_idx"] = celltypes
    
    # If celltypes is a list make sure its same length as adata)
    elif isinstance(celltypes, list):
        
        # Make sure celltype is in our list of available celltypes
        for celltype in celltypes:
            assert celltype in ucdmetadata['celltypes_all'], "Celltype name mismatch"
        
        # Ensure same length as metadata
        assert len(celltypes) == len(adata), "CellTypes must be length of anndata"
        
        # Index each celltype in list
        celltypes = [ucdmetadata['celltypes_all'].index(c) for c in celltypes]
        
        # Set explain index
        adata.obs["explain_idx"] = celltypes
        
    # If celltypes is a dictionary, make sure we also pass a groupby column
    elif isinstance(celltypes, dict):
        assert groupby is not None, "if dict used for celltypes, groupby must be passed"
        
        # If we are using a dict, we may be subsetting our anndata object
        # Get group indices
        group_indices = pd.RangeIndex(len(adata)).groupby(adata.obs[groupby])

        # Get random values for group indices according to size
        group_indices_subset = {
            k : list(np.random.choice(v.values, 
                min(group_n if not group_frac else int(len(v) * group_frac), 
                    len(v)), 
                replace = False)) for k,v in group_indices.items()}

        # Get the final collapsed group indices
        collapsed_group_indices = list(collapse([v for k,v in \
                              group_indices_subset.items()]))
        
        # Subset anndata, copy as we are adding a column
        adata = adata[collapsed_group_indices].copy()
        
        # Convert dict celltype names to indices
        top_celltypes_groups_idxs = {
            k : ucdmetadata['celltypes_all'].index(v) for k,v in \
                             celltypes.items()}
        
        # Assign explain index that guides which celltype to explain for given sample
        adata.obs["explain_idx"] = adata.obs[groupby].map(top_celltypes_groups_idxs)
    
    else:
        raise Exception("Celltypes input type not compatible")
        
    return adata