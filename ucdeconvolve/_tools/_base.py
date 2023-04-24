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

from .. import _utils as ucdutils
from .._settings import settings
from .._api import _api as ucdapi
from .._api import _apiutils as ucdapiutils
from .._api import _statewatchers as ucdstatewatchers

def deconvolve_base(
    data : Union[anndata.AnnData, pd.DataFrame],
    token : Optional[str] = None,
    split : bool = True,
    sort : bool = True,
    propagate : bool = True,
    return_results : bool = False,
    key_added : str = 'ucdbase',
    use_raw : Union[bool, Tuple[bool, bool]] = True,
    verbosity : Optional[int] = None,
    ) -> Optional[anndata.AnnData]:
    """\
    
    UniCell Deconvolve: Base
    
    Predicts cell type fractions for provided transcriptomic data.
    
    Params
    -------
    data
        Transcriptomic data (obs x genes) to predict cell type fractions. Can
        be either a dataframe or annotated dataset. Note that in any case
        data will be converted to an annotated datset object before proceeding.
    token
        UCDeconvolve API access token. If None, defaults to settings parameter.
    split
        Whether or not to split underlying data into
        three categories, primary, cancer cell_line.
        Helps with interpretability downstream, default
        is True.  
    sort
        Sort columns of results by mean predictions. Default True.
    propagate
        Whether or not to perform belief propagation and
        pass predictions up a cell-type heiarchy. helpful
        in interpreting overall deconvolution results.
        default is True.
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
        with ucdutils.timer(ucdlogger, "Starting UCDeconvolveBASE Run.", "Run Complete", logging.INFO):
            
            # Wrap the script in a try/except to catch keyboard interrupts
            try:

                # Ensure inputs are correct type
                assert isinstance(data, (anndata.AnnData)), "data must either be a annotated dataset"

                # Assign vars
                adata_mixture = data.copy()

                # Keep references to originals
                adata_mixture_orig = data
                
                # Drop all obs data from mixture
                columns_to_drop = list(adata_mixture.obs.columns)
                adata_mixture.obs = adata_mixture.obs.drop(columns = adata_mixture.obs.columns)

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
                                    progress = verbosity <= logging.INFO, minmax = False, zscore = False)

                # Package both objects into mudata object
                mdata = mudata.MuData({"mixture" : adata_mixture})

                # Format params 
                params = {"propagate" : propagate, "split" : split, "sort" : sort}

                # Get signed URL
                metadata, url = ucdapi.start_run_request("base", params = params)

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
                    params['method'] = 'base'
                    
                    adata_mixture_orig = \
                        ucdapiutils.download_and_attach_results(
                            adata_mixture_orig,
                            results['download_url'],
                            metadata['runtype'],
                            key_added,
                            params)

                # Remove and clean up run
                ucdapi.remove_run(metadata['run_id'], token)

                # Return original data with attached results
                # if return is set or if input was dataframe then return either way
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
                
        