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
import gzip
import base64

from .. import _utils as ucdutils
from .._settings import settings
from .._api import _api as ucdapi
from .._api import _apiutils as ucdapiutils
from .._api import _statewatchers as ucdstatewatchers

def deconvolve_select(
    data : Union[anndata.AnnData, pd.DataFrame],
    reference : Union[anndata.AnnData, pd.DataFrame, List[str], str],
    token : Optional[str] = None,
    reference_key : str = 'celltype',
    ignore_categories : Optional[Iterable[str]] = None,
    method : str = 'both',
    return_results : bool = False,
    key_added : str = 'ucdselect',
    use_raw : Union[bool, Tuple[bool, bool]] = True,
    verbosity : Optional[int] = None
    ) -> Optional[anndata.AnnData]:
    """\
    
    UniCell Deconvolve: Select
    
    Predicts cell type fractions for provided transcriptomic data
    using a user-specified reference. Leverages transfer learning
    from base UCD model embeddings.
    
    Params
    -------
    data
        Transcriptomic data (obs x genes) to predict cell type fractions. Can
        be either a dataframe or annotated dataset. Note that in any case
        data will be converted to an annotated datset object before proceeding.
    reference
        Transcriptomic data (obs x genes) to be used as a reference. Can be
        either a dataframe or annotated dataset. Note that if a dataframe is
        passed, row indices should correspond to categories for reference. If a list
        of strings is passed, these strings should correspond to reference profiles
        from the unicell cell type registry as any other names will throw an error. 
        If a string alone is passed, we look for a pre-built reference in the ucd backend.
        
        Currently valid prebuilt references include:
            allen-mouse-cortex : Mouse whole-brain cortex (44 cell types)
            enge2017-human-pancreas : Human pancreas (6 cell types)
            lee-human-pbmc-covid : Human PBMC (24 cell types)
            
    token
        UCDeconvolve API access token. If None, defaults to settings parameter.
    reference_key
        The key in reference.obs or index if reference is a dataframe to use to perform
        the grouping operation.
    method
        The method used for building a reference matrix. Must be one of two strings,
        either "embeddings" or "features". If "embeddings", the UCD base model is queried
        to return an embedding vector to represent celltype mixtures, and is
        used to generated representations for transfer learning. If "features", model
        defaults to using features in the reference matrix, similar to other available
        methods. Reccomended to use "both" in all cases.
    return_results
        Whether or not to return the predictions dict from the function,
        default to false as all data is written to anndata object either
        passed in, or created when passing in a dataframe, which will
        in that case be returned by default. Also returns the underlying
        anndata if it is a view as copying can destroy context internally. 
    ignore_categories
        Categories in 'reference.obs['reference_key']' to ignore. Default is None.
    use_raw
        Use counts in 'adata.raw'. Default True, as by convention in 
        single cell analysis, log1p scaled counts before HVG filter are kept here. Note
        that if a tuple is passed, it will selectively apply use_raw to DATA and then REF
        in that order.
    verbosity
        Logging verbosity, if None defaults to settings value.
        
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
        
    # If reference is dataframe make annotated dataset
    # Set classes reference_key to be the index of the original dataframe
    if isinstance(reference, pd.DataFrame):
        reference = anndata.AnnData(reference)
        reference.obs['classes'] = list(reference.obs_names)
        reference_key = "classes"
    
    # Wrap the entire lifecycle in a logging context
    with ucdutils.log(verbosity) as ucdlogger:
        
        # Wrap function in a timer
        with ucdutils.timer(ucdlogger, "Starting UCDeconvolveSELECT Run.", "Run Complete", logging.INFO):
            
            # Wrap the script in a try/except to catch keyboard interrupts
            try:
                
                # Check if reference is a list of strings, if so, attempt to get a custom reference object
                if isinstance(reference, list):
                    reference = ucdapi.build_custom_reference(reference)
                    reference_key = "classes"
                elif isinstance(reference, str):
                    # If instance is a string, collect prebuilt reference
                    reference = ucdapi.get_prebuilt_reference(reference)
                    reference_key = "classes"

                # Ensure inputs are correct type
                assert isinstance(data, (anndata.AnnData)),\
                    "data must either be a annotated dataset"
                assert isinstance(reference, (anndata.AnnData, pd.DataFrame)),\
                    "reference must either be a dataframe or annotated dataset"
                assert (method == 'embeddings') or (method == 'features') or (method == 'both'),\
                    "method must either be given as 'embeddings' or 'features'"

                # Assign vars as copies of originals that we will be modifying
                adata_mixture = data.copy()
                adata_reference = reference.copy()

                # Keep references to originals
                adata_mixture_orig = data
                adata_reference_orig = reference
                
                # Drop all obs data from mixture
                columns_to_drop = list(adata_mixture.obs.columns)
                adata_mixture.obs = adata_mixture.obs.drop(columns = adata_mixture.obs.columns)
                
                # Drop all columns minus reference key column from reference if it exists
                columns_to_drop = list(set(adata_reference.obs.columns) - (set([reference_key]) if reference_key else set()))
                adata_reference.obs = adata_reference.obs.drop(columns = adata_mixture.obs.columns)

                # Use raw from adata if specified, keep original 'adata' as reference to append to and return
                use_raw_mixture, use_raw_reference = (use_raw, use_raw) if isinstance(use_raw, bool) else use_raw
                adata_mixture = ucdutils.try_get_raw(adata_mixture, use_raw = use_raw_mixture)
                adata_reference = ucdutils.try_get_raw(adata_reference, use_raw = use_raw_reference)

                # Ensure scaled data is not provided
                assert adata_mixture.X.min() >= 0, "Detected expression < 0, make sure counts are provided"
                assert adata_reference.X.min() >= 0, "Detected expression < 0, make sure counts are provided"

                # Standardize gene names
                adata_reference.var_names = ucdutils.match_to_gene(adata_reference.var_names.tolist())
                adata_mixture.var_names = ucdutils.match_to_gene(adata_mixture.var_names.tolist())

                # Drop NA's and make copy for implicit view
                adata_reference = adata_reference[:,
                    ~adata_reference.var_names.isna() & ~adata_reference.var_names.duplicated()].copy()
                adata_mixture = adata_mixture[:,
                    ~adata_mixture.var_names.isna() & ~adata_mixture.var_names.duplicated()].copy()
                
                # Establish the celltype reference column
                adata_reference.obs['classes'] = adata_reference.obs[reference_key].astype(str) if \
                    isinstance(reference, anndata.AnnData) else reference.index.values.astype(str)

                # Eliminate unwanted categories from reference
                ignore_categories = ignore_categories if ignore_categories else (None,)
                ignore_categories = [ignore_categories] if isinstance(ignore_categories, str) else ignore_categories
                adata_reference = adata_reference[~adata_reference.obs[reference_key].isin(tuple(ignore_categories))]

                # Collapse the reference based on the grouping parameter
                adata_reference = ucdutils.collapse_anndata(adata_reference, "classes", op = np.mean)
                classes = adata_reference.obs['classes'].tolist()
                
                # Format inputs for streaming to UCD API, note that for the embeddings model
                # we need to preprocess / scale data before uploading as the model endpoint
                # does not perform these functions server-side.
                adata_mixture = ucdutils.get_preprocessed_anndata(adata_mixture, desc = "Mix", 
                                                                  progress = verbosity <= logging.INFO,
                                                                 minmax = False, zscore = False)
                adata_reference = ucdutils.get_preprocessed_anndata(adata_reference, desc = "Ref",
                                                                   progress = verbosity <= logging.INFO,
                                                                   minmax = False, zscore = False)
            
                # Package both objects into mudata object
                mdata = mudata.MuData({"mixture" : adata_mixture, "reference" : adata_reference})

                # Find common genes if features are being used in addition to embeddings
                if (method == 'features') or (method == 'both'):
                    common_genes = ucdutils.find_common_genes(mdata['mixture'], mdata['reference'])
                    mdata.uns['common_gene_indices'] =  base64.urlsafe_b64encode(gzip.compress(common_genes.tobytes())).decode()

                # Format params 
                params = {"method" : method }

                # Get signed URL
                metadata, url = ucdapi.start_run_request("select", params = params)

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
                    params['method'] = 'select'
                    
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