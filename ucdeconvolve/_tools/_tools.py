####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List, Dict
import anndata
import pandas as pd
import numpy as np
import scipy
import logging
import time
import math
import progressbar

from concurrent.futures import ThreadPoolExecutor
from functools import partial

from .. import _data as ucddata
from .. import _utils as ucdutils
from .. import _io as ucdio
from .. import _postprocessing as ucdpp

from .._exceptions import InvalidTokenException
from .._exceptions import PacketSizeException

from .._metadata import __version__ as ucdversion


def deconvolve(
    data : Union[anndata.AnnData, pd.DataFrame],
    token : str,
    split : bool = True,
    normalize_split_totals : bool = True,
    sort : bool = True,
    propagate : bool = True,
    return_results : bool = False,
    key_added : str = ucddata.metadata['default_key'],
    batchsize : int = 256,
    n_threads : int = 4,
    use_raw : bool = True,
    showprogress : bool = True,
    verbosity : int = logging.WARNING
    ) -> Optional[Union[anndata.AnnData, Dict[str, pd.DataFrame]]]:
    """\
    UniCell Deconvolve
    
    Predicts cell type fractions for provided transcriptomic data.
    
    Params
    -------
    data
        Transcriptomic data (obs x genes) to predict cell type fractions. Can
        be either a dataframe or annotated dataset. Note that in any case
        data will be converted to an annotated datset object before proceeding.
    token
        The API token to use for authorization. To sign up and recieve a token,
        go to https://forms.gle/NVrcEzn7Brva44QH6 and complete the form.
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
    batchsize
        Size of batches to process and send for predictions. Ideally should be
        divisible by 4 as batches are eventually split and sent as 4 prediction
        requests across separate threads for added performance. Maximum size
        is limited data bandwidth limits on prediction service, we reccomend
        no more than 256 sparse samples (i.e. cells) per batch.
    n_threads
        Number of threads to split prediction packets across. Improves speed of
        predictions. Can result in 504 timeout errors under high user load, 
        if so it is reccomended to reduce this value. 
    use_raw
        Use counts in 'adata.raw'. Default True, as by convention in 
        single cell analysis, log1p scaled counts before HVG filter are kept here.
    showprogress
        Show progressbar during predictions, default True.
    verbosity
        Level of verbosity for function information. Default is WARNING,
        set to 'logging.DEBUG' for more detailed information.
    """
    
    # Wrap the entire lifecycle in a logging context
    with ucdutils.log(verbosity) as ucdlogger:
        
        ucdlogger.info("Starting UniCell Deconvolve Run.")

        # Ensure data is correct type
        assert isinstance(data, (anndata.AnnData, pd.DataFrame)), "data must either be a dataframe or annotated dataset"
        
        # If data is a dataframe, convert to annotated dataset
        adata = anndata.AnnData(data) if isinstance(data, pd.DataFrame) else data
                
        # Record view status, will need to know if we have to return the eventual copy made as a result
        # of running this function.
        adata_is_view = adata.is_view
        
        # Use raw from adata if specified, keep original 'adata' as reference to append to and return
        adata_use = adata.raw.to_adata() if use_raw and adata.raw else adata

        # Generate the gene remapping index features to organize data columns
        features = ucdutils.remap_genes(ucdutils.match_to_gene(adata_use.var_names.tolist()))

        ###################################
        # Delayed Data Preprocessing Pipeline
        ###################################
        
        # Cap batchsize & thread counts
        batchsize = min(256, batchsize)
        effective_threads = min(4, n_threads)

        # Iterate adata in batches
        dataset = ucdutils.chunk(adata_use.X, batchsize // effective_threads)
        
        # calculate effective batch size
        batchsize_effective = min(batchsize, max(adata_use.shape[0], adata_use.shape[0] // 4))
        
        # Log dataset metrics if desired
        ucdlogger.debug(f"Using batchsize: {batchsize_effective} on dataset of shape: {adata_use.shape}")
        ucdlogger.debug(f"Data will be sent in {math.ceil(adata_use.shape[0] / batchsize_effective)} packet(s)," + \
                            f" each split into 4 sub-buffers of at most ~{batchsize_effective // 4} samples.")

        # Convert batch to dense, cast to float32, and preprocess for predictions
        dataset = map(lambda b : b.toarray() if  scipy.sparse.issparse(b) else b, dataset)
        dataset = map(lambda b : b.astype(np.float32, copy = False), dataset)
        dataset = map(lambda b : ucdutils.preprocess_expression(b, features), dataset)

        # Return batch to sparse CSR representation for efficient transfer
        dataset = map(lambda b : scipy.sparse.csr_matrix(b), dataset)

        # Encode data transfer packet
        dataset = map(lambda b : ucdutils.write_sparse_packet(b), dataset)

        # Create batch of buffered packets to send over HTTPS and compress
        dataset = ucdutils.buffer(dataset, 4)

        # Format packets list into instances dict required by prediction endpoint
        dataset = map(lambda b : {"instances" : b}, dataset)

        # Dump dict into json bytestring and perform final compression
        dataset = map(lambda b : ucdutils.json_dump_and_compress(b), dataset)

        # Report status
        ucdlogger.debug("Data pipeline ready.")
        
        ###################################
        # Prediction Process
        ###################################

        # Create a progressbar if showprogress is true otherwise passthrough
        pbar = progressbar.progressbar if showprogress else lambda x, max_value : x

        # Create a partial predict function and pass token information 
        predict = partial(ucdio.post_online_prediction, token = token)

        # Calculate number of iterations for progressbar
        n_iter = adata.shape[0] // batchsize
        
        # Log starting time for predictions for performance analysis
        ucdlogger.info("Initiating data stream for prediction.")
        start = time.time()
        
        # Send prediction requests using a threadpool
        ucdlogger.debug(f"Streaming using {effective_threads} thread(s).")
        with ThreadPoolExecutor(effective_threads) as executor:
            predictions = list(pbar(executor.map(predict, dataset), max_value = n_iter))
        
        # Get total time to complete streaming, rate, and hold as dictionary
        elapsed = round(time.time() - start,3)
        rate = round(adata.shape[0] / elapsed,3)
        runinfo = {'stream_time' : elapsed, 'stream_rate' : rate, 'split' : split}
        
        # Report end of data streaming and performance metrics
        ucdlogger.info("Data streaming complete.")
        ucdlogger.debug(f"Streaming time: {elapsed} (sec)")
        ucdlogger.debug(f"Streaming rate: {rate} (samples / sec)")

        # Concatenate predictions back together
        predictions = np.concatenate(predictions)
        
        # Postprocess predictions
        ucdlogger.info("Postprocessing predictions.")
        predictions = ucdpp.postprocess_predictions(predictions, split, normalize_split_totals, sort, propagate)
    
        # Append results to anndata with run metadata
        ucdlogger.info("Writing results to anndata object.")
        adata = ucdpp.append_predictions_to_anndata(predictions, adata, additional_runinfo = runinfo)
        
        # Return predictions
        ucdlogger.info("UniCell Deconvolve Run Complete.")
        
        # If the initial data was a dataframe, or a view of a truncated 
        # anndata, we need to return the annotated dataset
        # object made from it containing predictions
        if isinstance(data, pd.DataFrame) or adata_is_view:
            if return_results:
                return predictions, adata
            else:
                return adata
        else:
            # Initial object was an anotated dataset
            if return_results:
                return predictions
            
            
def read_results(
    adata : anndata.AnnData,
    category : Optional[str] = None,
    key : str = ucddata.metadata['default_key']
) -> pd.DataFrame:
    """\
    Read deconvolution results from an annotated dataset
    and return a dataframe.
    
    Params
    -------
    adata
        Annotated dataset with deconvolution results stored.
    category
        Which split category to use, defaults to 'all' if no
        split was made, or 'primary' if a split is detected from
        the run specificed by 'key'.
    key
        Key for deconvolution results, default key is 'ucd_results'
    
    """
    
    # Get the default split categry to show if none is provided
    default_category = 'all' if not adata.uns[key]['runinfo']['split'] else 'primary'
    category = default_category if not category else category
    
    # Get headers
    headers = adata.uns[key]['headers'][category]
    
    # Get data
    data = adata.obsm[f"{key}_{category}"]
    
    # Return dataframe
    return pd.DataFrame(data, columns = headers, index = adata.obs_names)