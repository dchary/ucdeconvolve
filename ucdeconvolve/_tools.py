####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List
import anndata
import pandas as pd
import numpy as np
import scipy
import logging
import time
import math

from progressbar import progressbar
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from scipy.sparse import issparse

from ._data import metadata
from . import _utils as ucdutils
from ._io import post_online_prediction
from ._exceptions import InvalidTokenException
from ._postprocessing import postprocess_predictions

def deconvolve(
    data : Union[anndata.AnnData, pd.DataFrame],
    token : str,
    split : bool = True,
    normalize_split_totals : bool = True,
    sort : bool = True,
    propagate : bool = True,
    batchsize : int = 256,
    n_threads : int = 4,
    use_raw : bool = True,
    showprogress : bool = True,
    verbosity : int = logging.WARNING
    ) -> Optional[Union[anndata.AnnData, pd.DataFrame]]:
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
        whether or not to perform belief propagation and
        pass predictions up a cell-type heiarchy. helpful
        in interpreting overall deconvolution results.
        default is True.
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

        # Use raw from adata if specified
        adata = adata.raw.to_adata() if use_raw and adata.raw else adata

        # Generate the gene remapping index features to organize data columns
        features = ucdutils.remap_genes(ucdutils.match_to_gene(adata.var_names.tolist()))

        ###################################
        # Delayed Data Preprocessing Pipeline
        ###################################
        
        # Cap batchsize to 256
        batchsize = min(256, batchsize)

        # Iterate adata in batches
        dataset = ucdutils.chunk(adata.X, batchsize // 4)
        
        # calculate effective batch size
        batchsize_effective = min(batchsize, max(adata.shape[0], adata.shape[0] // batchsize))
        
        # Log dataset metrics if desired
        ucdlogger.debug(f"Using batchsize: {batchsize_effective} on dataset of shape: {adata.shape}")
        ucdlogger.debug(f"Data will be sent in {math.ceil(adata.shape[0] / batchsize_effective)} packet(s)," + \
                            f" each split into 4 sub-buffers of at most ~{batchsize_effective // 4} cells.")

        # Convert batch to dense, cast to float32, and preprocess for predictions
        dataset = map(lambda b : b.toarray() if issparse(b) else b, dataset)
        dataset = map(lambda b : b.astype(np.float32, copy = False), dataset)
        dataset = map(lambda b : ucdutils.preprocess_expression(b, features), dataset)

        # Return batch to sparse CSR representation for efficient transfer
        dataset = map(lambda b : scipy.sparse.csr_matrix(b), dataset)

        # Encode data transfer packet
        dataset = map(lambda b : ucdutils.write_sparse_packet(b), dataset)

        # Create batch of 4 buffered packets to send over HTTPS and compress
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
        pbar = progressbar if showprogress else lambda x, max_value : x

        # Create a partial predict function and pass token information 
        predict = partial(post_online_prediction, token = token)

        # Calculate number of iterations for progressbar
        n_iter = adata.shape[0] // batchsize
        
        # Log starting time for predictions for performance analysis
        ucdlogger.info("Initiating data stream for prediction.")
        start = time.time()
        
        # Send prediction requests using a threadpool
        with ThreadPoolExecutor(min(4, n_threads)) as executor:
            predictions = list(pbar(executor.map(predict, dataset), max_value = n_iter))
        
        # Get total time to complete streaming
        elapsed = time.time() - start
        
        # Report end of data streaming and performance metrics
        ucdlogger.info("Data streaming complete.")
        ucdlogger.debug(f"Streaming time: {round(elapsed,3)} (sec)")
        ucdlogger.debug(f"Streaming rate: {round(adata.shape[0] / elapsed,3)} (samples / sec)")

        # Concatenate predictions back together
        predictions = np.concatenate(predictions)
        
        # Postprocess predictions
        ucdlogger.info("Postprocessing predictions.")
        predictions = postprocess_predictions(predictions, split, normalize_split_totals, sort, propagate)

        # Return predictions
        ucdlogger.info("UniCell Deconvolve Run Complete.")
        return predictions