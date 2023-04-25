####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

# Support functions for interacting with API
from typing import Optional, Dict

from .._data import metadata as ucdmetadata
from .._settings import settings
from .. import _utils as ucdutils
from .._postprocessing import append_predictions_to_anndata

import pandas as pd
import urllib
import requests
import os
import json
import random
import hashlib
import time
import base64
import gzip
import h5py
import scipy
import warnings
import tempfile
import requests
import mudata
import anndata

def _make_url(base_url , *components, **params) -> str:
    """\
    Make URL
    
    Use urllib to construct url from componenets
    
    """
    
    # Set up the base url
    url = base_url
    
    # Attach components to path
    for component in components:
        url = '{}/{}'.format(url, component)
    
    # Add optional params
    if params:
        url = '{}?{}'.format(url, urllib.urlencode(params))
    
    # Yield URL
    return url


def send_request(
    endpoint : str,
    payload : Optional[str] = None,
    token : Optional[str] = None,
    method : str = 'POST'
    ) -> requests.Response:
    """\
    
    Build and send an authorized API request
    
    Params
    -------
    endpoint
        Target API endpoint
    payload
        Target payload
    token
        User access token, if none infer from credentials
    method
        API method, default POST
    
    Returns
    -------
    response : requests.Response
        A response object
    """

    # If token not passed, retrieve from credentials object
    token = token if token else settings.token
    
    # Use a context manager
    with requests.Session() as session:
        
        # Create headers
        session.headers.update({"Content-Type" : "application/json"})
        session.headers.update({"x-api-key" : ucdmetadata['api_gateway_key']})
        session.headers.update({"Authorization" : f"Token {token}"})
        
        # Combine the base gateway URL with the target request endpoint
        url = _make_url(ucdmetadata['api_gateway_url'], endpoint)
        
        # Send request
        response = session.request(method, url, json = payload)
        
        # Get response
        return response
    
    
def process_response(
    response : requests.Response,
    as_json : bool = True,
    ignore_errors : bool = False,
    ) -> Optional[str]:
    """\
    
    Default response processor
    
    Params
    -------
    response
        A valid http response
    as_json
        Flag to return valid response as a json
        object, default = True/
    Returns
    -------
    content : Optional[str]
        A content string that may contain JSON data
        If error returns None
    
    """
    
    # Initialize logger
    with ucdutils.log(settings.verbosity) as ucdlogger:

        # Get decoded content from response
        content = response.content.decode()
        
        try:
            
            # Check for errors in response
            response.raise_for_status()

            # Convert to json if desired
            content = json.loads(content) if as_json else content
            
            # Return data object
            return content
            
        except Exception as e:
            
            # Log error
            if not ignore_errors: ucdlogger.error(content)
            
            # Return token as None
            return content
        
        
def get_random_string(length=12,
                      allowed_chars='abcdefghijklmnopqrstuvwxyz'
                                    'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789',
                     seed :str = "SEED"):
    """
    Returns a securely generated random string.
    The default length of 12 with the a-z, A-Z, 0-9 character set returns
    a 71-bit value. log_2((26+26+10)^12) =~ 71 bits
    
    Borrowed from Django API
    """
    
    # Generate a random seed each time
    random.seed(
        hashlib.sha256(
            ("%s%s%s" % (
                random.getstate(),
                time.time(),seed)).encode('utf-8')
        ).digest())
    return ''.join(random.choice(allowed_chars) for i in range(length))


def upload_mudata_to_cloud_storage_url(
    mdata : mudata.MuData,
    url : str,
    metadata : dict,
) -> str:
    """\
    Upload mudata object to cloud storage signed_url
    
    Params
    -------
    mdata
        Mudata object to upload
    url
        Upload path
        
    Returns
    -------
    response : requests.Response
        A response object
        
    """
    
    # Check object is mudata
    assert isinstance(mdata, mudata.MuData), f"Object {mdata} must be mudata object"

    with tempfile.TemporaryDirectory() as tmpdirname:
        
        # Create tempdir and filename combination
        tmpfilename = os.path.join(tmpdirname, f"{get_random_string()}.h5mu")
    
        # Write out packaged anndatas to tmp, 
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            mdata.write_h5mu(tmpfilename, compression = "gzip")
        
        # Check file size
        filesize = os.path.getsize(tmpfilename)
        
        # If filesize will exceed limits, raise an error
        max_size = 100000000 if metadata['runtype'] == 'explain' else 500000000
        assert filesize < max_size, f"File for '{metadata['runtype']}' mode must be <500MB,"+\
            " for larger queries please contact ucdeconvolve@gmail.com to enable support"
        
        #Temporarily open the file
        with open(tmpfilename, 'rb') as data:

            headers =   {"content_type":"application/octet-stream"}
            headers.update({f"x-goog-meta-{k}" : str(v) for k,v in metadata.items()})
        
            # Create request
            response = requests.put(url, data = data, headers = headers)
            
    return response



def download_and_attach_results(
    adata : anndata.AnnData,
    download_url : str,
    runtype : str,
    key_added : str,
    additional_runinfo : Optional[dict] = dict()
    ) -> str:
    """\
    
    Download Result
    
    Downloads results at attaches them to the anndata object
    depending on the type of run.
    
    Params
    -------
    adata
        Annotated dataset to save results to
    download_url
        A valid signed URL to download results content to
    runtype
        The type of run which changes how we save data.
    key_added
        Key to use
        
    Returns
    -------
    adata : anndata.AnnData
        An annotated dataset with saved results
    
    """
    
    # Make a GET request to the signed URL to download the file
    response = requests.get(download_url)
    
    # Write the downloaded content to a temporary file
    with tempfile.NamedTemporaryFile() as temp_file:
        temp_file.write(response.content)
        
        # Save function for base
        if runtype == 'base':
            
            # Get keys from file
            with h5py.File(temp_file.name) as hfile:
                keys = list(hfile.keys())

            # Iterate each results key and load
            preds = dict()
            for key in keys:
                preds[key] = pd.read_hdf(temp_file.name, key = key)

            # attach to anndata
            append_predictions_to_anndata(preds, adata, 
                key_added = key_added, 
                additional_runinfo = additional_runinfo)
            
        if runtype == 'select':

            # Get keys from file
            with h5py.File(temp_file.name) as hfile:
                keys = list(hfile.keys())

            # Iterate each results key and load
            preds = dict()
            for key in keys:
                preds[key] = pd.read_hdf(temp_file.name, key = key)

            # attach to anndata
            append_predictions_to_anndata(preds, adata, 
                key_added = key_added, 
                additional_runinfo = additional_runinfo)
            
            
        # Save function for explain
        if runtype == 'explain':
            
            # Load explanations file as sparse CSR matrix
            explanations = _read_csr_matrix_from_h5(temp_file.name)
            
            # Don't attach here, return directly?
            #adata.obsm[key_added] = explanations
            
            # attach runinfo
            adata.uns[key_added] = dict(runinfo=additional_runinfo)
            
            # We need to fix this later
            return explanations
            
    # return anndata
    return adata
    
    
def _read_csr_matrix_from_h5(
    filepath : str
    ) -> scipy.sparse.csr_matrix:
    """\
    
    Read CSR from H5 File
    
    Params
    -------
    filepath : str
        A valid filepath
        
    Returns
    -------
    mtx : scipy.sparse.csr_matrix
        A csr sparse matrix
        
    """
    
    with h5py.File(filepath, 'r') as hf:
        data = hf['data'][:]
        indices = hf['indices'][:]
        indptr = hf['indptr'][:]
        shape = hf['shape'][:]
        return scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)


    
def waitfor(func, null_signal = "", poll_interval : float = 0.5, timeout : int = 60):
    """\
    
    WaitFor
    
    Wrapper around http async functions that waits for some non-null
    signal to be released by a function before continuing
    
    Params
    -------
    func
        The function returning some signal
    null_signal
        The signal from function specifying it is not ready yet, if the return
        is anything but this, then release the result
    poll_interval
        Time to wait between polling functiion
    timeout
        Time before raising an exception, default is 1 minute
        
    Returns
    -------
    val
        Return value from func after non-null signal.
    
    """
    
    # Mark initial start time
    start_time = time.time()
    
    # Infinite loop
    while True:
        
        # Get result from function
        result = func()
        
        # If result is not null signal, break loops
        if result != null_signal: 
            break
            
        # Otherwise, sleep for the polling interval
        time.sleep(poll_interval)
        
        # If time exceeds timeout length, break anyway
        if (time.time() - start_time) > timeout:
            break
            
    # Return the result
    return result