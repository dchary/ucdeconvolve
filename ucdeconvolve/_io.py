####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List
from ._data import metadata
import numpy as np
import base64
import gzip
import logging

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests.packages.urllib3.util.retry import MaxRetryError

from . import _utils as ucdutils
from ._exceptions import InvalidTokenException

def post_online_prediction(
    packet : bytes,
    token : str,
    max_retries : int = 6,
    retry_backoff_factor : float = 1.0,
    timeout : int = 60
) -> np.ndarray:
    """\
    
    Post Online Prediction
    
    Send formatted data packet for generating deconvolution
    predictions.
    
    Params
    -------
    data
        The data packet as a bytestring
    token
        API token
    max_retries
        Maximum number of times to retry connection if there is some type
        of server 500 class error. Can occur when server is under load or
        new instances are being spawned during autoscaling.
    retry_backoff_factor
        Factor multiplier for how long to wait before each retry.
    timeout
        Maximum time in seconds to wait for response before exiting.
    
    Returns
    -------
    response_data : np.ndarray
        A numpy array of raw deconvolution predictions
        
    """
    
    # Get logging object
    ucdlogger = logging.getLogger("UCD")
    
    # Create a request session context
    with requests.Session() as session:

        # Setup three maximum retries with a 0.25s backoff factor
        retry_strategy = Retry(total = max_retries, 
                               backoff_factor = retry_backoff_factor,
                               status_forcelist = [500, 502, 503, 504, 429],
                               allowed_methods = ["HEAD","GET","POST","OPTIONS"])
        
        # Setup a retry adapter and add it to session
        session.mount("https://", HTTPAdapter(max_retries=retry_strategy))

        # Update session headers
        headers = {"Content-Type": "application/gzip", "x-api-key" : token}

        # Create the actual request
        request = requests.Request('POST', metadata['gateway_url'], data = packet, headers = headers)

        try:
            
            # Attempt to get a response
            response = session.send(request.prepare(), timeout = timeout)
            
            # Check status code of response
            if response.status_code == 200:

                # If status is OK [200] then process the response as a valid prediction
                # Begin by reading the response as a json
                response_data = response.json()

                # Decompress predictions into raw bytes
                predictions = ucdutils.decompress_b64_gzip(response_data['predictions'])

                # Convert raw predictions bytes into a numpy array and reshape
                predictions = np.frombuffer(predictions, dtype = np.float16)
                
                # Reshape predictions according to shape attribute from response data
                predictions = predictions.reshape(response_data['shape'])

                # Return the predictions
                return predictions
            
            elif response.status_code == 503:
                # Report status that service is busy
                ucdlogger.error(f"Failed status code {response.status_code}," + \
                                "service in high-demend, consider reducing threads or try again later. returning empty predictions.")
                
            elif response.status_code == 400:
                raise InvalidTokenException("Please pass a valid API token to use this service.")
            else:
                # Report failed response code
                ucdlogger.error(f"Failed status code {response.status_code}, returning empty predictions.")
                
        except Exception as e:
            raise e
            
            # If the error was an invalid token exception, raise the error above the packet
            # level as we want to terminate the entire procese.
            if isinstance(e, InvalidTokenException):
                raise e
                
            # Log general error
            ucdlogger.error("Fatal error sending packet, returning empty predictions.")
            

        # If there is an error in generating predicitions that is not
        # recoverable, yield an empty array of negative 1s.
        packet_shapes = ucdutils.read_sparse_packet_attr(packet, attr = 'shape')
        
        # We only need to take the batch dimension from the packet shape, the returned value
        # is going to be the prediction width of 842
        return np.concatenate([-1 * np.ones((tuple(x)[0], 842), np.float16) for x in packet_shapes])
    