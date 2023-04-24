####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

# Functions for interacting with api
from typing import Optional, Dict, List

from .._data import metadata as ucdmetadata
from .._settings import settings
from .. import _utils as ucdutils
from . import _apiutils as apiutils

from requests.status_codes import codes as status

from getpass import getpass
import base64
import gzip
import anndata
import numpy as np
import pandas as pd
import tempfile
import requests
import logging
import anndata
import json
import os

def register(
    username : Optional[str] = None,
    password : Optional[str] = None,
    firstname : Optional[str] = None,
    lastname : Optional[str] = None,
    email : Optional[str] = None,
    institution : Optional[str] = None,
    dynamic : bool = True
    ) -> None:
    """\
    
    Registers a New User
    
    Params
    -------
    username
        Username for new account
    password
        Password for new account
    firstname
        First name of new user
    lastname
        Last name of new user
    email
        Valid email address of new user. Note that an email will be sent
        for account activation.
    institution
        The insitution, academic or private, the user is affiliated with.
    dynamic
        Whether or not to prompt for inputs dynamically, default is True.
        
    Returns
    -------
    Either nothing or waits for user to complete.
    
    """


    # If all parameters passed skip dynamic
    if username and password and firstname and lastname and email and institution:
        dynamic = False

    # Params for payload
    params = dict()

    # If dynamic is true, build request from user inputs, otherwise look for
    # inputs directly to function
    params['username'] = username if username and not dynamic else input("Username: ")
    params['password'] = password if password and not dynamic  else getpass("Password: ")

    # Ensure password is correct
    if dynamic:
        password2 = getpass("Confirm Password: ")

        match = params['password'] == password2

        while not match:
            password2 = getpass("Passwords do not Match | Try Again Password (type 'e' to exit): ")
            match = params['password'] == password2
            if password2 == "e": return

    # Get other parameters for registration
    params['first_name'] = firstname if firstname and not dynamic else input("First Name: ")
    params['last_name'] = lastname if lastname and not dynamic else input("Last Name: ")
    params['email'] = email if email and not dynamic else input("Email Address: ")

    # Ensure email is correct
    if dynamic:
        email2 = input("Confirm Email: ")

        match =  params['email'] == email2
        while not match:
            email2 = input("Emails do not match | Try Again Email (type 'e' to exit): ")
            match = params['email'] == email2
            if email2 == "e": return

    params['institution'] = institution if institution and not dynamic else input("Institution: ")

    # Send request and parse response
    response = apiutils.send_request("/api/register_user", params)

    # Get response
    if not dynamic:
        return apiutils.process_response(response)
    else:
        # Tell user to go check email
        print(f"Please check inbox at dmcharytonowicz@gmail.com and copy/paste activation code into prompt below.")
        
        # If we are dynamic, we can wait for user prompt to get the activation code
        activation_code = getpass("Activation Code: ")
        
        # Process activation
        return activate(activation_code)


def activate(
    code : Optional[str] = None,
    ) -> None:
    """\
    
    Activate User Account
    
    Activates account with an acitvation code recieved via email
    
    Params
    -------
    code
        Activation code emailed to user upon registration
    
    Returns
    -------
    None
    
    """


    if not code:
        code = getpass("Activation Code:")

    # Create params object
    params = {"activation_code" : code }

    # Send request and get response 
    response = apiutils.send_request("/api/activate_user", params)

    # Get response data
    data =  apiutils.process_response(response)

    # Assign key to client settings, None if not exist
    settings.token = data['token'] if data else None        
        

def authenticate(
    token : str
    ) -> None:
    """\
    
    Authenticate
    
    Updates user access token credentials
    
    Params
    -------
    token
        Valid user token
        
    Returns
    -------
    None
    
    """
    # Initialize logger
    with ucdutils.log(settings.verbosity) as ucdlogger:
        
        # Check if token is valid
        response = apiutils.send_request("api/validate_token", None, token)
        
        # If is valid set and notify
        if response.status_code == status.OK:

            # Set token 
            settings.token = token
            ucdlogger.info("Updated valid user access token.")
        
        else:
            ucdlogger.error("Invalid user access token.")
        


def get_task_state(
    run_id : str,
    runtype : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    Get Task state
    
    When a new run is submitted, it queues as a task and will be initiated
    when resources are available to process it. In most cases this will be
    instataneous. This function polls the task state.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    runtype
        The runtype being called.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """
    
    # Create payload
    payload = {"run_id" : run_id, "runtype" : runtype}

    # Send request and process response
    response = apiutils.send_request("api/get_task_status", payload, token)
    data = apiutils.process_response(response)
    
    return data

def get_run_exists(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    Get Run Exists
    
    Check if a run has started and exists, otherwise we are still
    waiting on the queue.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """
    
    # Create payload
    payload = {"run_id" : run_id}

    # Send request and process response
    response = apiutils.send_request("api/get_run_exists", payload, token)
    data = apiutils.process_response(response)
    
    return data
        
    
    
def get_run_state(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    Get Run Status
    
    Get the status of the current run.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """

    # Send request and process response
    response = apiutils.send_request("api/get_run_state", {"run_id" : run_id}, token)
    data = apiutils.process_response(response)

    return data['state']


def get_run_result(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    Get Run Result
    
    Get results link to results from a complete run.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """

    # Send request and process response 
    response = apiutils.send_request("api/get_run_result", {"run_id" : run_id}, token)
    data = apiutils.process_response(response)
    
    return data

def kill_run(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    Kill Run
    
    Send a kill signal to an ongoing run to free compute resources and allow for
    start of a new run.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """

    # Send request and process response 
    response = apiutils.send_request("api/kill_run", {"run_id" : run_id}, token)
    data = apiutils.process_response(response)
    
    return data

def list_run_states(
    token : Optional[str] = None,
    ) -> Dict:
    """\
    
    List Run States
    
    Lists status of any saved run states for user.
    
    Params
    -------
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """

    # Send request and process response 
    response = apiutils.send_request("api/list_run_states", None, token)
    data = apiutils.process_response(response)
    
    # Convert to dataframe
    data = pd.DataFrame(data).T
    data.index.name = "run_id"
    
    return data


def get_run_progress(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\ dd
    Get Run Progress
    
    Get progressbar associated with an active run if available.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.

    """

    # Send request
    response = apiutils.send_request("api/get_run_progress", {"run_id" : run_id}, token)
    
    # Process response, note that we ignore logging errors here
    # as get_run_result tends to be called too early before the progressbar
    # metadata is loaded. We wrap this in a waitfor callback to hold until
    # the result is available.
    data = apiutils.process_response(response, ignore_errors = True)

    return data


def remove_run(
    run_id : str,
    token : Optional[str] = None,
    ) -> Dict:
    """\
    Remove Run
    
    Cleans up run for user after completion.
    
    Params
    -------
    run_id
        The run_id for the run in question.
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    data : Dict
        A dictionary response object.
    
    """

    # Send request
    response = apiutils.send_request("api/remove_run", {"run_id" : run_id}, token)
    
    # Process response
    data = apiutils.process_response(response)

    return data['message']

def build_custom_reference(
    targets : List[str],
    token : Optional[str] = None,
    ) -> anndata.AnnData:
    """\
    Build Custom Reference
    
    Create a reference from UCD base dataset.
    
    Params
    -------
    targets
        A list of strings representing cell types from the UCD database.
        
    token
        Optional access token, if not passed default token is retrieved from settings.
    
    Returns
    -------
    adata : anndata.AnnData
        Annotated dataset object containing the custom reference.
    
    """

    # Send request
    response = apiutils.send_request("api/build_custom_reference", {"targets" : targets}, token)
    
    # Process response
    data = apiutils.process_response(response)
    
    # Extract results
    exp = data['exp']
    celltypes = data['targets']
    
    # Convert exp to array
    exp = np.clip(np.frombuffer(gzip.decompress(base64.urlsafe_b64decode(exp)), dtype = np.float16), 0, None)
    
    # Reshape
    exp = exp.reshape((len(celltypes), 28867))
    
    # Convert to pandas dataframe
    exp = pd.DataFrame(exp, columns = ucdmetadata['target_genes'])
    
    # Convert to annotated dataset
    adata = anndata.AnnData(exp)
    
    # Add celltypes as classes
    adata.obs['classes'] = celltypes
    
    return adata

def list_prebuilt_references(
    token : Optional[str] = None,
    ) -> anndata.AnnData:
    """\
    List Prebuilt References
    
    Lists available pre-built references
    
    Params
    -------
    token
        Optional access token, if not passed default token is retrieved from settings.
        
    Returns
    -------
    adata : anndata.AnnData
        Annotated dataset object containing the prebuilt reference.
    
    """
    
    # Send request
    response = apiutils.send_request("api/list_prebuilt_references", {}, token)
    
    # Process response
    data = apiutils.process_response(response)
    
    # Return response
    return data["references"]
    

def get_prebuilt_reference(
    reference : str,
    token : Optional[str] = None,
    cache : bool = True,
    ) -> anndata.AnnData:
    """\
    Get Prebuilt Reference
    
    Downloads a pre-made reference dataset.
    
    Params
    -------
    reference
        String name of prebuilt reference to get.
    token
        Optional access token, if not passed default token is retrieved from settings.
    cache
        If true, save a loaded file into cachedir and try to reload it before downloading again.
        
    Returns
    -------
    adata : anndata.AnnData
        Annotated dataset object containing the prebuilt reference.
    
    """
        
    # Check cachedir if file exists
    cachefile = os.path.join(settings.cachedir, f"{reference}.h5ad")
    if os.path.exists(cachefile):
        return anndata.read_h5ad(cachefile)

    # Send request
    response = apiutils.send_request("api/get_prebuilt_reference", {"reference" : reference}, token)
    
    # Process response
    data = apiutils.process_response(response)
    
    # Get download url
    download_url = data['download_url']
    
    # Make a GET request to the signed URL to download the file
    response = requests.get(download_url)
    
    # If we enable caching, we save the ifle into cachedir location under its reference name
    if cache:
        
        # Make the cache if it doesn't exist already
        os.makedirs(os.path.dirname(cachefile), exist_ok = True)
        
        # Open cachefile path for writing
        with open(cachefile, mode = "wb") as cache_file:
            
            # Write file into cache
            cache_file.write(response.content)
            
        # Read file from cache
        adata = anndata.read_h5ad(cachefile)
            
    else:
        # Write the downloaded content to a temporary file
        with tempfile.NamedTemporaryFile(suffix = ".h5ad") as temp_file:

            # Write file to temp
            temp_file.write(response.content)

            # Load the file from temp
            adata = anndata.read_h5ad(temp_file.name)
        
    # Return the dataset loaded from temp or cache
    return adata
    

def start_run_request(
    runtype : str,
    params : Optional[Dict] = None,
    token : Optional[str] = None
    ) -> None:
    """\
    
    Starts a run request by getting an upload URL with
    assigned metadata headers
    
    Params
    -------
    runtype
        The runtype to request. Can be one of
        'base', 'select', or 'explain'.
    params
        Params to include in metadata headers, specific
        to each method.
    token
        User access token, if none get from settings.
    
    Returns
    -------
    metadata : Dict
        Metadata to be passed back with upload
    url : str
        Signed upload URL
        
    """
    
    # Validate runtype input
    assert runtype in ('base', 'select', 'explain'), "Runtype invalid"
    
    # Format payload
    payload = {
        "runtype" : runtype,
        "params" : params
    }

    # Send request
    response = apiutils.send_request("api/start_run", payload, token)

    # Process response
    data = apiutils.process_response(response)

    # Split response
    metadata = data['metadata']
    url = data['url']
    
    # Return split response
    return metadata, url