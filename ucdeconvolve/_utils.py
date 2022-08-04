####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List, Iterable, Dict
from ._data import metadata
import pandas as pd
import numpy as np
import re
import anndata
from scipy.sparse import issparse, spmatrix
import gzip
import base64
import scipy
import json
import logging
from contextlib import contextmanager

def chunk(
    obj : Iterable,
    chunksize : int = None
) -> Iterable:
    """\
    Iterable chunking
    
    Chunks an iterable (obj) into blocks of chunksize, 
    if the length of the last yielded value of this
    generaotr is smaller than the chunksize, then
    it will return a smaller chunk.
    
    Params
    -------
    obj
        Iterable-like object that can be chunked and
        has a length
    chunksize
        size of chunks to return
        
    Returns
    -------
    chunk : Iterable
        chunk of the underlying iterable
        
    """
    
    # Get chunksize, if not default to length of object
    chunksize = chunksize if chunksize else len(obj)
    
    # Return block of iterable object in chunks
    for i in range(0, len(obj) if not issparse(obj) else obj.shape[0], chunksize):
        yield obj[i : i + chunksize]

        
def buffer(
    iterable : Iterable,
    size : int
    ) -> List[object]:
    """\
    Buffer
    
    Creates a buffer that collects items from an iterable until a
    predefined size and then releases the items as a list of given
    size.
    
    Params
    -------
    iterable
        Iterable to accumulate frmo
    size
        Size of buffer to fill before release
    
    Returns
    -------
    buffer : List[object]
        A list of objects from the iterable of size.
    
    """
    
    # List buffer to hold elements
    buffer = []
    
    # Iterate items in iterable
    for item in iterable:
        
        # Append item to buffer
        buffer.append(item)
        
        # If buffer is full, return it
        # and then reset the buffer
        if len(buffer) == size:
            yield buffer
            buffer = []
            
    # If we are done iterating, return anything
    # left in the buffer even if it is not full
    if len(buffer) > 0:
        yield buffer


def match_to_gene(
    query : Union[List[str], str],
) -> List[str]:
    """
    Match to Gene
    
    Matches a string to known aliases for gene symbols using 
    an exhaustive dictionary of gene data.
    
    Params
    -------
    query
        query for gene symbols, can be a list or a single gene
    
    Returns
    -------
    gene : List[str]
        genes matched to dictionary, if not found returns None
        in the form of a list
        
    
    """
    
    # If query is not a list, make it one
    query = query if isinstance(query, (list, pd.Series)) else [query]    
    
    # Holds our matched genes
    genes = []
    
    # Iterate each gene
    for s in query:
        
        # Check first if the gene is even 
        # something remotely close to a valid type
        if isinstance(s, (str, int, float)):
            
            # Clean up name and standardize
            s = str(s)
            s = s.upper()
            s = re.sub("NM_[0-9]+\.[0-9]+","KRAKOW", s) # Removes transcript names

            # Attempt to match gene including delimiter split substring
            gene = metadata['gene_symbol_set'].get(re.split("\.|\_|\|", s)[0], None)

            # Tries the other side of the delimiter split if we fail on the first side
            if not gene:
                split = re.split("\.|\_|\|", s)[-1]
                
                # If the otherside of the split is too small then it failed
                # and return none
                if len(split) == 1:
                    gene = None
                else:
                    gene = metadata['gene_symbol_set'].get(split, None)
        else:
            gene = None
            
        # Append the potentially identified gene symbol to list
        genes.append(gene)
        
    # Return the final set of genes with None if it was not matched
    return genes


def remap_genes(
    features : List[str],
) -> np.ndarray:
    """\
    Remap Genes
    
    Remaps gene list into indices from UniCell reference
    gene list. Required to standardize input data.
    
    Params
    -------
    features
        list of genes as provided in dataset
    
    Returns
    -------
    indices
        array of indices which maps each gene position
        from source data to standard input reference.
    
    
    """
    
    # Convert features to a pandas series
    features = pd.Series(features)
    
    # Generate a gene dictionary and gene index for our positions with respect to unicell genes
    # This will be used to map our data to the ppi matrix gene order
    gd = {v:k for k,v in dict(enumerate(metadata['target_genes'])).items()}
    
    # Find each gene in our target set within the gene dictionary and get its column index
    gi = features.map(gd)
    
    # Assign any duplicates to be -1 (i.e. they go unmapped)
    gi.loc[gi.duplicated()] = -1
    
    # Get numpy array values
    gi = gi.values
    
    # Get any NaN values (unmapped features) and assign to -1 (unmapped)
    gi[np.isnan(gi)] = -1
    
    # Standardize column vector dtype
    gi = gi.astype(np.int32)
    
    # Return indices that determine which column position each input col from the source data
    # should recieve
    return gi


def reindex_array_by_columns(
    x : np.ndarray,
    y : Union[np.ndarray, List[int]],
    w : int
) -> np.ndarray:
    """\
    Reindex columns in x based on an index segments defined in y into an array of width w.
    
    Params
    -------
    x
        array with columns to be remapped
    y
        definition of index segments, of same size as x, as to where each column
        belongs in a new array of same length, but width w.
    w
        width of target resulting array.
        
    Returns
    -------
    arr : np.ndarray
        resulting data with remapped columns
    
    """
    
    # Create a buffer array to hold results based on the length of x and width given by w
    # Note that we add an additional column to the buffer to hold our null (-1) assignments
    # which we will discard before returning.
    buffer = np.zeros((x.shape[0], w + 1))
    
    # Put the slices into their positions in the buffer according to y
    np.put_along_axis(buffer, np.expand_dims(np.asarray(y), 0), x, axis=1)
    
    # Return the buffer without the additional column
    return buffer[:, :-1]


def scale_to_min_max(
    x : np.ndarray,
    t_min_value : float = 0.0,
    t_max_value : float = 1.0,
    axis : Optional[int] = None
) -> np.ndarray:
    """\
    Scale to min max. Scales array either globally or along
    specified axis to min-max.
    
    Params
    -------
    x
        an np.ndarray to scale
    t_min_value
        target minimum value for scaling
    t_max_value
        target maximum value for scaling
    axis
        axis along which to perform scaling
        
    Returns
    -------
    x : np.ndarray
        rescaled array
        
        
    """
    
    # Keepdims if we are using an axis
    kd = True if axis else False
    
    # Get the actual min and max along each axis
    r_min_value = np.min(x, axis = axis, keepdims = kd)
    r_max_value = np.max(x, axis = axis, keepdims = kd)

    # Normalize data and get the rescaled values
    return (((x - r_min_value) / (r_max_value - r_min_value)) \
                 * (t_max_value - t_min_value)) + t_min_value
    

def scale_to_z_score(
    x : np.ndarray,
    axis : Optional[int] = None
) -> np.ndarray:
    """\
    Peform z-score normalization on array.
    
    Params
    -------
    x
        an np.ndarray to normalize
    axis
        axis along which to perform scaling
        
    Returns
    -------
    x : np.ndarray
        rescaled array
        
        
    """
    
    # Keepdims if we are using an axis
    kd = True if axis else False
    
    # Get mean and std deviation
    mu = np.mean(x, axis = axis, keepdims = kd)
    std = np.std(x, axis = axis, keepdims = kd)
    
    # Return z-score
    return (x - mu) / std


def normalize_total(
    x : np.ndarray,
    total : float = 1e4,
    axis : Optional[int] = None
) -> np.ndarray:
    """\
    Peform total normalization on array.
    
    Params
    -------
    x
        an np.ndarray to normalize.
    total
        Number to normalize to.
    axis
        axis along which to perform normalization.
        
    Returns
    -------
    x : np.ndarray
        rescaled array
        
        
    """
    
    # Keepdims if we are using an axis
    kd = True if axis else False
    
    # Return normalized total values
    return x * (total / np.sum(x, axis, keepdims = kd))


def preprocess_expression(
    exp : np.ndarray,
    feature_index : Union[Iterable[str], np.ndarray],
    sum_total : float = 1e4,
    log_detect_threshold : float = 80.0,
    clip_threshold : float = 10.0
) -> np.ndarray:
    """\
    
    Preprocesses expression data for input into deconvolution pipeline.
    
    Params
    -------
    exp
        Chunk of unnormalized expression data.
    feature_index
        index of remapped column names.
    sum_total
        total counts to normalize data to.
    log_detect_threshold
        maximum value to use to detect if genes are log-normalized already.
    clip_threshold
        threshold for clipping z-scores
        
    Returns
    -------
    exp : np.ndarray
        normalized expression data
        
    """
    
    # Start by reindexing data columns to match required shape by model
    exp = reindex_array_by_columns(exp, feature_index, len(metadata['target_genes']))
    
    # Make sure that each row is not log normalized already.
    exp = np.expm1(exp) if exp.max() < log_detect_threshold else exp
    
    # Normalize row counts
    exp = normalize_total(exp, total = sum_total, axis = 1)
    
    # Log normalize and add pseudocount log10(x + 1)
    exp = np.log1p(exp)
    
    # Center data and normalize variance with z-score
    exp = scale_to_z_score(exp, axis = 1)
    
    # Clip z-score range 
    exp = np.clip(exp, -clip_threshold, clip_threshold)
    
    # Rescale data to minmax from 0 - 1
    exp = scale_to_min_max(exp, 0.0, 1.0, axis = 1)
    
    # Ensure no nans or infinities
    exp = np.nan_to_num(exp, nan = 0.0, posinf = 0.0, neginf = 0.0)
    
    # Return normalized transformed data
    return exp


def encode_arr1d_as_str(
    arr : np.ndarray,
    dtype : Optional[np.dtype] = None,
    compresslevel : Optional[int] = 9,
    b64_altchars : Optional[bytes] = b'-_',
    b64_rep : str = 'ascii'
    ) -> str:
    """\
    
    Encode a numpy array as a compressed web-safe string
    that can be sent as a json string value attribute.
    
    Params
    -------
    arr
        A 1-D numpy array to pass, is validated as a numpy array before continuing
    dtype
        The datatype to transform to, if None then keeps original type
    compresslevel
        Level of gzip compression to use. If None do not gzip
    b64_altchars
        Alternate characters for base64 encoding of compressed string. Default is
        web-safe chars compatible with tensorflow base64 decoding. 
    b64_rep
        String representation to use for base64 encoding, default is ascii.
        
    Returns
    -------
    rep : str
        An encoded compressed string representation of the underlying data that can
        be sent as a json string value.
        
    """
    
    # Validate the array
    assert isinstance(arr, np.ndarray), f"Array {arr} is not a numpy array"
    assert arr.ndim == 1, f"Array {arr} has {arr.ndim} dimensions but must only have 1"
    
    # Convert the data type to desired dtype and then to a bytes representation
    asbytes = arr.astype(dtype if dtype else arr.dtype).tobytes()
    
    # Compress the data if compresslevel is provided otherwise just pass through
    compressed = gzip.compress(asbytes, compresslevel = compresslevel) if compresslevel else asbytes
    
    # Peform web-safe base64 encoding
    encoded = base64.b64encode(compressed, altchars = b64_altchars)
    
    # Get string representation of base64 encoding as desired type
    rep = encoded.decode(b64_rep)
    
    return rep    


def decode_arr1d_from_str(
    bytestring : str,
    dtype : Optional[np.dtype] = None,
    b64_altchars : Optional[bytes] = b'-_',
    b64_rep : str = 'ascii'
    ) -> str:
    """\
    
    Decodes a a json string value attribute back into
    a numpy array from a gzipped b64 encoding.
    
    Params
    -------
    arr
        A 1-D numpy array to pass, is validated as a numpy array before continuing
    dtype
        The datatype to transform to, if None then keeps original type
    b64_altchars
        Alternate characters for base64 decoding of compressed string. Default is
        web-safe chars.
    b64_rep
        String representation to use for base64 decoding, default is ascii.
        
    Returns
    -------
    arr : np.ndarray
        Decoded underlying array
        
    """
    

    # Convert the data type to desired dtype and then to a bytes representation
    asbytes = arr.astype(dtype if dtype else arr.dtype).tobytes()
    
    # Compress the data if compresslevel is provided otherwise just pass through
    compressed = gzip.compress(asbytes, compresslevel = compresslevel) if compresslevel else asbytes
    
    # Peform web-safe base64 encoding
    encoded = base64.b64encode(compressed, altchars = b64_altchars)
    
    # Get string representation of base64 encoding as desired type
    rep = encoded.decode(b64_rep)
    
    return rep    


        
def write_sparse_packet(
    arr : Union[scipy.sparse.csr_matrix],
    dtype_indices = np.uint16,
    dtype_indptr = np.int32,
    dtype_data = np.float16,
    dtype_shape = np.uint16,
    compresslevel = 9
    ) -> Dict:
    """\
    Convert a data array into a row-major sparse compressed data packet
    that can be sent over a web connection.
    
    Params
    -------
    arr
        Underlying data array, must be a csr_sparse_matrix
    dtype_indices
        Datatype to cast indices, default uint16
    dtype_indptr
        Datatype to cast indptr, default int32
    dtype_data
        Datatype to cast expression data, default is float16
    dtype_shape
        Datatype to cast shape, default uint16
    compresslevel
        GZIP compression level, default is 9 (maximum)
    
    Returns
    -------
    
    
    """

    # If the data matrix is not sparse convert to CSR matrix 
    # (note that a dense matrix may take up more space using this approach)
    assert isinstance(arr, scipy.sparse.csr_matrix), "Arr must be a CSR matrix"

    # Convert data into compressed base64 encoded strings
    indices = encode_arr1d_as_str(arr.indices, dtype_indices, compresslevel)
    indptr = encode_arr1d_as_str(arr.indptr, dtype_indptr, compresslevel)
    data = encode_arr1d_as_str(arr.data, dtype_data, compresslevel)
    shape = encode_arr1d_as_str(np.asarray(arr.shape), dtype_shape, compresslevel)

    # Return json dictionary with compressed data
    return {"data" : data, "indices" : indices, "indptr" : indptr, "shape" : shape}

def decompress_b64_gzip(
    data : str,
    b64_altchars : Optional[bytes] = b'-_',
    ) -> bytes:
    """\
    Decompress data that was previously
    b64 encoded and then gzipped.
    
    Params
    -------
    data
        The data package as an underlying b64 encoded string
    b64_altchars
        Alternate characters for base64 decoding of compressed string. Default is
        web-safe chars.
    
    Returns
    -------
    data : bytes
        Data as a decoded bytestring
    
    """
    
    return gzip.decompress(base64.b64decode(data, altchars = b64_altchars))
    
    

def read_sparse_packet_attr(
    packet : bytes,
    attr : str = 'shape',
    dtype = np.uint16,
    ) -> Dict:
    """\
    Decompress a a row-major sparse compressed data packet.
    Basically invert the work done in write_sparse_packet
    in order to read an attribute, mainly shape.
    
    Params
    -------
    packet
        A bytestring to decompress
    attr
        Attribute to read, default is 'shape', useful to know
        what this packet was actually holding.
    dtype
        Datatype to read attribute as, default uint16
        
    Returns
    -------
    packet : dict
        A dictionary of our uncompressed packet
        
    """
    
    # Decompress the entire packet
    packet = json.loads(gzip.decompress(packet))
    
    # Function to decompress specific arrays
    decompress = lambda x, d : np.frombuffer(
                                        gzip.decompress(
                                            base64.b64decode(x, altchars = b'-_')), dtype = d)
    
    # Decompress desired attribute array from packet instances and return list
    return [decompress(x[attr], dtype) for x in packet['instances']]


def json_dump_and_compress(
    data : Dict,
    compresslevel : int = 9
    ) -> bytes:
    """\
    
    Dump a dictionary containing data into a json bytestring and
    compress it with gzip. 
    
    Params
    ------
    data
        A dictionary with data
    compresslevel
        Level of gzip compression, default is 9
    
    Returns
    ------
    data : bytes
        A compressed bytestring of data
    
    """
    
    return gzip.compress(json.dumps(data).encode(), compresslevel = compresslevel)
    
@contextmanager
def log(
    level : int
) -> None:
    """\
    Logging Context Manager
    
    Allows UCD to spawn a temporary logging session during each run in order to
    have control over writing status updates.
    
    Params
    -------
    level
        The level to set the logger to. Will temporarily set the root logging level
        to this, and then reset it back to its original upon exiting the context.
        
    """
    
    # Get the default logger
    logger = logging.getLogger("UCD")
    
    # Get the current level and set logger to it
    current_level = logger.getEffectiveLevel()
    logger.setLevel(level)

    # Create a console handler to print information
    handler = logging.StreamHandler()

    # Create a formatter to format logging output
    formatter = logging.Formatter('%(asctime)s|[%(name)s]|%(levelname)s: %(message)s')

    # Add the formatter to the handler
    handler.setFormatter(formatter)

    # Add the handler to the logger
    logger.addHandler(handler)
    
    # prevent logger propagation
    logger.propagate = False
    
    try:
        yield logger
        
    finally:
        
        # Set the level back to what it once was
        logger.setLevel(current_level)
        
        # Remove the handler we added
        logger.removeHandler(handler)
        
