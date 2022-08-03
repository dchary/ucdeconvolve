####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

import bz2
import pickle
import os
from typing import Optional

def __update_metadata_dict(
    key : object,
    value : Optional[object] = None,
    remove : bool = False,
) -> None:
    """\
    
    Update metadata dictionary. Should only be used for development purposes.
    
    Params
    -------
    key
        key to add to metadata
    value
        value corresponding to key
    remove
        if true then remove the key and corresponding value
    Returns
    -------
    None
    
    """
    
    
    # Load metadata files
    with bz2.BZ2File(os.path.join(os.path.dirname(__file__), "data/metadata.pkl"), 'rb') as f:
        metadata = pickle.load(f)
    
    
    if not remove:
        # update value
        metadata[key] = value
    else:
        # If remove is true then remove the key if it exists
        metadata.pop(key, None)
    
    # Write updated metadata file
    with bz2.BZ2File(os.path.join(os.path.dirname(__file__), "data/metadata.pkl"), 'wb') as f:
        pickle.dump(metadata, f)
        
# Load metadata files
with bz2.BZ2File(os.path.join(os.path.dirname(__file__), "data/metadata.pkl"), 'rb') as f:
    metadata = pickle.load(f)