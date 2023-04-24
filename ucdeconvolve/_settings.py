####################################################################################################
# # Copyright (C) 2023-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

import logging


# Config Object
class UCDeconvolveConfig:
    """\
    Config manager for ucdeconvolve.
    
    Params
    -------
    verbosity
        Default verbosity for loggers, baseline is info
    token
        UCDeconovlve API token, default is None, can be set with 'ucd.api.authenticate'
        or diretly with 'ucd.settings.credentials'
    cachedir
        Location to store cached datafiles such as premade references for ucdselect.
    """
    
    def __init__(
        self,
        verbosity : int = logging.DEBUG,
        token : str = None,
        cachedir : str = "./cache/ucdeconvolve/"
        ) -> object:
        
        self.verbosity = verbosity
        self.token = token
        self.cachedir = cachedir
        
        
        
settings = UCDeconvolveConfig()