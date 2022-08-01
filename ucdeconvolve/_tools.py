####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Union, Optional, Tuple, List
import anndata
import pandas as pd

from ._data import metadata

def deconvolve(
    data : Union[anndata.AnnData, pd.DataFrame],
    token : str,
    use_raw : bool = True,
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
    use_raw
        Use counts in 'adata.raw'. Default True, as by convention in 
        single cell analysis, log1p scaled counts before HVG filter are kept here.
    
    """
    
    
    
    
    
    return adata