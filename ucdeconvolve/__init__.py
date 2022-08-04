####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

"""UniCell Deconvolve - Cell Type Deconvolution For Transcriptomic Data."""

# Load metadata for package
from ._metadata import __version__, __author__, __email__
from ._metadata import __date__, __institution__, __laboratory__

from . import _tools as tl
from . import _plotting as pl