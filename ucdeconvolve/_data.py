####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

import bz2, pickle, os

# Load metadata files
with bz2.BZ2File(os.path.join(os.path.dirname(__file__), "data/metadata.pkl"), 'rb') as f:
    metadata = pickle.load(f)