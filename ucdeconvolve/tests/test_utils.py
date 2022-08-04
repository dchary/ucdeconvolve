####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

"""\

UniCell Deconvolve Tests

"""

import pytest
import os
import sys
import pandas as pd
import numpy as np

# Bring package onto the path
sys.path.append(os.path.abspath(os.path.join('..')))

# Import ucd utils for testing
import ucdeconvolve._utils as ucd_utils


# Test chunk function
@pytest.mark.parametrize("test_input,expected", [([0,1,2,3,4,5,6,7,8,9,10], [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10]])])
def test_chunk(test_input, expected):
    assert list(ucd_utils.chunk(test_input, 3)) == expected

    
# Test buffer function
@pytest.mark.parametrize("test_input,expected", [(list(range(9)), [[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
                                                (list(range(10)), [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9]])])
def test_buffer(test_input, expected):
    assert list(ucd_utils.buffer(test_input, 3)) == expected

    
# Test match_to_gene function
@pytest.mark.parametrize("test_input,expected", [(['SFTPC','krAS','CD103'], ['SFTPC', 'KRAS', 'ITGAE']),
                                                (['HGNC:24086','oinsgoinv','epcam'], ['A1CF', None, 'EPCAM']),
                                                (pd.Series(['HGNC:24086','oinsgoinv','epcam']), ['A1CF', None, 'EPCAM'])])
def test_match_to_gene(test_input, expected):
    assert ucd_utils.match_to_gene(test_input) == expected
    
    
# Test remap_genes function
@pytest.mark.parametrize("test_input,expected", [(['SFTPC','KRAS','ITGAE'], [22267, 12080, 11389])])
def test_remap_genes(test_input, expected):
    assert list(ucd_utils.remap_genes(test_input)) == expected

    
# Test reindex_array_by_columns function
@pytest.mark.parametrize("test_input,expected", [(np.array([[9, 2], [5, 5]]), np.array([[9., 0., 2., 0., 0.], [5., 0., 5., 0., 0.]]))])
def test_reindex_array_by_columns(test_input, expected):
    assert np.all(ucd_utils.reindex_array_by_columns(test_input, [0,2], 5) == expected)
    

# Test scale_to_min_max function
@pytest.mark.parametrize("test_input,expected", [(np.array([[10,5], [2,3]]), np.array([[1.0, 0.0], [0.0,1.0]]))])
def test_scale_to_min_max(test_input, expected):
    assert np.all(ucd_utils.scale_to_min_max(test_input, axis = 1) == expected)

    
@pytest.mark.parametrize("test_input,expected", [(np.array([-0.98058068,  1.37281295, -0.39223227]), 0.0)])
def test_scale_to_z_score(test_input, expected):
    assert ucd_utils.scale_to_z_score(test_input).mean() == expected

    
@pytest.mark.parametrize("test_input_axis,expected", [(0, 10000.0), (1, 10000.0)])
def test_normalize_total(test_input_axis, expected):
    
    # Create some random array to pass as input
    a = np.random.randint(0, 100, (32,512))
    
    assert ucd_utils.normalize_total(a, axis = test_input_axis).sum(test_input_axis).mean() == expected