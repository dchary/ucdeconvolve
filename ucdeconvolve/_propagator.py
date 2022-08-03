####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################


"""

======================
UniCell: Deconvolve - Belief Propagator
======================

Description: Propagating deconvolution predictions across a cell heiarchy.

"""

from sklearn.base import BaseEstimator, TransformerMixin, RegressorMixin
import pandas as pd
import numpy as np
import networkx as nx
from typing import List, Union, Optional
import sys, os
from numba import njit
import warnings
import pickle
import logging

from ._data import metadata
from . import _utils as ucdutils

@njit
def pass_each_cell(
    data : np.ndarray,
    A : np.ndarray,
    I : np.ndarray
) -> np.ndarray:
    """
    Pass Each Cell
    
    Perform message passing on each individual cell (row) of
    a raw dataframe. Uses numba to pre-compile for faster
    execution.
    
    Params
    -------
    data
        the raw data for message passing to occur
    A
        adjacency matrix of the underlying node graph.
    I
        indexes mapping each column (cell type) of the raw 'data'
        to a node in the graph defined by A.
    
    Returns
    ------
    buffer : np.ndarray
        buffer contained the message passed data.
        
    """
    
    # Output buffer of extended data
    buffer = np.zeros((data.shape[0], A.shape[0]), dtype = np.float32)
    
    for i in np.arange(data.shape[1]):
        if I[i] == -1: continue
        message_pass(buffer, data, A, I[i], data[:,i])
    
    return buffer

@njit
def message_pass(
    B : np.ndarray,
    D : np.ndarray,
    A : np.ndarray,
    pos : int,
    value : float):
    """\
    Message Pass
    
    Cell-level recursive message passing utility.
    
    Params
    -------
    B
        buffer storing message passed values.
    D
        original data matrix.
    A
        adjacency matrix to traverse
    pos
        current position (node) we are at
    value
        current predicted value at the current node
    
    Returns
    ------
    None
        passes predictions along buffer reference.
        
    """
        
    # Update buffer column with the latest value
    B[:,pos] = B[:,pos] + value
    
    next_stops = np.argwhere(A[pos]).flatten()
    
    if len(next_stops) == 0: return
    
    for next_stop in next_stops:
        message_pass(B, D, A, next_stop, (1.0 / len(next_stops)) * value)


class UCDeconvolvePredPropagator(TransformerMixin, BaseEstimator):
    """Unicell Deconvolve Prediction Belief Progataion
    
    This module takes the raw 840 unit vector predictions and expands
    it to the entire cell type annotation tree for easier interpretation.
    
    """
    
    def __init__(
        self,
        sort_output : bool = True,
        round_output : bool = True,
        remove_empty : bool = True
        ):
        """\
        Initalizes a UCDeconvolvePredictionPropagator object.
        
        Params
        -------
        sort_output
            whether or not to sort final output, defaults to True.
        round_output
            whether or not to round the final predictions. default True
        remove_empty
            remove any empty / zero predictions from the final result.
            default is True.
            
        Returns
        -------
        object
            A propagator object
            
        """
        
        # Load networkX graph directly from compressed dictionary
        self._G = pickle.loads(ucdutils.decompress_b64_gzip(metadata['cellgraph']))
        
        # Save adjaceency matrix and node list
        self._adj = nx.adjacency_matrix(self._G).todense().T
        self._nodes = list(self._G.nodes)
        self.is_fitted = False
        self.sort_output = sort_output
        self.round_output = round_output
        self.remove_empty = remove_empty
        
    def fit(
        self,
        X : Union[np.ndarray, pd.DataFrame],
        y : np.ndarray = None,
        feature_names : List[str] = None):
        """
        Fit dataset
        
        Fit here basically looks at X column names, or the fature name dict
        if X is not a datafarme and creates a matching index I to allow us
        to know which column of the node graph corresponds to the current
        output preictors.
        
        Params
        -------
        X
            array-like to fit data to
        y
            None, not used
        feature_names
            if X is not a dataframe, feature_names are the column names.
            
        Returns
        ------
        self
            itself
            
        """
        
        # Get feature names from appropriate source
        if isinstance(X, pd.DataFrame):
            feature_names = X.columns
        else:
            assert feature_names is not None, "Must pass feature names!"

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            # Warning-causing lines of code here
            self._I = pd.Series(list(feature_names))\
                        .str.replace('b-all','ball')\
                        .str.replace('splenic red pulp macrophagec','splenic red pulp macrophage')\
                        .str.replace('kidney loop of henle thick ascending limb epit...', 
                                     'kidney loop of henle thick ascending limb epithelial cell')\
                        .str.strip().str.lower()\
                        .map(lambda x : dict(zip(pd.Series(self._nodes).str.strip().str.lower().tolist(), 
                                                 range(len(self._nodes)))).get(x, -1)).astype(int).values

        self.is_fitted = True
        
        return self
    
    def transform(
        self,
        X : Union[np.ndarray, pd.DataFrame],
        y : np.ndarray = None,
    ) -> pd.DataFrame:
        """
        Transform / Message Pass
        
        Updates X by extending it. Pass probabilities of different levels
        of cell type prediction up a heierchecal tree structure.
        
        Params
        -------
        X
            initial predictions we want to expand
        y
            not used
        
        Returns
        -------
        X_extended : pd.DataFrame
            extended predictions
        
        """
        # Get logging object
        ucdlogger = logging.getLogger("UCD")
        ucdlogger.debug("Running belief propagation on predictions.")

        # Get extended predictions
        X_extended = pass_each_cell(X.values if isinstance(X, pd.DataFrame) else X, self._adj, self._I)
        
        # Format dataframe
        X_extended = pd.DataFrame(X_extended, columns = self._nodes)
        
        # Round off if desired
        if self.sort_output:
            X_extended = X_extended[X_extended.mean(0).sort_values(ascending = False).index.tolist()]
        
        # Remove empty columns
        if self.remove_empty:
            X_extended = X_extended[X_extended.columns[X_extended.mean(0).gt(0)]]
        
        # Remap indices if they exist
        if isinstance(X, pd.DataFrame):
            X_extended.index = X.index
        
        return X_extended
