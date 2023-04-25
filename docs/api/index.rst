****
API Overview
****

.. module:: ucdeconvolve

Import the ucdeconvolve package using ``import ucdeconvolve as ucd``. 
The package contains four main modules, described in detail below.

.. note::
   Authenticate a new python session using ``ucd.api.authenticate``

**API:** ``api``

API functions allow for user registration, account activation and authentication
for service invocation. 

.. autosummary::
   api.register
   api.activate
   api.authenticate

**Tools:** ``tl``

Tools module contains the three primary prediction functions of ucdeconvolve. 

.. autosummary::
   tl.base
   tl.explain
   tl.select

**Plotting:** ``pl``

Plotting functions for embedding and spatial are designed to interface as wrappers
around scanpy functions such as ``sc.pl.embedding`` and ``sc.pl.spatial`` with
additional functionality to enable construciton of plots similar to those in 
the ucdeconolve paper.

.. autosummary::
   pl.embedding
   pl.spatial
   pl.base_clustermap
   pl.explain_boxplot
   pl.explain_clustermap

**Utilities:** ``utils``

Utilities module contains useful functions for interfacing with results of
deconvolution functions and preparing prediction queries.

.. autosummary::
   utils.read_results
   utils.assign_top_celltypes
   utils.get_base_celltypes
   utils.get_prebuilt_reference
   utils.list_prebuilt_references


.. toctree::
   :maxdepth: 5
   :hidden:

   api/index
   tools/index
   plotting/index
   utilities/index

**Compatability:** ``compat``

The compatability module allows for users of earlier builds of ucdeconvolve
who have existing workflows to continue leveraging legacy code with minimal
required changes. 

.. warning::
   This module will be removed in the near future.

.. autosummary::
   compat.deconvolve
   compat.read_results

