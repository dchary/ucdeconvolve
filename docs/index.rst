|DocStatus| |Version| |RepoSize| |LastCommit| |CommitActivity| |Language|

UniCell Deconvolve: Cell Type Deconvolution For Transcriptomic Data
==============================================

.. |DocStatus| image:: https://readthedocs.org/projects/ucdeconvolve/badge/?version=latest
   :target: https://github.com/dchary/ucdeconvolve
.. |Version| image:: https://img.shields.io/github/v/release/dchary/ucdeconvolve?include_prereleases
   :target: https://github.com/dchary/ucdeconvolve
.. |RepoSize| image:: https://img.shields.io/github/repo-size/dchary/ucdeconvolve
   :target: https://github.com/dchary/ucdeconvolve
.. |LastCommit| image:: https://img.shields.io/github/last-commit/dchary/ucdeconvolve?style=flat-square
   :target: https://github.com/dchary/ucdeconvolve
.. |CommitActivity| image:: https://img.shields.io/github/commit-activity/w/dchary/ucdeconvolve?style=flat-square
   :target: https://github.com/dchary/ucdeconvolve
.. |Language| image:: https://img.shields.io/github/languages/top/dchary/ucdeconvolve
   :target: https://github.com/dchary/ucdeconvolve

UniCell Deconvolve (UCD) is a pre-trained deep learning model that provides context-free estimations of cell type fractions from whole transcriptome expression data for bulk, single-cell and spatial transcriptomics data. The model is trained on the world's largest fully-integrated scRNA-Seq training database, comprising 28M+ single cells spanning 840+ cell types from 899 studies to date. Extensive benchmarking shows UCD favors comperably when compared with reference-based deconvolution tools, without the need for pretraining. UCD demonstrates strong multi-task performance across a range of deconvolution challenges spanning several transcriptomic data modalities, disease types, and tissues.

* Install via ``pip install ucdeconvolve``.
* Discuss development on `GitHub <https://github.com/dchary/ucdeconvolve>`_.

.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials
   api
