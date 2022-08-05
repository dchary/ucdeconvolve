# UniCell Deconvolve: Cell Type Deconvolution For Transcriptomic Data
![image](https://user-images.githubusercontent.com/7418190/148014006-1c522756-05b2-4e6f-9ab5-dec627405b57.png)
========
<p align="center">
<img alt="Doc Status" src="https://readthedocs.org/projects/ucdeconvolve/badge/?version=latest">
<img alt="Version" src="https://img.shields.io/github/v/release/dchary/ucdeconvolve?include_prereleases">
<img alt="Repo Size" src="https://img.shields.io/github/repo-size/dchary/ucdeconvolve">
<img alt="Last Commit" src="https://img.shields.io/github/last-commit/dchary/ucdeconvolve?style=flat-square">
<img alt="Commit Activity" src="https://img.shields.io/github/commit-activity/w/dchary/ucdeconvolve?style=flat-square">
<img alt="Language" src="https://img.shields.io/github/languages/top/dchary/ucdeconvolve">
</p>

<p>
<figure>
  <img src="https://user-images.githubusercontent.com/7418190/182937919-511c30e1-5865-4af7-9ddf-1c32c1b4c700.png"/>
</figure>
<p align="center"><strong><small><i>UniCell Deconvolve applied to 10X Genomics Visium Gene Expression Slide of Breast Adenocarcinoma Sample</i></small></strong></p>
</p>

## Background
<p align="justify">
The amount of publically available  high-dimensional transcriptomic data, whether bulk-RNA, single-cell, or spatial, has increased exponentially in recent years. Although available for reanalysis, published data is often used in isolation to augment novel analyses. Particularly, the problem of cell type deconvolution, either from bulk or spatial transcriptomic datasets, has been addressed by numerous methods through the use of publicly available dataset as cell type specific references. The choice of reference profile however is not always readily apparent or available, and a mismatch between reference and actual cell types may potentially confound study results.
</p>
<p align="justify">
<strong><em>UniCell Deconvolve (UCD) is a pre-trained deep learning model that provides context-free estimations of cell type fractions</em></strong> from whole transcriptome expression data for bulk, single-cell and spatial transcriptomics data. The model is trained on the world's largest fully-integrated scRNA-Seq training database, comprising 28M+ single cells spanning 840+ cell types from 899 studies to date. Extensive benchmarking shows UCD favors comperably when compared with reference-based deconvolution tools, without the need for pretraining. UCD demonstrates strong multi-task performance across a range of deconvolution challenges spanning several transcriptomic data modalities, disease types, and tissues.
</p>
<p>
<figure>
  <img src="https://user-images.githubusercontent.com/7418190/182940544-45e4c99c-5277-4b51-bbde-0d08afe7e493.png"/>
</figure>
<p align="center"><strong><small><i>Nested rectangles visualizing cell type distribution heiarchy for 28 million single cells comprising the UCD Database to-date</i></small></strong></p>
</p>
</p>

## API Access
<p align="justify">
The UCD package offers the ability to directly integrate UCD predictions into any transcriptomics data analysis pipeline in the form of a web-based API. The package available here provides a secure and scalable connection to the latest pre-trained UCD model, built on top of Google Cloud Platform, which serves deconvolution requests. In order to access the current alpha build of UCD, we ask users to sign up for an early-access API key <a href="https://forms.gle/fhjRev977WAphQf58">here</a>. Please allow up to 24 hours to recieve a response.
</p>

Includes preprocessing and visualization capabilities. Designed to interface with the annotated dataset and scanpy workflows.

## Installation
UniCell Deconvolve can be installed from pyPI:
```bash
pip install ucdeconvolve
```

## Documentation
<p align="justify">
Full documentation with supporting tutorials is available <a href="https://ucdeconvolve.readthedocs.io/">here</a>.
</p>

## Quick Start Guide
<p>
<figure>
<img align="right" src="https://user-images.githubusercontent.com/7418190/182963249-6fde0191-7954-4d06-82cd-6111a1992f87.png">
</figure>

To demonstrate the functionality of UCD, we will perform a cell type deconvolution of a spatial gene expression section of the <a href="https://www.10xgenomics.com/resources/datasets/human-lymph-node-1-standard-1-0-0">human lymph node</a>, made available by 10X Genomics. We will utilize scanpy to quickly load the dataset, and then pass it into ucdeconvolve to obtain cell type predictions.
</p>

### 1. Begin by loading required packages

```python
import ucdeconvolve as ucd
import scanpy as sc
import logging
```
### 2. Load the human lymph node dataset

```python
adata = sc.datasets.visium_sge("V1_Human_Lymph_Node")
```
```
AnnData object with n_obs × n_vars = 4035 × 36601
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial'
    obsm: 'spatial'
```
### 3. Obtain cell type fraction predictions
<p align="justify">
Replace <i>API_KEY</i> with your personal API key. As data is sent for predictions in batches, you can enable a progressbar to show the status of your prediction request. To enable more information about your run, verbosity can be increased using the python logging flag <em>logging.DEBUG</em> which will output the status of each stage of the prediction pipeline as well as post-run metrics. By default logging levels are set to <em>logging.WARNING</em> so leaving this parameter out will result in the function only showing a progressbar.
</p>

```python
ucd.tl.deconvolve(adata, API_KEY, showprogress = True, verbosity = logging.DEBUG)
```
Example Console Output:

```
2022-08-04 20:53:14,436|[UCD]|INFO: Starting UniCell Deconvolve Run.
2022-08-04 20:53:14,622|[UCD]|DEBUG: Using batchsize: 256 on dataset of shape: (4035, 36601)
2022-08-04 20:53:14,625|[UCD]|DEBUG: Data will be sent in 16 packet(s), each split into 4 sub-buffers of ~64 samples.
2022-08-04 20:53:14,629|[UCD]|DEBUG: Data pipeline ready.
2022-08-04 20:53:14,634|[UCD]|INFO: Initiating data stream for prediction.
100% (15 of 15) |########################| Elapsed Time: 0:00:36 Time:  0:00:36
2022-08-04 20:54:16,370|[UCD]|INFO: Data streaming complete.
2022-08-04 20:54:16,372|[UCD]|DEBUG: Streaming time: 61.734 (sec)
2022-08-04 20:54:16,374|[UCD]|DEBUG: Streaming rate: 65.361 (samples / sec)
2022-08-04 20:54:16,379|[UCD]|INFO: Postprocessing predictions.
2022-08-04 20:54:16,429|[UCD]|DEBUG: Splitting predictions.
2022-08-04 20:54:16,673|[UCD]|DEBUG: Running belief propagation on predictions.
2022-08-04 20:54:20,085|[UCD]|DEBUG: Running belief propagation on predictions.
2022-08-04 20:54:21,131|[UCD]|DEBUG: Running belief propagation on predictions.
2022-08-04 20:54:23,387|[UCD]|INFO: Writing results to anndata object.
2022-08-04 20:54:23,390|[UCD]|INFO: UniCell Deconvolve Run Complete.
```
### 5. Reading and Visualizing Results

<p align="justify">
We can print our adata object to see what new information has been added to it. UCD appends the results of each deconvolution run into <i>'adata.obsm'</i> along with column names (i.e. celltypes) and run information into <i>'adata.uns'</i> under the default results stem <i>'ucd_results'</i>. Depending on whether or not the <i>split</i> parameter was set to True or False, you will either see a single new entry into <i>'adata.obsm'</i> or three entries. By default, <i>split = True</i> so predictions will be split into primary (non-malignat), cell lines, and primary cancer (malignant).
</p>

```
AnnData object with n_obs × n_vars = 4035 × 36601
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial', 'ucd_results'
    obsm: 'spatial', 'ucd_results_cancer', 'ucd_results_lines', 'ucd_results_primary'
```

<p align="justify">
We can visualize our results by using one of the built-in plotting functions in UCD, which wrap scanpy's plotting API.
</p>

```python
ucd.pl.deconvolve(adata, basis = "spatial", color = "germinal center b cell")
```

<div align="center">
  <img src="https://user-images.githubusercontent.com/7418190/182959702-9f9806bf-c699-414c-868a-523c2ef88ea6.png">
</div>
<p align="center"><strong><small><i>Predicted germinal center b cell distribution across lymph node section</i></small></strong></p>
</p>

