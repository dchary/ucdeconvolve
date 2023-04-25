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

### Conda (Recommended)
We recommend installing ucdeconvolve in a virtual environment using tools such as conda or miniconda. We suggest the following installation:

```bash
conda create -n ucdenv python=3.8 pytables jupyter jupyterlab
conda activate ucdenv
pip install ucdeconvolve
```

### PIP
UniCell Deconvolve can be installed from pyPI into an existing python workspace. The *pytables* package is required and may need to be installed separately using a package manage such as conda before installing ucdeconvolve. For detailed installation instruction see documentation. 

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

### 1. Create a New Account

#### Register
<p align="justify">
Load the ucdeconvolve package and run the "ucd.api.register()" command as shown below.
Follow the instructions by inputting the required information at each step.
</p>

```python
ucd.api.register()
```

#### Activate
<p align="justify">
Upon completion of the initial registration form, you will recieve an email at
the address specified with an activation code. Copy the code and paste it back
into the waiting input prompt in order to activate your account or paste the
activation code into the function "ucd.api.activate(code)"
</p>

```python
ucd.api.activate(code)
```

#### Authenticate
<p align="justify">
Upon completion of activation, you will recieve an emial with your user acess 
token. This token will be automatically appended to your current python instance
if you are running ucd.api.register, otherwise you can always authenticate a new
python instance with a valid api token using the function "ucd.api.authenticate"
</p>

```python
ucd.api.authenticate(token)
```

### 2. Load Required Packages

```python
import ucdeconvolve as ucd
import scanpy as sc
```

### 3. Load the human lymph node dataset

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
### 4. Run UCDBase to predict cell type fractions

```python
ucd.tl.base(adata)
```
Example Console Output:

```
2023-04-25 16:27:40,012|[UCD]|INFO: Starting UCDeconvolveBASE Run. | Timer Started.
Preprocessing Dataset | 100% (16 of 16) || Elapsed Time: 0:00:02 Time:  0:00:02
2023-04-25 16:27:43,509|[UCD]|INFO: Uploading Data | Timer Started.
2023-04-25 16:27:49,367|[UCD]|INFO: Upload Complete | Elapsed Time: 5.857 (s)
Waiting For Submission : UNKNOWN | Queue Size : 0 | \ |#| 2 Elapsed Time: 0:00:03
Waiting For Completion | 100% (4035 of 4035) || Elapsed Time: 0:00:45 Time:  0:00:45
2023-04-25 16:28:42,073|[UCD]|INFO: Download Results | Timer Started.
2023-04-25 16:28:42,817|[UCD]|INFO: Download Complete | Elapsed Time: 0.743 (s)
2023-04-25 16:28:43,466|[UCD]|INFO: Run Complete | Elapsed Time: 63.453 (s)
```

### 5. Reading and Visualizing Results

<p align="justify">
We can print our adata object to see what new information has been added to it. UCD appends the results of each deconvolution run into <i>'adata.obsm'</i> along with column names (i.e. celltypes) and run information into <i>'adata.uns'</i> under the default results stem <i>'ucdbase'</i>. Depending on whether or not the <i>split</i> parameter was set to True or False, you will either see a single new entry into <i>'adata.obsm'</i> or three entries. By default, <i>split = True</i> so predictions will be split into primary (non-malignat), cell lines, and primary cancer (malignant).
</p>

```
AnnData object with n_obs × n_vars = 4035 × 36601
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial', 'ucdbase'
    obsm: 'spatial', 'ucdbase_cancer', 'ucdbase_lines', 'ucdbase_primary', 'ucdbase_raw'
```

<p align="justify">
We can visualize our results by using one of the built-in plotting functions in UCD, which wrap scanpy's plotting API.
</p>

```python
ucd.pl.spatial(adata, color = "germinal center b cell")
```

<div align="center">
  <img src="https://user-images.githubusercontent.com/7418190/182959702-9f9806bf-c699-414c-868a-523c2ef88ea6.png">
</div>
<p align="center"><strong><small><i>Predicted germinal center b cell distribution across lymph node section</i></small></strong></p>
</p>

