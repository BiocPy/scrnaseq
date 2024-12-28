<!-- These are examples of badges you might want to add to your README:
     please update the URLs accordingly

[![Built Status](https://api.cirrus-ci.com/github/<USER>/scrnaseq.svg?branch=main)](https://cirrus-ci.com/github/<USER>/scrnaseq)
[![ReadTheDocs](https://readthedocs.org/projects/scrnaseq/badge/?version=latest)](https://scrnaseq.readthedocs.io/en/stable/)
[![Coveralls](https://img.shields.io/coveralls/github/<USER>/scrnaseq/main.svg)](https://coveralls.io/r/<USER>/scrnaseq)
[![PyPI-Server](https://img.shields.io/pypi/v/scrnaseq.svg)](https://pypi.org/project/scrnaseq/)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/scrnaseq.svg)](https://anaconda.org/conda-forge/scrnaseq)
[![Monthly Downloads](https://pepy.tech/badge/scrnaseq/month)](https://pepy.tech/project/scrnaseq)
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter)](https://twitter.com/scrnaseq)
-->

[![PyPI-Server](https://img.shields.io/pypi/v/scrnaseq.svg)](https://pypi.org/project/scrnaseq/)

# scrnaseq

The `scrnaseq` package provides convenient access to several publicly available single-cell datasets in the form of [SingleCellExperiment](https://github.com/biocpy/singlecellexperiment) objects. Users can obtain a `SingleCellExperiment` and transform it into analysis-ready representations for immediate use.

To enable discovery, each dataset is decorated with metadata such as the study title/abstract, the species, the number of cells, etc. Users can also contribute their own published datasets to enable re-use by the wider Bioconductor/BiocPy community.

**Also check out the R version of this library [here@scRNAseq](https://bioconductor.org/packages/devel/data/experiment/html/scRNAseq.html) published to Bioconductor.**

## Installation

Package is published to [PyPI](https://pypi.org/project/scrnaseq/),

```shell
pip install scrnaseq
```

## Find Datasets

The `list_datasets()` function will display all available datasets along with their metadata. This can be used to discover interesting datasets for further analysis.

```python
import scrnaseq
datasets = scrnaseq.list_datasets()
print(datasets[["name", "version"]].head(3))

## output
# |    | name                  | version    |
# |---:|:----------------------|:-----------|
# |  0 | romanov-brain-2017    | 2023-12-19 |
# |  1 | campbell-brain-2017   | 2023-12-14 |
# |  2 | zhong-prefrontal-2018 | 2023-12-22 |
```

This returns a pandas `DataFrame` to easily filter and download datasets of interest.

Users can also search on the metadata text using the `search_datasets()` function. This accepts both simple text queries as well as more complicated expressions involving boolean operations.

```python
# Find all datasets involving pancreas.
res = scrnaseq.search_datasets("pancreas")

# Find all mm10 datasets involving pancreas or neurons.
from gypsum_client import define_text_query

res = scrnaseq.search_datasets(
     define_text_query("GRCm38", field="genome")
     & (
          define_text_query("neuro%", partial=True)
          | define_text_query("pancrea%", partial=True)
     )
)
print(res[["name", "version"]].head(3))

## output
# |    | name                    | version    |
# |---:|:------------------------|:-----------|
# |  0 | romanov-brain-2017      | 2023-12-19 |
# |  1 | campbell-brain-2017     | 2023-12-14 |
# |  2 | fletcher-olfactory-2019 | 2023-12-21 |
```

Search results are not guaranteed to be reproducible - more datasets may be added over time, and existing datasets may be updated with new versions. Once a dataset of interest is identified, users should explicitly list the name and version of the dataset in their scripts to ensure reproducibility.

## Fetch Datasets

The `fetch_dataset()` function will download a particular dataset, as `SingleCellExperiment`:

```python
sce = scrnaseq.fetch_dataset("zeisel-brain-2015", "2023-12-14")
print(sce)

## output
# class: SingleCellExperiment
# dimensions: (20006, 3005)
# assays(1): ['counts']
# row_data columns(1): ['featureType']
# row_names(20006): ['Tspan12', 'Tshz1', 'Fnbp1l', ..., 'mt-Rnr2', 'mt-Rnr1', 'mt-Nd4l']
# column_data columns(9): ['tissue', 'group #', 'total mRNA mol', 'well', 'sex', 'age', 'diameter', 'level1class', 'level2class']
# column_names(3005): ['1772071015_C02', '1772071017_G12', '1772071017_A05', ..., '1772063068_D01', '1772066098_A12', '1772058148_F03']
# main_experiment_name: gene
# reduced_dims(0): []
# alternative_experiments(2): ['repeat', 'ERCC']
# row_pairs(0): []
# column_pairs(0): []
# metadata(0):
```

For studies that generate multiple datasets, the dataset of interest must be explicitly requested via the `path` argument:

```python
sce = scrnaseq.fetch_dataset("baron-pancreas-2016", "2023-12-14", path="human")
print(sce)

## output
# class: SingleCellExperiment
# dimensions: (20125, 8569)
# assays(1): ['counts']
# row_data columns(0): []
# row_names(20125): ['A1BG', 'A1CF', 'A2M', ..., 'ZZEF1', 'ZZZ3', 'pk']
# column_data columns(2): ['donor', 'label']
# column_names(8569): ['human1_lib1.final_cell_0001', 'human1_lib1.final_cell_0002', 'human1_lib1.final_cell_0003', ..., 'human4_lib3.final_cell_0699', 'human4_lib3.final_cell_0700', 'human4_lib3.final_cell_0701']
# main_experiment_name:
# reduced_dims(0): []
# alternative_experiments(0): []
# row_pairs(0): []
# column_pairs(0): []
# metadata(0):
```

By default, array data is loaded as a file-backed `DelayedArray` from the [HDF5Array](https://github.com/BiocPy/HDF5Array) package. Setting `realize_assays=True` and/or `realize_reduced_dims=True` will coerce file-backed arrays to numpy or scipy sparse (csr/csc) objects.

```python
sce = scrnaseq.fetch_dataset("baron-pancreas-2016", "2023-12-14", path="human", realize_assays=True)
print(sce)

## output
# class: SingleCellExperiment
# dimensions: (20125, 8569)
# assays(1): ['counts']
# row_data columns(0): []
# row_names(20125): ['A1BG', 'A1CF', 'A2M', ..., 'ZZEF1', 'ZZZ3', 'pk']
# column_data columns(2): ['donor', 'label']
# column_names(8569): ['human1_lib1.final_cell_0001', 'human1_lib1.final_cell_0002', 'human1_lib1.final_cell_0003', ..., 'human4_lib3.final_cell_0699', 'human4_lib3.final_cell_0700', 'human4_lib3.final_cell_0701']
# main_experiment_name:
# reduced_dims(0): []
# alternative_experiments(0): []
# row_pairs(0): []
# column_pairs(0): []
# metadata(0):
```

Users can also fetch the metadata associated with each dataset:

```python
meta = scrnaseq.fetch_metadata("zeisel-brain-2015", "2023-12-14")
print(meta.keys())

## output
# dict_keys(['title', 'description', 'taxonomy_id', 'genome', 'sources', 'maintainer_name', 'maintainer_email', 'bioconductor_version', 'applications'])
```


## Adding New Datasets

Want to contribute your own dataset to this package? It's easy! Just follow these simple steps:

1. Format your dataset as a `SummarizedExperiment` or `SingleCellExperiment`. Let's mock a dataset:

     ```python
     import numpy as np
     from singlecellexperiment import SingleCellExperiment
     from biocframe import BiocFrame

     mat = np.random.poisson(1, (100, 10))
     row_names = [f"GENE_{i}" for i in range(mat.shape[0])]
     col_names = list("ABCDEFGHIJ")
     sce = SingleCellExperiment(
          assays={"counts": mat},
          row_data=BiocFrame(row_names=row_names),
          column_data=BiocFrame(row_names=col_names),
     )
     ```

     Or if you want to save an `AnnData` object, easily convert that into a `SingleCellExperiment` object:

     ```python
     from singlecellexperiment import SingleCellExperiment

     sce = SingleCellExperiment.from_anndata(YOUR_ANNDATA_OBJ)
     ```

2. Assemble the metadata for your dataset. This should be a dictionary as specified in the [Bioconductor metadata schema](https://github.com/ArtifactDB/bioconductor-metadata-index). Check out some examples from `fetch_metadata()`  Note that the `application.takane` property will be automatically added later, and so can be omitted from the list that you create.

     ```python
     meta = {
          "title": "My dataset forked from ziesel brain",
          "description": "This is a copy of the ziesel",
          "taxonomy_id": ["10090"],  # NCBI ID
          "genome": ["GRCh38"],  # genome build
          "sources": [{"provider": "GEO", "id": "GSE12345"}],
          "maintainer_name": "Shizuka Mogami",
          "maintainer_email": "mogami.shizuka@765pro.com",
     }
     ```

3. Save your `SummarizedExperiment` or `SingleCellExperiment` object to disk with `save_dataset()`. This saves the dataset into a "staging directory" using language-agnostic file formats - check out the [ArtifactDB](https://github.com/artifactdb) framework for more details. In more complex cases involving multiple datasets, users may save each dataset into a subdirectory of the staging directory.

     ```python
     import tempfile
     from scrnaseq import save_dataset

     # replace tmp with a staging directory
     staging_dir = tempfile.mkdtemp()
     save_dataset(sce, staging_dir, meta)
     ```

     You can check that everything was correctly saved by reloading the on-disk data for inspection:

     ```python
     import dolomite_base as dl

     dl.read_object(staging_dir)
     ```

4. Open a [pull request (PR)](https://github.com/BiocPy/scRNAseq/pulls) for the addition of a new dataset. You will need to provide a few things here:
   - The name of your dataset. This typically follows the format of `{NAME}-{SYSTEM}-{YEAR}`, where `NAME` is the last name of the first author of the study, `SYSTEM` is the biological system (e.g., tissue, cell types) being studied, and `YEAR` is the year of publication for the dataset.
   - The version of your dataset. This is usually just the current date or whenever you started putting together the dataset for upload. The exact date doesn't really matter as long as we can establish a timeline for later versions.
   - An Python file containing the code used to assemble the dataset. This should be added to the [`scripts/`](https://github.com/BiocPy/scRNAseq/tree/master/scripts) directory of this package, in order to provide some record of how the dataset was created.

5. Wait for us to grant temporary upload permissions to your GitHub account.

6. Upload your staging directory to [**gypsum** backend](https://github.com/ArtifactDB/gypsum-worker) with `upload_dataset()`. On the first call to this function, it will automatically prompt you to log into GitHub so that the backend can authenticate you. If you are on a system without browser access (e.g., most computing clusters), a [token](https://github.com/settings/tokens) can be manually supplied via `set_access_token()`.

     ```python
     from scrnaseq import upload_dataset

     upload_dataset(staging_dir, "my_dataset_name", "my_version")
     ```

     You can check that everything was successfully uploaded by calling `fetch_dataset()` with the same name and version:

     ```python
     from scrnaseq import upload_dataset

     fetch_dataset("my_dataset_name", "my_version")
     ```

     If you realized you made a mistake, no worries. Use the following call to clear the erroneous dataset, and try again:

     ```python
     from gypsum_client import reject_probation

     reject_probation("scRNAseq", "my_dataset_name", "my_version")
     ```

7. Comment on the PR to notify us that the dataset has finished uploading and you're happy with it. We'll review it and make sure everything's in order. If some fixes are required, we'll just clear the dataset so that you can upload a new version with the necessary changes. Otherwise, we'll approve the dataset. Note that once a version of a dataset is approved, no further changes can be made to that version; you'll have to upload a new version if you want to modify something.

<!-- pyscaffold-notes -->

## Note

The tests for upload are skipped. To run them, include a GitHub token as an environment variable `gh_token`.

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
